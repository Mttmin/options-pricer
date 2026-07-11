use crate::ctmc::{price_american_put_heston, CTMCResult};
use nalgebra::{SMatrix, SVector};
use num_complex::Complex64;
use options::{Greeks, Options};
use std::cmp::Ordering;
use std::collections::BTreeMap;

pub const HESTON_DIM: usize = 5;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HestonParameters {
	pub v0: f64,
	pub v_bar: f64,
	pub kappa: f64,
	pub sigma: f64,
	pub rho: f64,
}

impl HestonParameters {
	pub fn to_vector(self) -> SVector<f64, HESTON_DIM> {
		SVector::from([self.v0, self.v_bar, self.kappa, self.sigma, self.rho])
	}

	pub fn from_vector(v: &SVector<f64, HESTON_DIM>) -> Self {
		Self {
			v0: v[0],
			v_bar: v[1],
			kappa: v[2],
			sigma: v[3],
			rho: v[4],
		}
	}

	pub fn feller_lhs(self) -> f64 {
		2.0 * self.kappa * self.v_bar - self.sigma * self.sigma
	}
}

#[derive(Debug, Clone, Copy)]
pub struct HestonBounds {
	pub v0: (f64, f64),
	pub v_bar: (f64, f64),
	pub kappa: (f64, f64),
	pub sigma: (f64, f64),
	pub rho: (f64, f64),
}

impl Default for HestonBounds {
	fn default() -> Self {
		Self {
			v0: (1e-5, 2.0),
			v_bar: (1e-5, 2.0),
			kappa: (1e-4, 30.0),
			sigma: (1e-4, 5.0),
			rho: (-0.999, 0.999),
		}
	}
}

impl HestonBounds {
	pub fn clamp(self, params: HestonParameters) -> HestonParameters {
		HestonParameters {
			v0: params.v0.clamp(self.v0.0, self.v0.1),
			v_bar: params.v_bar.clamp(self.v_bar.0, self.v_bar.1),
			kappa: params.kappa.clamp(self.kappa.0, self.kappa.1),
			sigma: params.sigma.clamp(self.sigma.0, self.sigma.1),
			rho: params.rho.clamp(self.rho.0, self.rho.1),
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum OptionSide {
	Call,
	Put,
}

#[derive(Debug, Clone)]
pub struct CalibrationInstrument {
	pub symbol: String,
	pub side: OptionSide,
	pub strike: f64,
	pub maturity_years: f64,
	pub market_price: f64,
	pub spot: f64,
	pub risk_free_rate: f64,
	pub dividend_yield: f64,
	pub implied_vol: Option<f64>,
	pub bid: Option<f64>,
	pub ask: Option<f64>,
	pub open_interest: Option<u64>,
	pub volume: Option<u64>,
}

impl CalibrationInstrument {
	pub fn moneyness(&self) -> f64 {
		self.strike / self.spot
	}

	pub fn spread(&self) -> Option<f64> {
		self.bid.zip(self.ask).map(|(b, a)| (a - b).max(0.0))
	}

	pub fn spread_ratio(&self) -> Option<f64> {
		let spread = self.spread()?;
		if self.market_price > 1e-12 {
			Some(spread / self.market_price)
		} else {
			None
		}
	}

	pub fn base_option(&self, vol: f64) -> Options {
		match self.side {
			OptionSide::Call => Options::new_call(
				self.strike,
				self.spot,
				vol,
				self.risk_free_rate,
				self.maturity_years,
				Some(self.dividend_yield),
			),
			OptionSide::Put => Options::new_put(
				self.strike,
				self.spot,
				vol,
				self.risk_free_rate,
				self.maturity_years,
				Some(self.dividend_yield),
			),
		}
	}
}

#[derive(Debug, Clone)]
pub struct CalibrationDataset {
	pub valuation_symbol: String,
	pub as_of_epoch_secs: i64,
	pub instruments: Vec<CalibrationInstrument>,
}

#[derive(Debug, Clone, Copy)]
pub struct CalibrationConfig {
	pub max_spread_ratio: f64,
	pub min_open_interest: u64,
	pub min_volume: u64,
	pub moneyness_min: f64,
	pub moneyness_max: f64,
	pub min_maturity_years: f64,
	pub spread_weight_floor: f64,
	pub spread_weight_scale: f64,
	pub oi_volume_weight_scale: f64,
	pub maturity_balance_strength: f64,
	pub vega_floor: f64,
	pub lambda_feller: f64,
	pub lambda_anchor: f64,
	pub anchor: Option<HestonParameters>,
	pub max_iters: usize,
	pub grad_epsilon: f64,
	pub step_tolerance: f64,
	pub objective_tolerance: f64,
	pub initial_damping: f64,
	pub damping_up: f64,
	pub damping_down: f64,
	pub bounds: HestonBounds,
}

impl Default for CalibrationConfig {
	fn default() -> Self {
		Self {
			max_spread_ratio: 0.25,
			min_open_interest: 50,
			min_volume: 10,
			moneyness_min: 0.7,
			moneyness_max: 1.3,
			min_maturity_years: 1.0 / 365.0,
			spread_weight_floor: 0.05,
			spread_weight_scale: 1.0,
			oi_volume_weight_scale: 1.0,
			maturity_balance_strength: 0.5,
			vega_floor: 1e-4,
			lambda_feller: 10.0,
			lambda_anchor: 0.0,
			anchor: None,
			max_iters: 75,
			grad_epsilon: 1e-4,
			step_tolerance: 1e-6,
			objective_tolerance: 1e-10,
			initial_damping: 1e-2,
			damping_up: 10.0,
			damping_down: 0.33,
			bounds: HestonBounds::default(),
		}
	}
}

#[derive(Debug, Clone)]
pub struct CalibrationReport {
	pub input_count: usize,
	pub kept_count: usize,
	pub dropped_count: usize,
	pub violations: Vec<String>,
}

#[derive(Debug, Clone, Copy)]
struct WeightingParams {
	pub spread_weight_floor: f64,
	pub spread_weight_scale: f64,
	pub oi_volume_weight_scale: f64,
	pub maturity_balance_strength: f64,
	pub vega_floor: f64,
}


impl From<&CalibrationConfig> for WeightingParams {
	fn from(cfg: &CalibrationConfig) -> Self {
		Self {
			spread_weight_floor: cfg.spread_weight_floor,
			spread_weight_scale: cfg.spread_weight_scale,
			oi_volume_weight_scale: cfg.oi_volume_weight_scale,
			maturity_balance_strength: cfg.maturity_balance_strength,
			vega_floor: cfg.vega_floor,
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjectiveKind {
	VegaWeightedPrice,
	PriceLeastSquares,
	ImpliedVolLeastSquares,
}

#[derive(Debug, Clone)]
pub struct CalibrationResult {
	pub parameters: HestonParameters,
	pub converged: bool,
	pub iterations: usize,
	pub objective: f64,
	pub rmse: f64,
	pub feller_lhs: f64,
	pub report: CalibrationReport,
}

#[derive(Debug)]
pub enum CalibrationError {
	EmptyDataset,
	MissingPricingSupport(String),
	NumericalFailure(String),
}

pub trait CalibrationPricer {
	fn price(
		&self,
		instrument: &CalibrationInstrument,
		params: HestonParameters,
	) -> Result<f64, CalibrationError>;
}

#[derive(Debug, Clone, Copy)]
pub struct CtmcPricingConfig {
	pub n_x: usize,
	pub m_v: usize,
	pub n_time: usize,
}

impl Default for CtmcPricingConfig {
	fn default() -> Self {
		Self {
			n_x: 80,
			m_v: 20,
			n_time: 50,
		}
	}
}

pub struct CtmcAmericanPutPricer {
	pub config: CtmcPricingConfig,
}

impl CtmcAmericanPutPricer {
	pub fn new(config: CtmcPricingConfig) -> Self {
		Self { config }
	}
}

impl CalibrationPricer for CtmcAmericanPutPricer {
	fn price(
		&self,
		instrument: &CalibrationInstrument,
		params: HestonParameters,
	) -> Result<f64, CalibrationError> {
		if instrument.side != OptionSide::Put {
			return Err(CalibrationError::MissingPricingSupport(
				"CTMC skeleton pricer currently supports puts only; provide a call pricer adapter for calls".to_string(),
			));
		}
		let result: CTMCResult = price_american_put_heston(
			instrument.strike,
			instrument.spot,
			params.v0,
			instrument.maturity_years,
			instrument.risk_free_rate,
			instrument.dividend_yield,
			params.kappa,
			params.v_bar,
			params.sigma,
			params.rho,
			self.config.n_x,
			self.config.m_v,
			self.config.n_time,
		);
		Ok(result.price)
	}
}

/// Heston characteristic function for X = ln(S_T/S_0) under Q.
///
/// Uses Formulation 2 (Lord & Kahl 2010) with g = (ξ−d)/(ξ+d), which keeps the
/// argument of the complex logarithm away from the branch cut for all u and τ.
fn heston_cf(u: f64, params: &HestonParameters, tau: f64, r: f64, q: f64) -> Complex64 {
	let v0 = params.v0;
	let v_bar = params.v_bar;
	let kappa = params.kappa;
	let sigma = params.sigma;
	let rho = params.rho;
	let sigma_sq = sigma * sigma;

	// ξ = κ − ρσiu  (real part κ, imaginary part −ρσu)
	let xi = Complex64::new(kappa, -rho * sigma * u);

	// d = sqrt(ξ² + σ²(u² + iu)) — treat u²+iu as Complex64::new(u²,u)
	let d = (xi * xi + Complex64::new(sigma_sq * u * u, sigma_sq * u)).sqrt();

	// Formulation 2: g = (ξ−d)/(ξ+d)
	let g = (xi - d) / (xi + d);

	let one = Complex64::new(1.0, 0.0);
	let exp_neg_dtau = (-d * tau).exp();
	let denom = one - g * exp_neg_dtau;

	// B(u,τ) = (ξ−d)(1−e^{−dτ}) / (σ²·denom)
	let b = (xi - d) * (one - exp_neg_dtau) / (sigma_sq * denom);

	// A(u,τ) = (κv̄/σ²)·[(ξ−d)τ − 2·ln(denom/(1−g))]
	let a = (kappa * v_bar / sigma_sq) * ((xi - d) * tau - 2.0 * (denom / (one - g)).ln());

	// φ(u) = exp(iu·(r−q)·τ + A + B·v₀)
	(Complex64::new(0.0, u) * (r - q) * tau + a + b * v0).exp()
}

/// Integral ∫_c^d e^x cos(ω(x−a)) dx used in COS payoff coefficients.
fn chi_k(omega: f64, c: f64, d: f64, a: f64) -> f64 {
	if omega == 0.0 {
		return d.exp() - c.exp();
	}
	let cos_d = (omega * (d - a)).cos();
	let sin_d = (omega * (d - a)).sin();
	let cos_c = (omega * (c - a)).cos();
	let sin_c = (omega * (c - a)).sin();
	(d.exp() * (cos_d + omega * sin_d) - c.exp() * (cos_c + omega * sin_c))
		/ (1.0 + omega * omega)
}

/// Integral ∫_c^d cos(ω(x−a)) dx used in COS payoff coefficients.
fn psi_k(omega: f64, c: f64, d: f64, a: f64) -> f64 {
	if omega == 0.0 {
		return d - c;
	}
	((omega * (d - a)).sin() - (omega * (c - a)).sin()) / omega
}

/// Price a European call or put under the Heston model using the COS method
/// (Fang & Oosterlee 2008). Accuracy is ~10⁻⁸ for n_terms = 128.
pub(crate) fn cos_european(
	side: OptionSide,
	s0: f64,
	strike: f64,
	r: f64,
	q: f64,
	tau: f64,
	params: &HestonParameters,
	n_terms: usize,
) -> f64 {
	let eps = (-params.kappa * tau).exp();
	let mean_reversion_factor = (1.0 - eps) / params.kappa;
	let integrated_var =
		params.v_bar * tau + (params.v0 - params.v_bar) * mean_reversion_factor;

	// First cumulant of log-return X = ln(S_T/S_0)
	let c1 = (r - q - 0.5 * params.v_bar) * tau
		+ 0.5 * (params.v0 - params.v_bar) * mean_reversion_factor;

	// Truncation window centred at 0, wide enough to cover BOTH the payoff kink
	// (y = 0) and the density of y = ln(S_T/K), centred near x0 + c1. A window
	// centred on x0 + c1 alone can exclude y = 0 for short-T deep-OTM cells,
	// which inverts the put/call integration bounds (c_lo > c_hi) below and
	// silently corrupts the price. Mirrors the Phase 2A-v2 fix in
	// deep-calibration `training data creation/v2/pricer.py`.
	let x0 = (s0 / strike).ln();
	let margin = 12.0 * integrated_var.abs().sqrt();
	let half = ((x0 + c1).abs() + margin).max(0.5);
	let a = -half;
	let b = half;
	let ba = b - a;
	let pi_ba = std::f64::consts::PI / ba;

	// Integration bounds for the payoff: calls from max(0,a) to b, puts from a to min(0,b)
	let (c_lo, c_hi) = match side {
		OptionSide::Call => (a.max(0.0), b),
		OptionSide::Put => (a, b.min(0.0)),
	};

	let mut price = 0.0;
	for k in 0..n_terms {
		let omega = k as f64 * pi_ba;
		let phi = heston_cf(omega, params, tau, r, q);
		// Phase shift: e^{iω(x₀−a)} aligns x₀ = ln(S₀/K) with the CF of ln(S_T/S₀)
		let phase = Complex64::from_polar(1.0, omega * (x0 - a));
		let re_phi = (phi * phase).re;

		let vk = 2.0 / ba * match side {
			OptionSide::Call => chi_k(omega, c_lo, c_hi, a) - psi_k(omega, c_lo, c_hi, a),
			OptionSide::Put => -chi_k(omega, c_lo, c_hi, a) + psi_k(omega, c_lo, c_hi, a),
		};

		let coeff = if k == 0 { 0.5 } else { 1.0 };
		price += coeff * re_phi * vk;
	}

	(strike * (-r * tau).exp() * price).max(0.0)
}

/// European option pricer for the Heston model using the COS Fourier method.
///
/// Use this inside `calibrate_heston` to fit Heston parameters to market data.
/// Once calibrated, pass the resulting `HestonParameters` to the CTMC engine
/// (`price_american_put_heston`) for American or path-dependent pricing.
pub struct HestonEuropeanPricer {
	pub n_cos_terms: usize,
}

impl Default for HestonEuropeanPricer {
	fn default() -> Self {
		Self { n_cos_terms: 128 }
	}
}

impl CalibrationPricer for HestonEuropeanPricer {
	fn price(
		&self,
		instrument: &CalibrationInstrument,
		params: HestonParameters,
	) -> Result<f64, CalibrationError> {
		let p = cos_european(
			instrument.side,
			instrument.spot,
			instrument.strike,
			instrument.risk_free_rate,
			instrument.dividend_yield,
			instrument.maturity_years,
			&params,
			self.n_cos_terms,
		);
		if p.is_finite() {
			Ok(p)
		} else {
			Err(CalibrationError::NumericalFailure(
				"COS price is non-finite".to_string(),
			))
		}
	}
}

pub fn clean_dataset(
	dataset: &CalibrationDataset,
	config: &CalibrationConfig,
) -> (CalibrationDataset, CalibrationReport) {
	let mut kept = Vec::with_capacity(dataset.instruments.len());
	for instrument in &dataset.instruments {
		if instrument.market_price <= 0.0 {
			continue;
		}
		if instrument.maturity_years < config.min_maturity_years {
			continue;
		}
		let m = instrument.moneyness();
		if m < config.moneyness_min || m > config.moneyness_max {
			continue;
		}
		if let Some(spread_ratio) = instrument.spread_ratio()
			&& spread_ratio > config.max_spread_ratio
		{
			continue;
		}
		if let Some(oi) = instrument.open_interest
			&& oi < config.min_open_interest
		{
			continue;
		}
		if let Some(vol) = instrument.volume
			&& vol < config.min_volume
		{
			continue;
		}
		kept.push(instrument.clone());
	}

	let filtered = CalibrationDataset {
		valuation_symbol: dataset.valuation_symbol.clone(),
		as_of_epoch_secs: dataset.as_of_epoch_secs,
		instruments: kept,
	};

	let violations = no_arbitrage_violations(&filtered);
	let report = CalibrationReport {
		input_count: dataset.instruments.len(),
		kept_count: filtered.instruments.len(),
		dropped_count: dataset.instruments.len().saturating_sub(filtered.instruments.len()),
		violations,
	};

	(filtered, report)
}

pub fn no_arbitrage_violations(dataset: &CalibrationDataset) -> Vec<String> {
	let mut violations = Vec::new();
	violations.extend(check_calendar_variance(dataset));
	violations.extend(check_strike_monotonicity(dataset));
	violations.extend(check_butterfly_convexity(dataset));
	violations
}

fn check_calendar_variance(dataset: &CalibrationDataset) -> Vec<String> {
	let mut violations = Vec::new();
	let mut by_key: BTreeMap<(OptionSide, i64), Vec<&CalibrationInstrument>> = BTreeMap::new();
	for instrument in &dataset.instruments {
		if let Some(iv) = instrument.implied_vol {
			let strike_key = (instrument.strike * 10_000.0).round() as i64;
			by_key
				.entry((instrument.side, strike_key))
				.or_default()
				.push(instrument);
			if iv <= 0.0 {
				violations.push(format!(
					"calendar_variance: non-positive implied vol for {} strike {:.4}, T={:.6}",
					instrument.symbol, instrument.strike, instrument.maturity_years
				));
			}
		}
	}

	for ((_side, strike_key), mut quotes) in by_key {
		quotes.sort_by(|a, b| {
			a.maturity_years
				.partial_cmp(&b.maturity_years)
				.unwrap_or(Ordering::Equal)
		});
		let mut last_w: Option<f64> = None;
		for quote in quotes {
			if let Some(iv) = quote.implied_vol {
				let w = iv * iv * quote.maturity_years;
				if let Some(prev) = last_w
					&& w + 1e-10 < prev
				{
					violations.push(format!(
						"calendar_variance: total variance decreases at strike key {} around T={:.6}",
						strike_key, quote.maturity_years
					));
				}
				last_w = Some(w);
			}
		}
	}
	violations
}

fn check_strike_monotonicity(dataset: &CalibrationDataset) -> Vec<String> {
	let mut violations = Vec::new();
	let mut by_slice: BTreeMap<(OptionSide, i64), Vec<&CalibrationInstrument>> = BTreeMap::new();
	for instrument in &dataset.instruments {
		let t_key = (instrument.maturity_years * 3650.0).round() as i64;
		by_slice
			.entry((instrument.side, t_key))
			.or_default()
			.push(instrument);
	}

	for ((side, t_key), mut quotes) in by_slice {
		quotes.sort_by(|a, b| a.strike.partial_cmp(&b.strike).unwrap_or(Ordering::Equal));
		let mut prev: Option<f64> = None;
		for quote in quotes {
			if let Some(last) = prev {
				match side {
					OptionSide::Call if quote.market_price > last + 1e-10 => {
						violations.push(format!(
							"strike_monotonicity: call slice not decreasing for maturity key {}",
							t_key
						));
						break;
					}
					OptionSide::Put if quote.market_price + 1e-10 < last => {
						violations.push(format!(
							"strike_monotonicity: put slice not increasing for maturity key {}",
							t_key
						));
						break;
					}
					_ => {}
				}
			}
			prev = Some(quote.market_price);
		}
	}
	violations
}

fn check_butterfly_convexity(dataset: &CalibrationDataset) -> Vec<String> {
	let mut violations = Vec::new();
	let mut by_slice: BTreeMap<(OptionSide, i64), Vec<&CalibrationInstrument>> = BTreeMap::new();
	for instrument in &dataset.instruments {
		let t_key = (instrument.maturity_years * 3650.0).round() as i64;
		by_slice
			.entry((instrument.side, t_key))
			.or_default()
			.push(instrument);
	}

	for ((side, t_key), mut quotes) in by_slice {
		quotes.sort_by(|a, b| a.strike.partial_cmp(&b.strike).unwrap_or(Ordering::Equal));
		if quotes.len() < 3 {
			continue;
		}
		for idx in 1..(quotes.len() - 1) {
			let k0 = quotes[idx - 1].strike;
			let k1 = quotes[idx].strike;
			let k2 = quotes[idx + 1].strike;
			let c0 = quotes[idx - 1].market_price;
			let c1 = quotes[idx].market_price;
			let c2 = quotes[idx + 1].market_price;
			let left = (c1 - c0) / (k1 - k0).max(1e-12);
			let right = (c2 - c1) / (k2 - k1).max(1e-12);
			if right + 1e-8 < left {
				violations.push(format!(
					"butterfly_convexity: {:?} slice violates convexity at maturity key {}",
					side, t_key
				));
				break;
			}
		}
	}
	violations
}

fn maturity_counts(dataset: &CalibrationDataset) -> BTreeMap<i64, usize> {
	let mut counts: BTreeMap<i64, usize> = BTreeMap::new();
	for quote in &dataset.instruments {
		let key = (quote.maturity_years * 3650.0).round() as i64;
		*counts.entry(key).or_insert(0) += 1;
	}
	counts
}

fn instrument_weight(
	instrument: &CalibrationInstrument,
	weights: WeightingParams,
	maturity_count: usize,
) -> f64 {
	let spread_component = instrument
		.spread_ratio()
		.map(|x| 1.0 / (1.0 + weights.spread_weight_scale * x))
		.unwrap_or(1.0)
		.max(weights.spread_weight_floor);

	let liq = (instrument.open_interest.unwrap_or(0) + instrument.volume.unwrap_or(0)) as f64;
	let liq_component = (1.0 + weights.oi_volume_weight_scale * liq).ln().max(1e-6);

	let maturity_component = (maturity_count as f64)
		.powf(-weights.maturity_balance_strength)
		.max(1e-6);

	spread_component * liq_component * maturity_component
}

struct PreparedCalibration<'a> {
	instruments: Vec<&'a CalibrationInstrument>,
	sqrt_weights: Vec<f64>,
}

fn prepare_calibration<'a>(
	dataset: &'a CalibrationDataset,
	weights: WeightingParams,
) -> Result<PreparedCalibration<'a>, CalibrationError> {
	if dataset.instruments.is_empty() {
		return Err(CalibrationError::EmptyDataset);
	}

	let maturity_count_map = maturity_counts(dataset);
	let mut instruments = Vec::with_capacity(dataset.instruments.len());
	let mut sqrt_weights = Vec::with_capacity(dataset.instruments.len());

	for instrument in &dataset.instruments {
		let key = (instrument.maturity_years * 3650.0).round() as i64;
		let maturity_count = *maturity_count_map.get(&key).unwrap_or(&1);
		instruments.push(instrument);
		sqrt_weights.push(instrument_weight(instrument, weights, maturity_count).sqrt());
	}

	Ok(PreparedCalibration {
		instruments,
		sqrt_weights,
	})
}

fn feller_penalty(params: HestonParameters, config: &CalibrationConfig) -> f64 {
	let violation = (params.sigma * params.sigma - 2.0 * params.kappa * params.v_bar).max(0.0);
	let mut penalty = config.lambda_feller * violation * violation;
	if let Some(anchor) = config.anchor {
		let dv = params.to_vector() - anchor.to_vector();
		penalty += config.lambda_anchor * dv.dot(&dv);
	}
	penalty
}

fn implied_vol_from_price(instrument: &CalibrationInstrument, target: f64) -> Option<f64> {
	if target <= 0.0 {
		return None;
	}
	let mut lo = 1e-4;
	let mut hi = 5.0;

	let price_with_vol = |vol: f64| instrument.base_option(vol).bs_pricing();
	let mut plo = price_with_vol(lo) - target;
	let mut phi = price_with_vol(hi) - target;

	if plo * phi > 0.0 {
		for _ in 0..8 {
			hi *= 1.5;
			phi = price_with_vol(hi) - target;
			if plo * phi <= 0.0 {
				break;
			}
		}
	}
	if plo * phi > 0.0 {
		return None;
	}

	for _ in 0..80 {
		let mid = 0.5 * (lo + hi);
		let pmid = price_with_vol(mid) - target;
		if pmid.abs() < 1e-8 {
			return Some(mid);
		}
		if plo * pmid <= 0.0 {
			hi = mid;
		} else {
			lo = mid;
			plo = pmid;
		}
	}
	Some(0.5 * (lo + hi))
}

fn residual_for_instrument(
	instrument: &CalibrationInstrument,
	model_price: f64,
	objective: ObjectiveKind,
	weights: WeightingParams,
	sqrt_weight: f64,
) -> Option<f64> {
	let raw = match objective {
		ObjectiveKind::PriceLeastSquares => model_price - instrument.market_price,
		ObjectiveKind::VegaWeightedPrice => {
			let vol = instrument.implied_vol.unwrap_or(0.20);
			let opt = instrument.base_option(vol);
			let vega = opt.vega(instrument.spot).max(weights.vega_floor);
			(model_price - instrument.market_price) / vega
		}
		ObjectiveKind::ImpliedVolLeastSquares => {
			let market_iv = instrument
				.implied_vol
				.or_else(|| implied_vol_from_price(instrument, instrument.market_price))?;
			let model_iv = implied_vol_from_price(instrument, model_price)?;
			model_iv - market_iv
		}
	};
	Some(sqrt_weight * raw)
}

fn objective_with_pricer<P: CalibrationPricer>(
	pricer: &P,
	prepared: &PreparedCalibration,
	params: HestonParameters,
	objective: ObjectiveKind,
	weights: WeightingParams,
	config: &CalibrationConfig,
) -> Result<(f64, Vec<f64>), CalibrationError> {
	if prepared.instruments.is_empty() {
		return Err(CalibrationError::EmptyDataset);
	}
	let mut residuals = Vec::with_capacity(prepared.instruments.len());

	for (instrument, sqrt_weight) in prepared.instruments.iter().zip(&prepared.sqrt_weights) {
		let model_price = pricer.price(instrument, params)?;
		if let Some(residual) = residual_for_instrument(
			instrument,
			model_price,
			objective,
			weights,
			*sqrt_weight,
		) {
			residuals.push(residual);
		}
	}

	if residuals.is_empty() {
		return Err(CalibrationError::NumericalFailure(
			"All residuals dropped after objective transform".to_string(),
		));
	}

	let objective_value = residuals.iter().map(|r| r * r).sum::<f64>()
		+ feller_penalty(params, config);
	Ok((objective_value, residuals))
}

fn normal_equations_from_finite_diff<P: CalibrationPricer>(
	pricer: &P,
	prepared: &PreparedCalibration,
	params: HestonParameters,
	base_residuals: &[f64],
	objective: ObjectiveKind,
	weights: WeightingParams,
	config: &CalibrationConfig,
) -> Result<(SMatrix<f64, HESTON_DIM, HESTON_DIM>, SVector<f64, HESTON_DIM>), CalibrationError> {
	let mut jt_j = SMatrix::<f64, HESTON_DIM, HESTON_DIM>::zeros();
	let mut jt_r = SVector::<f64, HESTON_DIM>::zeros();
	let base_vec = params.to_vector();
	let mut jac_cols: Vec<Vec<f64>> = Vec::with_capacity(HESTON_DIM);

	for j in 0..HESTON_DIM {
		let mut bumped = base_vec;
		let eps_j = config.grad_epsilon.max(1e-8 * base_vec[j].abs());
		bumped[j] += eps_j;
		let bumped_params = config.bounds.clamp(HestonParameters::from_vector(&bumped));
		let (_, residuals) = objective_with_pricer(
			pricer,
			prepared,
			bumped_params,
			objective,
			weights,
			config,
		)?;

		if residuals.len() != base_residuals.len() {
			return Err(CalibrationError::NumericalFailure(
				"Residual dimensions changed during finite differences".to_string(),
			));
		}

		let col: Vec<f64> = residuals
			.iter()
			.zip(base_residuals.iter())
			.map(|(&r1, &r0)| (r1 - r0) / eps_j)
			.collect();
		jac_cols.push(col);
	}

	for i in 0..HESTON_DIM {
		for j in 0..HESTON_DIM {
			let mut acc = 0.0;
			for k in 0..base_residuals.len() {
				acc += jac_cols[i][k] * jac_cols[j][k];
			}
			jt_j[(i, j)] = acc;
		}
		let mut rhs = 0.0;
		for k in 0..base_residuals.len() {
			rhs += jac_cols[i][k] * base_residuals[k];
		}
		jt_r[i] = rhs;
	}
	Ok((jt_j, jt_r))
}

pub fn calibrate_heston<P: CalibrationPricer>(
	pricer: &P,
	dataset: &CalibrationDataset,
	initial: HestonParameters,
	objective: ObjectiveKind,
	config: &CalibrationConfig,
	report: CalibrationReport,
) -> Result<CalibrationResult, CalibrationError> {
	let weights = WeightingParams::from(config);
	let prepared = prepare_calibration(dataset, weights)?;

	let mut params = config.bounds.clamp(initial);
	let mut damping = config.initial_damping;
	let mut iter_count = 0usize;
	let mut converged = false;

	let (mut objective_value, mut residuals) = objective_with_pricer(
		pricer,
		&prepared,
		params,
		objective,
		weights,
		config,
	)?;

	for iter in 0..config.max_iters {
		iter_count = iter + 1;
		let (mut h, jt_r) = normal_equations_from_finite_diff(
			pricer,
			&prepared,
			params,
			&residuals,
			objective,
			weights,
			config,
		)?;

		for i in 0..HESTON_DIM {
			h[(i, i)] += damping;
		}

		let step = h.lu().solve(&(-jt_r)).ok_or_else(|| {
			CalibrationError::NumericalFailure("Unable to solve damped normal equations".to_string())
		})?;

		let step_norm = step.norm();
		if step_norm < config.step_tolerance {
			converged = true;
			break;
		}

		let candidate = config
			.bounds
			.clamp(HestonParameters::from_vector(&(params.to_vector() + step)));
		let (candidate_obj, candidate_residuals) = objective_with_pricer(
			pricer,
			&prepared,
			candidate,
			objective,
			weights,
			config,
		)?;

		let accepted = candidate_obj + config.objective_tolerance < objective_value;
		if accepted {
			params = candidate;
			objective_value = candidate_obj;
			residuals = candidate_residuals;
			damping = (damping * config.damping_down).max(1e-10);
		} else {
			damping = (damping * config.damping_up).min(1e10);
		}
	}

	let rmse_price = (residuals.iter().map(|r| r * r).sum::<f64>() / residuals.len() as f64).sqrt();
	Ok(CalibrationResult {
		parameters: params,
		converged,
		iterations: iter_count,
		objective: objective_value,
		rmse: rmse_price,
		feller_lhs: params.feller_lhs(),
		report,
	})
}

#[cfg(test)]
mod tests {
	use super::*;

	struct ToyPricer;
	impl CalibrationPricer for ToyPricer {
		fn price(&self, ins: &CalibrationInstrument, p: HestonParameters) -> Result<f64, CalibrationError> {
			Ok(0.5 * ins.strike + 2.0 * p.v0 + 3.0 * p.v_bar + 0.1 * p.kappa + 0.4 * p.sigma - 0.2 * p.rho)
		}
	}

	fn toy_dataset() -> CalibrationDataset {
		let instruments = [90.0, 95.0, 100.0, 105.0, 110.0_f64].iter().enumerate().map(|(i, &k)| {
			CalibrationInstrument {
				symbol: format!("OPT-{i}"),
				side: OptionSide::Put, strike: k, maturity_years: 0.5,
				market_price: 0.0, spot: 100.0, risk_free_rate: 0.03, dividend_yield: 0.01,
				implied_vol: Some(0.2), bid: Some(1.0), ask: Some(1.2),
				open_interest: Some(500), volume: Some(250),
			}
		}).collect();
		CalibrationDataset { valuation_symbol: "TEST".to_string(), as_of_epoch_secs: 0, instruments }
	}

	#[test]
	fn clean_dataset_drops_wide_spread_and_wings() {
		let mut ds = toy_dataset();
		ds.instruments[0].ask = Some(9.0);
		ds.instruments[4].strike = 180.0;
		let (filtered, report) = clean_dataset(&ds, &CalibrationConfig::default());
		assert!(filtered.instruments.len() < ds.instruments.len());
		assert_eq!(report.kept_count, filtered.instruments.len());
	}

	#[test]
	fn calibrate_heston_lm_smoke() {
		let pricer = ToyPricer;
		let true_p = HestonParameters { v0: 0.06, v_bar: 0.05, kappa: 1.5, sigma: 0.4, rho: -0.3 };
		let mut ds = toy_dataset();
		for q in &mut ds.instruments { q.market_price = pricer.price(q, true_p).unwrap(); }
		let result = calibrate_heston(
			&pricer, &ds,
			HestonParameters { v0: 0.04, v_bar: 0.04, kappa: 2.0, sigma: 0.3, rho: 0.0 },
			ObjectiveKind::PriceLeastSquares,
			&CalibrationConfig { max_iters: 30, ..CalibrationConfig::default() },
			CalibrationReport { input_count: ds.instruments.len(), kept_count: ds.instruments.len(), dropped_count: 0, violations: vec![] },
		).unwrap();
		assert!(result.objective >= 0.0);
		assert!(result.parameters.v0 > 0.0);
	}

	#[test]
	fn cos_pricer_sanity() {
		use options::Options;
		let p = HestonParameters { v0: 0.04, v_bar: 0.04, kappa: 2.0, sigma: 0.3, rho: -0.5 };
		// CF normalization: φ(0) = 1
		let phi0 = heston_cf(0.0, &p, 1.0, 0.05, 0.02);
		assert!((phi0.re - 1.0).abs() < 1e-12);
		assert!(phi0.im.abs() < 1e-12);
		// Put-call parity: C - P = S·e^{-qT} - K·e^{-rT}
		let (s0, k, r, q, tau) = (100.0, 100.0, 0.05, 0.02, 0.5);
		let call = cos_european(OptionSide::Call, s0, k, r, q, tau, &p, 256);
		let put  = cos_european(OptionSide::Put,  s0, k, r, q, tau, &p, 256);
		assert!((call - put - (s0 * (-q * tau).exp() - k * (-r * tau).exp())).abs() < 1e-6);
		// BS limit: sigma->0, v0=v_bar => matches Black-Scholes within 0.5%
		let bs_params = HestonParameters { v0: 0.04, v_bar: 0.04, kappa: 10.0, sigma: 1e-6, rho: 0.0 };
		let heston_c = cos_european(OptionSide::Call, 100.0, 100.0, 0.05, 0.0, 1.0, &bs_params, 256);
		let bs_c = Options::new_call(100.0, 100.0, 0.2, 0.05, 1.0, None).bs_pricing();
		assert!((heston_c - bs_c).abs() / bs_c < 0.005);
	}

	#[test]
	fn calibrate_heston_with_fourier_pricer() {
		let pricer = HestonEuropeanPricer::default();
		let true_p = HestonParameters { v0: 0.04, v_bar: 0.04, kappa: 2.0, sigma: 0.4, rho: -0.6 };
		let mut instruments = Vec::new();
		for &strike in &[90.0, 95.0, 100.0, 105.0, 110.0_f64] {
			for &tau in &[0.25, 1.0_f64] {
				for &side in &[OptionSide::Call, OptionSide::Put] {
					let mut ins = CalibrationInstrument {
						symbol: format!("{side:?}-{strike}"), side, strike, maturity_years: tau,
						market_price: 0.0, spot: 100.0, risk_free_rate: 0.05, dividend_yield: 0.01,
						implied_vol: Some(0.2), bid: Some(1.0), ask: Some(1.2),
						open_interest: Some(500), volume: Some(200),
					};
					ins.market_price = pricer.price(&ins, true_p).unwrap();
					instruments.push(ins);
				}
			}
		}
		let ds = CalibrationDataset { valuation_symbol: "TEST".to_string(), as_of_epoch_secs: 0, instruments };
		let result = calibrate_heston(
			&pricer, &ds,
			HestonParameters { v0: 0.06, v_bar: 0.05, kappa: 1.5, sigma: 0.3, rho: -0.4 },
			ObjectiveKind::VegaWeightedPrice,
			&CalibrationConfig { max_iters: 50, ..CalibrationConfig::default() },
			CalibrationReport { input_count: ds.instruments.len(), kept_count: ds.instruments.len(), dropped_count: 0, violations: vec![] },
		).unwrap();
		assert!(result.objective < 1e-4, "residual: {}", result.objective);
	}

	/// Full end-to-end pipeline:
	/// 1. Calibrate Heston parameters from a synthetic European surface.
	/// 2. Feed calibrated parameters into the CTMC American put pricer.
	/// 3. Verify the price is positive, finite, and close to the true-parameter price.
	///
	/// The American ≥ European parity at coarse grid sizes is not asserted here
	/// because ATM put early-exercise premium (q=0.02 < r=0.05) can be smaller
	/// than CTMC grid discretization error at nm=900. See
	/// `american_exceeds_european_high_res` for that assertion at high resolution.
	#[test]
	#[ignore]
	fn full_pipeline_calibrate_then_ctmc_american_put() {
		use crate::ctmc::price_american_put_heston;

		let true_params = HestonParameters {
			v0: 0.04,
			v_bar: 0.04,
			kappa: 2.0,
			sigma: 0.35,
			rho: -0.65,
		};
		let (spot, r, q) = (100.0_f64, 0.05, 0.02);
		let pricer = HestonEuropeanPricer::default();

		// Synthetic surface: 6 strikes × 3 maturities × both sides
		let mut instruments = Vec::new();
		for &strike in &[85.0, 90.0, 95.0, 100.0, 105.0, 110.0] {
			for &tau in &[0.25, 0.5, 1.0] {
				for &side in &[OptionSide::Call, OptionSide::Put] {
					let mut ins = CalibrationInstrument {
						symbol: format!("{side:?}-{strike}-{tau}"),
						side,
						strike,
						maturity_years: tau,
						market_price: 0.0,
						spot,
						risk_free_rate: r,
						dividend_yield: q,
						implied_vol: Some(0.2),
						bid: Some(1.0),
						ask: Some(1.2),
						open_interest: Some(500),
						volume: Some(200),
					};
					ins.market_price = pricer.price(&ins, true_params).unwrap();
					instruments.push(ins);
				}
			}
		}
		let ds = CalibrationDataset {
			valuation_symbol: "PIPELINE-TEST".to_string(),
			as_of_epoch_secs: 0,
			instruments,
		};

		// Step 1 — calibrate Heston parameters
		let init = HestonParameters {
			v0: 0.05,
			v_bar: 0.05,
			kappa: 1.5,
			sigma: 0.3,
			rho: -0.5,
		};
		let cal = calibrate_heston(
			&pricer,
			&ds,
			init,
			ObjectiveKind::VegaWeightedPrice,
			&CalibrationConfig { max_iters: 75, ..CalibrationConfig::default() },
			CalibrationReport {
				input_count: ds.instruments.len(),
				kept_count: ds.instruments.len(),
				dropped_count: 0,
				violations: vec![],
			},
		)
		.unwrap();
		assert!(cal.objective < 1e-3, "calibration residual: {}", cal.objective);
		let p = cal.parameters;

		// Step 2 — price American put via CTMC with calibrated parameters
		let (strike, maturity) = (100.0_f64, 0.5_f64);
		let ctmc_result = price_american_put_heston(
			strike, spot, p.v0, maturity, r, q,
			p.kappa, p.v_bar, p.sigma, p.rho,
			80, 20, 50,
		);
		let american_put = ctmc_result.price;

		// Reference European put via COS at true parameters
		let european_ref =
			cos_european(OptionSide::Put, spot, strike, r, q, maturity, &true_params, 256);

		assert!(american_put > 0.0, "American put must be positive, got {american_put}");
		assert!(american_put.is_finite(), "American put must be finite");
		// American put should be within 10% of the European reference (calibration + grid noise)
		let rel_diff = (american_put - european_ref).abs() / european_ref;
		assert!(rel_diff < 0.10, "American {american_put:.4} vs European ref {european_ref:.4}: {:.1}%", rel_diff * 100.0);
	}

	/// Verify American put ≥ European put (no-arbitrage) for the same Heston parameters
	/// at a resolution where CTMC discretization error is well below the early-exercise premium.
	/// Uses a deep ITM put (S=80, K=100) where the premium is large and unambiguous.
	#[test]
	#[ignore]
	fn american_exceeds_european_high_res() {
		use crate::ctmc::price_american_put_heston;

		let params = HestonParameters {
			v0: 0.04,
			v_bar: 0.04,
			kappa: 2.0,
			sigma: 0.35,
			rho: -0.65,
		};
		// Deep ITM put: S=80, K=100, r=0.05, q=0.02 → large early exercise premium
		let (spot, strike, r, q, maturity) = (80.0_f64, 100.0_f64, 0.05, 0.02, 0.5);

		let european_put =
			cos_european(OptionSide::Put, spot, strike, r, q, maturity, &params, 256);
		let ctmc_result = price_american_put_heston(
			strike, spot, params.v0, maturity, r, q,
			params.kappa, params.v_bar, params.sigma, params.rho,
			120, 25, 75,
		);
		let american_put = ctmc_result.price;

		assert!(american_put > european_put,
			"American {american_put:.4} should exceed European {european_put:.4} for deep ITM put"
		);
	}
}
