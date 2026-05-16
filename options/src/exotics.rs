use crate::{BlackScholesPrice, Call, MonteCarloParameters, Options, optionspreads::{OptionSpreads, Position}};
use statrs::distribution::{ContinuousCDF, Normal};

pub enum ExoticOptions {
    ConvertibleBond(ConvertibleBond),
    AsianOption(AsianOption),
    CliquetOption(CliquetOption),
}

#[derive(Debug, Clone, Copy)]
pub struct ConvertibleBond {
    // Bond parameters
    pub face_value: f64,
    pub coupon_rate: f64,
    pub maturity: f64,
    pub payment_frequency: u32,
    pub credit_spread: f64,

    // Conversion option parameters
    pub conversion_price: f64,
    pub stock_price: f64,
    pub volatility: f64,
    pub time_to_maturity: f64,
    pub risk_free_rate: f64,
    pub dividend_yield: Option<f64>,
}

impl ConvertibleBond {
    fn npv(&self) -> f64 {
        // NPV calculation logic
        let periods = (self.maturity * self.payment_frequency as f64) as u32;
        let coupon_payment = self.face_value * self.coupon_rate / self.payment_frequency as f64;
        let mut npv = 0.0;
        let mut discounter = 1.0;
        for _ in 1..=periods {
            discounter *=
                1.0 + (self.risk_free_rate + self.credit_spread) / self.payment_frequency as f64;
            npv += coupon_payment / discounter;
        }
        npv += self.face_value / discounter;
        npv
    }
    fn conversion_option_price(&self) -> f64 {
        // Use Black-Scholes to price the conversion option
        let underlying_call = Call {
            strike_price: self.conversion_price,
            spot_price: self.stock_price,
            volatility: self.volatility,
            risk_free_rate: self.risk_free_rate,
            time_to_maturity: self.time_to_maturity,
            dividend_yield: self.dividend_yield,
        };
        underlying_call.bs_pricing()
    }
}

impl BlackScholesPrice for ConvertibleBond {
    /// Calculate the total price of the convertible bond using Black-Scholes for the conversion option and NPV for the bond component
    fn bs_pricing(&self) -> f64 {
        self.npv() + self.face_value / self.conversion_price * self.conversion_option_price()
    }
}

pub struct ChooserOption {
    underlying: OptionSpreads,
}

impl ChooserOption{
    pub fn new(
        spot_price: f64,
        strike_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
        choose_time: f64,
    ) -> Self {
        let call: Options = Options::new_call(strike_price, spot_price, volatility, risk_free_rate, time_to_maturity, dividend_yield);
        let put: Options = Options::new_put(strike_price*(-(risk_free_rate-dividend_yield.unwrap_or(0.0))*(time_to_maturity - choose_time)).exp(), spot_price, volatility, risk_free_rate, choose_time, dividend_yield);
        let portfolio: Vec<Position> = vec![
            Position {
                option: call,
                weight: 1.0,
            },
            Position {
                option: put,
                weight: ((-dividend_yield.unwrap_or(0.0)*(time_to_maturity-choose_time)).exp()) as f32,
            },
        ];
        let underlying = OptionSpreads::new(portfolio);
        Self { underlying }
    }
}

impl BlackScholesPrice for ChooserOption {
    fn bs_pricing(&self) -> f64 {
        self.underlying.bs_pricing()
    }
}

// =====================================================================
// Asian options: path-dependent payoff aggregated over the life of the
// option. Two pooling modes:
//   - Average: arithmetic mean of monitored spots (Asian).
//   - Max: extremum of monitored spots — path max for call,
//          path min for put (fixed-strike lookback).
//
// Pricing: closed-form analytics (Kemna-Vorst geometric Asian /
// Conze-Viswanathan lookback) live here; path-MC pricing lives in
// the `numerical-methods` crate.
// =====================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsianKind {
    Call,
    Put,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsianPooling {
    Average,
    Max,
}

#[derive(Debug, Clone, Copy)]
pub struct AsianOption {
    pub kind: AsianKind,
    pub pooling: AsianPooling,
    pub strike_price: f64,
    pub spot_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    pub num_observations: u32,
}

impl AsianOption {
    pub fn new(
        kind: AsianKind,
        pooling: AsianPooling,
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
        num_observations: u32,
    ) -> Self {
        Self {
            kind,
            pooling,
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
            num_observations: num_observations.max(1),
        }
    }

    pub fn dividend(&self) -> f64 {
        self.dividend_yield.unwrap_or(0.0)
    }

    pub fn closed_form_price(&self) -> f64 {
        match self.pooling {
            AsianPooling::Average => kemna_vorst_geometric_price(self),
            AsianPooling::Max => conze_lookback_fixed_strike_price(self),
        }
    }
}

impl BlackScholesPrice for AsianOption {
    fn bs_pricing(&self) -> f64 {
        self.closed_form_price()
    }
}

impl MonteCarloParameters for AsianOption {
    fn spot_price(&self) -> f64 { self.spot_price }
    fn volatility(&self) -> f64 { self.volatility }
    fn risk_free_rate(&self) -> f64 { self.risk_free_rate }
    fn time_to_maturity(&self) -> f64 { self.time_to_maturity }
}

/// Kemna-Vorst closed-form for a discretely-monitored geometric-average
/// Asian option. Acts as a tight surrogate for the arithmetic Asian and
/// as a control variate in MC. Monitoring dates t_i = iT/N, i=1..N.
pub fn kemna_vorst_geometric_price(opt: &AsianOption) -> f64 {
    let n = opt.num_observations as f64;
    let t = opt.time_to_maturity;
    let r = opt.risk_free_rate;
    let q = opt.dividend();
    let sigma = opt.volatility;
    let s0 = opt.spot_price;
    let k = opt.strike_price;

    // Moments of ln(G) under risk-neutral measure for equally-spaced obs.
    let drift = r - q - 0.5 * sigma * sigma;
    let mu_lng = s0.ln() + drift * (n + 1.0) * t / (2.0 * n);
    let var_lng = sigma * sigma * t * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
    let sigma_avg = (var_lng / t).sqrt();
    let sd_lng = var_lng.sqrt();

    // E[G] under risk-neutral, then BS-style with F = E[G].
    let forward = (mu_lng + 0.5 * var_lng).exp();

    let d1 = ((forward / k).ln() + 0.5 * sigma_avg * sigma_avg * t) / (sigma_avg * t.sqrt());
    let d2 = d1 - sd_lng;

    let std_norm = Normal::new(0.0, 1.0).unwrap();
    let disc = (-r * t).exp();

    match opt.kind {
        AsianKind::Call => disc * (forward * std_norm.cdf(d1) - k * std_norm.cdf(d2)),
        AsianKind::Put  => disc * (k * std_norm.cdf(-d2) - forward * std_norm.cdf(-d1)),
    }
}

/// Conze-Viswanathan / Goldman-Sosin-Gatto fixed-strike lookback option
/// under GBM with continuous monitoring. At-issue pricing assumes the
/// running extremum equals the spot. For calls the path-max is used;
/// for puts the path-min. b = r - q (cost of carry); clamped away from
/// zero to avoid σ²/(2b) blow-up.
pub fn conze_lookback_fixed_strike_price(opt: &AsianOption) -> f64 {
    let s = opt.spot_price;
    let k = opt.strike_price;
    let t = opt.time_to_maturity;
    let r = opt.risk_free_rate;
    let q = opt.dividend();
    let sigma = opt.volatility;
    let std_norm = Normal::new(0.0, 1.0).unwrap();

    let mut b = r - q;
    if b.abs() < 1e-6 { b = if b >= 0.0 { 1e-6 } else { -1e-6 }; }
    let sqrt_t = t.sqrt();
    let sigma_sqrt_t = sigma * sqrt_t;
    let disc_r = (-r * t).exp();
    let disc_q = (-q * t).exp();

    match opt.kind {
        AsianKind::Call => {
            // At issue the running max equals S, so effective strike = max(K, S).
            let m = k.max(s);
            let d1 = ((s / m).ln() + (b + 0.5 * sigma * sigma) * t) / sigma_sqrt_t;
            let d2 = d1 - sigma_sqrt_t;
            let extra_arg = d1 - 2.0 * b * sqrt_t / sigma;
            let power = (s / m).powf(-2.0 * b / (sigma * sigma));
            let core = s * disc_q * std_norm.cdf(d1) - m * disc_r * std_norm.cdf(d2);
            let lookback_term = s * disc_r * (sigma * sigma / (2.0 * b))
                * (-power * std_norm.cdf(extra_arg) + (b * t).exp() * std_norm.cdf(d1));
            // If K < S the contract already has intrinsic on the M-K = S-K piece.
            let intrinsic = disc_r * (m - k);
            core + lookback_term + intrinsic
        }
        AsianKind::Put => {
            // At issue running min equals S, so effective strike = min(K, S).
            let m = k.min(s);
            let d1 = ((s / m).ln() + (b + 0.5 * sigma * sigma) * t) / sigma_sqrt_t;
            let d2 = d1 - sigma_sqrt_t;
            let extra_arg = -d1 + 2.0 * b * sqrt_t / sigma;
            let power = (s / m).powf(-2.0 * b / (sigma * sigma));
            let core = m * disc_r * std_norm.cdf(-d2) - s * disc_q * std_norm.cdf(-d1);
            let lookback_term = s * disc_r * (sigma * sigma / (2.0 * b))
                * (power * std_norm.cdf(extra_arg) - (b * t).exp() * std_norm.cdf(-d1));
            let intrinsic = disc_r * (k - m);
            core + lookback_term + intrinsic
        }
    }
}

// =====================================================================
// Cliquet (ratchet) options: a strip of consecutive forward-start
// options. The underlying is monitored at N reset dates; each period
// contributes a return r_i = S_{t_{i+1}}/S_{t_i} - 1. Per-period returns
// are clamped to [local_floor, local_cap], summed, and the sum is
// clamped to [global_floor, global_cap]. Calls use +r_i; puts use -r_i.
//
// Pricing: the analytic closed form here values the UNCAPPED strip with
// each leg floored at zero as a sum of forward-start ATM Black-Scholes
// options (Rubinstein, 1991), in fractional-return units. It is a
// benchmark only; the path-MC in the numerical-methods crate handles
// arbitrary caps/floors exactly.
// =====================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CliquetKind {
    Call,
    Put,
}

#[derive(Debug, Clone, Copy)]
pub struct CliquetOption {
    pub kind: CliquetKind,
    pub spot_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    pub num_resets: u32,
    pub local_cap: Option<f64>,
    pub local_floor: Option<f64>,
    pub global_cap: Option<f64>,
    pub global_floor: Option<f64>,
}

impl CliquetOption {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        kind: CliquetKind,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
        num_resets: u32,
        local_cap: Option<f64>,
        local_floor: Option<f64>,
        global_cap: Option<f64>,
        global_floor: Option<f64>,
    ) -> Self {
        Self {
            kind,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
            num_resets: num_resets.max(1),
            local_cap,
            local_floor,
            global_cap,
            global_floor,
        }
    }

    pub fn dividend(&self) -> f64 {
        self.dividend_yield.unwrap_or(0.0)
    }

    /// Analytic benchmark: value of the UNCAPPED strip with each leg
    /// floored at zero (sum of forward-start ATM Black-Scholes options).
    /// Caps/floors are NOT reflected here — use `monte_carlo_cliquet`
    /// for the exact capped/floored price.
    pub fn closed_form_price(&self) -> f64 {
        forward_start_strip_price(self)
    }
}

impl BlackScholesPrice for CliquetOption {
    fn bs_pricing(&self) -> f64 {
        self.closed_form_price()
    }
}

impl MonteCarloParameters for CliquetOption {
    fn spot_price(&self) -> f64 { self.spot_price }
    fn volatility(&self) -> f64 { self.volatility }
    fn risk_free_rate(&self) -> f64 { self.risk_free_rate }
    fn time_to_maturity(&self) -> f64 { self.time_to_maturity }
}

/// Sum of N consecutive forward-start ATM options (Rubinstein, 1991),
/// in fractional-return units, settled at maturity to match the MC
/// payoff: r_i = S_{t_{i+1}}/S_{t_i} - 1 floored at zero, accumulated
/// and discounted once by e^{-rT}.
///
/// For equally spaced resets each period has length tau = T/N. By the
/// scale invariance of GBM the per-period return is i.i.d.; the
/// undiscounted expectation of one floored leg (S = K = 1) is the
/// forward-style BS amount with forward F = e^{(r-q)·tau}:
///
///   d1 = ((r - q + 0.5σ²)·tau) / (σ·sqrt(tau)),  d2 = d1 - σ·sqrt(tau)
///   leg = F·N(d1) - N(d2)            (call)
///         N(-d2) - F·N(-d1)          (put)
///   price = N · e^{-r·T} · leg       (fractional units; S0-independent)
pub fn forward_start_strip_price(opt: &CliquetOption) -> f64 {
    let n = opt.num_resets as f64;
    let t = opt.time_to_maturity;
    let r = opt.risk_free_rate;
    let q = opt.dividend();
    let sigma = opt.volatility;
    let tau = t / n;

    let std_norm = Normal::new(0.0, 1.0).unwrap();
    let sqrt_tau = tau.sqrt();
    let d1 = ((r - q + 0.5 * sigma * sigma) * tau) / (sigma * sqrt_tau);
    let d2 = d1 - sigma * sqrt_tau;
    let fwd = ((r - q) * tau).exp();

    let leg = match opt.kind {
        CliquetKind::Call => fwd * std_norm.cdf(d1) - std_norm.cdf(d2),
        CliquetKind::Put => std_norm.cdf(-d2) - fwd * std_norm.cdf(-d1),
    };

    n * (-r * t).exp() * leg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convertible_bond_pricing() {
        let cb = ConvertibleBond {
            face_value: 1000.0,
            coupon_rate: 0.05,
            maturity: 5.0,
            payment_frequency: 2,
            credit_spread: 0.02,
            risk_free_rate: 0.03,
            conversion_price: 50.0,
            stock_price: 55.0,
            volatility: 0.2,
            time_to_maturity: 5.0,
            dividend_yield: None,
        };
        let price = cb.bs_pricing();
        println!("Convertible Bond Price: {:.4}", price);
        assert!((price - 1_318.0).abs() < 1e-1); // expected value
    }

    #[test]
    fn test_kemna_vorst_below_vanilla() {
        // Geometric Asian must be cheaper than vanilla European (vol averaging
        // reduces effective variance).
        let s = 100.0;
        let k = 100.0;
        let sigma = 0.3;
        let r = 0.05;
        let t = 1.0;
        let asian = AsianOption::new(
            AsianKind::Call, AsianPooling::Average,
            k, s, sigma, r, t, None, 252,
        );
        let asian_price = asian.bs_pricing();
        let vanilla = Options::new_call(k, s, sigma, r, t, None).bs_pricing();
        assert!(asian_price > 0.0 && asian_price < vanilla);
    }

    #[test]
    fn test_kemna_vorst_put_call_parity() {
        // C - P = e^{-rT}·(F - K) where F is the forward of the geometric mean.
        let s = 100.0;
        let k = 100.0;
        let sigma = 0.2;
        let r = 0.05;
        let t = 1.0;
        let n = 50u32;
        let call = AsianOption::new(
            AsianKind::Call, AsianPooling::Average, k, s, sigma, r, t, None, n,
        ).bs_pricing();
        let put = AsianOption::new(
            AsianKind::Put, AsianPooling::Average, k, s, sigma, r, t, None, n,
        ).bs_pricing();
        // Recompute F from same moments used inside Kemna-Vorst.
        let nf = n as f64;
        let drift = r - 0.5 * sigma * sigma;
        let mu = s.ln() + drift * (nf + 1.0) * t / (2.0 * nf);
        let var = sigma * sigma * t * (nf + 1.0) * (2.0 * nf + 1.0) / (6.0 * nf * nf);
        let f = (mu + 0.5 * var).exp();
        let lhs = call - put;
        let rhs = (-r * t).exp() * (f - k);
        assert!((lhs - rhs).abs() < 1e-8, "parity violation: {} vs {}", lhs, rhs);
    }

    #[test]
    fn test_lookback_above_vanilla() {
        // Fixed-strike lookback call ≥ vanilla European call (path max ≥ S_T).
        let s = 100.0;
        let k = 100.0;
        let sigma = 0.3;
        let r = 0.05;
        let t = 1.0;
        let lookback = AsianOption::new(
            AsianKind::Call, AsianPooling::Max,
            k, s, sigma, r, t, None, 252,
        ).bs_pricing();
        let vanilla = Options::new_call(k, s, sigma, r, t, None).bs_pricing();
        assert!(lookback > vanilla);
    }

    #[test]
    fn test_lookback_put_above_vanilla() {
        let s = 100.0;
        let k = 100.0;
        let sigma = 0.3;
        let r = 0.05;
        let t = 1.0;
        let lookback = AsianOption::new(
            AsianKind::Put, AsianPooling::Max,
            k, s, sigma, r, t, None, 252,
        ).bs_pricing();
        let vanilla = Options::new_put(k, s, sigma, r, t, None).bs_pricing();
        assert!(lookback > vanilla);
    }

    #[test]
    fn test_cliquet_strip_positive_and_finite() {
        let c = CliquetOption::new(
            CliquetKind::Call, 100.0, 0.3, 0.05, 1.0, None, 4,
            None, Some(0.0), None, None,
        );
        let p = c.bs_pricing();
        assert!(p > 0.0 && p.is_finite());
    }

    #[test]
    fn test_cliquet_single_reset_matches_vanilla_atm() {
        // N=1, t_0=0 → V_0 = bs_unit(tau=T) = fractional ATM BS call
        // = vanilla ATM BS call value / S0 (scale invariance).
        let s = 100.0;
        let c = CliquetOption::new(
            CliquetKind::Call, s, 0.2, 0.05, 1.0, None, 1,
            None, Some(0.0), None, None,
        );
        let vanilla = Options::new_call(100.0, s, 0.2, 0.05, 1.0, None).bs_pricing();
        assert!((c.bs_pricing() - vanilla / s).abs() < 1e-9);
    }

    #[test]
    fn test_cliquet_put_below_call_zero_div() {
        let call = CliquetOption::new(
            CliquetKind::Call, 100.0, 0.25, 0.05, 1.0, None, 4,
            None, Some(0.0), None, None,
        ).bs_pricing();
        let put = CliquetOption::new(
            CliquetKind::Put, 100.0, 0.25, 0.05, 1.0, None, 4,
            None, Some(0.0), None, None,
        ).bs_pricing();
        assert!(call > put && put > 0.0);
    }
}
