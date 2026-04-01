//! Deep calibration via the BatesSurrogate ONNX model.
//!
//! Loads the exported `bates_surrogate.onnx` (+ companion `.onnx.data`) and
//! calibrates Heston parameters from a market IV surface using Nelder-Mead
//! through the frozen surrogate (Option A: jump params stripped, pure Heston
//! output).
//!
//! # Grid
//! The surrogate was trained on the `high41x14` preset:
//!   - NK = 49 strikes,  log-moneyness ∈ [−0.80, +0.40]
//!   - NT = 14 maturities (weeks → years, see `MATURITIES`)
//!   - N_FLAT = 686  (row-major: strike-first, then maturity)
//!
//! # Parameter order (normalised [0,1])
//! kappa, theta, sigma_v, rho, v0  (pure Heston)
//!
//! # Mapping to HestonParameters
//! kappa→kappa, theta→v_bar, sigma_v→sigma, rho→rho, v0→v0.
//!
//! Market rates r (risk-free) and q (dividend yield) are applied separately,
//! not part of the calibration search. They should be known from market data.

use ort::{
    session::{builder::GraphOptimizationLevel, Session},
    value::Tensor,
};

use crate::slv::{CalibrationDataset, HestonParameters};

// ---------------------------------------------------------------------------
// Model path resolution
// ---------------------------------------------------------------------------

/// Resolve the ONNX model path, checking locations in priority order:
///
/// 1. `<workspace>/weights/bates_surrogate.onnx`  — drop weights here for
///    local use (symlinks to the training project work fine).
/// 2. `/home/mttmin/coding/deep-calibration/model/bates_surrogate.onnx`
///    — fallback to the training project location on the dev machine.
///
/// Returns `None` if neither path exists.
pub fn default_onnx_path() -> Option<String> {
    // Path 1: workspace weights/ directory (sibling of this crate)
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../weights/bates_surrogate.onnx");
    if std::path::Path::new(workspace).exists() {
        return Some(workspace.to_string());
    }
    // Path 2: training project on the dev machine
    let dev = "/home/mttmin/coding/deep-calibration/model/bates_surrogate.onnx";
    if std::path::Path::new(dev).exists() {
        return Some(dev.to_string());
    }
    None
}

// ---------------------------------------------------------------------------
// Grid constants  (must match heston_datagen.py high41x14 preset)
// ---------------------------------------------------------------------------
//
// Systematic IV upward bias note:
//   The surrogate is a Bates model where μ_j ≤ 0 (downward jumps only).
//   Unconstrained calibration inflates Heston params (especially v0) and uses
//   negative jumps as an offset, so extracting Heston-only always overshoots.
//   Fix: pin λ_j = 0 so the surrogate degenerates to pure Heston during search.

pub const NK: usize = 49;
pub const NT: usize = 14;
pub const N_FLAT: usize = NK * NT; // 686
pub const N_PARAMS: usize = 5;     // Pure Heston: kappa, theta, sigma, rho, v0

/// Maturity grid in years (14 points).
pub const MATURITIES: [f32; NT] = [
    1.0 / 52.0,
    2.0 / 52.0,
    3.0 / 52.0,
    1.0 / 12.0,
    2.0 / 12.0,
    3.0 / 12.0,
    4.0 / 12.0,
    6.0 / 12.0,
    9.0 / 12.0,
    1.0,
    1.25,
    1.5,
    2.0,
    3.0,
];

/// Log-moneyness grid: linspace(−0.80, +0.40, 49).
pub fn log_moneyness_grid() -> [f32; NK] {
    let mut g = [0f32; NK];
    let lo = -0.80_f32;
    let hi = 0.40_f32;
    for i in 0..NK {
        g[i] = lo + (hi - lo) * i as f32 / (NK - 1) as f32;
    }
    g
}

// ---------------------------------------------------------------------------
// Parameter bounds  (physical units, matches calibrate.py PARAM_BOUNDS_PHYSICAL)
// ---------------------------------------------------------------------------

/// Physical lower bounds: [kappa, theta, sigma_v, rho, v0]
/// Tightened to prevent degenerate solutions where sigma maxes out.
const PARAM_LO: [f32; N_PARAMS] = [0.10, 0.02, 0.10, -0.85, 0.02];
/// Physical upper bounds: [kappa, theta, sigma_v, rho, v0]
const PARAM_HI: [f32; N_PARAMS] = [5.00, 0.15, 1.00, 0.05, 0.15];

fn denormalise(norm: &[f32; N_PARAMS]) -> [f32; N_PARAMS] {
    let mut out = [0f32; N_PARAMS];
    for i in 0..N_PARAMS {
        out[i] = norm[i] * (PARAM_HI[i] - PARAM_LO[i]) + PARAM_LO[i];
    }
    out
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Controls the calibration search over pure Heston parameters.
///
/// All 5 Heston parameters (kappa, theta, sigma, rho, v0) are freely optimized.
/// Market rates (r, q) are applied separately during pricing, not calibrated.
#[derive(Debug, Clone)]
pub struct CalibrationConfig {
    /// Gaussian ATM weighting: each strike gets weight = exp(−lm² / (2σ²)).
    /// Near-ATM options have tighter bid-ask spreads and better price discovery.
    /// `None` = uniform weights. Typical value: 0.30.
    pub moneyness_weight_sigma: Option<f32>,

    // Deprecated fields kept for backward compatibility (no-op)
    #[deprecated(note = "pure Heston model, pin_jumps not used")]
    pub pin_jumps: bool,
    #[deprecated(note = "r is applied separately, not calibrated")]
    pub risk_free_rate: Option<f32>,
    #[deprecated(note = "q is applied separately, not calibrated")]
    pub dividend_yield: Option<f32>,
}

impl Default for CalibrationConfig {
    fn default() -> Self {
        Self {
            moneyness_weight_sigma: Some(0.30),
            #[allow(deprecated)]
            pin_jumps: false,
            #[allow(deprecated)]
            risk_free_rate: None,
            #[allow(deprecated)]
            dividend_yield: None,
        }
    }
}

/// Result of one deep-calibration run.
#[derive(Debug, Clone)]
pub struct DeepCalibrationResult {
    /// Calibrated Heston parameters extracted from the surrogate output.
    pub heston: HestonParameters,
    /// Full normalised parameter vector [0,1]^7 at the optimum.
    pub bates_norm: [f32; N_PARAMS],
    /// IV RMSE over valid (masked) market cells, in IV units (not bps).
    pub ivrmse: f32,
    /// Total surrogate forward-pass evaluations across all restarts.
    pub n_evals: usize,
}

/// Wraps the ORT session for repeated calibration calls.
pub struct BatesCalibrator {
    session: Session,
}

impl BatesCalibrator {
    /// Load from an ONNX file.
    ///
    /// The `bates_surrogate.onnx.data` external-weights file must live in the
    /// same directory as `onnx_path`.
    pub fn new(onnx_path: &str) -> ort::Result<Self> {
        let session = Session::builder()?
            .with_optimization_level(GraphOptimizationLevel::Level3)?
            .commit_from_file(onnx_path)?;
        Ok(Self { session })
    }

    /// Calibrate Heston parameters to a market IV surface.
    ///
    /// * `iv_market`  – `[f32; 686]` flat IV surface (NK×NT, strike-major).
    ///   Unobserved cells should be `0.0` and masked out via `mask`.
    /// * `mask`       – `[bool; 686]` validity; `true` = genuine market quote.
    /// * `n_restarts` – number of independent Nelder-Mead restarts (≥1).
    pub fn calibrate(
        &mut self,
        iv_market: &[f32; N_FLAT],
        mask: &[bool; N_FLAT],
        n_restarts: usize,
        config: &CalibrationConfig,
    ) -> ort::Result<DeepCalibrationResult> {
        let weights = build_weights(iv_market, mask, config.moneyness_weight_sigma);
        let seeds = make_seeds(n_restarts.max(1), config);

        let mut best: Option<(f32, [f32; N_PARAMS], usize)> = None;
        let mut total_evals = 0usize;

        for seed in &seeds {
            let (theta, loss, evals) = self.nelder_mead(seed, iv_market, &weights, config)?;
            total_evals += evals;
            let rmse = loss.sqrt();
            if best.as_ref().map_or(true, |b| rmse < b.0) {
                best = Some((rmse, theta, total_evals));
            }
        }

        let (rmse, mut bates_norm, n_evals) = best.unwrap();
        // Enforce pins in the stored result so bates_norm is consistent
        // with what was actually passed to the surrogate.
        pin_params(&mut bates_norm, config);
        let phys = denormalise(&bates_norm);
        Ok(DeepCalibrationResult {
            heston: bates_to_heston(&phys),
            bates_norm,
            ivrmse: rmse,
            n_evals,
        })
    }

    // -----------------------------------------------------------------------
    // Internals
    // -----------------------------------------------------------------------

    fn forward(&mut self, theta_norm: &[f32; N_PARAMS], config: &CalibrationConfig) -> ort::Result<[f32; N_FLAT]> {
        let mut input = *theta_norm;
        pin_params(&mut input, config);
        let tensor = Tensor::from_array(([1usize, N_PARAMS], input.to_vec()))?;
        let outputs = self.session.run(ort::inputs!["parameters" => tensor])?;
        let (_shape, slice) = outputs["iv_surface"].try_extract_tensor::<f32>()?;
        let mut out = [0f32; N_FLAT];
        out.copy_from_slice(&slice[..N_FLAT]);
        Ok(out)
    }

    fn objective(
        &mut self,
        theta: &[f32; N_PARAMS],
        iv_market: &[f32; N_FLAT],
        weights: &[f32; N_FLAT],
        config: &CalibrationConfig,
    ) -> ort::Result<f32> {
        let iv_pred = self.forward(theta, config)?;
        let mut num = 0f32;
        let mut den = 0f32;
        for i in 0..N_FLAT {
            let w = weights[i];
            if w > 0.0 {
                let e = iv_pred[i] - iv_market[i];
                num += w * e * e;
                den += w;
            }
        }
        let iv_loss = if den > 0.0 { num / den } else { f32::MAX };

        // L2 regularization: penalise Heston params (0-4) far from center [0.5].
        // This disambiguates parameters when IV surface alone is underdetermined.
        // Rate/dividend (5-6) are not regularised as they have clear market meaning.
        const REG_WEIGHT: f32 = 0.005; // Weak regularization; weight by tuning this
        let mut reg = 0.0_f32;
        for i in 0..5 {
            let d = theta[i] - 0.5;
            reg += d * d;
        }
        let total_loss = iv_loss + REG_WEIGHT * reg;
        Ok(total_loss)
    }

    /// Nelder-Mead minimiser in [0,1]^7.
    ///
    /// Box constraints are enforced by clamping every vertex to [0,1] after
    /// each geometric operation.
    fn nelder_mead(
        &mut self,
        init: &[f32; N_PARAMS],
        iv_market: &[f32; N_FLAT],
        weights: &[f32; N_FLAT],
        config: &CalibrationConfig,
    ) -> ort::Result<([f32; N_PARAMS], f32, usize)> {
        const MAX_ITERS: usize = 600;
        const TOL: f32 = 1e-7;
        // Standard NM coefficients
        const ALPHA: f32 = 1.0; // reflection
        const GAMMA: f32 = 2.0; // expansion
        const RHO: f32 = 0.5; // contraction
        const SIGMA: f32 = 0.5; // shrink

        let n = N_PARAMS;

        // Build initial simplex: init + n vertices perturbed by ±0.05.
        // Pins are applied to every vertex so the search stays in the
        // Heston (or configured) subspace from the first iteration.
        let mut simplex: Vec<[f32; N_PARAMS]> = Vec::with_capacity(n + 1);
        let mut init_pinned = *init;
        pin_params(&mut init_pinned, config);
        simplex.push(init_pinned);
        for i in 0..n {
            let mut v = init_pinned;
            v[i] = if v[i] < 0.5 {
                (v[i] + 0.05).min(1.0)
            } else {
                (v[i] - 0.05).max(0.0)
            };
            pin_params(&mut v, config);
            simplex.push(v);
        }

        let mut fvals: Vec<f32> = simplex
            .iter()
            .map(|v| self.objective(v, iv_market, weights, config))
            .collect::<ort::Result<_>>()?;
        let mut n_evals = n + 1;

        for _ in 0..MAX_ITERS {
            // Sort ascending by function value
            let mut order: Vec<usize> = (0..=n).collect();
            order.sort_unstable_by(|&a, &b| fvals[a].partial_cmp(&fvals[b]).unwrap());
            let sorted_simplex: Vec<_> = order.iter().map(|&i| simplex[i]).collect();
            let sorted_fvals: Vec<_> = order.iter().map(|&i| fvals[i]).collect();
            simplex = sorted_simplex;
            fvals = sorted_fvals;

            // Convergence check
            if fvals[n] - fvals[0] < TOL {
                break;
            }

            // Centroid of all but the worst
            let c = centroid(&simplex[..n]);

            // Reflection: x_r = c + α*(c − x_worst)
            let x_r = lerp_clamp(&c, &simplex[n], -ALPHA);
            let f_r = self.objective(&x_r, iv_market, weights, config)?;
            n_evals += 1;

            if f_r < fvals[0] {
                // Better than best: try expansion
                let x_e = lerp_clamp(&c, &simplex[n], -GAMMA);
                let f_e = self.objective(&x_e, iv_market, weights, config)?;
                n_evals += 1;
                if f_e < f_r {
                    simplex[n] = x_e;
                    fvals[n] = f_e;
                } else {
                    simplex[n] = x_r;
                    fvals[n] = f_r;
                }
            } else if f_r < fvals[n - 1] {
                // Better than second-worst: accept reflection
                simplex[n] = x_r;
                fvals[n] = f_r;
            } else {
                // Contraction
                let x_c = lerp_clamp(&c, &simplex[n], RHO);
                let f_c = self.objective(&x_c, iv_market, weights, config)?;
                n_evals += 1;
                if f_c < fvals[n] {
                    simplex[n] = x_c;
                    fvals[n] = f_c;
                } else {
                    // Shrink all toward best
                    for i in 1..=n {
                        simplex[i] = lerp_clamp(&simplex[0], &simplex[i], SIGMA);
                        fvals[i] = self.objective(&simplex[i], iv_market, weights, config)?;
                        n_evals += 1;
                    }
                }
            }
        }

        Ok((simplex[0], fvals[0], n_evals))
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Extract `HestonParameters` from the 7-dim physical parameter vector.
/// Indices: kappa=0, theta=1, sigma_v=2, rho=3, v0=4.
fn bates_to_heston(phys: &[f32; N_PARAMS]) -> HestonParameters {
    HestonParameters {
        kappa: phys[0] as f64,
        v_bar: phys[1] as f64,
        sigma: phys[2] as f64,
        rho: phys[3] as f64,
        v0: phys[4] as f64,
    }
}

/// Build per-cell weights: valid mask × optional Gaussian ATM down-weighting.
///
/// When `moneyness_sigma` is `Some(s)`, weight = exp(−lm² / (2s²)) so that
/// near-ATM cells (tighter bid-ask, better price discovery) dominate the loss.
fn build_weights(
    iv_market: &[f32; N_FLAT],
    mask: &[bool; N_FLAT],
    moneyness_sigma: Option<f32>,
) -> [f32; N_FLAT] {
    let lm_grid = log_moneyness_grid();
    let mut w = [0f32; N_FLAT];
    for ik in 0..NK {
        let lm = lm_grid[ik];
        let atm_w = match moneyness_sigma {
            Some(s) if s > 0.0 => (-0.5 * (lm / s).powi(2)).exp(),
            _ => 1.0,
        };
        for it in 0..NT {
            let i = ik * NT + it;
            if mask[i] && iv_market[i] > 0.0 {
                w[i] = atm_w;
            }
        }
    }
    w
}

/// Apply `CalibrationConfig` constraints to a normalised [0,1] parameter vector.
///
/// With pure Heston (5 parameters), there are no pins — all Heston parameters
/// are free to vary during calibration. Market rates (r, q) are applied
/// separately during pricing, not optimised.
fn pin_params(_theta: &mut [f32; N_PARAMS], _config: &CalibrationConfig) {
    // No-op for pure Heston. Kept for API compatibility.
    // If needed in future, can constrain specific Heston parameters here.
}

/// Initial starting points for multi-start, respecting config pins.
fn make_seeds(n: usize, config: &CalibrationConfig) -> Vec<[f32; N_PARAMS]> {
    let mut seeds: Vec<[f32; N_PARAMS]> = Vec::with_capacity(n);
    // Centroid of the Heston subspace
    let mut s0 = [0.5; N_PARAMS];
    pin_params(&mut s0, config);
    seeds.push(s0);
    // Equity-like: fast mean-reversion, moderate vol-of-vol, steep skew
    if n >= 2 {
        let mut s = [0.7, 0.3, 0.3, 0.1, 0.3];
        pin_params(&mut s, config);
        seeds.push(s);
    }
    // High vol-of-vol regime
    if n >= 3 {
        let mut s = [0.3, 0.5, 0.7, 0.2, 0.5];
        pin_params(&mut s, config);
        seeds.push(s);
    }
    // Fill remaining with evenly-spaced scalar values
    for k in seeds.len()..n {
        let t = (k as f32 + 1.0) / (n as f32 + 1.0);
        let mut s = [t; N_PARAMS];
        pin_params(&mut s, config);
        seeds.push(s);
    }
    seeds
}

/// Centroid of a slice of vertices.
fn centroid(verts: &[[f32; N_PARAMS]]) -> [f32; N_PARAMS] {
    let mut c = [0f32; N_PARAMS];
    let n = verts.len() as f32;
    for v in verts {
        for i in 0..N_PARAMS {
            c[i] += v[i];
        }
    }
    for i in 0..N_PARAMS {
        c[i] /= n;
    }
    c
}

/// `a + t * (b − a)`, clamped element-wise to [0, 1].
///
/// This covers all Nelder-Mead moves:
///   reflection  t = −α  →  a + α*(a − b)
///   expansion   t = −γ
///   contraction t = +ρ  →  a + ρ*(b − a)
///   shrink      t = +σ  (call with a=best, b=vertex)
fn lerp_clamp(a: &[f32; N_PARAMS], b: &[f32; N_PARAMS], t: f32) -> [f32; N_PARAMS] {
    let mut out = [0f32; N_PARAMS];
    for i in 0..N_PARAMS {
        out[i] = (a[i] + t * (b[i] - a[i])).clamp(0.0, 1.0);
    }
    out
}

// ---------------------------------------------------------------------------
// Market surface builder
// ---------------------------------------------------------------------------

/// Map a `CalibrationDataset` onto the fixed 49×14 ONNX grid.
///
/// Each instrument with a known `implied_vol` is assigned to its nearest
/// (strike, maturity) grid cell via log-moneyness = ln(K/S).  Multiple
/// instruments that fall on the same cell have their IVs averaged.  Cells
/// with no market quotes are left as 0.0 with `mask = false`.
///
/// Returns `(surface, mask)` where the flat index is `ik * NT + it`
/// (strike-major, matching the Python training layout).
pub fn build_iv_surface(dataset: &CalibrationDataset) -> ([f32; N_FLAT], [bool; N_FLAT]) {
    let lm_grid = log_moneyness_grid();
    let mut sum   = [0f32; N_FLAT];
    let mut count = [0u32; N_FLAT];

    for inst in &dataset.instruments {
        let iv = match inst.implied_vol {
            Some(v) if v > 0.0 => v as f32,
            _ => continue,
        };
        if inst.spot <= 0.0 || inst.strike <= 0.0 {
            continue;
        }
        let lm = (inst.strike / inst.spot).ln() as f32;
        let t  = inst.maturity_years as f32;

        let ik = lm_grid
            .iter()
            .enumerate()
            .min_by(|a, b| (a.1 - lm).abs().partial_cmp(&(b.1 - lm).abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        let it = MATURITIES
            .iter()
            .enumerate()
            .min_by(|a, b| (a.1 - t).abs().partial_cmp(&(b.1 - t).abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        let idx = ik * NT + it;
        sum[idx]   += iv;
        count[idx] += 1;
    }

    let mut surface = [0f32; N_FLAT];
    let mut mask    = [false; N_FLAT];
    for i in 0..N_FLAT {
        if count[i] > 0 {
            surface[i] = sum[i] / count[i] as f32;
            mask[i]    = true;
        }
    }
    (surface, mask)
}

/// Calibrate Heston parameters to a `CalibrationDataset` using the deep surrogate.
///
/// Builds the IV surface from `dataset`, then runs `BatesCalibrator::calibrate`.
pub fn calibrate_heston_from_dataset(
    calibrator: &mut BatesCalibrator,
    dataset: &CalibrationDataset,
    n_restarts: usize,
    config: &CalibrationConfig,
) -> ort::Result<DeepCalibrationResult> {
    let (surface, mask) = build_iv_surface(dataset);
    calibrator.calibrate(&surface, &mask, n_restarts, config)
}

/// End-to-end pipeline: load the ONNX surrogate, calibrate Heston parameters
/// from market data, and price an American option with the CTMC algorithm.
///
///   1. Load `bates_surrogate.onnx` from `onnx_path`.
///   2. Build IV surface from `dataset` and run Nelder-Mead calibration.
///   3. Price `option` with `price_american_option_heston`.
///
/// `n_x`, `m_v`, `n_time` control CTMC grid resolution (typical: 80, 20, 50).
pub fn price_american_option_deep(
    onnx_path: &str,
    dataset: &CalibrationDataset,
    option: options::Options,
    n_restarts: usize,
    n_x: usize,
    m_v: usize,
    n_time: usize,
    config: &CalibrationConfig,
) -> ort::Result<f64> {
    let mut calibrator = BatesCalibrator::new(onnx_path)?;
    let result = calibrate_heston_from_dataset(&mut calibrator, dataset, n_restarts, config)?;
    let h = result.heston;
    let ctmc = crate::ctmc::price_american_option_heston(
        option, h.v0, h.kappa, h.v_bar, h.sigma, h.rho, n_x, m_v, n_time,
    );
    Ok(ctmc.price)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::slv::{CalibrationInstrument, OptionSide};

    fn make_instrument(strike: f64, maturity: f64, iv: f64, spot: f64) -> CalibrationInstrument {
        CalibrationInstrument {
            symbol: "TEST".to_string(),
            side: OptionSide::Put,
            strike,
            maturity_years: maturity,
            market_price: 0.0,
            spot,
            risk_free_rate: 0.03,
            dividend_yield: 0.0,
            implied_vol: Some(iv),
            bid: None,
            ask: None,
            open_interest: None,
            volume: None,
        }
    }

    #[test]
    fn build_iv_surface_maps_instruments_to_grid() {
        let spot = 100.0_f64;
        let instruments = vec![
            make_instrument(100.0, 0.5,  0.20, spot),
            make_instrument(90.0,  0.25, 0.25, spot),
            make_instrument(110.0, 1.0,  0.18, spot),
        ];
        let dataset = CalibrationDataset {
            valuation_symbol: "TEST".to_string(),
            as_of_epoch_secs: 0,
            instruments,
        };

        let (surface, mask) = build_iv_surface(&dataset);

        let n_valid = mask.iter().filter(|&&m| m).count();
        assert!(n_valid >= 3, "expected >=3 valid cells, got {n_valid}");
        for i in 0..N_FLAT {
            if mask[i] {
                assert!(surface[i] > 0.0, "masked cell {i} has non-positive IV");
            } else {
                assert_eq!(surface[i], 0.0, "unmasked cell {i} should be 0.0");
            }
        }
    }

    #[test]
    fn build_iv_surface_averages_duplicate_cells() {
        let spot = 100.0_f64;
        let instruments = vec![
            make_instrument(100.0, 0.5, 0.20, spot),
            make_instrument(100.0, 0.5, 0.30, spot),
        ];
        let dataset = CalibrationDataset {
            valuation_symbol: "TEST".to_string(),
            as_of_epoch_secs: 0,
            instruments,
        };

        let (surface, mask) = build_iv_surface(&dataset);

        let n_valid = mask.iter().filter(|&&m| m).count();
        assert_eq!(n_valid, 1, "two instruments on same cell should produce 1 masked cell");

        let iv = surface.iter().zip(mask.iter()).find(|(_, m)| **m).map(|(v, _)| *v).unwrap();
        assert!((iv - 0.25).abs() < 1e-5, "expected averaged IV ~0.25, got {iv}");
    }

    /// Comparison test: deep-calibration pipeline vs direct CTMC with true params.
    ///
    /// Scenario: S=K=100, T=1yr, r=5%, q=2%.
    /// True Heston params: kappa=2.0, theta=0.04, sigma=0.35, rho=-0.65, v0=0.04.
    /// Synthetic IV surface is built by approximating the Heston skew via the ATM
    /// level (sqrt(theta)) and a linear skew slope from rho/sigma.
    ///
    /// Prints a table comparing:
    ///   1. CTMC (true params)     — ground truth
    ///   2. CTMC (Fourier LM cal)  — from full_pipeline_calibrate_then_ctmc_american_put
    ///   3. CTMC (deep cal)        — this pipeline
    ///   4. Black-Scholes European — lower bound reference
    #[test]
    #[ignore]
    fn price_american_option_deep_smoke() {
        use crate::ctmc::price_american_put_heston;

        let onnx_path = default_onnx_path().expect("ONNX model not found — run from workspace root");

        // Shared scenario
        let (spot, strike, r, q, tau) = (100.0_f64, 100.0_f64, 0.05_f64, 0.02_f64, 1.0_f64);
        // True Heston parameters (same as full_pipeline_calibrate_then_ctmc_american_put)
        let (kappa, theta, sigma_v, rho_h, v0) = (2.0_f64, 0.04_f64, 0.35_f64, -0.65_f64, 0.04_f64);
        let (n_x, m_v, n_time) = (100, 25, 75);

        // 1. Ground truth: CTMC with true Heston params
        let true_result = price_american_put_heston(
            strike, spot, v0, tau, r, q, kappa, theta, sigma_v, rho_h, n_x, m_v, n_time,
        );

        // 2. Deep calibration
        //    Build a synthetic surface using ATM IV ≈ sqrt(theta) and a
        //    first-order skew approximation: IV(m) ≈ atm_iv - skew * log_moneyness
        //    where skew ≈ -rho * sigma_v / (2 * sqrt(theta)) (Heston skew expansion).
        let atm_iv = theta.sqrt();
        let skew   = -rho_h * sigma_v / (2.0 * atm_iv);
        let moneyness_grid = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20_f64];
        let maturity_grid  = [0.25, 0.5, 1.0, 2.0_f64];
        let instruments: Vec<_> = moneyness_grid.iter().flat_map(|&m| {
            maturity_grid.iter().map(move |&t| {
                let lm = m.ln();
                let iv = (atm_iv - skew * lm).max(0.05);
                make_instrument(m * spot, t, iv, spot)
            })
        }).collect();
        let dataset = CalibrationDataset {
            valuation_symbol: "SYN".to_string(),
            as_of_epoch_secs: 0,
            instruments,
        };
        // Run deep calibration and also capture calibrated params
        let mut calibrator = BatesCalibrator::new(&onnx_path).expect("failed to load ONNX");
        let cal_result = calibrate_heston_from_dataset(&mut calibrator, &dataset, 3, &CalibrationConfig::default())
            .expect("deep calibration failed");
        let h = cal_result.heston;
        let deep_result = crate::ctmc::price_american_put_heston(
            strike, spot, h.v0, tau, r, q, h.kappa, h.v_bar, h.sigma, h.rho,
            n_x, m_v, n_time,
        );

        // 3. Black-Scholes European (lower bound, using ATM IV)
        let bs_put = options::Options::new_put(strike, spot, atm_iv, r, tau, Some(q));
        let bs_price = bs_put.bs_pricing();

        println!("\n=== American Put Pricing Comparison ===");
        println!("Scenario: S={spot}, K={strike}, T={tau}yr, r={:.1}%, q={:.1}%", r*100.0, q*100.0);
        println!();
        println!("True Heston params:");
        println!("  kappa={kappa}  theta={theta}  sigma={sigma_v}  rho={rho_h}  v0={v0}");
        println!();
        println!("Deep-calibrated Heston params (IVRMSE={:.4}):", cal_result.ivrmse);
        println!("  kappa={:.4}  theta={:.4}  sigma={:.4}  rho={:.4}  v0={:.4}",
            h.kappa, h.v_bar, h.sigma, h.rho, h.v0);
        println!();
        println!("{:<32} {:>10}", "Method", "Price");
        println!("{}", "-".repeat(44));
        println!("{:<32} {:>10.4}", "CTMC (true params)",      true_result.price);
        println!("{:<32} {:>10.4}", "CTMC (deep cal)",          deep_result.price);
        println!("{:<32} {:>10.4}", "BS European (ATM vol lb)", bs_price);
        println!("{}", "-".repeat(44));
        println!("{:<32} {:>10.4}", "Deep cal vs true", deep_result.price - true_result.price);

        assert!(deep_result.price > 0.0);
        assert!(deep_result.price < spot);
    }
}
