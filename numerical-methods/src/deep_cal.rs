//! Deep calibration via the BatesSurrogate ONNX model.
//!
//! Loads the exported `bates_surrogate.onnx` (+ companion `.onnx.data`) and
//! calibrates Heston parameters from a market IV surface using Levenberg-Marquardt
//! through the frozen surrogate (pure Heston, 5 parameters).
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
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../weights/bates_surrogate.onnx");
    if std::path::Path::new(workspace).exists() {
        return Some(workspace.to_string());
    }
    let dev = "/home/mttmin/coding/deep-calibration/model/bates_surrogate.onnx";
    if std::path::Path::new(dev).exists() {
        return Some(dev.to_string());
    }
    None
}

// ---------------------------------------------------------------------------
// Grid constants  (must match heston_datagen.py high41x14 preset)
// ---------------------------------------------------------------------------

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
// Parameter bounds  (must match PARAM_BOUNDS in heston_datagen.py exactly)
// ---------------------------------------------------------------------------
//
// The ONNX surrogate's [0,1] normalised space is defined by these bounds.
// Training data was normalised: theta_norm = (theta - LO) / (HI - LO).
// Inference must use the same bounds; any mismatch maps every calibrated
// parameter to the wrong physical value.
//
// Source of truth: PARAM_BOUNDS in training data creation/heston_datagen.py.

/// Physical lower bounds: [kappa, theta, sigma_v, rho, v0]
const PARAM_LO: [f32; N_PARAMS] = [0.30, 0.02, 0.05, -0.98, 0.02];
/// Physical upper bounds: [kappa, theta, sigma_v, rho, v0]
const PARAM_HI: [f32; N_PARAMS] = [8.00, 0.12, 1.20,  0.10, 0.12];

fn normalise(phys: &[f32; N_PARAMS]) -> [f32; N_PARAMS] {
    let mut out = [0f32; N_PARAMS];
    for i in 0..N_PARAMS {
        out[i] = (phys[i] - PARAM_LO[i]) / (PARAM_HI[i] - PARAM_LO[i]);
    }
    out
}

fn denormalise(norm: &[f32; N_PARAMS]) -> [f32; N_PARAMS] {
    let mut out = [0f32; N_PARAMS];
    for i in 0..N_PARAMS {
        out[i] = norm[i] * (PARAM_HI[i] - PARAM_LO[i]) + PARAM_LO[i];
    }
    out
}

/// Market-realistic Heston prior in normalised [0,1] space.
/// Physical centre: κ=3.0, θ=0.04, σ_v=0.40, ρ=−0.70, v₀=0.04.
/// This is a gentle attractor during calibration; the IV loss dominates
/// when enough market quotes are present.
fn market_prior_norm() -> [f32; N_PARAMS] {
    let phys = [3.0_f32, 0.04, 0.40, -0.70, 0.04];
    normalise(&phys)
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Controls the calibration search over pure Heston parameters.
#[derive(Debug, Clone)]
pub struct CalibrationConfig {
    /// Gaussian ATM weighting: each strike gets weight = exp(−lm² / (2σ²)).
    /// Near-ATM options have tighter bid-ask spreads and better price discovery.
    /// `None` = uniform weights. Typical value: 0.30.
    pub moneyness_weight_sigma: Option<f32>,

    /// L2 regularisation weight toward the market-realistic prior.
    /// The prior is κ=3.0, θ=0.04, σ_v=0.40, ρ=−0.70, v₀=0.04 in physical space.
    /// Prevents degenerate solutions when the surface has few/sparse quotes.
    /// Set to 0.0 to disable. Default: 0.001.
    pub reg_weight: f32,

    /// Weight residuals by Black-Scholes vega × confidence.
    ///
    /// Vega ∝ N′(d₁) · √T, where d₁ uses the observed IV at each cell.
    /// This correctly down-weights deep OTM and very short-dated contracts
    /// (large bid-ask spread relative to vega), and up-weights long-dated
    /// near-ATM contracts that carry the most calibration signal.
    ///
    /// When `false`, falls back to the Gaussian ATM weighting controlled by
    /// `moneyness_weight_sigma`. Default: `true`.
    pub use_vega_weights: bool,

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
            reg_weight: 0.001,
            use_vega_weights: true,
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
    /// Full normalised parameter vector [0,1]^5 at the optimum.
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
    ///
    /// # Panics
    /// Panics if the ONNX model's input is not `[batch, N_PARAMS]`.
    /// This catches model/code mismatches (e.g. model retrained with different
    /// n_params) at startup rather than mid-calibration.
    pub fn new(onnx_path: &str) -> ort::Result<Self> {
        let mut session = Session::builder()?
            .with_optimization_level(GraphOptimizationLevel::Level3)?
            .commit_from_file(onnx_path)?;

        // Shape sanity: run a 1×N_PARAMS forward pass. If the ONNX model was
        // exported with a different n_params, ORT will return a shape-mismatch
        // error here rather than silently producing wrong outputs mid-calibration.
        let dummy = Tensor::from_array(([1usize, N_PARAMS], vec![0.5f32; N_PARAMS]))?;
        session.run(ort::inputs!["parameters" => dummy])?;

        Ok(Self { session })
    }

    // -----------------------------------------------------------------------
    // Public calibration API
    // -----------------------------------------------------------------------

    /// Calibrate Heston parameters given pre-computed weights.
    ///
    /// Weights already incorporate mask, moneyness Gaussian down-weighting,
    /// and confidence (from IDW interpolation). Zero-weight cells are excluded
    /// from the loss.
    ///
    /// * `iv_market` – `[f32; 686]` flat IV surface (NK×NT, strike-major).
    /// * `weights`   – `[f32; 686]` per-cell weights; 0.0 = excluded.
    /// * `n_restarts` – number of independent LM restarts (≥1).
    pub fn calibrate_with_weights(
        &mut self,
        iv_market: &[f32; N_FLAT],
        weights: &[f32; N_FLAT],
        n_restarts: usize,
        config: &CalibrationConfig,
    ) -> ort::Result<DeepCalibrationResult> {
        let seeds = make_seeds(n_restarts.max(1));

        let mut best: Option<(f32, [f32; N_PARAMS], usize)> = None;
        let mut total_evals = 0usize;

        for seed in &seeds {
            let (theta, iv_loss, evals) =
                self.levenberg_marquardt(seed, iv_market, weights, config)?;
            total_evals += evals;
            let rmse = iv_loss.sqrt();
            if best.as_ref().map_or(true, |b| rmse < b.0) {
                best = Some((rmse, theta, total_evals));
            }
        }

        let (rmse, bates_norm, n_evals) = best.unwrap();
        let phys = denormalise(&bates_norm);
        Ok(DeepCalibrationResult {
            heston: bates_to_heston(&phys),
            bates_norm,
            ivrmse: rmse,
            n_evals,
        })
    }

    /// Backward-compatible wrapper: builds weights from mask (uniform confidence = 1.0).
    ///
    /// Prefer `calibrate_with_weights` when confidence information is available
    /// (e.g. from `build_iv_surface`).
    pub fn calibrate(
        &mut self,
        iv_market: &[f32; N_FLAT],
        mask: &[bool; N_FLAT],
        n_restarts: usize,
        config: &CalibrationConfig,
    ) -> ort::Result<DeepCalibrationResult> {
        let confidence = [1.0f32; N_FLAT];
        let weights = build_weights(iv_market, mask, &confidence, config.moneyness_weight_sigma);
        self.calibrate_with_weights(iv_market, &weights, n_restarts, config)
    }

    // -----------------------------------------------------------------------
    // Internals
    // -----------------------------------------------------------------------

    fn forward(&mut self, theta_norm: &[f32; N_PARAMS]) -> ort::Result<[f32; N_FLAT]> {
        let tensor = Tensor::from_array(([1usize, N_PARAMS], theta_norm.to_vec()))?;
        let outputs = self.session.run(ort::inputs!["parameters" => tensor])?;
        let (_shape, slice) = outputs["iv_surface"].try_extract_tensor::<f32>()?;
        let mut out = [0f32; N_FLAT];
        out.copy_from_slice(&slice[..N_FLAT]);
        Ok(out)
    }

    /// Levenberg-Marquardt calibration in normalised [0,1]^5 space.
    ///
    /// Uses finite-difference Jacobians (N_PARAMS forward passes per iteration).
    /// Regularises toward the market-realistic prior with weight `config.reg_weight`.
    ///
    /// Damping schedule (matches slv.rs):
    ///   initial = 1e-2, up_factor = 10.0, down_factor = 0.33.
    ///
    /// Returns `(theta_norm, weighted_iv_mse_no_reg, n_evals)`.
    fn levenberg_marquardt(
        &mut self,
        init: &[f32; N_PARAMS],
        iv_market: &[f32; N_FLAT],
        weights: &[f32; N_FLAT],
        config: &CalibrationConfig,
    ) -> ort::Result<([f32; N_PARAMS], f32, usize)> {
        const MAX_ITERS: usize = 50;
        const EPS: f32 = 1e-4;      // finite-difference step in normalised space
        const DAMP_INIT: f32 = 1e-2;
        const DAMP_UP: f32 = 10.0;
        const DAMP_DOWN: f32 = 0.33;
        const STEP_TOL: f32 = 1e-6;

        let prior = market_prior_norm();
        let reg_w = config.reg_weight;

        let clamp = |v: [f32; N_PARAMS]| -> [f32; N_PARAMS] {
            let mut c = v;
            for x in &mut c {
                *x = x.clamp(0.0, 1.0);
            }
            c
        };

        let mut theta = clamp(*init);
        let mut iv_pred = self.forward(&theta)?;
        let mut cost = lm_cost(&iv_pred, iv_market, weights, &theta, &prior, reg_w);
        let mut damping = DAMP_INIT;
        let mut n_evals = 1usize;

        // Collect indices of valid (non-zero weight) cells for fast inner loops.
        let valid: Vec<usize> = (0..N_FLAT).filter(|&i| weights[i] > 0.0).collect();

        for _ in 0..MAX_ITERS {
            if valid.is_empty() {
                break;
            }

            // Build Jacobian columns via forward differences.
            // jac[j][k] = d(iv_pred[valid[k]]) / d(theta[j])
            let n_valid = valid.len();
            let mut jac = vec![0f32; N_PARAMS * n_valid];

            for j in 0..N_PARAMS {
                let mut th_b = theta;
                let step = if theta[j] + EPS <= 1.0 { EPS } else { -EPS };
                th_b[j] = (theta[j] + step).clamp(0.0, 1.0);
                let eps_actual = th_b[j] - theta[j];
                if eps_actual.abs() < 1e-10 {
                    continue;
                }
                let iv_b = self.forward(&th_b)?;
                n_evals += 1;
                for (k, &i) in valid.iter().enumerate() {
                    jac[j * n_valid + k] = (iv_b[i] - iv_pred[i]) / eps_actual;
                }
            }

            // Build H = J^T W J + reg_w*I + damping*I
            // and g = J^T W r + reg_w*(theta - prior)
            // where W[i] = weights[i] and r[i] = iv_pred[i] - iv_market[i].
            let mut h = [[0f32; N_PARAMS]; N_PARAMS];
            let mut g = [0f32; N_PARAMS];

            for a in 0..N_PARAMS {
                for b in a..N_PARAMS {
                    let mut s = 0f32;
                    for (k, &i) in valid.iter().enumerate() {
                        s += weights[i] * jac[a * n_valid + k] * jac[b * n_valid + k];
                    }
                    h[a][b] = s;
                    h[b][a] = s;
                }
                let mut rhs = 0f32;
                for (k, &i) in valid.iter().enumerate() {
                    rhs += weights[i] * jac[a * n_valid + k] * (iv_pred[i] - iv_market[i]);
                }
                // Prior regularisation: pulls toward market_prior_norm().
                rhs += reg_w * (theta[a] - prior[a]);
                g[a] = rhs;
                h[a][a] += reg_w + damping;
            }

            // Solve H * delta = -g  (5×5 linear system).
            let delta = match solve_5x5(&h, &g.map(|x| -x)) {
                Some(d) => d,
                None => break,
            };

            let step_norm: f32 = delta.iter().map(|d| d * d).sum::<f32>().sqrt();
            if step_norm < STEP_TOL {
                break;
            }

            let candidate = clamp({
                let mut c = theta;
                for i in 0..N_PARAMS {
                    c[i] += delta[i];
                }
                c
            });

            let iv_cand = self.forward(&candidate)?;
            n_evals += 1;
            let cost_cand = lm_cost(&iv_cand, iv_market, weights, &candidate, &prior, reg_w);

            if cost_cand < cost {
                theta = candidate;
                iv_pred = iv_cand;
                cost = cost_cand;
                damping = (damping * DAMP_DOWN).max(1e-10);
            } else {
                damping = (damping * DAMP_UP).min(1e8);
            }
        }

        // Return unregularised weighted IV MSE for cross-restart RMSE comparison.
        let iv_loss = weighted_iv_mse(&iv_pred, iv_market, weights);
        Ok((theta, iv_loss, n_evals))
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Extract `HestonParameters` from the 5-dim physical parameter vector.
fn bates_to_heston(phys: &[f32; N_PARAMS]) -> HestonParameters {
    HestonParameters {
        kappa: phys[0] as f64,
        v_bar: phys[1] as f64,
        sigma: phys[2] as f64,
        rho: phys[3] as f64,
        v0: phys[4] as f64,
    }
}

/// LM total cost: weighted IV MSE + prior regularisation.
fn lm_cost(
    iv_pred: &[f32; N_FLAT],
    iv_market: &[f32; N_FLAT],
    weights: &[f32; N_FLAT],
    theta: &[f32; N_PARAMS],
    prior: &[f32; N_PARAMS],
    reg_weight: f32,
) -> f32 {
    let iv_loss = weighted_iv_mse(iv_pred, iv_market, weights);
    let mut reg = 0f32;
    for i in 0..N_PARAMS {
        let d = theta[i] - prior[i];
        reg += d * d;
    }
    iv_loss + reg_weight * reg
}

/// Weighted mean-squared IV error (excluding zero-weight cells).
fn weighted_iv_mse(
    iv_pred: &[f32; N_FLAT],
    iv_market: &[f32; N_FLAT],
    weights: &[f32; N_FLAT],
) -> f32 {
    let mut num = 0f32;
    let mut den = 0f32;
    for i in 0..N_FLAT {
        if weights[i] > 0.0 {
            let e = iv_pred[i] - iv_market[i];
            num += weights[i] * e * e;
            den += weights[i];
        }
    }
    if den > 0.0 { num / den } else { f32::MAX }
}

/// Gaussian elimination with partial pivoting: Ax = b for A ∈ R^{5×5}.
/// Returns `None` if the system is numerically singular.
fn solve_5x5(a: &[[f32; N_PARAMS]; N_PARAMS], b: &[f32; N_PARAMS]) -> Option<[f32; N_PARAMS]> {
    // Augmented matrix [A | b], shape (5, 6)
    let mut m = [[0f32; N_PARAMS + 1]; N_PARAMS];
    for i in 0..N_PARAMS {
        m[i][..N_PARAMS].copy_from_slice(&a[i]);
        m[i][N_PARAMS] = b[i];
    }
    for col in 0..N_PARAMS {
        // Partial pivot
        let mut max_row = col;
        for row in (col + 1)..N_PARAMS {
            if m[row][col].abs() > m[max_row][col].abs() {
                max_row = row;
            }
        }
        m.swap(col, max_row);
        let pivot = m[col][col];
        if pivot.abs() < 1e-10 {
            return None;
        }
        for j in col..=N_PARAMS {
            m[col][j] /= pivot;
        }
        for row in 0..N_PARAMS {
            if row != col {
                let factor = m[row][col];
                for j in col..=N_PARAMS {
                    m[row][j] -= factor * m[col][j];
                }
            }
        }
    }
    let mut x = [0f32; N_PARAMS];
    for i in 0..N_PARAMS {
        x[i] = m[i][N_PARAMS];
    }
    Some(x)
}

/// Build per-cell weights using Black-Scholes vega × confidence.
///
/// Vega proxy: N′(d₁) · √T, where d₁ = (−lm + ½σ²T) / (σ√T) uses the
/// observed IV at each grid cell as σ.  This correctly down-weights deep
/// OTM short-dated contracts (tiny vega, large relative bid-ask) and
/// up-weights long-dated near-ATM contracts that anchor the smile shape.
///
/// `confidence[i]` ∈ (0, 1] further scales each cell by its IDW reliability.
pub fn build_vega_weights(
    iv_surface: &[f32; N_FLAT],
    mask: &[bool; N_FLAT],
    confidence: &[f32; N_FLAT],
) -> [f32; N_FLAT] {
    const INV_SQRT_2PI: f32 = 0.398_942_28; // 1 / sqrt(2π)
    let lm_grid = log_moneyness_grid();
    let mut w = [0f32; N_FLAT];
    for ik in 0..NK {
        let lm = lm_grid[ik];
        for it in 0..NT {
            let i = ik * NT + it;
            if !mask[i] || iv_surface[i] <= 0.0 {
                continue;
            }
            let sigma = iv_surface[i];
            let t = MATURITIES[it];
            let sigma_sqrt_t = sigma * t.sqrt();
            // BS d₁ (r = q = 0 approximation — only the relative weighting matters)
            let d1 = (-lm + 0.5 * sigma * sigma * t) / (sigma_sqrt_t + 1e-8);
            // N′(d₁) · √T: large for near-ATM, near-zero for deep OTM
            let vega_proxy = INV_SQRT_2PI * (-0.5 * d1 * d1).exp() * t.sqrt();
            w[i] = vega_proxy * confidence[i].max(0.0);
        }
    }
    w
}

/// Build per-cell weights: mask × confidence × Gaussian ATM down-weighting.
///
/// Fallback when `use_vega_weights = false`. Applies exp(−lm² / (2σ²)) when
/// `moneyness_sigma` is `Some(s)`, otherwise uniform within the mask.
pub fn build_weights(
    iv_market: &[f32; N_FLAT],
    mask: &[bool; N_FLAT],
    confidence: &[f32; N_FLAT],
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
                w[i] = atm_w * confidence[i].max(0.0);
            }
        }
    }
    w
}

/// Multi-start seeds in normalised [0,1]^5.
fn make_seeds(n: usize) -> Vec<[f32; N_PARAMS]> {
    let mut seeds: Vec<[f32; N_PARAMS]> = Vec::with_capacity(n);
    // Centroid of parameter space
    seeds.push([0.5; N_PARAMS]);
    // Near-market prior in normalised space (κ=3, θ=0.04, σ_v=0.4, ρ=−0.7, v₀=0.04)
    if n >= 2 {
        seeds.push(market_prior_norm());
    }
    // Equity-like: fast mean-reversion, moderate vol-of-vol, steep skew
    // (normalised kappa≈0.7 → physical ≈5.7, rho≈0.1 → physical ≈−0.87)
    if n >= 3 {
        seeds.push([0.7, 0.3, 0.3, 0.1, 0.3]);
    }
    // High vol-of-vol regime
    if n >= 4 {
        seeds.push([0.3, 0.5, 0.7, 0.2, 0.5]);
    }
    // Fill remaining with evenly-spaced scalar values
    for k in seeds.len()..n {
        let t = (k as f32 + 1.0) / (n as f32 + 1.0);
        seeds.push([t; N_PARAMS]);
    }
    seeds
}

// ---------------------------------------------------------------------------
// Market surface builder
// ---------------------------------------------------------------------------

/// Map a `CalibrationDataset` onto the fixed 49×14 ONNX grid using
/// inverse-distance-weighted (IDW) interpolation.
///
/// For each grid cell, the K=4 nearest observed quotes in (log-moneyness, T)
/// space are combined via inverse-distance weighting. Cells near real quotes
/// receive `confidence ≈ 1.0`; extrapolated cells receive as low as `0.1`.
///
/// Returns `(surface, mask, confidence)` where:
///   - `surface[ik * NT + it]` is the interpolated IV (always finite for a
///     non-empty dataset),
///   - `mask[i]` is `true` for every filled cell (the entire grid when at
///     least one quote is present),
///   - `confidence[i]` ∈ [0.1, 1.0] encodes proximity to real market quotes.
///
/// Flat index is `ik * NT + it` (strike-major, matching Python training layout).
pub fn build_iv_surface(
    dataset: &CalibrationDataset,
) -> ([f32; N_FLAT], [bool; N_FLAT], [f32; N_FLAT]) {
    let lm_grid = log_moneyness_grid();

    // Collect valid observed quotes: (log_moneyness, maturity, iv)
    let quotes: Vec<(f32, f32, f32)> = dataset
        .instruments
        .iter()
        .filter_map(|inst| {
            let iv = inst.implied_vol?;
            if iv <= 0.0 || inst.spot <= 0.0 || inst.strike <= 0.0 {
                return None;
            }
            let lm = (inst.strike / inst.spot).ln() as f32;
            let t = inst.maturity_years as f32;
            Some((lm, t, iv as f32))
        })
        .collect();

    if quotes.is_empty() {
        return ([0f32; N_FLAT], [false; N_FLAT], [0f32; N_FLAT]);
    }

    // Normalisation scales for the distance metric (grid extent in each axis)
    const LM_SCALE: f32 = 1.20; // range of log-moneyness grid
    const T_SCALE: f32 = 3.00;  // approximate range of maturity grid
    const K_NEAREST: usize = 4;

    // Confidence thresholds:
    //   distance < NEAR → confidence = 1.0 (cell has a real nearby quote)
    //   distance > FAR  → confidence = MIN_CONF (fully extrapolated)
    //   in between      → linear decay
    const NEAR_THRESH: f32 = 0.15;
    const FAR_THRESH: f32 = 0.50;
    const MIN_CONF: f32 = 0.10;

    let mut surface = [0f32; N_FLAT];
    let mut mask = [false; N_FLAT];
    let mut confidence = [0f32; N_FLAT];

    for ik in 0..NK {
        let lm = lm_grid[ik];
        for it in 0..NT {
            let t = MATURITIES[it];

            // Normalised distance to each observed quote
            let mut dists: Vec<(f32, f32)> = quotes
                .iter()
                .map(|&(qlm, qt, qiv)| {
                    let dlm = (qlm - lm) / LM_SCALE;
                    let dt = (qt - t) / T_SCALE;
                    let d = (dlm * dlm + dt * dt).sqrt();
                    (d, qiv)
                })
                .collect();
            dists.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

            let min_dist = dists[0].0;

            let iv = if min_dist < 1e-6 {
                // Exact hit(s): average all zero-distance quotes.
                // This handles duplicate instruments at the same (lm, T).
                let n_exact = dists.iter().take_while(|&&(d, _)| d < 1e-6).count();
                dists[..n_exact].iter().map(|&(_, iv)| iv).sum::<f32>() / n_exact as f32
            } else {
                // Inverse-distance weighted average of K nearest quotes.
                let k = K_NEAREST.min(dists.len());
                let mut num = 0f32;
                let mut den = 0f32;
                for &(d, qiv) in &dists[..k] {
                    let w = 1.0 / d;
                    num += w * qiv;
                    den += w;
                }
                num / den
            };

            // Confidence based on normalised distance to the nearest real quote.
            let conf = if min_dist < NEAR_THRESH {
                1.0
            } else if min_dist > FAR_THRESH {
                MIN_CONF
            } else {
                let t_c = (min_dist - NEAR_THRESH) / (FAR_THRESH - NEAR_THRESH);
                1.0 - (1.0 - MIN_CONF) * t_c
            };

            let idx = ik * NT + it;
            surface[idx] = iv;
            mask[idx] = true;
            confidence[idx] = conf;
        }
    }

    (surface, mask, confidence)
}

/// Calibrate Heston parameters to a `CalibrationDataset` using the deep surrogate.
///
/// Builds the IV surface from `dataset` using IDW interpolation with confidence
/// weighting, then runs multi-start LM calibration.
pub fn calibrate_heston_from_dataset(
    calibrator: &mut BatesCalibrator,
    dataset: &CalibrationDataset,
    n_restarts: usize,
    config: &CalibrationConfig,
) -> ort::Result<DeepCalibrationResult> {
    let (surface, mask, confidence) = build_iv_surface(dataset);
    let weights = if config.use_vega_weights {
        build_vega_weights(&surface, &mask, &confidence)
    } else {
        build_weights(&surface, &mask, &confidence, config.moneyness_weight_sigma)
    };
    calibrator.calibrate_with_weights(&surface, &weights, n_restarts, config)
}

/// End-to-end pipeline: load the ONNX surrogate, calibrate Heston parameters
/// from market data, and price an American option with the CTMC algorithm.
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

        let (surface, mask, confidence) = build_iv_surface(&dataset);

        // IDW fills all 686 cells (one non-empty quote set → entire grid interpolated)
        let n_valid = mask.iter().filter(|&&m| m).count();
        assert_eq!(n_valid, N_FLAT, "IDW should fill all cells; got {n_valid}");

        for i in 0..N_FLAT {
            assert!(surface[i] > 0.0, "IDW cell {i} has non-positive IV");
            assert!(
                confidence[i] >= 0.1 && confidence[i] <= 1.0,
                "confidence out of [0.1, 1.0] at cell {i}: {}",
                confidence[i]
            );
        }
    }

    #[test]
    fn build_iv_surface_averages_duplicate_cells() {
        // Two instruments at the same (K, T) → IDW exact-hit averages them.
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

        let (surface, mask, _confidence) = build_iv_surface(&dataset);

        // All cells filled via IDW with only one unique location
        let n_valid = mask.iter().filter(|&&m| m).count();
        assert_eq!(n_valid, N_FLAT, "IDW fills all cells; got {n_valid}");

        // The cell closest to (lm=0, T=0.5) should have IV ≈ 0.25 (average)
        let lm_grid = log_moneyness_grid();
        let ik = lm_grid
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.abs().partial_cmp(&b.1.abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        let it = MATURITIES
            .iter()
            .enumerate()
            .min_by(|a, b| (a.1 - 0.5_f32).abs().partial_cmp(&(b.1 - 0.5_f32).abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        let cell_iv = surface[ik * NT + it];
        assert!(
            (cell_iv - 0.25).abs() < 1e-4,
            "expected IV ≈ 0.25 at ATM 6M cell, got {cell_iv}"
        );
    }

    /// Synthetic round-trip test: calibrate known Heston parameters from a
    /// COS-generated IV surface and verify recovery within 15% relative error.
    ///
    /// This test loads the ONNX model — run with:
    ///   cargo test --release -- --ignored synthetic_round_trip_recovers_known_params
    ///
    /// If this test fails after all other fixes, the surrogate itself is the
    /// bottleneck and retraining is required.
    #[test]
    #[ignore]
    fn synthetic_round_trip_recovers_known_params() {
        use crate::slv::{cos_european, HestonParameters, OptionSide};
        use options::black_scholes::calculate_iv;

        let onnx_path =
            default_onnx_path().expect("ONNX model not found — run from workspace root");
        let mut calibrator = BatesCalibrator::new(&onnx_path).expect("failed to load ONNX");

        // True Heston parameters — comfortably inside training bounds.
        let true_phys = [2.0_f64, 0.04, 0.35, -0.65, 0.04];
        let true_heston = HestonParameters {
            kappa: true_phys[0],
            v_bar: true_phys[1],
            sigma: true_phys[2],
            rho:   true_phys[3],
            v0:    true_phys[4],
        };

        // Scenario: S=100, r=3%, q=1%
        let spot = 100.0_f64;
        let r = 0.03_f64;
        let q = 0.01_f64;
        let n_cos = 256;

        // Generate IV surface using Heston COS pricer for each grid cell.
        let lm_grid = log_moneyness_grid();
        let mut iv_market = [0f32; N_FLAT];
        let confidence = [1.0f32; N_FLAT]; // All cells are "real" quotes
        let mut mask = [true; N_FLAT];

        for ik in 0..NK {
            let lm = lm_grid[ik] as f64;
            let strike = spot * lm.exp();
            for it in 0..NT {
                let tau = MATURITIES[it] as f64;
                let idx = ik * NT + it;

                // Price a put option (OTM for negative lm, ITM otherwise)
                let price = cos_european(
                    OptionSide::Put,
                    spot,
                    strike,
                    r,
                    q,
                    tau,
                    &true_heston,
                    n_cos,
                );

                // Extract implied vol via Newton-Raphson bisection
                let iv = calculate_iv(price, spot, strike, tau, r, q, false);

                if iv > 0.005 && iv < 5.0 {
                    iv_market[idx] = iv as f32;
                } else {
                    // Fallback: ATM vol approximation for problematic cells
                    iv_market[idx] = (true_heston.v_bar as f32).sqrt();
                    mask[idx] = false; // Don't use as calibration target
                }
            }
        }

        let weights = build_weights(&iv_market, &mask, &confidence, Some(0.30));
        let result = calibrator
            .calibrate_with_weights(&iv_market, &weights, 5, &CalibrationConfig::default())
            .expect("calibration failed");

        let h = result.heston;
        println!("\n=== Synthetic Round-Trip Test ===");
        println!(
            "True:      κ={:.4}  θ={:.4}  σ_v={:.4}  ρ={:.4}  v₀={:.4}",
            true_phys[0], true_phys[1], true_phys[2], true_phys[3], true_phys[4]
        );
        println!(
            "Recovered: κ={:.4}  θ={:.4}  σ_v={:.4}  ρ={:.4}  v₀={:.4}",
            h.kappa, h.v_bar, h.sigma, h.rho, h.v0
        );
        println!("IVRMSE: {:.4} ({:.1} bps)", result.ivrmse, result.ivrmse * 10000.0);

        assert!(
            (h.kappa - true_phys[0]).abs() / true_phys[0] < 0.15,
            "kappa: true={:.4} got={:.4}",
            true_phys[0], h.kappa
        );
        assert!(
            (h.v_bar - true_phys[1]).abs() / true_phys[1] < 0.15,
            "theta: true={:.4} got={:.4}",
            true_phys[1], h.v_bar
        );
        assert!(
            (h.sigma - true_phys[2]).abs() / true_phys[2] < 0.15,
            "sigma_v: true={:.4} got={:.4}",
            true_phys[2], h.sigma
        );
        assert!(
            (h.rho - true_phys[3]).abs() / true_phys[3].abs() < 0.15,
            "rho: true={:.4} got={:.4}",
            true_phys[3], h.rho
        );
        assert!(
            (h.v0 - true_phys[4]).abs() / true_phys[4] < 0.15,
            "v0: true={:.4} got={:.4}",
            true_phys[4], h.v0
        );
    }

    /// Real-data smoke test: calibrate a synthetic option chain and verify
    /// the recovered parameters fall in market-realistic ranges.
    ///
    /// Exercises the full pipeline: build_iv_surface → IDW interpolation →
    /// LM calibration. Does NOT require a specific real dataset.
    ///
    /// Run with:
    ///   cargo test --release -- --ignored real_data_parameters_in_market_range
    #[test]
    #[ignore]
    fn real_data_parameters_in_market_range() {
        let onnx_path =
            default_onnx_path().expect("ONNX model not found — run from workspace root");
        let mut calibrator = BatesCalibrator::new(&onnx_path).expect("failed to load ONNX");

        // Synthetic SPY-like option chain: ~30 instruments across 4 expiries
        let spot = 450.0_f64;
        let atm_iv = 0.18_f64;
        let skew = 0.15_f64; // ≈ −∂IV/∂log-moneyness

        let moneyness = [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15_f64];
        let expiries = [1.0 / 12.0, 3.0 / 12.0, 6.0 / 12.0, 1.0_f64];

        let instruments: Vec<CalibrationInstrument> = moneyness
            .iter()
            .flat_map(|&m| {
                expiries.iter().map(move |&t| {
                    let lm = m.ln();
                    let iv = (atm_iv - skew * lm).clamp(0.05, 1.0);
                    CalibrationInstrument {
                        symbol: "SPY".to_string(),
                        side: if m < 1.0 { OptionSide::Put } else { OptionSide::Call },
                        strike: m * spot,
                        maturity_years: t,
                        market_price: 0.0,
                        spot,
                        risk_free_rate: 0.05,
                        dividend_yield: 0.015,
                        implied_vol: Some(iv),
                        bid: None,
                        ask: None,
                        open_interest: None,
                        volume: None,
                    }
                })
            })
            .collect();

        let dataset = CalibrationDataset {
            valuation_symbol: "SPY".to_string(),
            as_of_epoch_secs: 0,
            instruments,
        };

        let result = calibrate_heston_from_dataset(&mut calibrator, &dataset, 5, &CalibrationConfig::default())
            .expect("calibration failed");

        let p = result.heston;
        println!("\n=== Real-Data Smoke Test (SPY-like) ===");
        println!(
            "κ={:.4}  θ={:.4}  σ_v={:.4}  ρ={:.4}  v₀={:.4}",
            p.kappa, p.v_bar, p.sigma, p.rho, p.v0
        );
        println!("IVRMSE: {:.4} ({:.1} bps)", result.ivrmse, result.ivrmse * 10000.0);

        // Assert market-realistic parameter ranges
        assert!(p.kappa > 0.5 && p.kappa < 10.0, "kappa={:.4} out of [0.5, 10]", p.kappa);
        assert!(p.v_bar > 0.01 && p.v_bar < 0.10, "theta={:.4} out of [0.01, 0.10]", p.v_bar);
        assert!(p.sigma > 0.10 && p.sigma < 1.20, "sigma_v={:.4} out of [0.10, 1.20]", p.sigma);
        assert!(p.rho < -0.30 && p.rho > -0.98, "rho={:.4} out of (-0.98, -0.30)", p.rho);
        assert!(p.v0 > 0.01 && p.v0 < 0.10, "v0={:.4} out of [0.01, 0.10]", p.v0);

        // Loose IVRMSE threshold — if fit is good it will be much lower
        let ivrmse_bps = result.ivrmse * 10000.0;
        assert!(
            ivrmse_bps < 200.0,
            "IVRMSE={:.1} bps exceeds 200 bps threshold",
            ivrmse_bps
        );
    }

    /// Comparison test: deep-calibration pipeline vs direct CTMC with true params.
    #[test]
    #[ignore]
    fn price_american_option_deep_smoke() {
        use crate::ctmc::price_american_put_heston;

        let onnx_path = default_onnx_path().expect("ONNX model not found — run from workspace root");

        let (spot, strike, r, q, tau) = (100.0_f64, 100.0_f64, 0.05_f64, 0.02_f64, 1.0_f64);
        let (kappa, theta, sigma_v, rho_h, v0) = (2.0_f64, 0.04_f64, 0.35_f64, -0.65_f64, 0.04_f64);
        let (n_x, m_v, n_time) = (100, 25, 75);

        // Ground truth: CTMC with true Heston params
        let true_result = price_american_put_heston(
            strike, spot, v0, tau, r, q, kappa, theta, sigma_v, rho_h, n_x, m_v, n_time,
        );

        // Build a synthetic IV surface via ATM + linear skew approximation
        let atm_iv = theta.sqrt();
        let skew = -rho_h * sigma_v / (2.0 * atm_iv);
        let moneyness_grid = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20_f64];
        let maturity_grid = [0.25, 0.5, 1.0, 2.0_f64];
        let instruments: Vec<_> = moneyness_grid
            .iter()
            .flat_map(|&m| {
                maturity_grid.iter().map(move |&t| {
                    let lm = m.ln();
                    let iv = (atm_iv - skew * lm).max(0.05);
                    CalibrationInstrument {
                        symbol: "SYN".to_string(),
                        side: OptionSide::Put,
                        strike: m * spot,
                        maturity_years: t,
                        market_price: 0.0,
                        spot,
                        risk_free_rate: r,
                        dividend_yield: q,
                        implied_vol: Some(iv),
                        bid: None,
                        ask: None,
                        open_interest: None,
                        volume: None,
                    }
                })
            })
            .collect();
        let dataset = CalibrationDataset {
            valuation_symbol: "SYN".to_string(),
            as_of_epoch_secs: 0,
            instruments,
        };

        let mut calibrator = BatesCalibrator::new(&onnx_path).expect("failed to load ONNX");
        let cal_result =
            calibrate_heston_from_dataset(&mut calibrator, &dataset, 3, &CalibrationConfig::default())
                .expect("deep calibration failed");
        let h = cal_result.heston;
        let deep_result = crate::ctmc::price_american_put_heston(
            strike, spot, h.v0, tau, r, q, h.kappa, h.v_bar, h.sigma, h.rho,
            n_x, m_v, n_time,
        );

        let bs_put = options::Options::new_put(strike, spot, atm_iv, r, tau, Some(q));
        let bs_price = bs_put.bs_pricing();

        println!("\n=== American Put Pricing Comparison ===");
        println!("Scenario: S={spot}, K={strike}, T={tau}yr, r={:.1}%, q={:.1}%", r * 100.0, q * 100.0);
        println!("True Heston params:");
        println!("  kappa={kappa}  theta={theta}  sigma={sigma_v}  rho={rho_h}  v0={v0}");
        println!("Deep-calibrated Heston params (IVRMSE={:.4}):", cal_result.ivrmse);
        println!(
            "  kappa={:.4}  theta={:.4}  sigma={:.4}  rho={:.4}  v0={:.4}",
            h.kappa, h.v_bar, h.sigma, h.rho, h.v0
        );
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
