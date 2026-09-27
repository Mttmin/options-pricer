//! Deep calibration of the Grasselli 4/2 model through an ONNX surrogate.
//!
//! Input `parameters` is f32 `[batch, 9]`:
//! kappa, theta, sigma_v, rho, v0, a, b, r_norm, q_norm.
//! Output `iv_surface` is f32 `[batch, 686]`, the same 49×14 grid as the
//! Heston surrogate. `q_norm` is pinned at 0; physical dividend yield is
//! folded into the carry `r − q`, matching the training generator.
//!
use crate::deep_cal::{self, N_FLAT};

use ort::session::Session;
use ort::session::builder::GraphOptimizationLevel;
use ort::value::Tensor;

use crate::slv::CalibrationDataset;

pub const N_PARAMS: usize = 7;
const N_MODEL_INPUTS: usize = N_PARAMS + 2;

/// Physical lower bounds: kappa, theta, sigma_v, rho, v0, a, b.
const PARAM_LO: [f32; N_PARAMS] = [0.30, 0.01, 0.10, -0.90, 0.02, 0.40, 0.00];
/// Physical upper bounds.
const PARAM_HI: [f32; N_PARAMS] = [10.0, 0.16, 1.50, -0.30, 0.20, 1.60, 0.04];
const R_LO: f32 = 0.00;
const R_HI: f32 = 0.05;

#[derive(Debug, Clone, Copy)]
pub struct FourTwoParameters {
    pub kappa: f64,
    pub v_bar: f64,
    pub sigma: f64,
    pub rho: f64,
    pub v0: f64,
    pub a: f64,
    pub b: f64,
}

#[derive(Debug, Clone)]
pub struct FourTwoCalibrationConfig {
    pub reg_weight: f32,
    /// Multipliers on `reg_weight`: kappa, theta, sigma, rho, v0, a, b.
    /// `b` is pulled harder so a Heston-like surface does not grow a spurious
    /// 3/2 loading. kappa and sigma stay heavy, as in the Heston calibrator.
    pub reg_per_param: [f32; N_PARAMS],
    pub use_vega_weights: bool,
}

impl Default for FourTwoCalibrationConfig {
    fn default() -> Self {
        Self {
            reg_weight: 3.0,
            reg_per_param: [1.5, 0.5, 1.5, 1.0, 0.3, 1.2, 1.8],
            use_vega_weights: true,
        }
    }
}

#[derive(Debug, Clone)]
pub struct FourTwoCalibrationResult {
    pub params: FourTwoParameters,
    pub norm: [f32; N_PARAMS],
    pub ivrmse: f32,
    pub n_evals: usize,
}

pub struct FourTwoCalibrator {
    session: Session,
}

pub fn default_onnx_path() -> Option<String> {
    let workspace = concat!(env!("CARGO_MANIFEST_DIR"), "/../weights/four_two_surrogate.onnx");
    if std::path::Path::new(workspace).exists() {
        return Some(workspace.to_string());
    }
    let dev = "/home/mttmin/coding/deep-calibration/model/four_two_surrogate.onnx";
    if std::path::Path::new(dev).exists() {
        return Some(dev.to_string());
    }
    None
}

fn normalise_r(carry: f32) -> f32 {
    ((carry - R_LO) / (R_HI - R_LO)).clamp(0.0, 1.0)
}

fn denormalise(norm: &[f32; N_PARAMS]) -> [f32; N_PARAMS] {
    let mut out = [0f32; N_PARAMS];
    for i in 0..N_PARAMS {
        out[i] = PARAM_LO[i] + norm[i].clamp(0.0, 1.0) * (PARAM_HI[i] - PARAM_LO[i]);
    }
    out
}

fn to_params(phys: &[f32; N_PARAMS]) -> FourTwoParameters {
    FourTwoParameters {
        kappa: phys[0] as f64,
        v_bar: phys[1] as f64,
        sigma: phys[2] as f64,
        rho: phys[3] as f64,
        v0: phys[4] as f64,
        a: phys[5] as f64,
        b: phys[6] as f64,
    }
}

fn market_prior_norm() -> [f32; N_PARAMS] {
    // kappa=3, theta=0.04, sigma=0.40, rho=-0.70, v0=0.04, a=1, b=0.004
    let phys = [3.0, 0.04, 0.40, -0.70, 0.04, 1.0, 0.004];
    let mut out = [0f32; N_PARAMS];
    for i in 0..N_PARAMS {
        out[i] = (phys[i] - PARAM_LO[i]) / (PARAM_HI[i] - PARAM_LO[i]);
    }
    out
}

fn dataset_rq(dataset: &CalibrationDataset) -> (f32, f32) {
    let Some(first) = dataset.instruments.first() else {
        return (0.0, 0.0);
    };
    (first.risk_free_rate as f32, first.dividend_yield as f32)
}

impl FourTwoCalibrator {
    pub fn new(onnx_path: &str) -> ort::Result<Self> {
        let mut session = Session::builder()?
            .with_optimization_level(GraphOptimizationLevel::Level3)?
            .commit_from_file(onnx_path)?;
        let dummy = Tensor::from_array(([1usize, N_MODEL_INPUTS], vec![0.5f32; N_MODEL_INPUTS]))?;
        let outputs = session.run(ort::inputs!["parameters" => dummy]).unwrap_or_else(|_| {
            panic!(
                "ONNX model rejected input shape [1, {N_MODEL_INPUTS}]. \
                 Re-export four_two_surrogate.onnx."
            )
        });
        let (_shape, slice) = outputs["iv_surface"]
            .try_extract_tensor::<f32>()
            .expect("ONNX output 'iv_surface' missing or not f32");
        assert_eq!(slice.len(), N_FLAT, "ONNX iv_surface length {}, expected {N_FLAT}", slice.len());
        drop(outputs);
        Ok(Self { session })
    }

    fn forward(&mut self, theta: &[f32; N_PARAMS], r_norm: f32, q_norm: f32) -> ort::Result<[f32; N_FLAT]> {
        let mut inputs = [0f32; N_MODEL_INPUTS];
        inputs[..N_PARAMS].copy_from_slice(theta);
        inputs[N_PARAMS] = r_norm;
        inputs[N_PARAMS + 1] = q_norm;
        let tensor = Tensor::from_array(([1usize, N_MODEL_INPUTS], inputs.to_vec()))?;
        let outputs = self.session.run(ort::inputs!["parameters" => tensor])?;
        let (_shape, slice) = outputs["iv_surface"].try_extract_tensor::<f32>()?;
        let mut out = [0f32; N_FLAT];
        out.copy_from_slice(&slice[..N_FLAT]);
        Ok(out)
    }

    pub fn calibrate_with_weights(
        &mut self,
        iv_market: &[f32; N_FLAT],
        weights: &[f32; N_FLAT],
        r_norm: f32,
        q_norm: f32,
        n_restarts: usize,
        config: &FourTwoCalibrationConfig,
    ) -> ort::Result<FourTwoCalibrationResult> {
        let seeds = make_seeds(n_restarts.max(1));
        let mut best: Option<(f32, [f32; N_PARAMS], usize)> = None;
        let mut total_evals = 0usize;
        for seed in &seeds {
            let (theta, iv_loss, evals) =
                self.levenberg_marquardt(seed, r_norm, q_norm, iv_market, weights, config)?;
            total_evals += evals;
            let rmse = iv_loss.sqrt();
            if best.as_ref().map_or(true, |b| rmse < b.0) {
                best = Some((rmse, theta, total_evals));
            }
        }
        let (rmse, norm, n_evals) = best.expect("at least one restart");
        Ok(FourTwoCalibrationResult {
            params: to_params(&denormalise(&norm)),
            norm,
            ivrmse: rmse,
            n_evals,
        })
    }

    fn levenberg_marquardt(
        &mut self,
        init: &[f32; N_PARAMS],
        r_norm: f32,
        q_norm: f32,
        iv_market: &[f32; N_FLAT],
        weights: &[f32; N_FLAT],
        config: &FourTwoCalibrationConfig,
    ) -> ort::Result<([f32; N_PARAMS], f32, usize)> {
        const MAX_ITERS: usize = 60;
        const EPS: [f32; N_PARAMS] = [8e-4, 1e-4, 1e-4, 1e-4, 1e-4, 2e-4, 4e-4];
        const DAMP_INIT: f32 = 1e-2;
        const DAMP_UP: f32 = 10.0;
        const DAMP_DOWN: f32 = 0.33;
        const STEP_TOL: f32 = 1e-6;

        let prior = market_prior_norm();
        let mut reg_per = [0f32; N_PARAMS];
        for i in 0..N_PARAMS {
            reg_per[i] = config.reg_weight * config.reg_per_param[i];
        }
        let clamp = |v: [f32; N_PARAMS]| -> [f32; N_PARAMS] {
            let mut c = v;
            for x in &mut c {
                *x = x.clamp(0.0, 1.0);
            }
            c
        };

        let mut theta = clamp(*init);
        let mut iv_pred = self.forward(&theta, r_norm, q_norm)?;
        let mut cost = lm_cost(&iv_pred, iv_market, weights, &theta, &prior, &reg_per);
        let mut damping = DAMP_INIT;
        let mut n_evals = 1usize;
        let valid: Vec<usize> = (0..N_FLAT).filter(|&i| weights[i] > 0.0).collect();

        for _ in 0..MAX_ITERS {
            if valid.is_empty() {
                break;
            }
            let n_valid = valid.len();
            let mut jac = vec![0f32; N_PARAMS * n_valid];
            for j in 0..N_PARAMS {
                let mut th_b = theta;
                let step = if theta[j] + EPS[j] <= 1.0 { EPS[j] } else { -EPS[j] };
                th_b[j] = (theta[j] + step).clamp(0.0, 1.0);
                let eps_actual = th_b[j] - theta[j];
                if eps_actual.abs() < 1e-10 {
                    continue;
                }
                let iv_b = self.forward(&th_b, r_norm, q_norm)?;
                n_evals += 1;
                for (k, &i) in valid.iter().enumerate() {
                    jac[j * n_valid + k] = (iv_b[i] - iv_pred[i]) / eps_actual;
                }
            }

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
                rhs += reg_per[a] * (theta[a] - prior[a]);
                g[a] = rhs;
                h[a][a] += reg_per[a] + damping;
            }

            let Some(delta) = solve_nxn(&h, &g.map(|x| -x)) else {
                break;
            };
            let step_norm = delta.iter().map(|d| d * d).sum::<f32>().sqrt();
            if step_norm < STEP_TOL {
                break;
            }
            let mut candidate = theta;
            for i in 0..N_PARAMS {
                candidate[i] += delta[i];
            }
            let candidate = clamp(candidate);
            let iv_cand = self.forward(&candidate, r_norm, q_norm)?;
            n_evals += 1;
            let cost_cand = lm_cost(&iv_cand, iv_market, weights, &candidate, &prior, &reg_per);
            if cost_cand < cost {
                theta = candidate;
                iv_pred = iv_cand;
                cost = cost_cand;
                damping = (damping * DAMP_DOWN).max(1e-10);
            } else {
                damping = (damping * DAMP_UP).min(1e8);
            }
        }
        Ok((theta, weighted_iv_mse(&iv_pred, iv_market, weights), n_evals))
    }
}

pub fn calibrate_four_two_from_dataset(
    calibrator: &mut FourTwoCalibrator,
    dataset: &CalibrationDataset,
    n_restarts: usize,
    config: &FourTwoCalibrationConfig,
) -> ort::Result<FourTwoCalibrationResult> {
    let (surface, mask, confidence) = deep_cal::build_iv_surface(dataset);
    let weights = if config.use_vega_weights {
        deep_cal::build_vega_weights(&surface, &mask, &confidence)
    } else {
        deep_cal::build_weights(&surface, &mask, &confidence, Some(0.30))
    };
    let (r, q) = dataset_rq(dataset);
    calibrator.calibrate_with_weights(&surface, &weights, normalise_r(r - q), 0.0, n_restarts, config)
}

fn lm_cost(
    iv_pred: &[f32; N_FLAT],
    iv_market: &[f32; N_FLAT],
    weights: &[f32; N_FLAT],
    theta: &[f32; N_PARAMS],
    prior: &[f32; N_PARAMS],
    reg_per: &[f32; N_PARAMS],
) -> f32 {
    let mut cost = weighted_iv_mse(iv_pred, iv_market, weights);
    for i in 0..N_PARAMS {
        let d = theta[i] - prior[i];
        cost += reg_per[i] * d * d;
    }
    cost
}

fn weighted_iv_mse(iv_pred: &[f32; N_FLAT], iv_market: &[f32; N_FLAT], weights: &[f32; N_FLAT]) -> f32 {
    let mut num = 0f32;
    let mut den = 0f32;
    for i in 0..N_FLAT {
        if weights[i] > 0.0 {
            let d = iv_pred[i] - iv_market[i];
            num += weights[i] * d * d;
            den += weights[i];
        }
    }
    if den == 0.0 { 0.0 } else { num / den }
}

fn solve_nxn(a: &[[f32; N_PARAMS]; N_PARAMS], b: &[f32; N_PARAMS]) -> Option<[f32; N_PARAMS]> {
    let mut m = [[0f32; N_PARAMS + 1]; N_PARAMS];
    for i in 0..N_PARAMS {
        m[i][..N_PARAMS].copy_from_slice(&a[i]);
        m[i][N_PARAMS] = b[i];
    }
    for col in 0..N_PARAMS {
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

fn make_seeds(n: usize) -> Vec<[f32; N_PARAMS]> {
    let prior = market_prior_norm();
    let mut heston_limit = prior;
    heston_limit[5] = (1.0 - PARAM_LO[5]) / (PARAM_HI[5] - PARAM_LO[5]);
    heston_limit[6] = 0.0;
    let mut low_k = prior;
    low_k[0] = (0.8 - PARAM_LO[0]) / (PARAM_HI[0] - PARAM_LO[0]);
    let mut high_k = prior;
    high_k[0] = (6.0 - PARAM_LO[0]) / (PARAM_HI[0] - PARAM_LO[0]);
    let mut high_b = prior;
    high_b[6] = 0.55;
    let mut low_a = prior;
    low_a[5] = 0.15;
    let mut steep = prior;
    steep[3] = 0.05;
    let catalog = [
        [0.5; N_PARAMS],
        prior,
        heston_limit,
        low_k,
        high_k,
        high_b,
        low_a,
        steep,
    ];
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        if i < catalog.len() {
            out.push(catalog[i]);
        } else {
            let t = i as f32 / n as f32;
            out.push([t, 0.45, 0.35, 0.30, 0.40, 0.50, 0.15]);
        }
    }
    out
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounds_round_trip_and_solver() {
        let phys = [3.0, 0.04, 0.40, -0.70, 0.04, 1.0, 0.004];
        let mut norm = [0f32; N_PARAMS];
        for i in 0..N_PARAMS {
            norm[i] = (phys[i] - PARAM_LO[i]) / (PARAM_HI[i] - PARAM_LO[i]);
        }
        let back = denormalise(&norm);
        for i in 0..N_PARAMS {
            assert!((back[i] - phys[i]).abs() < 1e-5, "param {i}");
        }
        let mut a = [[0f32; N_PARAMS]; N_PARAMS];
        let mut b = [0f32; N_PARAMS];
        for i in 0..N_PARAMS {
            a[i][i] = 2.0;
            b[i] = 2.0 * (i as f32 + 1.0);
        }
        let x = solve_nxn(&a, &b).unwrap();
        for i in 0..N_PARAMS {
            assert!((x[i] - (i as f32 + 1.0)).abs() < 1e-5);
        }
    }

    #[test]
    fn onnx_round_trip_recovers_parameters() {
        let Some(path) = default_onnx_path() else {
            return;
        };
        let mut cal = FourTwoCalibrator::new(&path).expect("load 4/2 onnx");
        let mut theta = market_prior_norm();
        theta[0] = 0.22;
        theta[3] = 0.15;
        theta[5] = 0.62;
        theta[6] = 0.28;
        let r_norm = 0.40;
        let iv = cal.forward(&theta, r_norm, 0.0).expect("forward");
        assert!(iv.iter().all(|x| x.is_finite() && (0.01..1.5).contains(x)));
        let weights = [1.0f32; N_FLAT];
        let got = cal
            .calibrate_with_weights(
                &iv,
                &weights,
                r_norm,
                0.0,
                4,
                &FourTwoCalibrationConfig {
                    reg_weight: 0.0,
                    ..FourTwoCalibrationConfig::default()
                },
            )
            .expect("calibrate");
        assert!(
            got.ivrmse < 0.002,
            "repriced surface IVRMSE {} exceeds 20 bps; got={:?} true={theta:?}",
            got.ivrmse,
            got.norm,
        );
    }
}
