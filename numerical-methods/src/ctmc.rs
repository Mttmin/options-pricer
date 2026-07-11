use nalgebra::DMatrix;
use options::Options;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use crate::expm::MatExpOperator;

// SLV model interface
trait SLVModel {
    fn omega(&self, s: f64, v: f64) -> f64;
    fn m(&self, v: f64) -> f64;
    fn gamma(&self, s: f64) -> f64;
    fn mu(&self, current_vol: f64) -> f64;
    fn sigma_v(&self, v: f64) -> f64;
    fn rho(&self) -> f64;
}

pub struct Heston {
    rho: f64,
    sigma: f64,
    r_f: f64,
    q: Option<f64>,
    eta: f64,
    ref_vol: f64,
}

impl Heston {
    pub fn new(rho: f64, sigma: f64, r_f: f64, q: Option<f64>, eta: f64, ref_vol: f64) -> Self {
        Self { rho, sigma, r_f, q, eta, ref_vol }
    }
    pub fn from_options(option: &Options) -> Self {
        // Placeholders: calibrate eta, rho, sigma, ref_vol to market data in practice.
        // ref_vol is the long-run VARIANCE θ (mu(v) = η(θ − v) operates in variance
        // units), so the option's implied vol must be squared here.
        Self {
            rho: 0.0,
            sigma: option.volatility(),
            r_f: option.risk_free_rate(),
            q: option.dividend_yield(),
            eta: 1.0,
            ref_vol: option.volatility().powi(2),
        }
    }
    fn heston_g_integral(&self, s: f64) -> f64 { s.ln() }
    fn heston_f_integral(&self, v: f64) -> f64 { v / self.sigma }
    fn heston_g_inverse(&self, x: f64) -> f64 { x.exp() }
    fn heston_f_inverse(&self, x: f64) -> f64 { x * self.sigma }
    fn heston_m_prime(&self, v: f64) -> f64 { 0.5 / v.sqrt() }
    fn heston_sigma_prime_v(&self, v: f64) -> f64 { 0.5 * self.sigma / v.sqrt() }
    fn heston_gamma_prime(&self, _s: f64) -> f64 { 1.0 }
    // Uniform dispatch names (used by the dispatch! macro in SLVModelVariant)
    fn g_integral(&self, s: f64) -> f64 { self.heston_g_integral(s) }
    fn f_integral(&self, v: f64) -> f64 { self.heston_f_integral(v) }
    fn g_inverse(&self, x: f64) -> f64  { self.heston_g_inverse(x) }
    fn f_inverse(&self, y: f64) -> f64  { self.heston_f_inverse(y) }
    fn m_prime(&self, v: f64) -> f64    { self.heston_m_prime(v) }
    fn sigma_prime_v(&self, v: f64) -> f64 { self.heston_sigma_prime_v(v) }
    fn gamma_prime(&self, s: f64) -> f64   { self.heston_gamma_prime(s) }
}

impl SLVModel for Heston {
    fn omega(&self, s: f64, _v: f64) -> f64 { s * (self.r_f - self.q.unwrap_or(0.0)) }
    fn m(&self, v: f64) -> f64 { v.sqrt() }
    fn gamma(&self, s: f64) -> f64 { s }
    fn mu(&self, v: f64) -> f64 { self.eta * (self.ref_vol - v) }
    fn sigma_v(&self, v: f64) -> f64 { self.sigma * v.sqrt() }
    fn rho(&self) -> f64 { self.rho }
}

// 4/2 model (Grasselli 2017): dS = (r-q)S dt + S[a√v + b/√v] dW^(1), dv = η(ϑ-v)dt + σ_v√v dW^(2)
pub struct FourTwo {
    pub rho: f64,
    pub sigma: f64,
    pub r_f: f64,
    pub q: Option<f64>,
    pub eta: f64,
    pub ref_vol: f64,
    pub a: f64,
    pub b: f64,
}

impl FourTwo {
    pub fn new(rho: f64, sigma: f64, r_f: f64, q: Option<f64>, eta: f64, ref_vol: f64, a: f64, b: f64) -> Self {
        Self { rho, sigma, r_f, q, eta, ref_vol, a, b }
    }
    fn g_integral(&self, s: f64) -> f64 { s.ln() }
    fn f_integral(&self, v: f64) -> f64 { (self.a * v + self.b * v.ln()) / self.sigma }
    fn g_inverse(&self, x: f64) -> f64 { x.exp() }
    fn f_inverse(&self, y: f64) -> f64 {
        // Newton: solve a*v + b*ln(v) = sigma*y
        let target = self.sigma * y;
        let mut v = (target / self.a).max(1e-6);
        for _ in 0..60 {
            let fv = self.a * v + self.b * v.ln() - target;
            let dfv = self.a + self.b / v;
            let dv = fv / dfv;
            v = (v - dv).max(1e-10);
            if dv.abs() < 1e-12 * v { break; }
        }
        v
    }
    fn m_prime(&self, v: f64) -> f64 { (self.a * v - self.b) / (2.0 * v.powf(1.5)) }
    fn sigma_prime_v(&self, v: f64) -> f64 { 0.5 * self.sigma / v.sqrt() }
    fn gamma_prime(&self, _s: f64) -> f64 { 1.0 }
}

impl SLVModel for FourTwo {
    fn omega(&self, s: f64, _v: f64) -> f64 { s * (self.r_f - self.q.unwrap_or(0.0)) }
    fn m(&self, v: f64) -> f64 { self.a * v.sqrt() + self.b / v.sqrt() }
    fn gamma(&self, s: f64) -> f64 { s }
    fn mu(&self, v: f64) -> f64 { self.eta * (self.ref_vol - v) }
    fn sigma_v(&self, v: f64) -> f64 { self.sigma * v.sqrt() }
    fn rho(&self) -> f64 { self.rho }
}

// α-Hyper (Da Fonseca and Martini 2016): dS = (r-q)S dt + S·exp(v) dW^(1), dv = (η−ϑ·exp(α·v))dt + σ_v dW^(2)
pub struct AlphaHyper {
    pub rho: f64,
    pub sigma: f64,
    pub r_f: f64,
    pub q: Option<f64>,
    pub eta: f64,
    pub theta_lr: f64,
    pub alpha: f64,
}

impl AlphaHyper {
    pub fn new(rho: f64, sigma: f64, r_f: f64, q: Option<f64>, eta: f64, theta_lr: f64, alpha: f64) -> Self {
        Self { rho, sigma, r_f, q, eta, theta_lr, alpha }
    }
    fn g_integral(&self, s: f64) -> f64 { s.ln() }
    fn f_integral(&self, v: f64) -> f64 { v.exp() / self.sigma }
    fn g_inverse(&self, x: f64) -> f64 { x.exp() }
    fn f_inverse(&self, y: f64) -> f64 { (self.sigma * y).ln() }
    fn m_prime(&self, v: f64) -> f64 { v.exp() }
    fn sigma_prime_v(&self, _v: f64) -> f64 { 0.0 }
    fn gamma_prime(&self, _s: f64) -> f64 { 1.0 }
}

impl SLVModel for AlphaHyper {
    fn omega(&self, s: f64, _v: f64) -> f64 { s * (self.r_f - self.q.unwrap_or(0.0)) }
    fn m(&self, v: f64) -> f64 { v.exp() }
    fn gamma(&self, s: f64) -> f64 { s }
    fn mu(&self, v: f64) -> f64 { self.eta - self.theta_lr * (self.alpha * v).exp() }
    fn sigma_v(&self, _v: f64) -> f64 { self.sigma }
    fn rho(&self) -> f64 { self.rho }
}

// SABR (Hagan et al. 2002): dS = (r-q)S dt + v·S^β dW^(1), dv = σ_v·v dW^(2)
pub struct SABR {
    pub rho: f64,
    pub sigma: f64,
    pub r_f: f64,
    pub q: Option<f64>,
    pub beta: f64,
}

impl SABR {
    pub fn new(rho: f64, sigma: f64, r_f: f64, q: Option<f64>, beta: f64) -> Self {
        assert!((0.0..=1.0).contains(&beta), "SABR beta must be in [0, 1]");
        Self { rho, sigma, r_f, q, beta }
    }
    fn g_integral(&self, s: f64) -> f64 {
        if (self.beta - 1.0).abs() < 1e-10 { s.ln() }
        else { s.powf(1.0 - self.beta) / (1.0 - self.beta) }
    }
    fn f_integral(&self, v: f64) -> f64 { v / self.sigma }
    fn g_inverse(&self, x: f64) -> f64 {
        if (self.beta - 1.0).abs() < 1e-10 { x.exp() }
        else { (x * (1.0 - self.beta)).powf(1.0 / (1.0 - self.beta)) }
    }
    fn f_inverse(&self, y: f64) -> f64 { y * self.sigma }
    fn m_prime(&self, _v: f64) -> f64 { 1.0 }
    fn sigma_prime_v(&self, _v: f64) -> f64 { self.sigma }
    fn gamma_prime(&self, s: f64) -> f64 { self.beta * s.powf(self.beta - 1.0) }
}

impl SLVModel for SABR {
    fn omega(&self, s: f64, _v: f64) -> f64 { s * (self.r_f - self.q.unwrap_or(0.0)) }
    fn m(&self, v: f64) -> f64 { v }
    fn gamma(&self, s: f64) -> f64 { s.powf(self.beta) }
    fn mu(&self, _v: f64) -> f64 { 0.0 }
    fn sigma_v(&self, v: f64) -> f64 { self.sigma * v }
    fn rho(&self) -> f64 { self.rho }
}

// Heston-SABR (Van der Stoep et al. 2014): dS = (r-q)S dt + S·v dW^(1), dv = η(ϑ-v)dt + σ_v·v dW^(2)
pub struct HestonSABR {
    pub rho: f64,
    pub sigma: f64,
    pub r_f: f64,
    pub q: Option<f64>,
    pub eta: f64,
    pub ref_vol: f64,
}

impl HestonSABR {
    pub fn new(rho: f64, sigma: f64, r_f: f64, q: Option<f64>, eta: f64, ref_vol: f64) -> Self {
        Self { rho, sigma, r_f, q, eta, ref_vol }
    }
    fn g_integral(&self, s: f64) -> f64 { s.ln() }
    fn f_integral(&self, v: f64) -> f64 { v / self.sigma }
    fn g_inverse(&self, x: f64) -> f64 { x.exp() }
    fn f_inverse(&self, y: f64) -> f64 { y * self.sigma }
    fn m_prime(&self, _v: f64) -> f64 { 1.0 }
    fn sigma_prime_v(&self, _v: f64) -> f64 { self.sigma }
    fn gamma_prime(&self, _s: f64) -> f64 { 1.0 }
}

impl SLVModel for HestonSABR {
    fn omega(&self, s: f64, _v: f64) -> f64 { s * (self.r_f - self.q.unwrap_or(0.0)) }
    fn m(&self, v: f64) -> f64 { v }
    fn gamma(&self, s: f64) -> f64 { s }
    fn mu(&self, v: f64) -> f64 { self.eta * (self.ref_vol - v) }
    fn sigma_v(&self, v: f64) -> f64 { self.sigma * v }
    fn rho(&self) -> f64 { self.rho }
}

// Quadratic SLV (Lipton 2002): dS = (r-q)S dt + √v·(aS²+bS+c) dW^(1), dv = η(ϑ-v)dt + σ_v√v dW^(2)
// Requires 4ac > b² so the quadratic is positive definite.
pub struct QuadraticSLV {
    pub rho: f64,
    pub sigma: f64,
    pub r_f: f64,
    pub q: Option<f64>,
    pub eta: f64,
    pub ref_vol: f64,
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl QuadraticSLV {
    pub fn new(rho: f64, sigma: f64, r_f: f64, q: Option<f64>, eta: f64, ref_vol: f64, a: f64, b: f64, c: f64) -> Self {
        assert!(4.0 * a * c > b * b, "QuadraticSLV requires 4ac > b²");
        Self { rho, sigma, r_f, q, eta, ref_vol, a, b, c }
    }
    fn disc(&self) -> f64 { (4.0 * self.a * self.c - self.b * self.b).sqrt() }
    fn g_integral(&self, s: f64) -> f64 {
        let d = self.disc();
        2.0 / d * ((2.0 * self.a * s + self.b) / d).atan()
    }
    fn f_integral(&self, v: f64) -> f64 { v / self.sigma }
    fn g_inverse(&self, x: f64) -> f64 {
        let d = self.disc();
        (d * (d * x / 2.0).tan() - self.b) / (2.0 * self.a)
    }
    fn f_inverse(&self, y: f64) -> f64 { y * self.sigma }
    fn m_prime(&self, v: f64) -> f64 { 0.5 / v.sqrt() }
    fn sigma_prime_v(&self, v: f64) -> f64 { 0.5 * self.sigma / v.sqrt() }
    fn gamma_prime(&self, s: f64) -> f64 { 2.0 * self.a * s + self.b }
}

impl SLVModel for QuadraticSLV {
    fn omega(&self, s: f64, _v: f64) -> f64 { s * (self.r_f - self.q.unwrap_or(0.0)) }
    fn m(&self, v: f64) -> f64 { v.sqrt() }
    fn gamma(&self, s: f64) -> f64 { self.a * s * s + self.b * s + self.c }
    fn mu(&self, v: f64) -> f64 { self.eta * (self.ref_vol - v) }
    fn sigma_v(&self, v: f64) -> f64 { self.sigma * v.sqrt() }
    fn rho(&self) -> f64 { self.rho }
}

macro_rules! dispatch {
    ($self:expr, $method:ident $(, $arg:expr)*) => {
        match $self {
            SLVModelVariant::Heston(m)      => m.$method($($arg),*),
            SLVModelVariant::FourTwo(m)     => m.$method($($arg),*),
            SLVModelVariant::AlphaHyper(m)  => m.$method($($arg),*),
            SLVModelVariant::SABR(m)        => m.$method($($arg),*),
            SLVModelVariant::HestonSABR(m)  => m.$method($($arg),*),
            SLVModelVariant::QuadraticSLV(m) => m.$method($($arg),*),
        }
    };
}

pub enum SLVModelVariant {
    Heston(Heston),
    FourTwo(FourTwo),
    AlphaHyper(AlphaHyper),
    SABR(SABR),
    HestonSABR(HestonSABR),
    QuadraticSLV(QuadraticSLV),
}

impl SLVModelVariant {
    pub fn from_options(option: &Options, model_type: SLVModelVariant) -> Self {
        match model_type {
            SLVModelVariant::Heston(_) => SLVModelVariant::Heston(Heston::from_options(option)),
            other => other,
        }
    }
    #[inline] pub fn g_integral(&self, s: f64) -> f64 { dispatch!(self, g_integral, s) }
    #[inline] pub fn f_integral(&self, v: f64) -> f64 { dispatch!(self, f_integral, v) }
    #[inline] pub fn g_inverse(&self, x: f64) -> f64  { dispatch!(self, g_inverse, x) }
    #[inline] pub fn f_inverse(&self, y: f64) -> f64  { dispatch!(self, f_inverse, y) }
    #[inline] pub fn m_prime(&self, v: f64) -> f64    { dispatch!(self, m_prime, v) }
    #[inline] pub fn sigma_prime_v(&self, v: f64) -> f64 { dispatch!(self, sigma_prime_v, v) }
    #[inline] pub fn gamma_prime(&self, s: f64) -> f64   { dispatch!(self, gamma_prime, s) }
    #[inline] pub fn rho(&self) -> f64                { dispatch!(self, rho) }
    #[inline] pub fn omega(&self, s: f64, v: f64) -> f64 { dispatch!(self, omega, s, v) }
}

impl SLVModel for SLVModelVariant {
    fn omega(&self, s: f64, v: f64) -> f64  { dispatch!(self, omega, s, v) }
    fn m(&self, v: f64) -> f64              { dispatch!(self, m, v) }
    fn gamma(&self, s: f64) -> f64          { dispatch!(self, gamma, s) }
    fn mu(&self, v: f64) -> f64             { dispatch!(self, mu, v) }
    fn sigma_v(&self, v: f64) -> f64        { dispatch!(self, sigma_v, v) }
    fn rho(&self) -> f64                    { dispatch!(self, rho) }
}

// Grid construction

// Uniform variance grid v ∈ [step, 4·v0] with n steps, avoiding v=0.
//
// DEPRECATED for Heston pricing: this grid ignores the long-run variance θ, so
// whenever θ > 4·v0 the CTMC chain can never reach its own mean-reversion level
// and long-dated options are underpriced (−41% at θ = 4·v0, −95% at θ = 8·v0 in
// the European CTMC-vs-COS benchmark). Prefer `heston_variance_grid`.
pub fn uniform_variance_grid(vol_0: f64, n_vol_steps: i32) -> Vec<f64> {
    let step = 4.0 * vol_0 / (n_vol_steps as f64);
    (1..=n_vol_steps).map(|i| i as f64 * step).collect()
}

/// CIR conditional moments of v_T | v_0 = v0: returns (mean, variance).
///
///   E[v_T]   = θ + (v0 − θ)·e^{−κT}
///   Var[v_T] = v0·σ_v²·e^{−κT}(1 − e^{−κT})/κ + θ·σ_v²(1 − e^{−κT})²/(2κ)
fn cir_conditional_moments(v0: f64, kappa: f64, theta_lr: f64, sigma_v: f64, t: f64) -> (f64, f64) {
    let e = (-kappa * t).exp();
    let one_minus_e = -(-kappa * t).exp_m1();
    // e·(1−e)/κ → t and (1−e)²/(2κ) → κt²/2 as κ→0; guard only the division.
    let k = kappa.max(1e-12);
    let mean = theta_lr + (v0 - theta_lr) * e;
    let var = v0 * sigma_v * sigma_v * e * one_minus_e / k
        + theta_lr * sigma_v * sigma_v * one_minus_e * one_minus_e / (2.0 * k);
    (mean, var)
}

/// Variance grid sized from the CIR conditional moments of v_T | v_0
/// (Cui, Kirkby & Nguyen 2018 recommend moment-matched variance grids).
///
/// Uniformly spaced over [lo, hi] with
///   hi = max(v0, E[v_T]) + L·sd(v_T)
///   lo = max(hi·1e-4, min(v0, E[v_T]) − L·sd(v_T)),  L = 4,
/// so the grid always reaches the mean-reversion level θ regardless of v0.
/// Empirically (European CTMC vs COS benchmark, 80×20×50 grid) this cuts the
/// θ = 4·v0 error from −41% to ~+2% and the θ = 8·v0 error from −95% to ~+2%.
pub fn heston_variance_grid(
    v0: f64, kappa: f64, theta_lr: f64, sigma_v: f64, t: f64, n_vol_steps: usize,
) -> Vec<f64> {
    let (mean_t, var_t) = cir_conditional_moments(v0, kappa, theta_lr, sigma_v, t);
    let sd = var_t.max(0.0).sqrt();
    const L: f64 = 4.0;
    let hi = v0.max(mean_t) + L * sd;
    let mut lo = (v0.min(mean_t) - L * sd).max(hi * 1e-4);
    if lo >= v0 {
        // Degenerate (sd ≈ 0 with mean ≥ v0): keep v0 strictly interior.
        lo = 0.5 * v0;
    }
    let n = n_vol_steps.max(2);
    (0..n)
        .map(|i| lo + (hi - lo) * i as f64 / (n - 1) as f64)
        .collect()
}

/// X-grid half-width from the moments of X_T = g(S_T) − ρ·f(v_T).
///
/// Covers (a) the diffusion of X (std ≈ √((1−ρ²)·E[∫v dt])) with an L_x = 7
/// margin, (b) the drift terms |r−q|·T + ½·E[∫v dt], and (c) the spread of the
/// payoff-kink position ln K − ρ·f(v_l) across all variance states. The
/// paper's fixed range [g(0.001·S), g(4·S)] spends most of its n_x points on
/// probabilistically dead regions; at n_x = 80 that dominated the pricing
/// error once the variance grid was widened to reach θ.
fn heston_x_half_width(
    v0: f64, kappa: f64, theta_lr: f64, sigma_v: f64, rho: f64,
    r: f64, q: f64, t: f64, v_lo: f64, v_hi: f64,
) -> f64 {
    let one_minus_e = -(-kappa * t).exp_m1();
    let k = kappa.max(1e-12);
    let integ_var = theta_lr * t + (v0 - theta_lr) * one_minus_e / k; // E[∫₀ᵀ v_s ds]
    let std_x = ((1.0 - rho * rho) * integ_var).max(0.0).sqrt();
    let kink_span = rho.abs() * (v_hi - v_lo) / sigma_v; // |ρ|·(f(v_hi) − f(v_lo))
    const L_X: f64 = 7.0;
    (kink_span + 0.5 * integ_var + (r - q).abs() * t + L_X * std_x).max(0.75)
}

#[allow(dead_code)]
fn x_space_grid(slv_model: &SLVModelVariant, m_steps: i32, spot_0: f64, vol_0: f64) -> (Vec<f64>, Vec<f64>) {
    let x_0 = slv_model.g_integral(spot_0) - slv_model.rho() * slv_model.f_integral(vol_0);
    if x_0 <= 0.0 {
        panic!("Invalid initial x space grid point: x_0 must be positive");
    }
    let step = x_0 * (4.0 - 1e-3) / (m_steps as f64);
    let grid = (1..=m_steps).map(|i| i as f64 * step).collect();
    (grid, vec![step; m_steps as usize])
}

// Generator assembly
pub struct VarianceGeneratorCOO {
    pub rows: Vec<usize>,
    pub cols: Vec<usize>,
    pub vals: Vec<f64>,
    pub num_states: usize,
}

fn build_variance_generator(slv_model: &SLVModelVariant, v_grid: &[f64]) -> VarianceGeneratorCOO {
    let m = v_grid.len();
    // Interior 3 entries + 2-state reflecting boundaries at each end.
    let capacity = 3 * m.saturating_sub(2) + 4;
    let mut rows = Vec::with_capacity(capacity);
    let mut cols = Vec::with_capacity(capacity);
    let mut vals = Vec::with_capacity(capacity);

    for l in 1..(m - 1) {
        let dv_up = v_grid[l + 1] - v_grid[l];
        let dv_dn = v_grid[l] - v_grid[l - 1];
        let dv_sum = dv_up + dv_dn;
        let sigma_sq = slv_model.sigma_v(v_grid[l]).powi(2);
        let mu = slv_model.mu(v_grid[l]);

        let lambda_up = ((sigma_sq + mu * dv_dn) / (dv_up * dv_sum)).max(0.0);
        let lambda_dn = ((sigma_sq - mu * dv_up) / (dv_dn * dv_sum)).max(0.0);

        rows.extend([l, l, l]);
        cols.extend([l + 1, l - 1, l]);
        vals.extend([lambda_up, lambda_dn, -(lambda_up + lambda_dn)]);
    }

    // Reflecting boundaries: reuse interior formula with dv_dn=dv_up at bottom
    // and dv_up=dv_dn at top (no-flux reflection). Absorbing endpoints (no rates)
    // leak probability mass and systematically bias prices.
    if m >= 2 {
        let dv = v_grid[1] - v_grid[0];
        let sigma_sq = slv_model.sigma_v(v_grid[0]).powi(2);
        let mu = slv_model.mu(v_grid[0]);
        let lambda_up = ((sigma_sq + mu * dv) / (2.0 * dv * dv)).max(0.0);
        rows.extend([0, 0]);
        cols.extend([1, 0]);
        vals.extend([lambda_up, -lambda_up]);

        let dv = v_grid[m - 1] - v_grid[m - 2];
        let sigma_sq = slv_model.sigma_v(v_grid[m - 1]).powi(2);
        let mu = slv_model.mu(v_grid[m - 1]);
        let lambda_dn = ((sigma_sq - mu * dv) / (2.0 * dv * dv)).max(0.0);
        rows.extend([m - 1, m - 1]);
        cols.extend([m - 2, m - 1]);
        vals.extend([lambda_dn, -lambda_dn]);
    }

    VarianceGeneratorCOO { rows, cols, vals, num_states: m }
}

pub struct VarianceGeneratorDense {
    pub lambda: DMatrix<f64>,
}

impl VarianceGeneratorDense {
    pub fn from_coo(coo: &VarianceGeneratorCOO) -> Self {
        let m = coo.num_states;
        let mut lambda = DMatrix::zeros(m, m);
        for ((&r, &c), &v) in coo.rows.iter().zip(&coo.cols).zip(&coo.vals) {
            lambda[(r, c)] = v;
        }
        Self { lambda }
    }
}

fn sigma_tilde_sq(slv_model: &SLVModelVariant, vol_grid: &[f64]) -> Vec<f64> {
    let rho = slv_model.rho();
    let rho_factor = 1.0 - rho * rho;
    vol_grid.iter().map(|&v| {
        let m = slv_model.m(v);
        m * m * rho_factor
    }).collect()
}

fn psi_v(slv_model: &SLVModelVariant, v: f64) -> f64 {
    let m = slv_model.m(v);
    let sigma = slv_model.sigma_v(v);
    slv_model.mu(v) * m / sigma
        + 0.5 * (sigma * slv_model.m_prime(v) - slv_model.sigma_prime_v(v) * m)
}

fn theta(slv_model: &SLVModelVariant, x_grid: &[f64], vol_grid: &[f64]) -> Vec<Vec<f64>> {
    let rho = slv_model.rho();
    vol_grid.iter().map(|&vol| {
        let f_vl  = slv_model.f_integral(vol);
        let psi   = psi_v(slv_model, vol);
        let m_sq  = slv_model.m(vol).powi(2);
        let rho_psi = rho * psi;
        x_grid.iter().map(|&x| {
            let s = slv_model.g_inverse(x + rho * f_vl);
            slv_model.omega(s, vol) / slv_model.gamma(s)
                - slv_model.gamma_prime(s) * m_sq * 0.5
                - rho_psi
        }).collect()
    }).collect()
}

// N×N spatial generator for variance state v_l.
// Boundary rows absorbing; off-diagonal rates clamped if grid is too coarse.
fn build_g_l_matrix_inner(x_spacings: &[f64], theta_col: &[f64], sigma_tilde_sq_l: f64, n: usize) -> DMatrix<f64> {
    let mut g_l = DMatrix::zeros(n, n);
    for i in 1..(n - 1) {
        let d_i = x_spacings[i];
        let d_im1 = x_spacings[i - 1];
        let d_sum = d_i + d_im1;
        let th = theta_col[i];
        let q_up = ((sigma_tilde_sq_l + th * d_im1) / (d_i * d_sum)).max(0.0);
        let q_dn = ((sigma_tilde_sq_l - th * d_i) / (d_im1 * d_sum)).max(0.0);
        g_l[(i, i + 1)] = q_up;
        g_l[(i, i - 1)] = q_dn;
        g_l[(i, i)] = -(q_up + q_dn);
    }
    g_l
}

// Combined NM×NM generator G = Λ⊗I_N + blockdiag(G_1,...,G_M)
fn build_combined_generator(variance_gen: &VarianceGeneratorDense, g_l_matrices: &[DMatrix<f64>], n: usize, m: usize) -> DMatrix<f64> {
    let nm = n * m;
    let mut g = DMatrix::zeros(nm, nm);

    for l in 0..m {
        for lp in 0..m {
            let lambda = variance_gen.lambda[(l, lp)];
            if lambda.abs() < 1e-30 { continue; }
            for i in 0..n {
                g[(l * n + i, lp * n + i)] += lambda;
            }
        }
    }

    for l in 0..m {
        let g_l = &g_l_matrices[l];
        debug_assert_eq!(g_l.nrows(), n);
        debug_assert_eq!(g_l.ncols(), n);
        for i in 0..n {
            for j in 0..n {
                let val = g_l[(i, j)];
                if val.abs() < 1e-30 { continue; }
                g[(l * n + i, l * n + j)] += val;
            }
        }
    }
    g
}

pub struct CombinedCTMC {
    pub generator: DMatrix<f64>,
    pub x_grid: Vec<f64>,
    pub x_spacings: Vec<f64>,
    pub v_grid: Vec<f64>,
    pub n: usize,
    pub m: usize,
}

impl CombinedCTMC {
    #[inline]
    pub fn phi(&self, i: usize, l: usize) -> usize { l * self.n + i }

    #[inline]
    pub fn phi_inv(&self, idx: usize) -> (usize, usize) { (idx % self.n, idx / self.n) }

    #[inline]
    pub fn dim(&self) -> usize { self.n * self.m }

    pub fn selector_vec(&self, i: usize, l: usize) -> Vec<f64> {
        let mut e = vec![0.0; self.dim()];
        e[self.phi(i, l)] = 1.0;
        e
    }

    pub fn max_row_sum_error(&self) -> f64 {
        let nm = self.dim();
        (0..nm)
            .map(|row| (0..nm).map(|col| self.generator[(row, col)]).sum::<f64>().abs())
            .fold(0.0_f64, f64::max)
    }

    pub fn check_off_diagonal_non_negative(&self) -> (usize, f64) {
        let nm = self.dim();
        let mut count = 0;
        let mut worst = 0.0_f64;
        for row in 0..nm {
            for col in 0..nm {
                if row != col {
                    let val = self.generator[(row, col)];
                    if val < -1e-14 {
                        count += 1;
                        worst = worst.min(val);
                    }
                }
            }
        }
        (count, worst)
    }
}

/// Build the combined NM×NM CTMC generator G = Λ⊗I_N + blockdiag(G_1,...,G_M).
pub fn assemble_generator(slv_model: &SLVModelVariant, x_grid: &[f64], x_spacings: &[f64], v_grid: &[f64]) -> CombinedCTMC {
    let n = x_grid.len();
    let m = v_grid.len();
    let theta_mat = theta(slv_model, x_grid, v_grid);
    let sigma_tilde = sigma_tilde_sq(slv_model, v_grid);

    let g_l_matrices: Vec<DMatrix<f64>> = (0..m)
        .into_par_iter()
        .map(|l| build_g_l_matrix_inner(x_spacings, &theta_mat[l], sigma_tilde[l], n))
        .collect();

    let var_coo = build_variance_generator(slv_model, v_grid);
    let var_dense = VarianceGeneratorDense::from_coo(&var_coo);
    let generator = build_combined_generator(&var_dense, &g_l_matrices, n, m);

    CombinedCTMC { generator, x_grid: x_grid.to_vec(), x_spacings: x_spacings.to_vec(), v_grid: v_grid.to_vec(), n, m }
}

// Algorithm 3.1: CTMC Integral Equation method for American options (Ma, Yang & Cui 2021)

pub struct CTMCPricerConfig {
    pub option: options::Options,
    pub n_time_steps: usize,
}

impl CTMCPricerConfig {
    #[inline] pub fn dt(&self) -> f64 { self.option.time_to_maturity() / self.n_time_steps as f64 }

    // Terminal exercise boundary in spot-space at maturity.
    // Put: min(K, rK/q) — exercised early only when spot is low.
    // Call: max(K, rK/q) — exercised early only when dividend yield exceeds risk-free rate; q=0 → +∞.
    pub fn terminal_boundary_s(&self) -> f64 {
        let k = self.option.strike_price();
        let r = self.option.risk_free_rate();
        let q = self.option.dividend_yield().unwrap_or(0.0);
        match &self.option {
            options::Options::Put(_) => {
                if q > 1e-15 { k.min(k * r / q) } else { k }
            }
            options::Options::Call(_) => {
                if q > 1e-15 { k.max(k * r / q) } else { f64::INFINITY }
            }
        }
    }
}

pub struct CTMCResult {
    pub price: f64,
    /// Full J_0 vector: prices at every (x_i, v_l) grid point
    pub j0: Vec<f64>,
    /// Early exercise boundary[k][l] in X-space, k=0..=n_time_steps, l=0..M
    pub boundary_x: Vec<Vec<f64>>,
    pub boundary_s: Vec<Vec<f64>>,
}

// H^(1): European option payoff vector.
// excess[l*n+i] = intrinsic value (K−S for put, S−K for call), pre-computed, unclamped.
#[inline]
fn build_h1(excess: &[f64]) -> Vec<f64> {
    excess.iter().map(|&e| e.max(0.0)).collect()
}

// H^(3)_k: early exercise integrand.
//   Put  (is_call=false): exercise region x_i ≤ b_k(v_l); coeff = rK − (rS − ω).
//   Call (is_call=true) : exercise region x_i ≥ b_k(v_l); coeff = (rS − ω) − rK.
// is_call is pre-computed outside the time-step loop so the branch is hoisted.
#[inline]
fn build_h3(h3_coeff: &[f64], x_grid: &[f64], boundary: &[f64], m: usize, n: usize, is_call: bool) -> Vec<f64> {
    (0..m)
        .into_par_iter()
        .flat_map_iter(|l| {
            let b = boundary[l];
            (0..n).map(move |i| {
                let in_region = if is_call { x_grid[i] >= b } else { x_grid[i] <= b };
                if in_region { h3_coeff[l * n + i] } else { 0.0 }
            })
        })
        .collect()
}

// Find b_k(v_l): X-space boundary where intrinsic = continuation for variance state l.
// f(i) = excess[base+i] − J_k[base+i] using SIGNED excess = K-S (not clamped), monotone
// decreasing in i: f ≈ 0 across the exercise plateau (deep ITM, V_k = K-S), then f < 0
// past the boundary (continuation, V_k > K-S), ultimately very negative at deep OTM.
// Boundary = rightmost i with f(i) ≥ −tol. Negative tolerance absorbs numerical noise
// from the explicit H^(3) time-step that can push J_k slightly above intrinsic in the
// exercise region (previous f ≥ 0 threshold spuriously collapsed the exercise set).
#[inline]
fn find_boundary_for_state(j_k: &[f64], excess: &[f64], x_grid: &[f64], l: usize, n: usize, strike: f64) -> f64 {
    let base = l * n;
    let f = |i: usize| excess[base + i] - j_k[base + i];
    let tol = 1e-6 * strike;

    // Entire grid is exercise region (rare: deep ITM throughout).
    if f(n - 1) >= -tol { return x_grid[n - 1]; }
    // No exercise region (even deep ITM has f < -tol — exercise boundary outside grid).
    if f(0) < -tol { return x_grid[0]; }
    // Invariant: f(lo) ≥ -tol (exercise plateau), f(hi) < -tol (continuation).
    let mut lo = 0;
    let mut hi = n - 1;
    while lo + 1 < hi {
        let mid = (lo + hi) / 2;
        if f(mid) >= -tol { lo = mid; } else { hi = mid; }
    }
    let f_lo = f(lo);
    let f_hi = f(lo + 1);
    if (f_lo - f_hi) > 1e-30 {
        let alpha = ((f_lo + tol) / (f_lo - f_hi)).clamp(0.0, 1.0);
        return x_grid[lo] + alpha * (x_grid[lo + 1] - x_grid[lo]);
    }
    x_grid[lo]
}

// Find b_k(v_l) for a call: X-space boundary where exercise region is the RIGHT side.
// f(i) = excess[base+i] − J_k[base+i] with signed excess = S−K; monotone increasing in i:
// f very negative at deep OTM, crosses 0 at the boundary, plateau f ≈ 0 across exercise.
// Boundary = leftmost i with f(i) ≥ −tol. Returns +∞ when no exercise region exists
// (q=0 case: h3 then zero everywhere, recovering European Heston).
#[inline]
fn find_boundary_for_state_call(j_k: &[f64], excess: &[f64], x_grid: &[f64], l: usize, n: usize, strike: f64) -> f64 {
    let base = l * n;
    let f = |i: usize| excess[base + i] - j_k[base + i];
    let tol = 1e-6 * strike;

    if f(n - 1) < -tol { return f64::INFINITY; } // no exercise anywhere in grid
    if f(0) >= -tol { return x_grid[0]; }         // entire grid is exercise region
    // Invariant: f(lo) < -tol (continuation), f(hi) ≥ -tol (exercise plateau).
    let mut lo = 0;
    let mut hi = n - 1;
    while lo + 1 < hi {
        let mid = (lo + hi) / 2;
        if f(mid) < -tol { lo = mid; } else { hi = mid; }
    }
    let f_lo = f(lo);
    let f_hi = f(lo + 1);
    if (f_hi - f_lo) > 1e-30 {
        let alpha = ((-tol - f_lo) / (f_hi - f_lo)).clamp(0.0, 1.0);
        return x_grid[lo] + alpha * (x_grid[lo + 1] - x_grid[lo]);
    }
    x_grid[lo + 1]
}

#[inline]
fn interpolate_value(j: &[f64], x_grid: &[f64], v_grid: &[f64], x: f64, v: f64, n: usize, m: usize) -> f64 {
    let ix = bracket_search(x_grid, x);
    let ix1 = (ix + 1).min(n - 1);
    let iv = bracket_search(v_grid, v);
    let iv1 = (iv + 1).min(m - 1);

    let wx = if ix < ix1 { ((x - x_grid[ix]) / (x_grid[ix1] - x_grid[ix])).clamp(0.0, 1.0) } else { 0.0 };
    let wv = if iv < iv1 { ((v - v_grid[iv]) / (v_grid[iv1] - v_grid[iv])).clamp(0.0, 1.0) } else { 0.0 };

    (1.0 - wx) * (1.0 - wv) * j[iv * n + ix]
        + wx * (1.0 - wv) * j[iv * n + ix1]
        + (1.0 - wx) * wv * j[iv1 * n + ix]
        + wx * wv * j[iv1 * n + ix1]
}

#[inline]
fn bracket_search(grid: &[f64], target: f64) -> usize {
    let n = grid.len();
    if n == 0 || target <= grid[0] { return 0; }
    if target >= grid[n - 1] { return n.saturating_sub(2); }
    let (mut lo, mut hi) = (0, n - 1);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if grid[mid] <= target { lo = mid; } else { hi = mid; }
    }
    lo
}

/// Run Algorithm 3.1 (Ma, Yang & Cui 2021) backward in time.
pub fn run_algorithm_31(
    ctmc: &CombinedCTMC,
    config: &CTMCPricerConfig,
    model: &SLVModelVariant,
    operator: &MatExpOperator,
    spot_0: f64,
    v0: f64,
) -> CTMCResult {
    let n = ctmc.n;
    let m = ctmc.m;
    let nm = n * m;
    let n_steps = config.n_time_steps;
    let dt = config.dt();
    let rho = model.rho();

    let f_vals: Vec<f64> = ctmc.v_grid.iter().map(|&v| model.f_integral(v)).collect();

    // Pre-compute per-grid tables — constant across all time steps.
    let s_grid: Vec<f64> = (0..m)
        .flat_map(|l| {
            let f_vl = f_vals[l];
            (0..n).map(move |i| model.g_inverse(ctmc.x_grid[i] + rho * f_vl))
        })
        .collect();
    let r = config.option.risk_free_rate();
    let k_strike = config.option.strike_price();
    // Hoist call/put dispatch out of all inner loops.
    let is_call  = matches!(config.option, options::Options::Call(_));
    let eep_sign = if is_call { -1.0_f64 } else { 1.0_f64 };

    // excess[l*n+i]: intrinsic value (unsigned), used for exercise boundary search.
    //   Put: K − S(x_i, v_l)   Call: S(x_i, v_l) − K
    let excess: Vec<f64> = if is_call {
        s_grid.iter().map(|&s| s - k_strike).collect()
    } else {
        s_grid.iter().map(|&s| k_strike - s).collect()
    };
    // h3_coeff[l*n+i]: early-exercise integrand, pre-computed.
    //   Put:  rK − (rS − ω(S, v_l))   Call: (rS − ω(S, v_l)) − rK  (unified via eep_sign)
    let s_sl: &[f64] = &s_grid;
    let h3_coeff: Vec<f64> = (0..m)
        .flat_map(|l| {
            let v_l = ctmc.v_grid[l];
            (0..n).map(move |i| {
                let s = s_sl[l * n + i];
                eep_sign * (r * k_strike - (r * s - model.omega(s, v_l)))
            })
        })
        .collect();

    let mut j = build_h1(&excess);
    let mut boundary_x: Vec<Vec<f64>> = vec![vec![0.0; m]; n_steps + 1];

    let b_term_g = model.g_integral(config.terminal_boundary_s());
    for l in 0..m {
        boundary_x[n_steps][l] = b_term_g - rho * f_vals[l];
    }

    let mut rhs = vec![0.0_f64; nm];
    for k in (0..n_steps).rev() {
        let h3 = build_h3(&h3_coeff, &ctmc.x_grid, &boundary_x[k + 1], m, n, is_call);
        rhs.iter_mut()
            .zip(j.iter())
            .zip(h3.iter())
            .for_each(|((r, &jv), &h)| *r = jv + dt * h);
        j = operator.apply(&rhs);
        // Enforce American constraint J_k ≥ h1 = max(intrinsic, 0) pointwise.
        // The paper's Algorithm 3.1 requires a self-consistent boundary at each step;
        // the explicit lag (h3 built from boundary_{k+1}) lets J_k drift below intrinsic
        // under numerical noise, which then classifies the entire plateau as "no exercise"
        // (h3 → 0) and cascades into a European-style propagation. Projection restores the
        // invariant and stabilizes the boundary search on the next iteration.
        j.iter_mut().zip(excess.iter()).for_each(|(v, &e)| {
            let payoff = e.max(0.0);
            if *v < payoff { *v = payoff; }
        });
        boundary_x[k] = (0..m)
            .into_par_iter()
            .map(|l| if is_call {
                find_boundary_for_state_call(&j, &excess, &ctmc.x_grid, l, n, k_strike)
            } else {
                find_boundary_for_state(&j, &excess, &ctmc.x_grid, l, n, k_strike)
            })
            .collect();
    }

    let boundary_s = boundary_x.iter()
        .map(|bx| bx.iter().enumerate().map(|(l, &b)| model.g_inverse(b + rho * f_vals[l])).collect())
        .collect();

    let x0 = model.g_integral(spot_0) - rho * model.f_integral(v0);
    let price = interpolate_value(&j, &ctmc.x_grid, &ctmc.v_grid, x0, v0, n, m);

    CTMCResult { price, j0: j, boundary_x, boundary_s }
}

/// Price an American option under the Heston stochastic volatility model.
///
/// Accepts any `Options` variant (put or call). Builds grids, assembles G, selects Dense Padé
/// or Krylov based on problem size, then runs Algorithm 3.1.
///
/// Grids are moment-matched: the variance grid covers the CIR conditional
/// moments of v_T (so it reaches θ even when θ ≫ v0), and the x-grid is
/// centred on x₀ with a half-width from the moments of X_T (see
/// `heston_variance_grid` / `heston_x_half_width`).
pub fn price_american_option_heston(
    option: options::Options,
    v0: f64,
    kappa: f64,
    theta_lr: f64,
    sigma_v: f64,
    heston_rho: f64,
    n_x: usize,
    m_v: usize,
    n_time: usize,
) -> CTMCResult {
    let spot = option.spot_price();
    let r    = option.risk_free_rate();
    let q    = option.dividend_yield().unwrap_or(0.0);
    let t    = option.time_to_maturity();
    let model = SLVModelVariant::Heston(Heston::new(heston_rho, sigma_v, r, Some(q), kappa, theta_lr));

    let v_grid = heston_variance_grid(v0, kappa, theta_lr, sigma_v, t, m_v);

    let f_v0 = model.f_integral(v0);
    let x0 = model.g_integral(spot) - heston_rho * f_v0;
    let half = heston_x_half_width(
        v0, kappa, theta_lr, sigma_v, heston_rho, r, q, t,
        v_grid[0], *v_grid.last().unwrap(),
    );
    let (x_min, x_max) = (x0 - half, x0 + half);
    let x_grid: Vec<f64> = (0..n_x).map(|i| x_min + (x_max - x_min) * i as f64 / (n_x - 1) as f64).collect();
    let x_spacings: Vec<f64> = x_grid.windows(2).map(|w| w[1] - w[0]).collect();

    let ctmc = assemble_generator(&model, &x_grid, &x_spacings, &v_grid);
    let nm = ctmc.dim();
    let dt = t / n_time as f64;
    let operator = MatExpOperator::new(&ctmc.generator, r, dt, MatExpOperator::suggest_strategy(nm));
    let config = CTMCPricerConfig { option, n_time_steps: n_time };

    run_algorithm_31(&ctmc, &config, &model, &operator, spot, v0)
}

/// Price an American put under Heston — convenience wrapper around `price_american_option_heston`.
pub fn price_american_put_heston(
    strike: f64, spot: f64, v0: f64, maturity: f64,
    r: f64, q: f64, kappa: f64, theta_lr: f64, sigma_v: f64, rho: f64,
    n_x: usize, m_v: usize, n_time: usize,
) -> CTMCResult {
    price_american_option_heston(
        options::Options::new_put(strike, spot, 0.0, r, maturity, Some(q)),
        v0, kappa, theta_lr, sigma_v, rho, n_x, m_v, n_time,
    )
}

/// Price an American call under Heston — convenience wrapper around `price_american_option_heston`.
pub fn price_american_call_heston(
    strike: f64, spot: f64, v0: f64, maturity: f64,
    r: f64, q: f64, kappa: f64, theta_lr: f64, sigma_v: f64, rho: f64,
    n_x: usize, m_v: usize, n_time: usize,
) -> CTMCResult {
    price_american_option_heston(
        options::Options::new_call(strike, spot, 0.0, r, maturity, Some(q)),
        v0, kappa, theta_lr, sigma_v, rho, n_x, m_v, n_time,
    )
}

/// Price an American option spread under Heston.
///
/// Each leg is priced independently via `price_american_option_heston` and combined with its weight.
pub fn price_american_spread_heston(
    spread: &options::optionspreads::OptionSpreads,
    v0: f64, kappa: f64, theta_lr: f64, sigma_v: f64, heston_rho: f64,
    n_x: usize, m_v: usize, n_time: usize,
) -> f64 {
    spread.components().iter().map(|pos| {
        let result = price_american_option_heston(
            pos.option, v0, kappa, theta_lr, sigma_v, heston_rho, n_x, m_v, n_time,
        );
        result.price * pos.weight as f64
    }).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn h1_payoff_non_negative() {
        let x_grid = vec![-1.0, 0.0, 1.0, 2.0, 3.0];
        let x_spacings = vec![1.0; 4];
        let v_grid = vec![0.02, 0.04, 0.06];
        let model = SLVModelVariant::Heston(Heston::new(-0.7, 0.3, 0.05, Some(0.02), 2.0, 0.04));
        let ctmc = assemble_generator(&model, &x_grid, &x_spacings, &v_grid);
        let strike = 100.0_f64;
        let rho = model.rho();
        let f_vals: Vec<f64> = ctmc.v_grid.iter().map(|&v| model.f_integral(v)).collect();
        let x_sl: &[f64] = &ctmc.x_grid;
        let model_ref = &model;
        let excess: Vec<f64> = (0..ctmc.m).flat_map(|l| {
            let f_vl = f_vals[l];
            (0..ctmc.n).map(move |i| strike - model_ref.g_inverse(x_sl[i] + rho * f_vl))
        }).collect();
        let h1 = build_h1(&excess);
        assert_eq!(h1.len(), 15);
        for &v in &h1 { assert!(v >= 0.0, "Negative payoff: {v}"); }
        assert!(h1[0] > h1[4]);
    }

    #[test]
    fn variance_grid_reaches_theta() {
        // θ = 8·v0: the old (0, 4·v0] grid capped the chain at 0.16·v0/…, far
        // below the mean-reversion level. The moment-matched grid must cover θ.
        let g = heston_variance_grid(0.02, 2.0, 0.16, 0.3, 1.0, 25);
        assert_eq!(g.len(), 25);
        assert!(g[0] > 0.0, "grid must stay strictly positive");
        assert!(g.windows(2).all(|w| w[1] > w[0]), "grid must be increasing");
        assert!(g[0] < 0.02, "v0 must be interior (lower bound below v0)");
        assert!(*g.last().unwrap() > 0.16, "grid must reach beyond theta");

        // Downward term structure: grid still brackets both v0 and θ.
        let g = heston_variance_grid(0.08, 2.0, 0.04, 0.3, 1.0, 25);
        assert!(g[0] < 0.04 && *g.last().unwrap() > 0.08);

        // Degenerate vol-of-vol: v0 stays strictly interior.
        let g = heston_variance_grid(0.04, 2.0, 0.04, 1e-12, 1.0, 25);
        assert!(g[0] < 0.04 && *g.last().unwrap() >= 0.04);
    }

    #[test]
    #[ignore]
    fn american_put_theta_above_v0_exceeds_european() {
        // Regression for the variance-grid truncation bug: with θ = 4·v0 the
        // old (0, 4·v0] grid priced this put ~40% below its own European lower
        // bound. American CTMC must dominate the European COS price.
        use crate::slv::{cos_european, HestonParameters, OptionSide};
        let p = HestonParameters { v0: 0.04, v_bar: 0.16, kappa: 2.0, sigma: 0.3, rho: -0.7 };
        let (s, k, r, q, tau) = (100.0_f64, 100.0_f64, 0.05, 0.02, 1.0);
        let eur = cos_european(OptionSide::Put, s, k, r, q, tau, &p, 512);
        let amer = price_american_put_heston(
            k, s, p.v0, tau, r, q, p.kappa, p.v_bar, p.sigma, p.rho, 120, 25, 75,
        ).price;
        println!("European COS: {eur:.4}, American CTMC: {amer:.4}");
        assert!(amer >= eur * 0.98,
            "American {amer:.4} must not fall below European lower bound {eur:.4}");
        assert!(amer <= eur * 1.20,
            "American {amer:.4} implausibly far above European {eur:.4}");
    }

    #[test]
    fn bracket_search_correctness() {
        let grid = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        assert_eq!(bracket_search(&grid, -1.0), 0);
        assert_eq!(bracket_search(&grid, 0.0), 0);
        assert_eq!(bracket_search(&grid, 0.5), 0);
        assert_eq!(bracket_search(&grid, 1.0), 1);
        assert_eq!(bracket_search(&grid, 2.7), 2);
        assert_eq!(bracket_search(&grid, 5.0), 3);
    }
    #[test]
    #[ignore]
    fn american_put_heston_smoke() {
        // Realistic Heston params: kappa=2, theta=0.04 (LR var), sigma_v=0.3 (vol-of-vol), rho=-0.7
        // ATM American put S=K=100, T=1yr, r=5%, q=2%; expected ~$7-11 from literature
        let result = price_american_put_heston(
            100.0, 100.0, 0.04, 1.0, 0.05, 0.02,
            2.0, 0.04, 0.3, -0.7,
            100, 25, 75,
        );
        println!("American put price: {:.6}", result.price);
        println!("Boundary at t=0, mid v-state: {:.4}", result.boundary_s[0][12]);
        assert!(result.price > 5.0, "Price below European BS: {}", result.price);
        assert!(result.price < 20.0, "Price unreasonably high: {}", result.price);
        assert!(result.boundary_s[0][12] < 100.0, "Boundary should be below strike");
    }

    #[test]
    fn call_terminal_boundary_above_put() {
        // q=3%, r=5%, K=100 → call boundary = max(100, 5/3*100) ≈ 166.7; put = min(100, 166.7) = 100
        let put_config = CTMCPricerConfig {
            option: options::Options::new_put(100.0, 100.0, 0.0, 0.05, 1.0, Some(0.03)),
            n_time_steps: 10,
        };
        let call_config = CTMCPricerConfig {
            option: options::Options::new_call(100.0, 100.0, 0.0, 0.05, 1.0, Some(0.03)),
            n_time_steps: 10,
        };
        assert!(call_config.terminal_boundary_s() > put_config.terminal_boundary_s());
        assert!((put_config.terminal_boundary_s() - 100.0).abs() < 1e-10);
        assert!((call_config.terminal_boundary_s() - 500.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn call_no_dividend_boundary_is_infinite() {
        let config = CTMCPricerConfig {
            option: options::Options::new_call(100.0, 100.0, 0.0, 0.05, 1.0, None),
            n_time_steps: 10,
        };
        assert_eq!(config.terminal_boundary_s(), f64::INFINITY);
    }

    #[test]
    #[ignore]
    fn american_call_heston_smoke() {
        // ATM call S=K=100, T=1yr, r=5%, q=3%; Heston params mirror put smoke test
        // With q>0 early exercise is possible; expected price roughly $8-14
        let result = price_american_call_heston(
            100.0, 100.0, 0.04, 1.0, 0.05, 0.03,
            2.0, 0.04, 0.3, -0.7,
            100, 25, 75,
        );
        println!("American call price: {:.6}", result.price);
        println!("Boundary at t=0, mid v-state: {:.4}", result.boundary_s[0][12]);
        assert!(result.price > 5.0, "Price too low: {}", result.price);
        assert!(result.price < 25.0, "Price unreasonably high: {}", result.price);
        assert!(result.boundary_s[0][12] > 100.0, "Call boundary should be above strike");
    }

    #[test]
    #[ignore]
    fn american_call_zero_dividend_no_eep() {
        // q=0: American call = European call. Terminal boundary is +∞; no early exercise premium.
        // Validate: (1) terminal boundary row is +∞, (2) price is in European call range.
        let result = price_american_call_heston(
            100.0, 100.0, 0.04, 1.0, 0.05, 0.0,
            2.0, 0.04, 0.3, -0.7,
            100, 25, 75,
        );
        println!("American call (q=0) price: {:.6}", result.price);
        // Terminal boundary (last row of boundary_s) should be +∞ for all variance states.
        let n_steps = result.boundary_s.len() - 1;
        for &b in &result.boundary_s[n_steps] {
            assert!(b.is_infinite(), "Terminal call boundary with q=0 should be infinite: {b}");
        }
        // Price should be consistent with European call (BS benchmark: ~$10-12 with vol≈0.2)
        assert!(result.price > 8.0 && result.price < 18.0,
            "q=0 call price should be near European: {}", result.price);
    }

    #[test]
    #[ignore]
    fn spread_straddle_price_reasonable() {
        use options::optionspreads::{Direction, OptionSpreads};
        // ATM straddle: long call + long put at same strike. Price > each individual leg.
        let put_price = price_american_put_heston(
            100.0, 100.0, 0.04, 1.0, 0.05, 0.02, 2.0, 0.04, 0.3, -0.7, 60, 15, 50,
        ).price;
        let call_price = price_american_call_heston(
            100.0, 100.0, 0.04, 1.0, 0.05, 0.02, 2.0, 0.04, 0.3, -0.7, 60, 15, 50,
        ).price;
        let straddle = OptionSpreads::new_straddle(
            Direction::LONG, 100.0, 100.0, 0.2, 0.05, 1.0, Some(0.02),
        );
        let spread_price = price_american_spread_heston(
            &straddle, 0.04, 2.0, 0.04, 0.3, -0.7, 60, 15, 50,
        );
        println!("Put: {put_price:.4}, Call: {call_price:.4}, Straddle spread: {spread_price:.4}");
        assert!(spread_price > put_price, "Straddle should exceed put alone");
        assert!(spread_price > call_price, "Straddle should exceed call alone");
        assert!(spread_price < put_price + call_price + 0.5, "Straddle close to sum of legs");
    }

    // Run with: cargo test -p numerical-methods --release -- --ignored ctmc_speed --nocapture
    // Sweeps (n_x, m_v, n_time) combinations and prints wall-clock time.
    // Target: find the largest configuration that runs in ~0.1-0.2 seconds.
    #[test]
    #[ignore]
    fn ctmc_speed_sweep() {
        use std::time::Instant;

        // (n_x, m_v, n_time): NM = n_x*m_v is the matrix dimension driving cost
        let configs: &[(usize, usize, usize)] = &[
            (20,  5,  50),   // NM=100   — very fast baseline
            (40,  5,  50),   // NM=200
            (40, 10,  50),   // NM=400   — Dense→Krylov threshold boundary
            (60, 10,  50),   // NM=600   — Krylov
            (80, 10,  50),   // NM=800
            (80, 15,  50),   // NM=1200
            (80, 20,  50),   // NM=1600
            (100, 15, 50),   // NM=1500
            (100, 20, 50),   // NM=2000
            (100, 25, 50),   // NM=2500  — "default" paper grid
            (120, 25, 50),   // NM=3000
            (150, 25, 50),   // NM=3750
            (100, 20, 100),  // NM=2000, 2× time steps
            (100, 25, 100),  // NM=2500, 2× time steps
        ];

        println!("\n{:<10} {:<6} {:<6} {:<8} {:<8} {:<12} {:<8}",
            "NM", "n_x", "m_v", "n_time", "strategy", "elapsed(ms)", "price");
        println!("{}", "-".repeat(62));

        for &(n_x, m_v, n_time) in configs {
            let nm = n_x * m_v;
            let strategy = match MatExpOperator::suggest_strategy(nm) {
                crate::expm::MatExpStrategy::DensePade => "Dense",
                crate::expm::MatExpStrategy::Krylov { .. } => "Krylov",
            };
            let t0 = Instant::now();
            let result = price_american_put_heston(
                100.0, 100.0, 0.04, 1.0, 0.05, 0.02,
                2.0, 0.04, 0.3, -0.7,
                n_x, m_v, n_time,
            );
            let ms = t0.elapsed().as_secs_f64() * 1e3;
            println!("{:<10} {:<6} {:<6} {:<8} {:<8} {:<12.1} {:<8.4}",
                nm, n_x, m_v, n_time, strategy, ms, result.price);
        }
    }
}
