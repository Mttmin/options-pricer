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
        // Placeholders: calibrate eta, rho, sigma, ref_vol to market data in practice
        Self {
            rho: 0.0,
            sigma: option.volatility(),
            r_f: option.risk_free_rate(),
            q: option.dividend_yield(),
            eta: 1.0,
            ref_vol: option.volatility(),
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
    pub fn g_integral(&self, s: f64) -> f64 { dispatch!(self, g_integral, s) }
    pub fn f_integral(&self, v: f64) -> f64 { dispatch!(self, f_integral, v) }
    pub fn g_inverse(&self, x: f64) -> f64  { dispatch!(self, g_inverse, x) }
    pub fn f_inverse(&self, y: f64) -> f64  { dispatch!(self, f_inverse, y) }
    pub fn m_prime(&self, v: f64) -> f64    { dispatch!(self, m_prime, v) }
    pub fn sigma_prime_v(&self, v: f64) -> f64 { dispatch!(self, sigma_prime_v, v) }
    pub fn gamma_prime(&self, s: f64) -> f64   { dispatch!(self, gamma_prime, s) }
    pub fn rho(&self) -> f64                { dispatch!(self, rho) }
    pub fn omega(&self, s: f64, v: f64) -> f64 { dispatch!(self, omega, s, v) }
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
// Keep size small to balance numerical stability vs. matrix-exponentiation cost.
pub fn uniform_variance_grid(vol_0: f64, n_vol_steps: i32) -> Vec<f64> {
    let step = 4.0 * vol_0 / (n_vol_steps as f64);
    (1..=n_vol_steps).map(|i| i as f64 * step).collect()
}

#[allow(dead_code)]
/// Sinh-stretched variance grid per Cui, Zhenyu 2018.
fn sinh_variance_grid(_vol: f64, _n_vol_steps: i32) -> Vec<f64> {
    todo!()
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
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

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
    vol_grid.iter().map(|&v| {
        let m = slv_model.m(v);
        let rho = slv_model.rho();
        m * m * (1.0 - rho * rho)
    }).collect()
}

fn psi_v(slv_model: &SLVModelVariant, v: f64) -> f64 {
    let m = slv_model.m(v);
    let sigma = slv_model.sigma_v(v);
    slv_model.mu(v) * m / sigma
        + 0.5 * (sigma * slv_model.m_prime(v) - slv_model.sigma_prime_v(v) * m)
}

fn theta(slv_model: &SLVModelVariant, x_grid: &[f64], vol_grid: &[f64]) -> DMatrix<f64> {
    let mut t = DMatrix::zeros(vol_grid.len(), x_grid.len());
    for (i, &vol) in vol_grid.iter().enumerate() {
        let psi = psi_v(slv_model, vol);
        for j in 0..x_grid.len() {
            let s = slv_model.g_inverse(x_grid[j] + slv_model.rho() * slv_model.f_integral(vol));
            t[(i, j)] = slv_model.omega(s, vol) / slv_model.gamma(s)
                - slv_model.gamma_prime(s) * slv_model.m(vol).powi(2) * 0.5
                - slv_model.rho() * psi;
        }
    }
    t
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
        .map(|l| {
            let theta_col: Vec<f64> = (0..n).map(|i| theta_mat[(l, i)]).collect();
            build_g_l_matrix_inner(x_spacings, &theta_col, sigma_tilde[l], n)
        })
        .collect();

    let var_coo = build_variance_generator(slv_model, v_grid);
    let var_dense = VarianceGeneratorDense::from_coo(&var_coo);
    let generator = build_combined_generator(&var_dense, &g_l_matrices, n, m);

    CombinedCTMC { generator, x_grid: x_grid.to_vec(), x_spacings: x_spacings.to_vec(), v_grid: v_grid.to_vec(), n, m }
}

// Algorithm 3.1: CTMC Integral Equation method for American options (Ma, Yang & Cui 2021)

pub struct CTMCPricerConfig {
    pub strike: f64,
    pub risk_free: f64,
    pub dividend: f64,
    pub maturity: f64,
    pub n_time_steps: usize,
}

impl CTMCPricerConfig {
    pub fn dt(&self) -> f64 { self.maturity / self.n_time_steps as f64 }

    // Terminal boundary in spot-space: min(K, rK/q). Limit is K when q=0.
    pub fn terminal_boundary_s(&self) -> f64 {
        if self.dividend > 1e-15 {
            self.strike.min(self.strike * self.risk_free / self.dividend)
        } else {
            self.strike
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

// H^(1): European put payoff vector.
// Element at φ(i,l) = l·N+i: max(K − g^{−1}(x_i + ρ·f(v_l)), 0)
fn build_h1(ctmc: &CombinedCTMC, config: &CTMCPricerConfig, model: &SLVModelVariant) -> Vec<f64> {
    let n = ctmc.n;
    let m = ctmc.m;
    let rho = model.rho();
    (0..m)
        .into_par_iter()
        .flat_map_iter(|l| {
            let f_vl = model.f_integral(ctmc.v_grid[l]);
            (0..n).map(move |i| {
                (config.strike - model.g_inverse(ctmc.x_grid[i] + rho * f_vl)).max(0.0)
            })
        })
        .collect()
}

// H^(3)_k: early exercise integrand. For nodes in the exercise region (x_i ≤ b_k(v_l)):
//   rK − φ(x_i, v_l)  where  φ = r·S − ω(S, v_l)
// For Heston: ω = (r−q)·S, so φ = q·S and the integrand is rK − q·S.
fn build_h3(ctmc: &CombinedCTMC, config: &CTMCPricerConfig, model: &SLVModelVariant, boundary: &[f64], f_vals: &[f64]) -> Vec<f64> {
    let n = ctmc.n;
    let m = ctmc.m;
    let r = config.risk_free;
    let k = config.strike;
    let rho = model.rho();
    (0..m)
        .into_par_iter()
        .flat_map_iter(|l| {
            let b = boundary[l];
            let f_vl = f_vals[l];
            let v_l = ctmc.v_grid[l];
            (0..n).map(move |i| {
                if ctmc.x_grid[i] <= b {
                    let s = model.g_inverse(ctmc.x_grid[i] + rho * f_vl);
                    r * k - (r * s - model.omega(s, v_l))
                } else {
                    0.0
                }
            })
        })
        .collect()
}

// Find b_k(v_l): X-space boundary where intrinsic = continuation for variance state l.
// Scans right-to-left for the first crossing of K − S_i − J_k[l·N+i] ≥ 0,
// then refines with linear interpolation.
fn find_boundary_for_state(j_k: &[f64], x_grid: &[f64], l: usize, n: usize, strike: f64, rho: f64, f_vl: f64, model: &SLVModelVariant) -> f64 {
    let base = l * n;
    let mut cross = 0;
    for i in (0..n).rev() {
        if (strike - model.g_inverse(x_grid[i] + rho * f_vl)) - j_k[base + i] >= 0.0 {
            cross = i;
            break;
        }
    }
    if cross < n - 1 {
        let f_lo = (strike - model.g_inverse(x_grid[cross] + rho * f_vl)) - j_k[base + cross];
        let f_hi = (strike - model.g_inverse(x_grid[cross + 1] + rho * f_vl)) - j_k[base + cross + 1];
        if (f_lo - f_hi).abs() > 1e-30 {
            let alpha = f_lo / (f_lo - f_hi);
            return x_grid[cross] + alpha * (x_grid[cross + 1] - x_grid[cross]);
        }
    }
    x_grid[cross]
}

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

    let mut j = build_h1(ctmc, config, model);
    let mut boundary_x: Vec<Vec<f64>> = vec![vec![0.0; m]; n_steps + 1];

    let b_term_g = model.g_integral(config.terminal_boundary_s());
    for l in 0..m {
        boundary_x[n_steps][l] = b_term_g - rho * f_vals[l];
    }

    for k in (0..n_steps).rev() {
        let h3 = build_h3(ctmc, config, model, &boundary_x[k + 1], &f_vals);
        let mut rhs = vec![0.0; nm];
        for idx in 0..nm {
            rhs[idx] = j[idx] + dt * h3[idx];
        }
        j = operator.apply(&rhs);
        boundary_x[k] = (0..m)
            .into_par_iter()
            .map(|l| find_boundary_for_state(&j, &ctmc.x_grid, l, n, config.strike, rho, f_vals[l], model))
            .collect();
    }

    let boundary_s = boundary_x.iter()
        .map(|bx| bx.iter().enumerate().map(|(l, &b)| model.g_inverse(b + rho * f_vals[l])).collect())
        .collect();

    let x0 = model.g_integral(spot_0) - rho * model.f_integral(v0);
    let price = interpolate_value(&j, &ctmc.x_grid, &ctmc.v_grid, x0, v0, n, m);

    CTMCResult { price, j0: j, boundary_x, boundary_s }
}

/// Price an American put under the Heston stochastic volatility model.
///
/// Builds grids, assembles G, selects Dense Padé or Krylov based on problem size,
/// then runs Algorithm 3.1. Grid range follows the paper: x ∈ [g(0.001S), g(4S)].
pub fn price_american_put_heston(
    strike: f64,
    spot: f64,
    v0: f64,
    maturity: f64,
    r: f64,
    q: f64,
    kappa: f64,
    theta_lr: f64,
    sigma_v: f64,
    rho: f64,
    n_x: usize,
    m_v: usize,
    n_time: usize,
) -> CTMCResult {
    let model = SLVModelVariant::Heston(Heston::new(rho, sigma_v, r, Some(q), kappa, theta_lr));

    let v_grid = uniform_variance_grid(v0, m_v as i32);

    let f_v0 = model.f_integral(v0);
    let x_min = model.g_integral(1e-3 * spot) - rho * f_v0;
    let x_max = model.g_integral(4.0 * spot) - rho * f_v0;
    let x_grid: Vec<f64> = (0..n_x).map(|i| x_min + (x_max - x_min) * i as f64 / (n_x - 1) as f64).collect();
    let x_spacings: Vec<f64> = x_grid.windows(2).map(|w| w[1] - w[0]).collect();

    let ctmc = assemble_generator(&model, &x_grid, &x_spacings, &v_grid);
    let nm = ctmc.dim();
    let dt = maturity / n_time as f64;
    let operator = MatExpOperator::new(&ctmc.generator, r, dt, MatExpOperator::suggest_strategy(nm));
    let config = CTMCPricerConfig { strike, risk_free: r, dividend: q, maturity, n_time_steps: n_time };

    run_algorithm_31(&ctmc, &config, &model, &operator, spot, v0)
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
        let config = CTMCPricerConfig { strike: 100.0, risk_free: 0.05, dividend: 0.02, maturity: 1.0, n_time_steps: 10 };
        let h1 = build_h1(&ctmc, &config, &model);
        assert_eq!(h1.len(), 15);
        for &v in &h1 { assert!(v >= 0.0, "Negative payoff: {v}"); }
        assert!(h1[0] > h1[4]);
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
            100, 25, 50,
        );
        println!("American put price: {:.6}", result.price);
        println!("Boundary at t=0, mid v-state: {:.4}", result.boundary_s[0][12]);
        assert!(result.price > 5.0, "Price below European BS: {}", result.price);
        assert!(result.price < 20.0, "Price unreasonably high: {}", result.price);
        assert!(result.boundary_s[0][12] < 100.0, "Boundary should be below strike");
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
            (40, 10,  50),   // NM=400
            (60, 10,  50),   // NM=600
            (80, 10,  50),   // NM=800
            (80, 15,  50),   // NM=1200
            (80, 20,  50),   // NM=1600  — crosses Dense→Krylov threshold at 1500
            (100, 15, 50),   // NM=1500  — threshold boundary
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
            let strategy = if nm <= 600 { "Dense" } else { "Krylov" };
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
