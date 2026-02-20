use nalgebra::DMatrix;
use options::Options;

trait SLVModel {
    fn omega(&self, s: f64, v: f64) -> f64;
    fn m(&self, v: f64) -> f64;
    fn gamma(&self, s: f64) -> f64;
    fn mu(&self, current_vol: f64) -> f64;
    fn sigma_v(&self, v: f64) -> f64;
    fn rho(&self) -> f64;
    // Derived: theta_X from eq (9), g, f, g_inv, etc.
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
        Self {
            rho,
            sigma,
            r_f,
            q,
            eta,
            ref_vol,
        }
    }
    pub fn from_options(option: &Options) -> Self {
        let eta = 1.0; // placeholder; calibrate to market data in practice
        let rho = 0.0; // placeholder; calibrate to market data in practice
        let sigma = option.volatility(); // placeholder; use implied vol surface in practice
        let ref_vol = option.volatility(); // placeholder; calibrate to historical vol in practice
        Self {
            rho,
            sigma,
            r_f: option.risk_free_rate(),
            q: option.dividend_yield(),
            eta,
            ref_vol,
        }
    }
    fn heston_g_integral(&self, s: f64) -> f64 {
        s.ln()
    }

    fn heston_f_integral(&self, volatility: f64) -> f64 {
        volatility / self.sigma
    }

    fn heston_g_inverse(&self, x: f64) -> f64 {
        x.exp()
    }

    fn heston_f_inverse(&self, x: f64) -> f64 {
        x * self.sigma
    }
    fn heston_m_prime(&self, v: f64) -> f64 {
        0.5 / v.sqrt()
    }
    fn heston_sigma_prime_v(&self, v: f64) -> f64 {
        0.5 * self.sigma / v.sqrt()
    }
    fn heston_gamma_prime(&self, _s: f64) -> f64 {
        1.0
    }
}

impl SLVModel for Heston {
    fn omega(&self, s: f64, _v: f64) -> f64 {
        s * (self.r_f - self.q.unwrap_or(0.0))
    }
    fn m(&self, v: f64) -> f64 {
        v.sqrt()
    }
    fn gamma(&self, s: f64) -> f64 {
        s
    }
    fn mu(&self, current_vol: f64) -> f64 {
        self.eta * (self.ref_vol - current_vol)
    }
    fn sigma_v(&self, v: f64) -> f64 {
        self.sigma * v.sqrt()
    }
    fn rho(&self) -> f64 {
        self.rho
    }
}

pub enum SLVModelVariant {
    Heston(Heston),
    // SABR(SABR),
}

impl SLVModelVariant {
    pub fn from_options(option: &Options, model_type: SLVModelVariant) -> Self {
        match model_type {
            SLVModelVariant::Heston(_) => SLVModelVariant::Heston(Heston::from_options(option)),
        }
    }
    pub fn g_integral(&self, s: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_g_integral(s),
        }
    }
    pub fn f_integral(&self, volatility: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_f_integral(volatility),
        }
    }
    pub fn g_inverse(&self, x: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_g_inverse(x),
        }
    }
    pub fn f_inverse(&self, x: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_f_inverse(x),
        }
    }
    pub fn m_prime(&self, v: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_m_prime(v),
        }
    }
    pub fn sigma_prime_v(&self, v: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_sigma_prime_v(v),
        }
    }
    pub fn gamma_prime(&self, s: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_gamma_prime(s),
        }
    }
}

impl SLVModel for SLVModelVariant {
    fn omega(&self, s: f64, v: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.omega(s, v),
        }
    }
    fn m(&self, v: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.m(v),
        }
    }
    fn gamma(&self, s: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.gamma(s),
        }
    }
    fn mu(&self, current_vol: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.mu(current_vol),
        }
    }
    fn sigma_v(&self, v: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.sigma_v(v),
        }
    }
    fn rho(&self) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.rho(),
        }
    }
}

fn sigma_tilde_sq(slv_model: &SLVModelVariant, vol_grid: &[f64]) -> Vec<f64> {
    vol_grid
        .iter()
        .map(|&vol| {
            let m = slv_model.m(vol);
            let rho = slv_model.rho();
            m * m * (1.0 - rho * rho)
        })
        .collect()
}

fn psi_v(slv_model: &SLVModelVariant, v: f64) -> f64 {
    let m = slv_model.m(v);
    let sigma = slv_model.sigma_v(v);
    slv_model.mu(v) * m / sigma
        + 0.5 * (sigma * slv_model.m_prime(v) - slv_model.sigma_prime_v(v) * m)
}

fn theta(slv_model: &SLVModelVariant, x_grid: &[f64], vol_grid: &[f64]) -> DMatrix<f64> {
    let mut theta = DMatrix::zeros(vol_grid.len(), x_grid.len());
    for (i, &vol) in vol_grid.iter().enumerate() {
        let psi = psi_v(slv_model, vol);
        for j in 0..x_grid.len() {
            let s = slv_model.g_inverse(x_grid[j] + slv_model.rho() * slv_model.f_integral(vol));
            theta[(i, j)] = slv_model.omega(s, vol) / slv_model.gamma(s)
                - slv_model.gamma_prime(s) * slv_model.m(vol).powi(2) * 0.5
                - slv_model.rho() * psi;
        }
    }
    theta
}

// Size should be small enough to avoid numerical instability but balances against
// the cost of matrix exponentiation. 50 steps is a reasonable baseline.
fn uniform_variance_grid(vol_0: f64, n_vol_steps: i32) -> Vec<f64> {
    let step = 4.0 * vol_0 / (n_vol_steps as f64);
    (1..=n_vol_steps).map(|i| i as f64 * step).collect()
}

/// Sinh-stretched variance grid per Cui, Zhenyu 2018.
fn sinh_variance_grid(_vol: f64, _n_vol_steps: i32) -> Vec<f64> {
    todo!()
}

fn x_space_grid(
    slv_model: &SLVModelVariant,
    m_steps: i32,
    spot_0: f64,
    vol_0: f64,
) -> (Vec<f64>, Vec<f64>) {
    let x_0 = slv_model.g_integral(spot_0) - slv_model.rho() * slv_model.f_integral(vol_0);
    if x_0 <= 0.0 {
        panic!("Invalid initial x space grid point: x_0 must be positive");
    }
    let step = x_0 * (4.0 - 1e-3) / (m_steps as f64);
    let grid = (1..=m_steps).map(|i| i as f64 * step).collect();
    (grid, vec![step; m_steps as usize])
}

// ────────────────────────────────────────────────────────────────
// Variance generator Λ (M×M)
// ────────────────────────────────────────────────────────────────

pub struct VarianceGeneratorCOO {
    pub rows: Vec<usize>,
    pub cols: Vec<usize>,
    pub vals: Vec<f64>,
    pub num_states: usize,
}

/// Discretize the variance SDE dv = μ(v)dt + σ_v(v)dW_v onto the grid v_grid
/// using the same upwind finite-difference stencil as the x-space G_l matrices.
/// Boundary states (l=0, l=M-1) are absorbing.
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
        let lambda_diag = -(lambda_up + lambda_dn);

        rows.push(l);
        cols.push(l + 1);
        vals.push(lambda_up);

        rows.push(l);
        cols.push(l - 1);
        vals.push(lambda_dn);

        rows.push(l);
        cols.push(l);
        vals.push(lambda_diag);
    }

    VarianceGeneratorCOO { rows, cols, vals, num_states: m }
}

// Dense wrapper for Λ (used in Kronecker product)

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

// Per-regime x-space generator G_l (N×N)

/// Builds the N×N spatial generator for variance state v_l.
/// Boundary rows (i=0, i=N-1) are absorbing
/// Off-diagonal rates are clamped to zero if the grid is too coarse
/// relative to the drift.
fn build_g_l_matrix_inner(
    x_spacings: &[f64],
    theta_col: &[f64],
    sigma_tilde_sq_l: f64,
    n: usize,
) -> DMatrix<f64> {
    let mut g_l = DMatrix::zeros(n, n);
    for i in 1..(n - 1) {
        let delta_i = x_spacings[i];
        let delta_im1 = x_spacings[i - 1];
        let delta_sum = delta_i + delta_im1;
        let theta_i = theta_col[i];

        let q_up = ((sigma_tilde_sq_l + theta_i * delta_im1) / (delta_i * delta_sum)).max(0.0);
        let q_dn = ((sigma_tilde_sq_l - theta_i * delta_i) / (delta_im1 * delta_sum)).max(0.0);

        g_l[(i, i + 1)] = q_up;
        g_l[(i, i - 1)] = q_dn;
        g_l[(i, i)] = -(q_up + q_dn);
    }
    g_l
}
/// Combined NM×NM generator  G = Λ⊗I_N + blockdiag(G_1,...,G_M)
fn build_combined_generator(
    variance_gen: &VarianceGeneratorDense,
    g_l_matrices: &[DMatrix<f64>],
    n: usize,
    m: usize,
) -> DMatrix<f64> {
    let nm = n * m;
    let mut g = DMatrix::zeros(nm, nm);

    // Kronecker product Λ⊗I_N: each scalar λ_{l,l'} contributes to the
    // diagonal of the (l,l') block.
    for l in 0..m {
        for lp in 0..m {
            let lambda = variance_gen.lambda[(l, lp)];
            if lambda.abs() < 1e-30 {
                continue;
            }
            for i in 0..n {
                g[(l * n + i, lp * n + i)] += lambda;
            }
        }
    }

    // Block diagonal: add each G_l into its (l,l) diagonal block.
    for l in 0..m {
        let g_l = &g_l_matrices[l];
        debug_assert_eq!(g_l.nrows(), n);
        debug_assert_eq!(g_l.ncols(), n);
        for i in 0..n {
            for j in 0..n {
                let val = g_l[(i, j)];
                if val.abs() < 1e-30 {
                    continue;
                }
                g[(l * n + i, l * n + j)] += val;
            }
        }
    }

    g
}

// CombinedCTMC: state, index mapping, and diagnostics

pub struct CombinedCTMC {
    /// NM×NM combined generator matrix G
    pub generator: DMatrix<f64>,
    /// X-space grid points, length N
    pub x_grid: Vec<f64>,
    /// X-space spacings δ_i = x_{i+1} - x_i, length N-1
    pub x_spacings: Vec<f64>,
    /// Variance grid points, length M
    pub v_grid: Vec<f64>,
    pub n: usize,
    pub m: usize,
}

impl CombinedCTMC {
    /// Index mapping φ(i, l) = l·N + i  (0-indexed).
    #[inline]
    pub fn phi(&self, i: usize, l: usize) -> usize {
        l * self.n + i
    }
    /// Inverse: flat index → (x-index, v-index).
    #[inline]
    pub fn phi_inv(&self, idx: usize) -> (usize, usize) {
        (idx % self.n, idx / self.n)
    }
    /// Total dimension NM.
    #[inline]
    pub fn dim(&self) -> usize {
        self.n * self.m
    }
    /// Unit vector e_{φ(i,l)}: all zeros except 1 at position phi(i, l).
    pub fn selector_vec(&self, i: usize, l: usize) -> Vec<f64> {
        let mut e = vec![0.0; self.dim()];
        e[self.phi(i, l)] = 1.0;
        e
    }
    /// Maximum absolute row sum; should be < 1e-12 for a valid generator.
    pub fn max_row_sum_error(&self) -> f64 {
        let nm = self.dim();
        (0..nm)
            .map(|row| (0..nm).map(|col| self.generator[(row, col)]).sum::<f64>().abs())
            .fold(0.0_f64, f64::max)
    }
    /// Returns (count of violations, worst negative value) for off-diagonal entries.
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

    // pub fn print_diagnostics(&self) {
    //     let nm = self.dim();
    //     let nnz = self.generator.iter().filter(|&&v| v.abs() > 1e-30).count();
    //     let sparsity = 1.0 - (nnz as f64 / (nm * nm) as f64);
    //     println!("=== Combined Generator G Diagnostics ===");
    //     println!("Dimensions: {}x{} (N={}, M={})", nm, nm, self.n, self.m);
    //     println!("Nonzeros: {} / {} ({:.1}% sparse)", nnz, nm * nm, sparsity * 100.0);
    //     println!("Row-sum error: {:.2e}", self.max_row_sum_error());
    //     let (violations, worst) = self.check_off_diagonal_non_negative();
    //     if violations > 0 {
    //         println!(
    //             "WARNING: {} negative off-diagonal entries (worst: {:.2e})",
    //             violations, worst
    //         );
    //     } else {
    //         println!("Off-diagonal non-negativity: OK");
    //     }
    //     let diag_min =
    //         (0..nm).map(|i| self.generator[(i, i)]).fold(f64::INFINITY, f64::min);
    //     let diag_max =
    //         (0..nm).map(|i| self.generator[(i, i)]).fold(f64::NEG_INFINITY, f64::max);
    //     println!("Diagonal range: [{:.4}, {:.4}]", diag_min, diag_max);
    // }
}

/// Public orchestrator
/// Build the combined NM×NM CTMC generator G from model and grids.
///
/// Computes θ and σ̃² internally, constructs each G_l, builds Λ from the
/// variance SDE, and assembles G = Λ⊗I_N + blockdiag(G_1,...,G_M).
pub fn assemble_generator(
    slv_model: &SLVModelVariant,
    x_grid: &[f64],
    x_spacings: &[f64],
    v_grid: &[f64],
) -> CombinedCTMC {
    let n = x_grid.len();
    let m = v_grid.len();

    let theta_mat = theta(slv_model, x_grid, v_grid);
    let sigma_tilde = sigma_tilde_sq(slv_model, v_grid);

    let g_l_matrices: Vec<DMatrix<f64>> = (0..m)
        .map(|l| {
            let theta_col: Vec<f64> = (0..n).map(|i| theta_mat[(l, i)]).collect();
            build_g_l_matrix_inner(x_spacings, &theta_col, sigma_tilde[l], n)
        })
        .collect();

    let var_coo = build_variance_generator(slv_model, v_grid);
    let var_dense = VarianceGeneratorDense::from_coo(&var_coo);

    let generator = build_combined_generator(&var_dense, &g_l_matrices, n, m);

    CombinedCTMC {
        generator,
        x_grid: x_grid.to_vec(),
        x_spacings: x_spacings.to_vec(),
        v_grid: v_grid.to_vec(),
        n,
        m,
    }
}
