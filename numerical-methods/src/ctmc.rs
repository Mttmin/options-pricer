use nalgebra::{DMatrix};
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
        let eta = 1.0; // This is a placeholder. In practice, eta should be calibrated to market data
        let rho = 0.0; // This is a placeholder. In practice, rho should be calibrated to market data
        let sigma = option.volatility(); // This is a placeholder. In practice, sigma should be calibrated to market data, potentially using the implied volatility surface
        let ref_vol = option.volatility(); // This is a placeholder. In practice, ref_vol is calibrated to the historical volatility of the underlying asset
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
        // Heston-specific g integral
        s.ln()
    }

    fn heston_f_integral(&self, volatility: f64) -> f64 {
        // Heston-specific f integral
        volatility / self.sigma
    }

    fn heston_g_inverse(&self, x: f64) -> f64 {
        // Heston-specific g inverse
        x.exp()
    }

    fn heston_f_inverse(&self, x: f64) -> f64 {
        // Heston-specific f inverse
        x * self.sigma
    }
    fn heston_m_prime(&self, v: f64) -> f64 {
        0.5 / v.sqrt()
    }
    fn heston_sigma_prime_v(&self, v: f64) -> f64 {
        0.5 * self.sigma / v.sqrt()
    }
    fn heston_gamma_prime(&self, s: f64) -> f64 {
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
            SLVModelVariant::Heston(h) => SLVModelVariant::Heston(h),
        }
    }
    pub fn g_integral(&self, s: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_g_integral(s),
            // SLVModelVariant::SABR(s) => sabr_g_integral(s, option),
        }
    }
    pub fn f_integral(&self, volatility: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_f_integral(volatility),
            // SLVModelVariant::SABR(s) => sabr_f_integral(s, volatility),
        }
    }
    pub fn g_inverse(&self, x: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_g_inverse(x),
            // SLVModelVariant::SABR(s) => sabr_g_inverse(s, x),
        }
    }
    pub fn f_inverse(&self, x: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_f_inverse(x),
            // SLVModelVariant::SABR(s) => sabr_f_inverse(s, x),
        }
    }
    pub fn m_prime(&self, v: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_m_prime(v),
            // SLVModelVariant::SABR(s) => sabr_m_prime(s, v),
        }
    }
    pub fn sigma_prime_v(&self, v: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_sigma_prime_v(v),
            // SLVModelVariant::SABR(s) => sabr_sigma_prime_v(s, v),
        }
    }
    pub fn gamma_prime(&self, s: f64) -> f64 {
        match self {
            SLVModelVariant::Heston(h) => h.heston_gamma_prime(s),
            // SLVModelVariant::SABR(s) => sabr_gamma_prime(s, s),
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

fn sigma_tilde_sq(slv_model: &SLVModelVariant, vol_grid: &Vec<f64>) -> Vec<f64> {
    let mut sigma_tilde = Vec::new();
    for &vol in vol_grid.iter() {
        let m: f64 = slv_model.m(vol);
        let rho = slv_model.rho();
        sigma_tilde.push(m * m * (1.0 - rho * rho));
    }
    sigma_tilde
}

fn psi_v(slv_model: &SLVModelVariant, v: f64) -> f64 {
    let m = slv_model.m(v);
    let sigma = slv_model.sigma_v(v);
    slv_model.mu(v) * m / sigma
        + 0.5 * (sigma * slv_model.m_prime(v) - slv_model.sigma_prime_v(v) * m)
}

fn theta(slv_model: &SLVModelVariant, x_grid: &Vec<f64>, vol_grid: &Vec<f64>) -> DMatrix<f64> {
    let mut theta = DMatrix::zeros(vol_grid.len(), x_grid.len());
    for (i, &vol) in vol_grid.iter().enumerate() {
        let psi = psi_v(slv_model, vol);
        for j in 0..x_grid.len() {
            let s = slv_model.g_inverse(x_grid[j] + slv_model.rho() * slv_model.f_integral(vol));
            theta[(i, j)] = slv_model.omega(s, vol) / slv_model.gamma(s)
                - slv_model.gamma_prime(s) * slv_model.m(vol).powf(2.0) * 0.5
                - slv_model.rho() * psi;
        }
    }
    theta
}

// size should be relatively small to avoid numerical instability, but this contributes to the exponentiation of the G matrix which is time consuming, so we need to find a balance.
// As a baseline we choose 50 steps
fn uniform_variance_grid(vol_0: f64, n_vol_steps: i32) -> Vec<f64> {
    let mut grid = Vec::new();
    let step = 4.0 * vol_0 / (n_vol_steps as f64);
    for i in 1..=n_vol_steps {
        grid.push(i as f64 * step);
    }
    grid
}

/// Follows a sinh stretching of the volatility grid based on advice from:
/// "A General Valuation Framework for SABR and Stochastic Local Volatility Models- Cui, Zhenyu,2018"
fn sinh_variance_grid(vol: f64, n_vol_steps: i32) -> Vec<f64> {
    todo!()
}

fn x_space_grid(
    slv_model: SLVModelVariant,
    n_x_steps: i32,
    spot_0: f64,
    vol_0: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut grid = Vec::new();
    let x_0 = slv_model.g_integral(spot_0) - slv_model.rho() * slv_model.f_integral(vol_0);
    if x_0 <= 0.0 {
        panic!("Invalid initial x space grid point: x_0 must be positive");
    }
    let step = x_0 * (4.0 - 10.0_f64.powf(-3.0)) / (n_x_steps as f64);
    for i in 1..=n_x_steps {
        grid.push(i as f64 * step);
    }
    (grid, vec![step; n_x_steps as usize])
}

// fn cmtc_solver<T>(option : Options, n_steps: i32, slv_model: &T) -> f64
//     where T: SLVModel {
//     let j_matrix: DMatrix<f64> = calc_h1();
//     let strike = option.strike_price();
//     let volatility = option.volatility();
//     let bn_v = slv_model.g_integral(strike) - slv_model.rho()*slv_model.f_integral(volatility);
//     let a_matrix: DMatrix<f64> = exp(G - option.risk_free_rate()*Identity);
//     for k in (n_steps-1)..0 {
//         j_matrix = A*(j_matrix + delta_h3);
//         for l in 1..=m_dim {
//             let vol = 0.0;
//             let bk = find_bk(option.strike_price(),j_matrix, vol);
//         }
//     }
//     (b_tilde, j_matrix)
//     todo!()
// }
