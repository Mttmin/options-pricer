pub mod binomial;
use nalgebra::{DMatrix, DVector};
use options::{MonteCarloParameters, Options, Payoff};
use rand::SeedableRng;
use rand::rng;
use rand_distr::{Distribution, Normal};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

// fn brownian_step(spot_t: f64, vol: f64, random_step: f64, risk_free: f64, delta_t: f64) -> f64 {
//     spot_t * ((risk_free - vol * vol / 2.0) * delta_t + vol * random_step * delta_t.sqrt()).exp()
// }

pub fn price_path(
    spot_0: f64,
    vol: f64,
    risk_free: f64,
    num_steps: u32,
    time_simulated: f64,
) -> Vec<f64> {
    let normal = Normal::new(0.0, 1.0).unwrap();
    let delta_t = time_simulated / num_steps as f64;
    let sqrt_dt = delta_t.sqrt();
    let mut rng = rng();
    let mut result = vec![0.0; num_steps as usize + 1];

    result[0] = spot_0;
    let mut spot_t = spot_0;

    for result_entry in result.iter_mut().skip(1) {
        let z = normal.sample(&mut rng);
        spot_t *= ((risk_free - vol * vol / 2.0) * delta_t + vol * z * sqrt_dt).exp();
        *result_entry = spot_t;
    }

    result
}

pub struct MonteCarloResult {
    pub price: f64,
    pub std_error: f64,
}

impl MonteCarloResult {
    pub fn confidence_interval_95(&self) -> (f64, f64) {
        (
            self.price - 1.96 * self.std_error,
            self.price + 1.96 * self.std_error,
        )
    }
}

struct SimulationConstants {
    drift: f64,
    diffusion: f64,
    discount: f64,
}

impl SimulationConstants {
    fn new(opt: &Options) -> Self {
        let t = opt.time_to_maturity();
        let r = opt.risk_free_rate();
        let sigma = opt.volatility();

        Self {
            drift: (r - 0.5 * sigma * sigma) * t,
            diffusion: sigma * t.sqrt(),
            discount: (-r * t).exp(),
        }
    }
}

/// Prices European options using Monte Carlo with antithetic variates variance reduction.
///
/// Uses compile-time monomorphization for zero-cost abstraction on simple call/puts
pub fn price_european<P: Payoff>(opt: &options::Options, num_simulations: u32) -> MonteCarloResult {
    let constants = SimulationConstants::new(opt);

    let (sum, sum_sq): (f64, f64) = (0..num_simulations)
        .into_par_iter()
        .map_init(
            || {
                let rng = Xoshiro256PlusPlus::from_rng(&mut rand::rng());
                let normal = Normal::new(0.0, 1.0).unwrap();
                (rng, normal)
            },
            |(rng, normal), _| {
                let z: f64 = normal.sample(rng);

                // Compute payoff for +z path
                let spot_t_plus =
                    opt.spot_price() * (constants.drift + constants.diffusion * z).exp();
                let payoff_plus = P::compute_static(spot_t_plus, opt.strike_price());

                // Compute payoff for -z path (antithetic variate)
                let spot_t_minus =
                    opt.spot_price() * (constants.drift + constants.diffusion * (-z)).exp();
                let payoff_minus = P::compute_static(spot_t_minus, opt.strike_price());

                // Average the two payoffs
                let avg_payoff = (payoff_plus + payoff_minus) / 2.0;
                (avg_payoff, avg_payoff * avg_payoff)
            },
        )
        .reduce(|| (0.0, 0.0), |(a, a2), (b, b2)| (a + b, a2 + b2));

    let n = num_simulations as f64;
    let mean_payoff = sum / n;
    let price = constants.discount * mean_payoff;
    let variance = (sum_sq / n) - (mean_payoff * mean_payoff);
    let std_error = constants.discount * (variance / n).sqrt();

    MonteCarloResult { price, std_error }
}

/// Generic Monte Carlo pricing for any derivative implementing Payoff and MonteCarloParameters.
///
/// For simple Call/Put options, prefer `price_european()` for better performance.
pub fn monte_carlo<T>(derivative: &T, num_simulations: u32) -> MonteCarloResult
where
    T: Payoff + MonteCarloParameters + Sync,
{
    let t = derivative.time_to_maturity();
    let r = derivative.risk_free_rate();
    let sigma = derivative.volatility();
    let spot = derivative.spot_price();

    let drift = (r - 0.5 * sigma * sigma) * t;
    let diffusion = sigma * t.sqrt();
    let discount = (-r * t).exp();

    let (sum, sum_sq): (f64, f64) = (0..num_simulations)
        .into_par_iter()
        .map_init(
            || {
                let rng = Xoshiro256PlusPlus::from_rng(&mut rand::rng());
                let normal = Normal::new(0.0, 1.0).unwrap();
                (rng, normal)
            },
            |(rng, normal), _| {
                let z: f64 = normal.sample(rng);

                // Compute payoff for +z path
                let spot_t_plus = spot * (drift + diffusion * z).exp();
                let payoff_plus = derivative.compute(spot_t_plus);

                // Compute payoff for -z path (antithetic variate)
                let spot_t_minus = spot * (drift + diffusion * (-z)).exp();
                let payoff_minus = derivative.compute(spot_t_minus);
                // Average the two payoffs
                let avg_payoff = (payoff_plus + payoff_minus) / 2.0;
                (avg_payoff, avg_payoff * avg_payoff)
            },
        )
        .reduce(|| (0.0, 0.0), |(a, a2), (b, b2)| (a + b, a2 + b2));

    let n = num_simulations as f64;
    let mean_payoff = sum / n;
    let price = discount * mean_payoff;
    let variance = (sum_sq / n) - (mean_payoff * mean_payoff);
    let std_error = discount * (variance / n).sqrt();

    MonteCarloResult { price, std_error }
}

/// Spatial grid with log-spacing
pub struct SpatialGrid {
    pub prices: Vec<f64>,
    pub num_nodes: usize,
}

impl SpatialGrid {
    pub fn new_log_spaced<T>(derivative: &T, num_nodes: usize, width_multiplier: f64) -> Self
    where
        T: Payoff + MonteCarloParameters,
    {
        let char_price = derivative.spot_price();
        let sigma = derivative.volatility();
        let t = derivative.time_to_maturity();

        let std_range = sigma * t.sqrt() * char_price;
        // For American puts, extend grid much further down (early exercise near S=0)
        // Use 0.1% of spot as minimum instead of 1% for better put pricing
        let min_price = (char_price - width_multiplier * std_range).max(char_price * 0.001);
        let max_price = char_price + width_multiplier * std_range;

        let log_min = min_price.ln();
        let log_max = max_price.ln();

        let prices: Vec<f64> = (0..num_nodes)
            .map(|i| {
                let log_s = log_min + (log_max - log_min) * (i as f64) / ((num_nodes - 1) as f64);
                log_s.exp()
            })
            .collect();

        Self { prices, num_nodes }
    }

    /// A_i = (S_{i+1} - S_{i-1})/2
    pub fn cell_width(&self, i: usize) -> f64 {
        if i == 0 || i == self.num_nodes - 1 {
            panic!("Cell width undefined at boundaries");
        }
        (self.prices[i + 1] - self.prices[i - 1]) / 2.0
    }
}

/// Time discretization
/// Time discretization
pub struct TimeGrid {
    pub times: Vec<f64>, // τ values going backwards from expiry
    pub num_steps: usize,
}

impl TimeGrid {
    pub fn new_uniform(num_steps: usize, time_to_maturity: f64) -> Self {
        let dt = time_to_maturity / num_steps as f64;
        let times: Vec<f64> = (0..=num_steps).map(|n| n as f64 * dt).collect();
        Self { times, num_steps }
    }

    pub fn new_empty(num_steps: usize) -> Self {
        let times: Vec<f64> = vec![0.0; num_steps + 1];
        Self { times, num_steps }
    }

    pub fn delta_tau(&self, n: usize) -> f64 {
        self.times[n + 1] - self.times[n]
    }
}

pub enum TimestepMode {
    Constant,
    Adaptive { dnorm: f64, scale_d: f64 }, // dnorm = target relative change, scale_d for normalization
}

/// Main PDE solver for American options
pub struct PenaltySolver<T>
where
    T: Payoff + MonteCarloParameters,
{
    pub derivative: T,
    pub spatial_grid: SpatialGrid,
    pub time_grid: TimeGrid,

    /// V[n][i] = option value at time τ_n, spatial node i
    pub values: Vec<Vec<f64>>,

    /// Payoff V*(S_i) at each node
    pub payoff: Vec<f64>,

    /// Precomputed coefficients: (coeff_left, coeff_right, gamma_left, gamma_right) for each interior node
    pub coefficients: Vec<(f64, f64, f64, f64)>,

    /// Parameters
    pub penalty_large: f64,
    pub tolerance: f64,
    pub timestep_mode: TimestepMode,

    /// Tracking metrics
    pub max_american_error: f64,
    pub total_iterations: usize,
}

impl<T> PenaltySolver<T>
where
    T: Payoff + MonteCarloParameters,
{
    pub fn new(
        derivative: T,
        num_spatial_nodes: usize,
        num_timesteps: usize,
        penalty_tolerance: f64,
        timestep_mode: TimestepMode,
    ) -> Self {
        // Use 10 std devs for puts to ensure grid extends close to S=0 (critical for American exercise)
        let spatial_grid = SpatialGrid::new_log_spaced(&derivative, num_spatial_nodes, 10.0);
        let time_grid = match timestep_mode {
            TimestepMode::Constant => {
                TimeGrid::new_uniform(num_timesteps, derivative.time_to_maturity())
            }
            TimestepMode::Adaptive { .. } => TimeGrid::new_empty(num_timesteps),
        };

        let values = vec![vec![0.0; num_spatial_nodes]; num_timesteps + 1];
        let payoff: Vec<f64> = spatial_grid
            .prices
            .iter()
            .map(|&s| derivative.compute(s))
            .collect();

        // Precompute spatial discretization coefficients for all interior nodes
        let r = derivative.risk_free_rate();
        let sigma = derivative.volatility();
        let coefficients: Vec<(f64, f64, f64, f64)> = (1..num_spatial_nodes - 1)
            .into_par_iter()
            .map(|i| {
                let s_i = spatial_grid.prices[i];
                let a_i = (spatial_grid.prices[i + 1] - spatial_grid.prices[i - 1]) / 2.0;
                let s_diff_left = (spatial_grid.prices[i] - spatial_grid.prices[i - 1]).abs();
                let s_diff_right = (spatial_grid.prices[i + 1] - spatial_grid.prices[i]).abs();

                let gamma_left = (sigma * sigma * s_i * s_i) / (2.0 * a_i * s_diff_left);
                let gamma_right = (sigma * sigma * s_i * s_i) / (2.0 * a_i * s_diff_right);

                let u_i = r * s_i;
                let beta_left = if u_i < 0.0 { -u_i / a_i } else { 0.0 };
                let beta_right = if u_i > 0.0 { u_i / a_i } else { 0.0 };

                (
                    gamma_left + beta_left,
                    gamma_right + beta_right,
                    gamma_left,
                    gamma_right,
                )
            })
            .collect();

        // from paper recommendations, going under 10^-4 is overkill as prices are in cents
        let penalty_large = 1.0 / penalty_tolerance;
        Self {
            derivative,
            spatial_grid,
            time_grid,
            values,
            payoff,
            coefficients,
            penalty_large,
            tolerance: penalty_tolerance,
            timestep_mode,
            max_american_error: 0.0,
            total_iterations: 0,
        }
    }
    /// Initialize terminal condition V(S, τ=0) = payoff
    pub fn initialize(&mut self) {
        self.values[0] = self.payoff.clone();

        // // Debug: print payoff at a few key points
        // let spot_idx = self.spatial_grid.prices.iter()
        //     .position(|&s| (s - self.derivative.spot_price()).abs() < 1.0)
        //     .unwrap_or(self.spatial_grid.num_nodes / 2);
        // println!("Initial payoff at S={:.2}: {:.4}",
        //          self.spatial_grid.prices[spot_idx],
        //          self.payoff[spot_idx]);
        // println!("Max payoff: {:.4}", self.payoff.iter().cloned().fold(f64::NEG_INFINITY, f64::max));

        // Initialize adaptive timestep grid if needed
        if let TimestepMode::Adaptive { .. } = self.timestep_mode {
            self.initialize_adaptive_timesteps();
        }
    }

    #[allow(dead_code)]
    /// Get precomputed coefficients for interior node i (i in 1..N-1)
    fn get_coefficients(&self, i: usize) -> (f64, f64, f64, f64) {
        // Coefficients are stored for interior nodes only (indices 1..N-1)
        // So coefficient for node i is at index i-1
        self.coefficients[i - 1]
    }

    /// Initialize first timestep for adaptive mode (equation 10.1 heuristic)
    fn initialize_adaptive_timesteps(&mut self) {
        self.time_grid.times[0] = 0.0;
        // Start with reasonable first timestep (target average)
        let initial_dt = self.derivative.time_to_maturity() / (self.time_grid.num_steps as f64);
        self.time_grid.times[1] = initial_dt;
    }

    /// Penalty iteration for one timestep (Algorithm 5.2 from paper) - Optimized
    fn penalty_iteration(&mut self, n: usize, theta: f64) -> usize {
        let dt = self.time_grid.delta_tau(n);
        let n_nodes = self.spatial_grid.num_nodes;

        // Build M matrix using precomputed coefficients (30% faster than recomputing)
        let r = self.derivative.risk_free_rate();
        let mut m = DMatrix::zeros(n_nodes, n_nodes);

        // Vectorized matrix construction using precomputed coefficients
        for i in 1..n_nodes - 1 {
            let (coeff_left, coeff_right, _, _) = self.coefficients[i - 1];
            m[(i, i - 1)] = -dt * coeff_left;
            m[(i, i)] = dt * (coeff_left + coeff_right + r);
            m[(i, i + 1)] = -dt * coeff_right;
        }

        let v_n = DVector::from_vec(self.values[n].clone());
        let v_star = DVector::from_vec(self.payoff.clone());
        let mut v_current = v_n.clone();

        // Precompute matrices that don't change across penalty iterations
        let identity = DMatrix::identity(n_nodes, n_nodes);
        let lhs_base = &identity + &m * (1.0 - theta);
        let rhs_base = (&identity - &m * theta) * &v_n;

        let mut iterations: usize = 0;
        let max_iterations: usize = 50;

        for _ in 0..max_iterations {
            // Build penalty diagonal (vectorized, enforces V >= payoff)
            let p_diag: DVector<f64> = DVector::from_fn(n_nodes, |i, _| {
                if v_current[i] < self.payoff[i] {
                    self.penalty_large
                } else {
                    0.0
                }
            });

            // (I + (1-θ)M + P) V^{n+1} = (I - θM) V^n + P*V*
            let mut lhs = lhs_base.clone();
            for i in 0..n_nodes {
                lhs[(i, i)] += p_diag[i];
            }
            let rhs = &rhs_base + p_diag.component_mul(&v_star);

            let v_new = lhs.lu().solve(&rhs).expect("failed to solve system");

            // Parallel convergence check (faster for large grids)
            let max_change = (0..n_nodes)
                .into_par_iter()
                .map(|i| {
                    let change = (v_new[i] - v_current[i]).abs();
                    change / v_new[i].abs().max(1.0)
                })
                .reduce(|| 0.0, f64::max);

            let constraint_switches = (0..n_nodes).any(|i| {
                let was_constrained = v_current[i] < self.payoff[i];
                let is_constrained = v_new[i] < self.payoff[i];
                was_constrained != is_constrained
            });

            v_current = v_new;
            iterations += 1;

            if max_change < self.tolerance && !constraint_switches {
                break;
            }
        }

        if iterations >= max_iterations {
            eprintln!(
                "Warning: Penalty iteration did not converge at timestep {}",
                n
            );
        }
        self.values[n + 1] = v_current.as_slice().to_vec();

        iterations
    }

    /// Adaptive timestep selector (equation 10.1)
    fn compute_next_timestep(&self, n: usize, dnorm: f64, scale_d: f64) -> f64 {
        let n_nodes = self.spatial_grid.num_nodes;
        let mut min_ratio = f64::INFINITY;
        for i in 0..n_nodes {
            let v_new = self.values[n + 1][i];
            let v_old = self.values[n][i];
            let change = (v_new - v_old).abs();
            let denom = v_new.abs().max(v_old.abs()).max(scale_d);
            if change > 1e-10 {
                let ratio = dnorm * denom / change;
                min_ratio = min_ratio.min(ratio);
            }
        }
        let dt_current = self.time_grid.delta_tau(n);
        (min_ratio * dt_current).min(dt_current * 2.0) // Cap growth at 2x
    }

    /// Run the full solver with Rannacher smoothing
    pub fn solve(&mut self) {
        self.initialize();

        println!(
            "Solving with {} spatial nodes, {} timesteps...",
            self.spatial_grid.num_nodes, self.time_grid.num_steps
        );

        // Rannacher smoothing: 2 fully implicit steps, then Crank-Nicolson (Section 7)
        let rannacher_steps = 2;

        for n in 0..self.time_grid.num_steps {
            // Rannacher smoothing: fully implicit for first 2 steps, then Crank-Nicolson
            let theta = if n < rannacher_steps { 1.0 } else { 0.5 };

            // Solve timestep
            let iters = self.penalty_iteration(n, theta);
            self.total_iterations += iters;

            // // Print progress every 10 timesteps
            // if n % 10 == 0 || n == self.time_grid.num_steps - 1 {
            //     let dt = self.time_grid.delta_tau(n);
            //     let tau = self.time_grid.times[n + 1];
            //     println!("  Step {}/{}: τ={:.6}, dt={:.6}, iters={}",
            //              n + 1, self.time_grid.num_steps, tau, dt, iters);
            // }

            // Update timestep for adaptive mode
            if let TimestepMode::Adaptive { dnorm, scale_d } = self.timestep_mode
                && n + 2 <= self.time_grid.num_steps
            {
                let tau_max = self.derivative.time_to_maturity();
                let tau_current = self.time_grid.times[n + 1];
                let remaining_time = tau_max - tau_current;
                let remaining_steps = self.time_grid.num_steps - (n + 1);

                // For last few steps, ensure we reach tau_max exactly
                if remaining_steps <= 3 {
                    self.time_grid.times[n + 2] =
                        tau_current + remaining_time / remaining_steps as f64;
                } else {
                    let dt_next = self.compute_next_timestep(n, dnorm, scale_d);
                    // Don't let timestep overshoot remaining time
                    let dt_safe = dt_next.min(remaining_time / remaining_steps as f64 * 2.0);
                    self.time_grid.times[n + 2] = (tau_current + dt_safe).min(tau_max);
                }
            }
        }

        println!("Solving complete!");
    }

    #[allow(dead_code)]
    /// Compute max American constraint error (equation 4.20)
    fn compute_american_error(&self) -> f64 {
        let mut max_error: f64 = 0.0;
        for n in 0..=self.time_grid.num_steps {
            for i in 0..self.spatial_grid.num_nodes {
                let violation = (self.payoff[i] - self.values[n][i]).max(0.0);
                let normalized = violation / self.payoff[i].max(1.0);
                max_error = max_error.max(normalized);
            }
        }

        max_error
    }

    /// Get option value at initial spot price
    pub fn option_value(&self) -> f64 {
        // Find bracketing nodes for spot price
        let spot = self.derivative.spot_price();

        // Find i such that prices[i] <= spot < prices[i+1]
        let mut bracket_idx = 0;
        for i in 0..self.spatial_grid.num_nodes - 1 {
            if self.spatial_grid.prices[i] <= spot && spot < self.spatial_grid.prices[i + 1] {
                bracket_idx = i;
                break;
            }
        }

        // Linear interpolation between bracketing nodes
        let s_low = self.spatial_grid.prices[bracket_idx];
        let s_high = self.spatial_grid.prices[bracket_idx + 1];
        let v_low = self.values[self.time_grid.num_steps][bracket_idx];
        let v_high = self.values[self.time_grid.num_steps][bracket_idx + 1];

        let weight = (spot - s_low) / (s_high - s_low);
        v_low + weight * (v_high - v_low)
    }

    /// Print solver diagnostics
    pub fn print_diagnostics(&mut self) {
        println!("\n=== Solver Diagnostics ===");
        println!(
            "Spatial grid: {} nodes from ${:.2} to ${:.2}",
            self.spatial_grid.num_nodes,
            self.spatial_grid.prices[0],
            self.spatial_grid.prices[self.spatial_grid.num_nodes - 1]
        );
        println!(
            "Time grid: {} steps, final τ = {:.6}",
            self.time_grid.num_steps, self.time_grid.times[self.time_grid.num_steps]
        );
        println!("Total iterations: {}", self.total_iterations);
        println!(
            "Avg iterations/step: {:.2}",
            self.total_iterations as f64 / self.time_grid.num_steps as f64
        );

        // Compute and print max American error
        let error = self.compute_american_error();
        println!("Max American constraint violation: {:.2e}", error);
        self.max_american_error = error;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use binomial::{BinomialParameters, ExerciseStyle, binomial_price};
    use options::{Call, Put};
    use std::time::Instant;

    /// Speed comparison test: Monte Carlo vs Price Path
    /// Run all speedtests: cargo test --release speedtest -- --ignored --nocapture
    #[test]
    #[ignore]
    fn speedtest_monte_carlo_vs_price_path() {
        // Create a test call option
        let call_option = options::Options::new_call(
            100.0, // strike
            100.0, // spot
            0.2,   // volatility (20%)
            0.05,  // risk-free rate (5%)
            1.0,   // time to maturity (1 year)
            None,  // no dividend
        );

        let num_simulations = 100_000_000;
        let num_paths = 10_000;
        let num_steps = 356;

        println!("\n=== Monte Carlo Speed Test ===");
        println!("Option parameters: K=100, S=100, sigma=0.2, r=0.05, T=1.0");

        // Black-Scholes analytical price (baseline)
        println!("\n0. Black-Scholes Analytical Price (instant)...");
        let start = Instant::now();
        let bs_price = call_option.bs_pricing();
        let duration_bs = start.elapsed();
        println!("   Price: ${:.4}", bs_price);
        println!("   Time: {:?}", duration_bs);

        // Test price_european function
        println!(
            "\n1. Testing price_european ({} simulations, monomorphized, parallel)...",
            num_simulations
        );
        let start = Instant::now();
        let result_european = price_european::<options::Call>(&call_option, num_simulations);
        let duration_european = start.elapsed();

        println!("   Price: ${:.4}", result_european.price);
        println!("   Std Error: ${:.4}", result_european.std_error);
        let (lower, upper) = result_european.confidence_interval_95();
        println!("   95% CI: [${:.4}, ${:.4}]", lower, upper);
        println!(
            "   Error vs BS: ${:.4} ({:.2}%)",
            (result_european.price - bs_price).abs(),
            ((result_european.price - bs_price).abs() / bs_price * 100.0)
        );
        println!("   Time: {:?}", duration_european);

        // Test price_path function
        println!(
            "\n2. Testing price_path ({} paths, {} steps each)...",
            num_paths, num_steps
        );
        let start = Instant::now();

        // Simulate multiple paths and compute average payoff
        let mut total_payoff = 0.0;
        let discount = (-call_option.risk_free_rate() * call_option.time_to_maturity()).exp();

        for _ in 0..num_paths {
            let path = price_path(
                call_option.spot_price(),
                call_option.volatility(),
                call_option.risk_free_rate(),
                num_steps,
                call_option.time_to_maturity(),
            );
            let final_price = path[path.len() - 1];
            let payoff = (final_price - call_option.strike_price()).max(0.0);
            total_payoff += payoff;
        }

        let avg_price = discount * total_payoff / num_paths as f64;
        let duration_path = start.elapsed();

        println!("   Price: ${:.4}", avg_price);
        println!(
            "   Error vs BS: ${:.4} ({:.2}%)",
            (avg_price - bs_price).abs(),
            ((avg_price - bs_price).abs() / bs_price * 100.0)
        );
        println!("   Time: {:?}", duration_path);

        // Comparison
        println!("\n=== Performance Comparison ===");
        println!("Black-Scholes:   {:?} (baseline)", duration_bs);
        println!(
            "\nprice_european is {:.2}x faster than price_path",
            duration_path.as_secs_f64() / duration_european.as_secs_f64()
        );

        println!(
            "\n=== Accuracy Comparison (vs Black-Scholes ${:.4}) ===",
            bs_price
        );
        println!(
            "price_european error: ${:.4}",
            (result_european.price - bs_price).abs()
        );
        // println!("price_path error:     ${:.4}", (avg_price - bs_price).abs());
    }

    /// Speed comparison test: Straddle vs Call
    /// Run all speedtests: cargo test --release speedtest -- --ignored --nocapture
    #[test]
    #[ignore]
    fn speedtest_straddle_vs_call() {
        use options::optionspreads::{Direction, OptionSpreads};

        let strike = 100.0;
        let spot = 100.0;
        let vol = 0.2;
        let rfr = 0.05;
        let ttm = 1.0;

        let straddle =
            OptionSpreads::new_straddle(Direction::LONG, strike, spot, vol, rfr, ttm, None);

        let call = options::Options::new_call(strike, spot, vol, rfr, ttm, None);

        println!("\n=== Straddle vs Call Performance Comparison ===");
        // println!("Parameters: K={}, S={}, sigma={}, r={}, T={}\n", strike, spot, vol, rfr, ttm);

        // // BS Pricing
        // let start = Instant::now();
        // let straddle_bs = options::BlackScholesPrice::bs_pricing(&straddle);
        // let straddle_bs_time = start.elapsed();

        // let start = Instant::now();
        // let call_bs = options::BlackScholesPrice::bs_pricing(&call);
        // let call_bs_time = start.elapsed();

        // println!("Black-Scholes:");
        // println!("  Straddle: ${:.4} in {:?}", straddle_bs, straddle_bs_time);
        // println!("  Call:     ${:.4} in {:?}\n", call_bs, call_bs_time);

        // Monte Carlo with 1B simulations
        println!("Monte Carlo (1B simulations):");
        let num_sims = 1_000_000_000;

        let start = Instant::now();
        let _straddle_mc = monte_carlo(&straddle, num_sims);
        let straddle_mc_time = start.elapsed();

        let start = Instant::now();
        let _call_mc = price_european::<options::Call>(&call, num_sims);
        let call_mc_time = start.elapsed();

        // println!("  Straddle: ${:.4} ±{:.4} in {:?}", straddle_mc.price, straddle_mc.std_error, straddle_mc_time);
        // println!("  Call:     ${:.4} ±{:.4} in {:?}\n", call_mc.price, call_mc.std_error, call_mc_time);

        println!("Performance Ratio:");
        println!(
            "  Straddle/Call: {:.2}x",
            straddle_mc_time.as_secs_f64() / call_mc_time.as_secs_f64()
        );
    }

    /// Test binomial tree convergence to Black-Scholes as tree depth increases
    /// Run with: cargo test --release -- --ignored --nocapture binomial_convergence
    #[test]
    #[ignore]
    fn binomial_convergence_to_black_scholes() {
        println!("\n=== Binomial Tree Convergence to Black-Scholes ===");

        // Test parameters
        let spot = 100.0;
        let strike = 100.0;
        let volatility = 0.2;
        let risk_free = 0.05;
        let time = 1.0;

        let call_option = Options::new_call(strike, spot, volatility, risk_free, time, None);
        let put_option = Options::new_put(strike, spot, volatility, risk_free, time, None);

        let call = Call::new(strike, spot, volatility, risk_free, time, None);
        let put = Put::new(strike, spot, volatility, risk_free, time, None);

        // Calculate Black-Scholes prices
        let bs_call = call_option.bs_pricing();
        let bs_put = put_option.bs_pricing();

        println!("\nBlack-Scholes prices:");
        println!("  Call: ${:.6}", bs_call);
        println!("  Put:  ${:.6}", bs_put);

        let depths = vec![10, 50, 100, 200, 500, 1000, 2000];

        println!(
            "\n{:>6} | {:>12} | {:>12} | {:>10} | {:>10}",
            "Depth", "Call Price", "Put Price", "Call Err%", "Put Err%"
        );
        println!("{}", "-".repeat(62));

        for &depth in &depths {
            let params = BinomialParameters::new(spot, depth, time, volatility, risk_free);

            let binomial_call = binomial_price(&call, &params, ExerciseStyle::European);
            let binomial_put = binomial_price(&put, &params, ExerciseStyle::European);

            let call_error = ((binomial_call - bs_call) / bs_call * 100.0).abs();
            let put_error = ((binomial_put - bs_put) / bs_put * 100.0).abs();

            println!(
                "{:>6} | ${:>11.6} | ${:>11.6} | {:>9.4}% | {:>9.4}%",
                depth, binomial_call, binomial_put, call_error, put_error
            );
        }

        // Test convergence assertion with high depth
        let params_high = BinomialParameters::new(spot, 2000, time, volatility, risk_free);
        let binomial_call_high = binomial_price(&call, &params_high, ExerciseStyle::European);
        let binomial_put_high = binomial_price(&put, &params_high, ExerciseStyle::European);

        let call_error_pct = ((binomial_call_high - bs_call) / bs_call * 100.0).abs();
        let put_error_pct = ((binomial_put_high - bs_put) / bs_put * 100.0).abs();

        println!("\nConvergence test (depth=2000):");
        println!("  Call error: {:.4}%", call_error_pct);
        println!("  Put error:  {:.4}%", put_error_pct);

        assert!(
            call_error_pct < 0.1,
            "Call price should converge within 0.1%"
        );
        assert!(put_error_pct < 0.1, "Put price should converge within 0.1%");
    }

    /// Speed test for binomial tree pricing
    /// Run with: cargo test --release -- --ignored --nocapture binomial_speed
    #[test]
    #[ignore]
    fn binomial_speed_test() {
        println!("\n=== Binomial Tree Speed Test ===");
        println!("Option parameters: K=100, S=100, sigma=0.2, r=0.05, T=1.0");

        let spot = 100.0;
        let strike = 100.0;
        let volatility = 0.2;
        let risk_free = 0.05;
        let time = 1.0;

        let call_option = Options::new_call(strike, spot, volatility, risk_free, time, None);
        let call = Call::new(strike, spot, volatility, risk_free, time, None);

        // Black-Scholes baseline
        println!("\n--- Black-Scholes (analytical) ---");
        let start = Instant::now();
        let bs_price = call_option.bs_pricing();
        let duration_bs = start.elapsed();
        println!("Price: ${:.6}", bs_price);
        println!("Time: {:?}", duration_bs);

        let depths = vec![10, 50, 100, 500, 1000, 2000, 5000];

        println!(
            "\n{:>6} | {:>12} | {:>12} | {:>10}",
            "Depth", "European", "American", "Time"
        );
        println!("{}", "-".repeat(45));

        for &depth in &depths {
            let params = BinomialParameters::new(spot, depth, time, volatility, risk_free);

            let start = Instant::now();
            let european_price = binomial_price(&call, &params, ExerciseStyle::European);
            let american_price = binomial_price(&call, &params, ExerciseStyle::American);
            let duration = start.elapsed();

            println!(
                "{:>6} | ${:>11.6} | ${:>11.6} | {:>10?}",
                depth, european_price, american_price, duration
            );
        }

        // Compare with Monte Carlo
        println!("\n--- Monte Carlo Comparison (1M simulations) ---");
        let start = Instant::now();
        let mc_result = price_european::<Call>(&call_option, 1_000_000);
        let duration_mc = start.elapsed();
        println!("Price: ${:.6}", mc_result.price);
        println!("Std Error: ${:.6}", mc_result.std_error);
        println!("Time: {:?}", duration_mc);

        // Summary
        let params_1000 = BinomialParameters::new(spot, 1000, time, volatility, risk_free);
        let binomial_time = std::time::Instant::now();
        let _binomial_price_result = binomial_price(&call, &params_1000, ExerciseStyle::European);
        let binomial_duration = binomial_time.elapsed();

        println!("\n=== Speed Comparison ===");
        println!("Black-Scholes: {:?}", duration_bs);
        println!("Binomial (1000 steps): {:?}", binomial_duration);
        println!("Monte Carlo (1M sims): {:?}", duration_mc);
    }

    /// Compare American vs European options using binomial tree
    /// Run with: cargo test --release -- --ignored --nocapture american_vs_european
    #[test]
    #[ignore]
    fn american_vs_european_comparison() {
        println!("\n=== American vs European Options Comparison ===");

        let depth = 1000;
        let volatility = 0.2;
        let risk_free = 0.05;
        let time = 1.0;

        // Test scenarios: ATM, ITM, OTM for both calls and puts
        let scenarios = vec![
            ("ATM", 100.0, 100.0),
            ("ITM", 100.0, 110.0), // For call: spot > strike
            ("OTM", 100.0, 90.0),  // For call: spot < strike
            ("Deep ITM", 100.0, 120.0),
            ("Deep OTM", 100.0, 80.0),
        ];

        println!("\n--- Call Options ---");
        println!(
            "{:>10} | {:>12} | {:>12} | {:>12} | {:>10}",
            "Scenario", "European", "American", "Premium", "Premium %"
        );
        println!("{}", "-".repeat(63));

        for (name, strike, spot) in &scenarios {
            let call = Call::new(*strike, *spot, volatility, risk_free, time, None);
            let params = BinomialParameters::new(*spot, depth, time, volatility, risk_free);

            let european = binomial_price(&call, &params, ExerciseStyle::European);
            let american = binomial_price(&call, &params, ExerciseStyle::American);
            let premium = american - european;
            let premium_pct = if european > 0.0 {
                premium / european * 100.0
            } else {
                0.0
            };

            println!(
                "{:>10} | ${:>11.6} | ${:>11.6} | ${:>11.6} | {:>9.4}%",
                name, european, american, premium, premium_pct
            );
        }

        println!("\n--- Put Options ---");
        println!(
            "{:>10} | {:>12} | {:>12} | {:>12} | {:>10}",
            "Scenario", "European", "American", "Premium", "Premium %"
        );
        println!("{}", "-".repeat(63));

        for (name, strike, spot) in &scenarios {
            let put = Put::new(*strike, *spot, volatility, risk_free, time, None);
            let params = BinomialParameters::new(*spot, depth, time, volatility, risk_free);

            let european = binomial_price(&put, &params, ExerciseStyle::European);
            let american = binomial_price(&put, &params, ExerciseStyle::American);
            let premium = american - european;
            let premium_pct = if european > 0.0 {
                premium / european * 100.0
            } else {
                0.0
            };

            println!(
                "{:>10} | ${:>11.6} | ${:>11.6} | ${:>11.6} | {:>9.4}%",
                name, european, american, premium, premium_pct
            );
        }

        // Verify that American options are always worth at least as much as European
        println!("\n--- Verification Tests ---");
        let test_cases = vec![
            (100.0, 100.0, "ATM Call"),
            (100.0, 80.0, "ITM Put"),
            (100.0, 120.0, "ITM Call"),
            (100.0, 90.0, "OTM Put"),
        ];

        for (strike, spot, description) in test_cases {
            let call = Call::new(strike, spot, volatility, risk_free, time, None);
            let put = Put::new(strike, spot, volatility, risk_free, time, None);
            let params = BinomialParameters::new(spot, depth, time, volatility, risk_free);

            let euro_call = binomial_price(&call, &params, ExerciseStyle::European);
            let amer_call = binomial_price(&call, &params, ExerciseStyle::American);
            let euro_put = binomial_price(&put, &params, ExerciseStyle::European);
            let amer_put = binomial_price(&put, &params, ExerciseStyle::American);

            println!(
                "{}: Call American >= European: {} (${:.6} >= ${:.6})",
                description,
                amer_call >= euro_call - 1e-6,
                amer_call,
                euro_call
            );
            println!(
                "{}: Put American >= European: {} (${:.6} >= ${:.6})",
                description,
                amer_put >= euro_put - 1e-6,
                amer_put,
                euro_put
            );

            assert!(
                amer_call >= euro_call - 1e-6,
                "American call should be worth at least as much as European"
            );
            assert!(
                amer_put >= euro_put - 1e-6,
                "American put should be worth at least as much as European"
            );
        }

        println!("\nAll verification tests passed!");
    }

    /// Basic unit test for binomial pricing accuracy
    #[test]
    fn test_binomial_pricing_accuracy() {
        let spot = 100.0;
        let strike = 100.0;
        let volatility = 0.2;
        let risk_free = 0.05;
        let time = 1.0;

        let call = Call::new(strike, spot, volatility, risk_free, time, None);
        let call_option = Options::new_call(strike, spot, volatility, risk_free, time, None);

        let params = BinomialParameters::new(spot, 1000, time, volatility, risk_free);
        let binomial_price_result = binomial_price(&call, &params, ExerciseStyle::European);
        let bs_price = call_option.bs_pricing();

        let error_pct = ((binomial_price_result - bs_price) / bs_price * 100.0).abs();

        assert!(
            error_pct < 0.5,
            "Binomial price (${:.4}) should be within 0.5% of Black-Scholes (${:.4}), got {:.2}% error",
            binomial_price_result,
            bs_price,
            error_pct
        );
    }

    /// Basic unit test for American option premium
    #[test]
    fn test_american_premium() {
        let spot = 80.0;
        let strike = 100.0;
        let volatility = 0.2;
        let risk_free = 0.05;
        let time = 1.0;

        // Deep ITM put should have significant early exercise value
        let put = Put::new(strike, spot, volatility, risk_free, time, None);
        let params = BinomialParameters::new(spot, 500, time, volatility, risk_free);

        let european = binomial_price(&put, &params, ExerciseStyle::European);
        let american = binomial_price(&put, &params, ExerciseStyle::American);

        assert!(
            american > european,
            "American put should be more valuable than European for deep ITM"
        );

        let premium = american - european;
        assert!(
            premium > 0.5,
            "Deep ITM American put should have significant premium (${:.4})",
            premium
        );
    }
}
