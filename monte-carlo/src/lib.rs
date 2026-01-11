pub mod binomial;
use options::MonteCarloParameters;
use rand::SeedableRng;
use rand_distr::{Normal, Distribution};
use rand::rng;
use options::{Options, Payoff};
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rand_xoshiro::Xoshiro256PlusPlus;

// fn brownian_step(spot_t: f64, vol: f64, random_step: f64, risk_free: f64, delta_t: f64) -> f64 {
//     spot_t * ((risk_free - vol * vol / 2.0) * delta_t + vol * random_step * delta_t.sqrt()).exp()
// }

pub fn price_path(spot_0: f64, vol: f64, risk_free: f64, num_steps: u32, time_simulated: f64) -> Vec<f64> {
    let normal = Normal::new(0.0, 1.0).unwrap();
    let delta_t = time_simulated / num_steps as f64;
    let sqrt_dt = delta_t.sqrt();
    let mut rng = rng();
    let mut result = vec![0.0; num_steps as usize + 1];
    
    result[0] = spot_0;
    let mut spot_t = spot_0;
    
    for i in 1..=num_steps as usize {
        let z = normal.sample(&mut rng);
        spot_t *= ((risk_free - vol * vol / 2.0) * delta_t + vol * z * sqrt_dt).exp();
        result[i] = spot_t;
    }
    
    result
}

pub struct MonteCarloResult {
    pub price: f64,
    pub std_error: f64,
}

impl MonteCarloResult {
    pub fn confidence_interval_95(&self) -> (f64, f64) {
        (self.price - 1.96 * self.std_error, self.price + 1.96 * self.std_error)
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

/// Prices European options using Monte Carlo with antithetic variates variance reduction..
pub fn price_european<P: Payoff>(opt: &Options, num_simulations: u32) -> MonteCarloResult {
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
                let spot_t_plus = opt.spot_price() * (constants.drift + constants.diffusion * z).exp();
                let payoff_plus = P::compute_static(spot_t_plus, opt.strike_price());

                // Compute payoff for -z path (antithetic variate)
                let spot_t_minus = opt.spot_price() * (constants.drift + constants.diffusion * (-z)).exp();
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


#[cfg(test)]
mod tests {
    use super::*;
    use options::{Call, Put};
    use std::time::Instant;
    use binomial::{BinomialParameters, binomial_price, ExerciseStyle};

    /// Speed comparison test
    /// Run with: cargo test --release -- --ignored --nocapture speed_test
    #[test]
    #[ignore]
    fn speed_test_monte_carlo_vs_price_path() {
        // Create a test call option
        let call_option = Options::new_call(
            100.0,  // strike
            100.0,  // spot
            0.2,    // volatility (20%)
            0.05,   // risk-free rate (5%)
            1.0,    // time to maturity (1 year)
            None,   // no dividend
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
        println!("\n1. Testing price_european ({} simulations, monomorphized, parallel)...", num_simulations);
        let start = Instant::now();
        let result_european = price_european::<Call>(&call_option, num_simulations);
        let duration_european = start.elapsed();

        println!("   Price: ${:.4}", result_european.price);
        println!("   Std Error: ${:.4}", result_european.std_error);
        let (lower, upper) = result_european.confidence_interval_95();
        println!("   95% CI: [${:.4}, ${:.4}]", lower, upper);
        println!("   Error vs BS: ${:.4} ({:.2}%)",
            (result_european.price - bs_price).abs(),
            ((result_european.price - bs_price).abs() / bs_price * 100.0));
        println!("   Time: {:?}", duration_european);

        // Test price_path function
        println!("\n2. Testing price_path ({} paths, {} steps each)...", num_paths, num_steps);
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
        println!("   Error vs BS: ${:.4} ({:.2}%)",
            (avg_price - bs_price).abs(),
            ((avg_price - bs_price).abs() / bs_price * 100.0));
        println!("   Time: {:?}", duration_path);

        // Comparison
        println!("\n=== Performance Comparison ===");
        println!("Black-Scholes:   {:?} (baseline)", duration_bs);
        println!("\nprice_european is {:.2}x faster than price_path",
            duration_path.as_secs_f64() / duration_european.as_secs_f64());

        println!("\n=== Accuracy Comparison (vs Black-Scholes ${:.4}) ===", bs_price);
        println!("price_european error: ${:.4}", (result_european.price - bs_price).abs());
        println!("price_path error:     ${:.4}", (avg_price - bs_price).abs());
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

        println!("\n{:>6} | {:>12} | {:>12} | {:>10} | {:>10}",
            "Depth", "Call Price", "Put Price", "Call Err%", "Put Err%");
        println!("{}", "-".repeat(62));

        for &depth in &depths {
            let params = BinomialParameters::new(spot, depth, time, volatility, risk_free);

            let binomial_call = binomial_price(&call, &params, ExerciseStyle::European);
            let binomial_put = binomial_price(&put, &params, ExerciseStyle::European);

            let call_error = ((binomial_call - bs_call) / bs_call * 100.0).abs();
            let put_error = ((binomial_put - bs_put) / bs_put * 100.0).abs();

            println!("{:>6} | ${:>11.6} | ${:>11.6} | {:>9.4}% | {:>9.4}%",
                depth, binomial_call, binomial_put, call_error, put_error);
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

        assert!(call_error_pct < 0.1, "Call price should converge within 0.1%");
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

        println!("\n{:>6} | {:>12} | {:>12} | {:>10}",
            "Depth", "European", "American", "Time");
        println!("{}", "-".repeat(45));

        for &depth in &depths {
            let params = BinomialParameters::new(spot, depth, time, volatility, risk_free);

            let start = Instant::now();
            let european_price = binomial_price(&call, &params, ExerciseStyle::European);
            let american_price = binomial_price(&call, &params, ExerciseStyle::American);
            let duration = start.elapsed();

            println!("{:>6} | ${:>11.6} | ${:>11.6} | {:>10?}",
                depth, european_price, american_price, duration);
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
            ("ITM", 100.0, 110.0),  // For call: spot > strike
            ("OTM", 100.0, 90.0),   // For call: spot < strike
            ("Deep ITM", 100.0, 120.0),
            ("Deep OTM", 100.0, 80.0),
        ];

        println!("\n--- Call Options ---");
        println!("{:>10} | {:>12} | {:>12} | {:>12} | {:>10}",
            "Scenario", "European", "American", "Premium", "Premium %");
        println!("{}", "-".repeat(63));

        for (name, strike, spot) in &scenarios {
            let call = Call::new(*strike, *spot, volatility, risk_free, time, None);
            let params = BinomialParameters::new(*spot, depth, time, volatility, risk_free);

            let european = binomial_price(&call, &params, ExerciseStyle::European);
            let american = binomial_price(&call, &params, ExerciseStyle::American);
            let premium = american - european;
            let premium_pct = if european > 0.0 { premium / european * 100.0 } else { 0.0 };

            println!("{:>10} | ${:>11.6} | ${:>11.6} | ${:>11.6} | {:>9.4}%",
                name, european, american, premium, premium_pct);
        }

        println!("\n--- Put Options ---");
        println!("{:>10} | {:>12} | {:>12} | {:>12} | {:>10}",
            "Scenario", "European", "American", "Premium", "Premium %");
        println!("{}", "-".repeat(63));

        for (name, strike, spot) in &scenarios {
            let put = Put::new(*strike, *spot, volatility, risk_free, time, None);
            let params = BinomialParameters::new(*spot, depth, time, volatility, risk_free);

            let european = binomial_price(&put, &params, ExerciseStyle::European);
            let american = binomial_price(&put, &params, ExerciseStyle::American);
            let premium = american - european;
            let premium_pct = if european > 0.0 { premium / european * 100.0 } else { 0.0 };

            println!("{:>10} | ${:>11.6} | ${:>11.6} | ${:>11.6} | {:>9.4}%",
                name, european, american, premium, premium_pct);
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

            println!("{}: Call American >= European: {} (${:.6} >= ${:.6})",
                description, amer_call >= euro_call - 1e-6, amer_call, euro_call);
            println!("{}: Put American >= European: {} (${:.6} >= ${:.6})",
                description, amer_put >= euro_put - 1e-6, amer_put, euro_put);

            assert!(amer_call >= euro_call - 1e-6,
                "American call should be worth at least as much as European");
            assert!(amer_put >= euro_put - 1e-6,
                "American put should be worth at least as much as European");
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

        assert!(error_pct < 0.5,
            "Binomial price (${:.4}) should be within 0.5% of Black-Scholes (${:.4}), got {:.2}% error",
            binomial_price_result, bs_price, error_pct);
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

        assert!(american > european,
            "American put should be more valuable than European for deep ITM");

        let premium = american - european;
        assert!(premium > 0.5,
            "Deep ITM American put should have significant premium (${:.4})", premium);
    }
}
