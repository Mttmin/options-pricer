use rand::SeedableRng;
use rand_distr::{Normal, Distribution};
use rand::rng;
use options::{Options, Payoff, MonteCarloParameters};
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
                let payoff_plus = derivative.compute(spot_t_plus, 0.0);

                // Compute payoff for -z path (antithetic variate)
                let spot_t_minus = spot * (drift + diffusion * (-z)).exp();
                let payoff_minus = derivative.compute(spot_t_minus, 0.0);

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
    use std::time::Instant;

    /// Speed comparison test: Monte Carlo vs Price Path
    /// Run all speedtests: cargo test --release speedtest -- --ignored --nocapture
    #[test]
    #[ignore]
    fn speedtest_monte_carlo_vs_price_path() {
        // Create a test call option
        let call_option = options::Options::new_call(
            100.0,  // strike
            100.0,  // spot
            0.2,    // volatility (20%)
            0.05,   // risk-free rate (5%)
            1.0,    // time to maturity (1 year)
            None,   // no dividend
        );

        let num_simulations = 100_000_000;  
        // let num_paths = 10_000;         
        // let num_steps = 356;           

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
        let result_european = price_european::<options::Call>(&call_option, num_simulations);
        let duration_european = start.elapsed();
        
        println!("   Price: ${:.4}", result_european.price);
        println!("   Std Error: ${:.4}", result_european.std_error);
        let (lower, upper) = result_european.confidence_interval_95();
        println!("   95% CI: [${:.4}, ${:.4}]", lower, upper);
        println!("   Error vs BS: ${:.4} ({:.2}%)", 
            (result_european.price - bs_price).abs(),
            ((result_european.price - bs_price).abs() / bs_price * 100.0));
        println!("   Time: {:?}", duration_european);
        
        // // Test price_path function
        // println!("\n2. Testing price_path ({} paths, {} steps each)...", num_paths, num_steps);
        // let start = Instant::now();
        
        // // Simulate multiple paths and compute average payoff
        // let mut total_payoff = 0.0;
        // let discount = (-call_option.risk_free_rate() * call_option.time_to_maturity()).exp();
        
        // for _ in 0..num_paths {
        //     let path = price_path(
        //         call_option.spot_price(),
        //         call_option.volatility(),
        //         call_option.risk_free_rate(),
        //         num_steps,
        //         call_option.time_to_maturity(),
        //     );
        //     let final_price = path[path.len() - 1];
        //     let payoff = (final_price - call_option.strike_price()).max(0.0);
        //     total_payoff += payoff;
        // }
        
        // let avg_price = discount * total_payoff / num_paths as f64;
        // let duration_path = start.elapsed();
        
        // println!("   Price: ${:.4}", avg_price);
        // println!("   Error vs BS: ${:.4} ({:.2}%)", 
        //     (avg_price - bs_price).abs(),
        //     ((avg_price - bs_price).abs() / bs_price * 100.0));
        // println!("   Time: {:?}", duration_path);
        
        // Comparison
        println!("\n=== Performance Comparison ===");
        println!("Black-Scholes:   {:?} (baseline)", duration_bs);
        println!("Monte Carlo is {:.0}x slower than Black-Scholes for a total of {:.3}s", duration_european.as_secs_f64()/duration_bs.as_secs_f64(), duration_european.as_secs_f64());
        // println!("\nprice_european is {:.2}x faster than price_path", 
        //     duration_path.as_secs_f64() / duration_european.as_secs_f64());
        
        println!("\n=== Accuracy Comparison (vs Black-Scholes ${:.4}) ===", bs_price);
        println!("price_european error: ${:.4}", (result_european.price - bs_price).abs());
        // println!("price_path error:     ${:.4}", (avg_price - bs_price).abs());
    }

    /// Speed comparison test: Straddle vs Call
    /// Run all speedtests: cargo test --release speedtest -- --ignored --nocapture
    #[test]
    #[ignore]
    fn speedtest_straddle_vs_call() {
        use options::optionspreads::{OptionSpreads, Direction};

        let strike = 100.0;
        let spot = 100.0;
        let vol = 0.2;
        let rfr = 0.05;
        let ttm = 1.0;

        let straddle = OptionSpreads::new_straddle(
            Direction::LONG,
            strike,
            spot,
            vol,
            rfr,
            ttm,
            None,
        );

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
        println!("  Straddle/Call: {:.2}x", straddle_mc_time.as_secs_f64() / call_mc_time.as_secs_f64());
    }
}
