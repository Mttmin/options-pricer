pub mod fetcher;
pub mod options_chain;
pub mod volatility; 
use options::{BlackScholesPrice, Options, black_scholes::black_scholes_price};
use options::exotics::{AsianKind, AsianOption, AsianPooling, ConvertibleBond};
use numerical_methods::{PenaltySolver, TimestepMode, monte_carlo_asian};

#[tokio::main]
async fn main() {
    println!("Black-Scholes Option Pricing CLI");
    // create a dummy call option
    let call_option = Options::new_call(
        100.0,
        105.0,
        0.2,
        0.05,
        1.0,
        None
    );
    let price = black_scholes_price(call_option);
    println!("Call Option Price: {:.4}", price);

    // create a dummy convertible bond
    let convertible_bond = ConvertibleBond {
        face_value: 100.0,
        coupon_rate: 0.05,
        maturity: 5.0,
        payment_frequency: 2,
        risk_free_rate: 0.03,
        credit_spread: 0.02,
        conversion_price: 50.0,
        stock_price: 55.0,
        volatility: 0.25,
        time_to_maturity: 5.0,
        dividend_yield: None
    };
    let cb_price = convertible_bond.bs_pricing();
    println!("Convertible Bond Price: {:.4}", cb_price);

    // // Demonstrate live market data fetching
    // println!("\n Live Market");

    // match fetcher::load_api_keys() {
    //     Ok(_config) => {
    //         println!("API keys loaded successfully!");

    //         let fetcher = fetcher::DataFetcher::new();

    //         let symbol = "AAPL";
    //         let lookback_days = 180;

    //         println!("Fetching market data for {}...", symbol);

    //         match fetcher.update_symbol(symbol, lookback_days).await {
    //             Ok(market_data) => {
    //                 println!("\nMarket Data for {}:", symbol);
    //                 println!("  Spot Price: ${:.2}", market_data.spot_price);
    //                 println!("  Volatility (annualized): {:.2}%", market_data.volatility * 100.0);
    //                 println!("  Last Updated: {}", market_data.last_price_update);
    //                 println!("  Corporate Action Detected: {}", market_data.corporate_action_detected);

    //                 // Price a real option using live data
    //                 let strike = market_data.spot_price * 1.05; // 5% OTM call
    //                 let time_to_expiry = 0.25; // 3 months
    //                 let risk_free_rate = 0.05;

    //                 let live_option = Options::new_call(
    //                     market_data.spot_price,
    //                     strike,
    //                     market_data.volatility,
    //                     risk_free_rate,
    //                     time_to_expiry,
    //                     None,
    //                 );

    //                 let live_price = black_scholes_price(live_option);
    //                 println!("\nLive Option Pricing:");
    //                 println!("  Symbol: {}", symbol);
    //                 println!("  Strike: ${:.2}", strike);
    //                 println!("  Time to Expiry: {} months", time_to_expiry * 12.0);
    //                 println!("  Option Price: ${:.4}", live_price);
    //             },
    //             Err(e) => {
    //                 eprintln!("Error fetching market data: {}", e);
    //                 eprintln!("Note: API rate limits may apply. Alpha Vantage has a limit of 5 calls/minute for free tier.");
    //             }
    //         }
    //         let vol = volatility::vix_volatility("AAPL", fetcher,120).await.unwrap();
    //         println!("{}", vol);}
    //     Err(e) => {
    //         eprintln!("Failed to load API keys: {}", e);
    //         eprintln!("To use live market data, ensure api_keys.json exists in the project root.");
    //     }
    // }
    let put = options::Put::new(100.0, 100.0, 0.2, 0.1, 0.25, None);

    let mut solver = PenaltySolver::new(
        put,
        200,     // Spatial nodes
        100,     // Max timesteps
        1e-6,
        TimestepMode::Adaptive { dnorm: 0.05, scale_d: 1.0 },
    );
    solver.initialize();
    solver.solve();
    println!("value of the quadratic solver {}", solver.option_value());

    // Asian options demo
    println!("\nAsian Options");

    let asian_avg_call = AsianOption::new(
        AsianKind::Call, AsianPooling::Average,
        100.0, 100.0, 0.3, 0.05, 1.0, None, 50,
    );
    let avg_cf = asian_avg_call.bs_pricing();
    let avg_mc = monte_carlo_asian(&asian_avg_call, 200_000, 50, true);
    println!(
        "  Avg-pool Call  K=100 S=100 sigma=0.30 r=0.05 T=1  N=50",
    );
    println!(
        "    Kemna-Vorst (geom CF): {:.4}",
        avg_cf,
    );
    println!(
        "    Arithmetic MC (CV):    {:.4}  +/- {:.4}  (95% CI)",
        avg_mc.price, 1.96 * avg_mc.std_error,
    );

    let asian_max_put = AsianOption::new(
        AsianKind::Put, AsianPooling::Max,
        100.0, 100.0, 0.3, 0.05, 1.0, None, 100,
    );
    let max_cf = asian_max_put.bs_pricing();
    let max_mc = monte_carlo_asian(&asian_max_put, 200_000, 200, false);
    println!(
        "  Max-pool Put   K=100 S=100 sigma=0.30 r=0.05 T=1  N=100",
    );
    println!("    Conze lookback CF:     {:.4}", max_cf);
    println!(
        "    Path MC:               {:.4}  +/- {:.4}  (95% CI)",
        max_mc.price, 1.96 * max_mc.std_error,
    );
}
