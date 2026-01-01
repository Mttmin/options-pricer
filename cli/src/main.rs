use black_scholes::black_scholes_price;
use options::{Call, Options};

fn main() {
    println!("Black-Scholes Option Pricing CLI");
    // create a dummy call option
    let call_option = Options::Call(Call {
        strike_price: 100.0,
        spot_price: 105.0,
        volatility: 0.2,
        risk_free_rate: 0.05,
        time_to_maturity: 1.0,
    });
    let price = black_scholes_price(call_option);
    println!("Call Option Price: {:.4}", price);
}
