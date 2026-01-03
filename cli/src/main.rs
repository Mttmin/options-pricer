use options::{black_scholes::black_scholes_price, Options, exotics::ConvertibleBond};

fn main() {
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
}
