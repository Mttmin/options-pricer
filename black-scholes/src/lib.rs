pub use options::{Call, Options, Put};
use statrs::distribution::{ContinuousCDF, Normal};

// calculate d1 for the Black-Scholes formula
fn d_plus(t: f64, r: f64, sigma: f64, spot: f64, strike: f64) -> f64 {
    let numerator = (spot / strike).ln() + (r + 0.5 * sigma * sigma) * t;
    let denominator = sigma * t.sqrt();
    numerator / denominator
}

// calculate d2 for the Black-Scholes formula
fn d_minus(t: f64, r: f64, sigma: f64, spot: f64, strike: f64) -> f64 {
    let numerator = (spot / strike).ln() + (r - 0.5 * sigma * sigma) * t;
    let denominator = sigma * t.sqrt();
    numerator / denominator
}
/// Calculate the Black-Scholes price for a given option, either Call or Put.
///
/// Needs the option parameters encapsulated in the Options enum.
pub fn black_scholes_price(option: Options) -> f64 {
    let std_norm = Normal::new(0.0, 1.0).unwrap();
    if let Options::Call(call) = option {
        let d1 = d_plus(
            call.time_to_maturity,
            call.risk_free_rate,
            call.volatility,
            call.spot_price,
            call.strike_price,
        );
        let d2 = d_minus(
            call.time_to_maturity,
            call.risk_free_rate,
            call.volatility,
            call.spot_price,
            call.strike_price,
        );
        let nd1: f64 = std_norm.cdf(d1);
        let nd2: f64 = std_norm.cdf(d2);
        call.spot_price * nd1
            - call.strike_price * (-call.risk_free_rate * call.time_to_maturity).exp() * nd2
    } else if let Options::Put(put) = option {
        let d1 = d_plus(
            put.time_to_maturity,
            put.risk_free_rate,
            put.volatility,
            put.spot_price,
            put.strike_price,
        );
        let d2 = d_minus(
            put.time_to_maturity,
            put.risk_free_rate,
            put.volatility,
            put.spot_price,
            put.strike_price,
        );
        let n_neg_d1 = std_norm.cdf(-d1);
        let n_neg_d2 = std_norm.cdf(-d2);
        put.strike_price * (-put.risk_free_rate * put.time_to_maturity).exp() * n_neg_d2
            - put.spot_price * n_neg_d1
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_black_scholes_call() {
        let call_option = Options::Call(Call {
            strike_price: 100.0,
            spot_price: 105.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
        });
        let price = black_scholes_price(call_option);
        assert!((price - 13.8579).abs() < 0.0001);
    }

    #[test]
    fn test_black_scholes_put() {
        let put_option = Options::Put(Put {
            strike_price: 100.0,
            spot_price: 95.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
        });
        let price = black_scholes_price(put_option);
        // Expected price calculated using Black-Scholes formula
        assert!((price - 7.6338).abs() < 0.0001);
    }

    #[test]
    fn test_at_the_money_call() {
        let atm_call = Options::Call(Call {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.25,
            risk_free_rate: 0.03,
            time_to_maturity: 0.5,
        });
        let price = black_scholes_price(atm_call);
        // ATM option should have significant time value
        assert!(price > 0.0);
        assert!((price - 7.7603).abs() < 0.0001);
    }

    #[test]
    fn test_at_the_money_put() {
        let atm_put = Options::Put(Put {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.25,
            risk_free_rate: 0.03,
            time_to_maturity: 0.5,
        });
        let price = black_scholes_price(atm_put);
        assert!(price > 0.0);
        assert!((price - 6.2715).abs() < 0.0001);
    }

    #[test]
    fn test_out_of_money_call() {
        let otm_call = Options::Call(Call {
            strike_price: 110.0,
            spot_price: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
        });
        let price = black_scholes_price(otm_call);
        // OTM call should be worth less than ATM
        assert!(price > 0.0);
        assert!((price - 6.0401).abs() < 0.0001);
    }

    #[test]
    fn test_out_of_money_put() {
        let otm_put = Options::Put(Put {
            strike_price: 90.0,
            spot_price: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
        });
        let price = black_scholes_price(otm_put);
        // OTM put should have only time value
        assert!(price > 0.0);
        assert!((price - 2.3101).abs() < 0.0001);
    }

    #[test]
    fn test_put_call_parity() {
        // Put-Call Parity: C - P = S - K*e^(-rt)
        let strike = 100.0;
        let spot = 105.0;
        let vol = 0.2;
        let rate = 0.05;
        let time = 1.0;

        let call_option = Options::Call(Call {
            strike_price: strike,
            spot_price: spot,
            volatility: vol,
            risk_free_rate: rate,
            time_to_maturity: time,
        });
        let call_price = black_scholes_price(call_option);

        let put_option = Options::Put(Put {
            strike_price: strike,
            spot_price: spot,
            volatility: vol,
            risk_free_rate: rate,
            time_to_maturity: time,
        });
        let put_price = black_scholes_price(put_option);

        let left_side = call_price - put_price;
        let right_side = spot - strike * (-rate * time).exp();

        // Put-call parity should hold within numerical precision
        assert!((left_side - right_side).abs() < 0.0001);
    }

    #[test]
    fn test_deep_in_money_call() {
        let deep_itm_call = Options::Call(Call {
            strike_price: 50.0,
            spot_price: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
        });
        let price = black_scholes_price(deep_itm_call);
        // Deep ITM call should be approximately S - K*e^(-rt)
        let intrinsic_value = 100.0 - 50.0 * (-0.05 * 1.0_f64).exp();
        assert!((price - intrinsic_value).abs() < 1.0);
        assert!((price - 52.4389).abs() < 0.0001);
    }

    #[test]
    fn test_low_volatility() {
        let low_vol_call = Options::Call(Call {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.01,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
        });
        let price = black_scholes_price(low_vol_call);
        // Low volatility ATM option should have less time value
        assert!(price > 0.0);
        assert!((price - 4.8771).abs() < 0.0001);
    }

    #[test]
    fn test_high_volatility() {
        let high_vol_call = Options::Call(Call {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.5,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
        });
        let price = black_scholes_price(high_vol_call);
        // High volatility should increase option value
        assert!(price > 10.0);
        assert!((price - 21.7926).abs() < 0.0001);
    }
}
