use crate::{Options};
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
    use crate::{Options, Call, Put};
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
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
        });
        let price = black_scholes_price(atm_call);
        assert!(price > 0.0);
    }
}
