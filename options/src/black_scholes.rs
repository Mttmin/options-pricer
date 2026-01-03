use crate::Options;
use statrs::distribution::{ContinuousCDF, Normal};

// calculate d1 for the Black-Scholes formula
pub fn d_plus(t: f64, r: f64, q: Option<f64>, sigma: f64, spot: f64, strike: f64) -> f64 {
    let numerator = (spot / strike).ln() + (r - q.unwrap_or(0.0) + 0.5 * sigma * sigma) * t;
    let denominator = sigma * t.sqrt();
    numerator / denominator
}

// calculate d2 for the Black-Scholes formula
pub fn d_minus(t: f64, r: f64, q: Option<f64>, sigma: f64, spot: f64, strike: f64) -> f64 {
    let numerator = (spot / strike).ln() + (r - q.unwrap_or(0.0) - 0.5 * sigma * sigma) * t;
    let denominator = sigma * t.sqrt();
    numerator / denominator
}
/// Calculate the Black-Scholes price for a given option, either Call or Put.
///
/// Needs the option parameters encapsulated in the Options enum, handles dividend yields if they are present
pub fn black_scholes_price(option: Options) -> f64 {
    let std_norm = Normal::new(0.0, 1.0).unwrap();
    if let Options::Call(call) = option {
        let d1 = d_plus(
            call.time_to_maturity,
            call.risk_free_rate,
            call.dividend_yield,
            call.volatility,
            call.spot_price,
            call.strike_price,
        );
        let d2 = d_minus(
            call.time_to_maturity,
            call.risk_free_rate,
            call.dividend_yield,
            call.volatility,
            call.spot_price,
            call.strike_price,
        );
        let nd1: f64 = std_norm.cdf(d1);
        let nd2: f64 = std_norm.cdf(d2);
        call.spot_price * (-call.dividend_yield.unwrap_or(0.0) * call.time_to_maturity).exp() * nd1
            - call.strike_price * (-call.risk_free_rate * call.time_to_maturity).exp() * nd2
    } else if let Options::Put(put) = option {
        let d1 = d_plus(
            put.time_to_maturity,
            put.risk_free_rate,
            put.dividend_yield,
            put.volatility,
            put.spot_price,
            put.strike_price,
        );
        let d2 = d_minus(
            put.time_to_maturity,
            put.risk_free_rate,
            put.dividend_yield,
            put.volatility,
            put.spot_price,
            put.strike_price,
        );
        let n_neg_d1 = std_norm.cdf(-d1);
        let n_neg_d2 = std_norm.cdf(-d2);
        put.strike_price * (-put.risk_free_rate * put.time_to_maturity).exp() * n_neg_d2
            - put.spot_price
                * (-put.dividend_yield.unwrap_or(0.0) * put.time_to_maturity).exp()
                * n_neg_d1
    } else {
        0.0
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Call, Options, Put};
    #[test]
    fn test_black_scholes_pricing() {
        // Setup test options - ITM Call (spot > strike)
        let call_itm = Options::Call(Call {
            strike_price: 100.0,
            spot_price: 105.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: None,
        });

        // ITM Put (spot < strike)
        let put_itm = Options::Put(Put {
            strike_price: 100.0,
            spot_price: 95.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: None,
        });

        // ATM Call
        let call_atm = Options::Call(Call {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: None,
        });

        // Call with dividend
        let call_div = Options::Call(Call {
            strike_price: 100.0,
            spot_price: 105.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: Some(0.02),
        });

        // Put with dividend
        let put_div = Options::Put(Put {
            strike_price: 100.0,
            spot_price: 95.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: Some(0.02),
        });

        // Calculate all prices
        let call_itm_price = black_scholes_price(call_itm);
        let put_itm_price = black_scholes_price(put_itm);
        let call_atm_price = black_scholes_price(call_atm);
        let call_div_price = black_scholes_price(call_div);
        let put_div_price = black_scholes_price(put_div);

        // Assert all prices
        assert!(
            (call_itm_price - 13.8579).abs() < 0.0001,
            "ITM Call price incorrect"
        );
        assert!(
            (put_itm_price - 7.6338).abs() < 0.0001,
            "ITM Put price incorrect"
        );
        assert!(call_atm_price > 0.0, "ATM Call price should be positive");
        assert!(
            (call_div_price - 12.3885).abs() < 0.001,
            "Call with dividend price incorrect"
        );
        assert!(
            (put_div_price - 8.5416).abs() < 0.001,
            "Put with dividend price incorrect"
        );
    }

    #[test]
    fn test_greeks() {
        // Setup test options
        let call = Call {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: None,
        };

        let put = Put {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: None,
        };

        // Calculate all Greeks
        let call_delta = call.delta(0.2, 100.0);
        let call_gamma = call.gamma(0.2, 100.0);
        let call_vega = call.vega(100.0);
        let call_theta = call.theta(0.2, 100.0);
        let call_rho = call.rho(0.2, 100.0, 0.05);

        let put_delta = put.delta(0.2, 100.0);
        let put_gamma = put.gamma(0.2, 100.0);
        let put_vega = put.vega(100.0);
        let put_theta = put.theta(0.2, 100.0);
        let put_rho = put.rho(0.2, 100.0, 0.05);

        // Assert all Greeks
        assert!(
            call_delta > 0.0 && call_delta < 1.0,
            "Call delta should be between 0 and 1"
        );
        assert!(
            (call_delta - 0.6368).abs() < 0.001,
            "Call delta value incorrect"
        );

        assert!(call_gamma > 0.0, "Call gamma should be positive");
        assert!(
            (call_gamma - 0.01876).abs() < 0.00001,
            "Call gamma value incorrect"
        );

        assert!(call_vega > 0.0, "Call vega should be positive");
        assert!(
            (call_vega - 37.524).abs() < 0.001,
            "Call vega value incorrect"
        );

        assert!(call_theta < 0.0, "Call theta should be negative");

        assert!(call_rho > 0.0, "Call rho should be positive");

        assert!(
            (put_delta + 0.3632).abs() < 0.001,
            "Put delta value incorrect"
        );

        assert!(put_gamma > 0.0, "Put gamma should be positive");
        assert!(
            (put_gamma - 0.01876).abs() < 0.00001,
            "Put gamma value incorrect"
        );

        assert!(put_vega > 0.0, "Put vega should be positive");
        assert!(
            (put_vega - 37.524).abs() < 0.001,
            "Put vega value incorrect"
        );

        assert!(put_theta < 0.0, "Put theta should be negative");

        assert!(put_rho < 0.0, "Put rho should be negative");
    }

    #[test]
    fn test_greeks_with_dividend() {
        let call = Call {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: Some(0.03),
        };
        let delta = call.delta(0.2, 100.0);
        println!("Call delta with dividend: {}", delta);
        assert!(delta > 0.0 && delta < 1.0, "Delta with div incorrect");

        let put = Put {
            strike_price: 100.0,
            spot_price: 100.0,
            volatility: 0.2,
            risk_free_rate: 0.05,
            time_to_maturity: 1.0,
            dividend_yield: Some(0.03),
        };
        let gamma = put.gamma(0.2, 100.0);
        assert!(gamma > 0.0, "Gamma with div incorrect");
        assert!(
            (gamma - 0.01897).abs() < 0.00001,
            "Gamma with div incorrect"
        );
    }
}
