use crate::Options;
use statrs::distribution::{Continuous, ContinuousCDF, Normal};

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

/// Approximate the price of an American option using Black's approximation method
/// 
/// If a dividend date is provided, the option price is calculated just before the last dividend payment and compare it to exercising at maturity
/// If a stock doesn't pay dividends, the option should in theory not be exercised early, so the European price is returned.
pub fn black_scholes_approx_american(option: Options, dividend_date: Option<f64>) -> f64 {
    // we use the Black approximation for American options
    let european_price = black_scholes_price(option);
    if let Some(div_date) = dividend_date {
        // If there's a dividend before maturity, calculate the black scholes price just before the dividend date
        let adjusted_option = match option {
            Options::Call(call) => Options::Call(crate::Call {
                strike_price: call.strike_price,
                spot_price: call.spot_price - (call.spot_price * 0.02), // assuming a 2% dividend for simplicity
                volatility: call.volatility,
                risk_free_rate: call.risk_free_rate,
                time_to_maturity: div_date,
                dividend_yield: call.dividend_yield,
            }),
            Options::Put(put) => Options::Put(crate::Put {
                strike_price: put.strike_price,
                spot_price: put.spot_price - (put.spot_price * 0.02), // assuming a 2% dividend for simplicity
                volatility: put.volatility,
                risk_free_rate: put.risk_free_rate,
                time_to_maturity: div_date,
                dividend_yield: put.dividend_yield,
            }),
        };
        let pre_div_price = black_scholes_price(adjusted_option);
        if pre_div_price > european_price {
            pre_div_price
        } else {
            european_price
        }
    } else {
        european_price
    }
}

/// Calculate Implied Volatility using Newton-Raphson method
pub fn calculate_iv(
    price: f64,
    s: f64, // spot price
    k: f64, // strike price
    t: f64, // time to maturity (years)
    r: f64, // risk-free rate
    q: f64, // dividend yield
    is_call: bool,
) -> f64 {
    if price <= 0.0 || t <= 0.0 || s <= 0.0 || k <= 0.0 {
        return 0.0;
    }

    let mut sigma = 0.5; // Initial guess
    let max_iterations = 50;
    let tolerance = 1e-4;
    let std_norm = Normal::new(0.0, 1.0).unwrap();

    for _ in 0..max_iterations {
        let d1 = d_plus(t, r, Some(q), sigma, s, k);
        let d2 = d_minus(t, r, Some(q), sigma, s, k);

        let (price_implied, vega) = if is_call {
            let nd1 = std_norm.cdf(d1);
            let nd2 = std_norm.cdf(d2);
            let p = s * (-q * t).exp() * nd1 - k * (-r * t).exp() * nd2;
            let v = s * (-q * t).exp() * std_norm.pdf(d1) * t.sqrt();
            (p, v)
        } else {
             let n_neg_d1 = std_norm.cdf(-d1);
             let n_neg_d2 = std_norm.cdf(-d2);
             let p = k * (-r * t).exp() * n_neg_d2 - s * (-q * t).exp() * n_neg_d1;
             let v = s * (-q * t).exp() * std_norm.pdf(d1) * t.sqrt();
             (p, v)
        };

        let diff = price_implied - price;

        if diff.abs() < tolerance {
            return sigma;
        }

        if vega.abs() < 1e-8 {
            break;
        }

        sigma = sigma - diff / vega;

        if sigma <= 0.0 {
            sigma = 1e-5;
        }
    }

    sigma
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Call, Greeks, Options, Put};
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

    #[test]
    fn test_calculate_iv_convergence() {
        let spot = 100.0;
        let strike = 100.0;
        let time = 1.0;
        let rate = 0.05;
        let div = 0.0;
        let target_vol = 0.25;

        // Calculate price using Black-Scholes
        let call = Options::new_call(strike, spot, target_vol, rate, time, Some(div));
        let price = crate::black_scholes::black_scholes_price(call);

        // Calculate IV from price
        let implied_vol = calculate_iv(price, spot, strike, time, rate, div, true);

        assert!((implied_vol - target_vol).abs() < 1e-4, "IV should converge to input vol. Got: {}, Expected: {}", implied_vol, target_vol);
    }
}
