pub mod black_scholes;
pub mod exotics;

use black_scholes::*;
use statrs::distribution::{Continuous, ContinuousCDF, Normal};

// Core option contract types shared across pricing engines and front-ends.
#[derive(Debug, Clone, Copy)]
pub enum Options {
    Call(Call),
    Put(Put),
}

impl Options {
    pub fn bs_pricing(&self) -> f64 {
        black_scholes_price(*self)
    }
    pub fn new_call(
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> Self {
        Options::Call(Call::new(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        ))
    }
    pub fn new_put(
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> Self {
        Options::Put(Put::new(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        ))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Call {
    strike_price: f64,
    spot_price: f64,
    volatility: f64,
    risk_free_rate: f64,
    time_to_maturity: f64,
    dividend_yield: Option<f64>,
}

impl Call {
    pub fn bs_pricing(&self) -> f64 {
        black_scholes_price(Options::Call(*self))
    }
    pub fn new(
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> Self {
        Call {
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        }
    }
    pub fn payout(&self, spot: f64) -> f64 {
        f64::max(0.0, spot - self.strike_price)
    }
    /// Calculates Delta (Δ) - the rate of change of option price with respect to spot price.
    /// For calls, delta ranges from 0 to 1. Higher values indicate deeper in-the-money positions
    /// 
    /// Formula: Δ = e^(-qT) * N(d₁)
    pub fn delta(&self, imply_vol: f64, spot_price: f64) -> f64 {
        let d1 = d_plus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let dividend_correction = self
            .dividend_yield
            .map_or(1.0, |yield_val| (-yield_val * self.time_to_maturity).exp());
        std_norm.cdf(d1) * dividend_correction
    }
    /// Calculates Theta (Θ) - time decay of option value.
    /// Typically negative for calls. Divide by 365 for daily theta
    /// 
    /// Formula: Θ = -[S*N'(d₁)*σ*e^(-qT)] / [2√T] + qS*N(d₁)*e^(-qT) - rK*e^(-rT)*N(d₂)
    pub fn theta(&self, imply_vol: f64, spot_price: f64) -> f64 {
        let d1 = d_plus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        let d2 = d_minus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let n_d1 = std_norm.pdf(d1);
        let n_d2 = std_norm.cdf(d2);
        let dividend_correction = self
            .dividend_yield
            .map_or(1.0, |yield_val| (-yield_val * self.time_to_maturity).exp());
        let dividend_npv = self.dividend_yield.map_or(0.0, |yield_val| {
            yield_val * self.spot_price * dividend_correction * std_norm.cdf(d1)
        });
        -(spot_price * n_d1 * imply_vol * dividend_correction)
            / (2.0 * self.time_to_maturity.sqrt())
            + dividend_npv
            - self.risk_free_rate
                * self.strike_price
                * (-self.risk_free_rate * self.time_to_maturity).exp()
                * n_d2
    }
    /// Calculates Gamma (Γ) - the rate of change of delta with respect to spot price.
    /// Always positive. Highest for at-the-money options near expiration
    /// 
    /// Formula: Γ = N'(d₁) * e^(-qT) / (S * σ * √T)
    pub fn gamma(&self, imply_vol: f64, spot_price: f64) -> f64 {
        let d1 = d_plus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let dividend_correction = self
            .dividend_yield
            .map_or(1.0, |yield_val| (-yield_val * self.time_to_maturity).exp());
        std_norm.pdf(d1) * dividend_correction
            / (spot_price * imply_vol * self.time_to_maturity.sqrt())
    }
    /// Calculates Vega (ν) - sensitivity to volatility changes.
    /// Always positive. Highest for at-the-money options
    /// 
    /// Formula: ν = S * N'(d₁) * √T * e^(-qT)
    pub fn vega(&self, spot_price: f64) -> f64 {
        let d1 = d_plus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            self.volatility,
            spot_price,
            self.strike_price,
        );
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let dividend_correction = self
            .dividend_yield
            .map_or(1.0, |yield_val| (-yield_val * self.time_to_maturity).exp());
        spot_price * std_norm.pdf(d1) * self.time_to_maturity.sqrt() * dividend_correction
    }
    /// Calculates Rho (ρ) - sensitivity to interest rate changes.
    /// Positive for calls. Larger for longer-dated and in-the-money options
    /// 
    /// Formula: ρ = K * T * e^(-rT) * N(d₂)
    pub fn rho(&self, imply_vol: f64, spot_price: f64, interest_rate: f64) -> f64 {
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let d2 = d_minus(
            self.time_to_maturity,
            interest_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        self.strike_price
            * self.time_to_maturity
            * std_norm.cdf(d2)
            * (-interest_rate * self.time_to_maturity).exp()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Put {
    strike_price: f64,
    spot_price: f64,
    volatility: f64,
    risk_free_rate: f64,
    time_to_maturity: f64,
    dividend_yield: Option<f64>,
}

impl Put {
    pub fn bs_pricing(&self) -> f64 {
        black_scholes_price(Options::Put(*self))
    }
    pub fn new(
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> Self {
        Put {
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        }
    }
    pub fn payout(&self, spot: f64) -> f64 {
        f64::max(0.0, self.strike_price - spot)
    }
    /// Calculates Delta (Δ) - the rate of change of option price with respect to spot price.
    /// For puts, delta ranges from -1 to 0. More negative values indicate deeper in-the-money positions
    /// 
    /// Formula: Δ = e^(-qT) * [N(d₁) - 1]
    pub fn delta(&self, imply_vol: f64, spot_price: f64) -> f64 {
        let d1 = d_plus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let dividend_correction = self
            .dividend_yield
            .map_or(1.0, |yield_val| (-yield_val * self.time_to_maturity).exp());
        (std_norm.cdf(d1) - 1.0) * dividend_correction
    }
    /// Calculates Theta (Θ) - time decay of option value.
    /// Can be positive or negative for puts depending on moneyness. Divide by 365 for daily theta
    /// 
    /// Formula: Θ = -[S*N'(d₁)*σ*e^(-qT)] / [2√T] - qS*N(-d₁)*e^(-qT) + rK*e^(-rT)*N(-d₂)
    pub fn theta(&self, imply_vol: f64, spot_price: f64) -> f64 {
        let d1 = d_plus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        let d2 = d_minus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let n_d1 = std_norm.pdf(d1);
        let n_d2 = std_norm.cdf(-d2);
        let dividend_correction = self
            .dividend_yield
            .map_or(1.0, |yield_val| (-yield_val * self.time_to_maturity).exp());
        let dividend_npv = self.dividend_yield.map_or(0.0, |yield_val| {
            yield_val * self.spot_price * dividend_correction * std_norm.cdf(-d1)
        });
        -(spot_price * n_d1 * imply_vol * dividend_correction)
            / (2.0 * self.time_to_maturity.sqrt())
            - dividend_npv
            + self.risk_free_rate
                * self.strike_price
                * (-self.risk_free_rate * self.time_to_maturity).exp()
                * n_d2
    }
    /// Calculates Gamma (Γ) - the rate of change of delta with respect to spot price.
    /// Always positive. Highest for at-the-money options near expiration
    /// 
    /// Formula: Γ = N'(d₁) * e^(-qT) / (S * σ * √T)
    pub fn gamma(&self, imply_vol: f64, spot_price: f64) -> f64 {
        let d1 = d_plus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let dividend_correction = self
            .dividend_yield
            .map_or(1.0, |yield_val| (-yield_val * self.time_to_maturity).exp());
        (std_norm.pdf(d1) * dividend_correction) / (spot_price * imply_vol * self.time_to_maturity.sqrt())
    }
    /// Calculates Vega (ν) - sensitivity to volatility changes.
    /// Always positive. Highest for at-the-money options
    /// 
    /// Formula: ν = S * N'(d₁) * √T * e^(-qT)
    pub fn vega(&self, spot_price: f64) -> f64 {
        let d1 = d_plus(
            self.time_to_maturity,
            self.risk_free_rate,
            self.dividend_yield,
            self.volatility,
            spot_price,
            self.strike_price,
        );
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let dividend_correction = self
            .dividend_yield
            .map_or(1.0, |yield_val| (-yield_val * self.time_to_maturity).exp());
        spot_price * std_norm.pdf(d1) * self.time_to_maturity.sqrt() * dividend_correction
    }
    /// Calculates Rho (ρ) - sensitivity to interest rate changes.
    /// Negative for puts. Larger absolute value for longer-dated and in-the-money options
    /// 
    /// Formula: ρ = -K * T * e^(-rT) * N(-d₂)
    pub fn rho(&self, imply_vol: f64, spot_price: f64, interest_rate: f64) -> f64 {
        let std_norm = Normal::new(0.0, 1.0).unwrap();
        let d2 = d_minus(
            self.time_to_maturity,
            interest_rate,
            self.dividend_yield,
            imply_vol,
            spot_price,
            self.strike_price,
        );
        -self.strike_price
            * self.time_to_maturity
            * std_norm.cdf(-d2)
            * (-interest_rate * self.time_to_maturity).exp()
    }
}
