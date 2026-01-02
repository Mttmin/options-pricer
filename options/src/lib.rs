pub mod black_scholes;
pub mod exotics;
use black_scholes::*;

// Core option contract types shared across pricing engines and frontends.
#[derive(Debug, Clone, Copy)]
pub enum Options {
    Call(Call),
    Put(Put),
}

impl Options {
    pub fn bs_pricing(&self) -> f64 {
        black_scholes_price(*self)
    }
    pub fn new_call(strike_price: f64, spot_price: f64, volatility: f64, risk_free_rate: f64, time_to_maturity: f64) -> Self {
        Options::Call(Call::new(strike_price, spot_price, volatility, risk_free_rate, time_to_maturity))
    }
    pub fn new_put(strike_price: f64, spot_price: f64, volatility: f64, risk_free_rate: f64, time_to_maturity: f64) -> Self {
        Options::Put(Put::new(strike_price, spot_price, volatility, risk_free_rate, time_to_maturity))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Call {
    strike_price: f64,
    spot_price: f64,
    volatility: f64,
    risk_free_rate: f64,
    time_to_maturity: f64,
}

impl Call {
    pub fn bs_pricing (&self) -> f64 {
        black_scholes_price(Options::Call(*self))
    }
    pub fn new(strike_price: f64, spot_price: f64, volatility: f64, risk_free_rate: f64, time_to_maturity: f64) -> Self {
        Call {
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Put {
    strike_price: f64,
    spot_price: f64,
    volatility: f64,
    risk_free_rate: f64,
    time_to_maturity: f64,
}

impl Put {
    pub fn bs_pricing (&self) -> f64 {
        black_scholes_price(Options::Put(*self))
    }
    pub fn new(strike_price: f64, spot_price: f64, volatility: f64, risk_free_rate: f64, time_to_maturity: f64) -> Self {
        Put {
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
        }
    }
}
