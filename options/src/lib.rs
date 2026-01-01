// Core option contract types shared across pricing engines and frontends.
#[derive(Debug, Clone, Copy)]
pub enum Options {
    Call(Call),
    Put(Put),
}

#[derive(Debug, Clone, Copy)]
pub struct Call {
    pub strike_price: f64,
    pub spot_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Put {
    pub strike_price: f64,
    pub spot_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
}
