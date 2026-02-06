use serde::Serialize;

#[derive(Serialize)]
pub struct MarketDataResponse {
    pub symbol: String,
    pub spot_price: f64,
    pub historical_volatility: f64,
    pub implied_volatility: Option<f64>,
    pub corporate_action_detected: bool,
    pub dividend_yield: Option<f64>,
}

#[derive(Serialize)]
pub struct VolatilityResponse {
    pub symbol: String,
    pub historical: f64,
    pub ema: Option<f64>,
    pub vix_correlated: Option<f64>,
    pub implied: Option<f64>,
}

#[derive(Serialize)]
pub struct PriceResponse {
    pub structure_type: String,
    pub pricing: PricingResult,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub greeks: Option<GreeksResult>,
}

#[derive(Serialize)]
pub struct PricingResult {
    pub black_scholes: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub monte_carlo: Option<MonteCarloResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binomial_european: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binomial_american: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bs_american_approx: Option<f64>,
}

#[derive(Serialize)]
pub struct MonteCarloResult {
    pub price: f64,
    pub std_error: f64,
    pub ci_lower: f64,
    pub ci_upper: f64,
}

#[derive(Serialize)]
pub struct GreeksResult {
    pub delta: f64,
    pub gamma: f64,
    pub theta: f64,
    pub vega: f64,
    pub rho: f64,
}
