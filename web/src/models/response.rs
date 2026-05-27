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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vix_fallback_used: Option<bool>,
    pub implied: Option<f64>,
}

#[derive(Serialize)]
pub struct PriceResponse {
    pub structure_type: String,
    pub pricing: PricingResult,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub greeks: Option<GreeksResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub payoff_curve: Option<PayoffCurve>,
}

#[derive(Serialize)]
pub struct PayoffCurve {
    pub spot_prices: Vec<f64>,
    pub payoffs: Vec<f64>,
}

#[derive(Serialize)]
pub struct PricingResult {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub black_scholes: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub monte_carlo: Option<MonteCarloResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binomial_european: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binomial_american: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bs_american_approx: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub penalty_solver: Option<PenaltySolverResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub timings: Option<PricingTimings>,
}

#[derive(Serialize)]
pub struct PricingTimings {
    pub total_ms: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub black_scholes_ms: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub monte_carlo_ms: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binomial_european_ms: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binomial_american_ms: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bs_american_approx_ms: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub penalty_solver_ms: Option<f64>,
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

#[derive(Serialize)]
pub struct PenaltySolverResult {
    pub price: f64,
    pub diagnostics: PenaltySolverDiagnostics,
}

#[derive(Serialize)]
pub struct HestonResult {
    pub kappa: f64,
    pub theta: f64,
    pub sigma: f64,
    pub rho: f64,
    pub v0: f64,
}

#[derive(Serialize)]
pub struct CtmcPriceResponse {
    pub symbol: String,
    pub spot_price: f64,
    pub price: f64,
    pub heston: HestonResult,
    pub ivrmse: f32,
    pub n_evals: usize,
    pub n_instruments: usize,
    pub timing_ms: f64,
}

#[derive(Serialize)]
pub struct PenaltySolverDiagnostics {
    pub total_iterations: usize,
    pub avg_iterations_per_step: f64,
    pub max_american_error: f64,
    pub spatial_nodes: usize,
    pub timesteps: usize,
    pub grid_min_price: f64,
    pub grid_max_price: f64,
}
