use serde::Deserialize;

#[derive(Deserialize)]
#[serde(tag = "structure_type")]
pub enum PriceRequest {
    #[serde(rename = "single")]
    Single(SingleOptionRequest),
    #[serde(rename = "spread")]
    Spread(SpreadRequest),
    #[serde(rename = "exotic")]
    Exotic(ExoticRequest),
}

#[derive(Deserialize)]
pub struct SingleOptionRequest {
    pub option_type: OptionTypeInput,
    pub direction: DirectionInput,
    pub strike_price: f64,
    pub spot_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    #[serde(default = "default_mc_simulations")]
    pub mc_simulations: u32,
    #[serde(default = "default_binomial_depth")]
    pub binomial_depth: u16,
}

#[derive(Deserialize)]
pub struct SpreadRequest {
    pub spread_type: SpreadTypeInput,
    pub direction: DirectionInput,
    pub spot_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    pub strikes: SpreadStrikes,
    #[serde(default)]
    pub short_term_maturity: Option<f64>,
    #[serde(default)]
    pub long_term_maturity: Option<f64>,
    #[serde(default = "default_mc_simulations")]
    pub mc_simulations: u32,
}

#[derive(Deserialize)]
pub struct SpreadStrikes {
    pub strike: Option<f64>,
    pub strike_call: Option<f64>,
    pub strike_put: Option<f64>,
    pub strike_low: Option<f64>,
    pub strike_high: Option<f64>,
    pub strike_medium: Option<f64>,
}

#[derive(Deserialize)]
pub struct ExoticRequest {
    pub exotic_type: ExoticTypeInput,
    #[serde(flatten)]
    pub params: serde_json::Value,
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum OptionTypeInput {
    Call,
    Put,
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum DirectionInput {
    Long,
    Short,
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum SpreadTypeInput {
    Straddle,
    Strangle,
    Strip,
    Strap,
    SyntheticStock,
    BullSpreadCall,
    BullSpreadPut,
    BearSpreadCall,
    BearSpreadPut,
    Butterfly,
    CalendarSpread,
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum ExoticTypeInput {
    ConvertibleBond,
    ChooserOption,
}

#[derive(Deserialize)]
pub struct ConvertibleBondRequest {
    pub face_value: f64,
    pub coupon_rate: f64,
    pub maturity: f64,
    pub payment_frequency: u32,
    pub credit_spread: f64,
    pub conversion_price: f64,
    pub stock_price: f64,
    pub volatility: f64,
    pub time_to_maturity: f64,
    pub risk_free_rate: f64,
    pub dividend_yield: Option<f64>,
}

#[derive(Deserialize)]
pub struct ChooserOptionRequest {
    pub spot_price: f64,
    pub strike_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    pub choose_time: f64,
}

#[derive(Deserialize)]
pub struct MarketQuery {
    #[serde(default = "default_lookback")]
    pub lookback_days: usize,
}

fn default_mc_simulations() -> u32 {
    100_000
}

fn default_binomial_depth() -> u16 {
    1000
}

fn default_lookback() -> usize {
    180
}
