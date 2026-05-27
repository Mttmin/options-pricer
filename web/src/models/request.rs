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
    #[serde(default)]
    pub symbol: Option<String>,
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
    #[serde(default)]
    pub symbol: Option<String>,
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
    AsianOption,
    CliquetOption,
    BarrierOption,
    AutocallableNote,
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum BarrierKindInput {
    Call,
    Put,
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum BarrierTypeInput {
    UpAndIn,
    UpAndOut,
    DownAndIn,
    DownAndOut,
}

#[derive(Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum AsianPoolingInput {
    Average,
    Max,
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
pub struct AsianOptionRequest {
    pub option_type: OptionTypeInput,
    pub pooling: AsianPoolingInput,
    pub spot_price: f64,
    pub strike_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    #[serde(default = "default_asian_observations")]
    pub num_observations: u32,
    #[serde(default = "default_mc_simulations")]
    pub mc_simulations: u32,
}

fn default_asian_observations() -> u32 {
    50
}

#[derive(Deserialize)]
pub struct CliquetOptionRequest {
    pub option_type: OptionTypeInput,
    pub spot_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    #[serde(default = "default_cliquet_resets")]
    pub num_resets: u32,
    #[serde(default)]
    pub local_cap: Option<f64>,
    #[serde(default)]
    pub local_floor: Option<f64>,
    #[serde(default)]
    pub global_cap: Option<f64>,
    #[serde(default)]
    pub global_floor: Option<f64>,
    #[serde(default = "default_mc_simulations")]
    pub mc_simulations: u32,
}

fn default_cliquet_resets() -> u32 {
    12
}

#[derive(Deserialize)]
pub struct BarrierOptionRequest {
    pub option_type: BarrierKindInput,
    pub barrier_type: BarrierTypeInput,
    pub spot_price: f64,
    pub strike_price: f64,
    pub barrier: f64,
    #[serde(default = "default_rebate")]
    pub rebate: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    #[serde(default = "default_barrier_steps")]
    pub num_steps: u32,
    #[serde(default = "default_mc_simulations")]
    pub mc_simulations: u32,
}

#[derive(Deserialize)]
pub struct AutocallableNoteRequest {
    pub spot_price: f64,
    pub volatility: f64,
    pub risk_free_rate: f64,
    pub time_to_maturity: f64,
    pub dividend_yield: Option<f64>,
    #[serde(default = "default_notional")]
    pub notional: f64,
    pub autocall_barrier: f64,
    pub coupon_rate: f64,
    pub protection_barrier: f64,
    #[serde(default = "default_autocall_observations")]
    pub num_observations: u32,
    #[serde(default)]
    pub memory_coupon: bool,
    #[serde(default = "default_mc_simulations")]
    pub mc_simulations: u32,
}

fn default_rebate() -> f64 {
    0.0
}

fn default_barrier_steps() -> u32 {
    252
}

fn default_notional() -> f64 {
    1000.0
}

fn default_autocall_observations() -> u32 {
    4
}

#[derive(Deserialize)]
pub struct MarketQuery {
    #[serde(default = "default_lookback")]
    pub lookback_days: usize,
}

#[derive(Deserialize)]
pub struct CtmcPriceRequest {
    pub symbol: String,
    pub option_type: OptionTypeInput,
    pub direction: DirectionInput,
    pub strike_price: f64,
    pub time_to_maturity: f64,
    pub risk_free_rate: f64,
    #[serde(default)]
    pub dividend_yield: Option<f64>,
    #[serde(default)]
    pub n_restarts: Option<usize>,
    #[serde(default)]
    pub n_x: Option<usize>,
    #[serde(default)]
    pub m_v: Option<usize>,
    #[serde(default)]
    pub n_time: Option<usize>,
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
