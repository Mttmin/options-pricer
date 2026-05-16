export type OptionType = "call" | "put";

export type Direction = "long" | "short";

export type ExerciseStyle = "european" | "american";

export type SpreadType =
  | "straddle"
  | "strangle"
  | "strip"
  | "strap"
  | "synthetic_stock"
  | "bull_spread_call"
  | "bull_spread_put"
  | "bear_spread_call"
  | "bear_spread_put"
  | "butterfly"
  | "calendar_spread";

export type ExoticType = "convertible_bond" | "chooser_option" | "asian_option" | "cliquet_option";

export type AsianPooling = "average" | "max";

export type StructureType = "single" | "spread" | "exotic";

export type SpreadStrikes = {
  strike?: number;
  strike_call?: number;
  strike_put?: number;
  strike_low?: number;
  strike_high?: number;
  strike_medium?: number;
};

export type SingleOptionRequest = {
  structure_type: "single";
  option_type: OptionType;
  direction: Direction;
  strike_price: number;
  spot_price: number;
  volatility: number;
  risk_free_rate: number;
  time_to_maturity: number;
  dividend_yield: number | null;
  mc_simulations?: number;
  binomial_depth?: number;
};

export type SpreadRequest = {
  structure_type: "spread";
  spread_type: SpreadType;
  direction: Direction;
  spot_price: number;
  volatility: number;
  risk_free_rate: number;
  time_to_maturity: number;
  dividend_yield: number | null;
  strikes: SpreadStrikes;
  short_term_maturity?: number;
  long_term_maturity?: number;
  mc_simulations?: number;
};

export type ConvertibleBondRequest = {
  structure_type: "exotic";
  exotic_type: "convertible_bond";
  face_value: number;
  coupon_rate: number;
  maturity: number;
  payment_frequency: number;
  credit_spread: number;
  conversion_price: number;
  stock_price: number;
  volatility: number;
  time_to_maturity: number;
  risk_free_rate: number;
  dividend_yield: number | null;
};

export type ChooserOptionRequest = {
  structure_type: "exotic";
  exotic_type: "chooser_option";
  spot_price: number;
  strike_price: number;
  volatility: number;
  risk_free_rate: number;
  time_to_maturity: number;
  dividend_yield: number | null;
  choose_time: number;
};

export type AsianOptionRequest = {
  structure_type: "exotic";
  exotic_type: "asian_option";
  option_type: OptionType;
  pooling: AsianPooling;
  spot_price: number;
  strike_price: number;
  volatility: number;
  risk_free_rate: number;
  time_to_maturity: number;
  dividend_yield: number | null;
  num_observations: number;
  mc_simulations?: number;
};

export type CliquetOptionRequest = {
  structure_type: "exotic";
  exotic_type: "cliquet_option";
  option_type: OptionType;
  spot_price: number;
  volatility: number;
  risk_free_rate: number;
  time_to_maturity: number;
  dividend_yield: number | null;
  num_resets: number;
  local_cap: number | null;
  local_floor: number | null;
  global_cap: number | null;
  global_floor: number | null;
  mc_simulations?: number;
};

export type PriceRequest =
  | SingleOptionRequest
  | SpreadRequest
  | ConvertibleBondRequest
  | ChooserOptionRequest
  | AsianOptionRequest
  | CliquetOptionRequest;

export type MonteCarloResult = {
  price: number;
  std_error: number;
  ci_lower: number;
  ci_upper: number;
};

export type PenaltySolverDiagnostics = {
  total_iterations: number;
  avg_iterations_per_step: number;
  max_american_error: number;
  spatial_nodes: number;
  timesteps: number;
  grid_min_price: number;
  grid_max_price: number;
};

export type PenaltySolverResult = {
  price: number;
  diagnostics: PenaltySolverDiagnostics;
};

export type PricingResult = {
  black_scholes: number;
  monte_carlo?: MonteCarloResult;
  binomial_european?: number;
  binomial_american?: number;
  bs_american_approx?: number;
  penalty_solver?: PenaltySolverResult;
  timings?: PricingTimings;
};

export type PricingTimings = {
  total_ms: number;
  black_scholes_ms?: number | null;
  monte_carlo_ms?: number | null;
  binomial_european_ms?: number | null;
  binomial_american_ms?: number | null;
  bs_american_approx_ms?: number | null;
  penalty_solver_ms?: number | null;
};

export type GreeksResult = {
  delta: number;
  gamma: number;
  theta: number;
  vega: number;
  rho: number;
};

export type PayoffCurve = {
  spot_prices: number[];
  payoffs: number[];
};

export type PriceResponse = {
  structure_type: StructureType;
  pricing: PricingResult;
  greeks?: GreeksResult;
  payoff_curve?: PayoffCurve;
};

export type MarketDataResponse = {
  symbol: string;
  spot_price: number;
  historical_volatility: number;
  implied_volatility: number | null;
  corporate_action_detected: boolean;
  dividend_yield: number | null;
};

export type SofrResponse = {
  rate: number;
  date: string;
};

export type VolatilityResponse = {
  symbol: string;
  historical: number;
  ema: number | null;
  vix_correlated: number | null;
  vix_fallback_used?: boolean | null;
  implied: number | null;
};

export type VolSource = "implied" | "historical" | "ema" | "vix_correlated";

export type SmilePoint = {
  strike: number;
  implied_volatility: number;
  option_type: "Call" | "Put";
  volume: number;
  open_interest: number;
  moneyness: number;
};

export type IVSmileData = {
  symbol: string;
  underlying_price: number;
  data_date: string;
  smiles_by_expiry: Record<string, SmilePoint[]>;
};

export type CtmcPriceRequest = {
  symbol: string;
  option_type: OptionType;
  direction: Direction;
  strike_price: number;
  time_to_maturity: number;
  risk_free_rate: number;
  dividend_yield?: number | null;
  n_restarts?: number;
  n_x?: number;
  m_v?: number;
  n_time?: number;
};

export type HestonResult = {
  kappa: number;
  theta: number;
  sigma: number;
  rho: number;
  v0: number;
};

export type CtmcPriceResponse = {
  symbol: string;
  spot_price: number;
  price: number;
  heston: HestonResult;
  ivrmse: number;
  n_evals: number;
  n_instruments: number;
  timing_ms: number;
};

export const SPREAD_LABELS: Record<SpreadType, string> = {
  straddle: "Straddle",
  strangle: "Strangle",
  strip: "Strip",
  strap: "Strap",
  synthetic_stock: "Synthetic Stock",
  bull_spread_call: "Bull Spread (Call)",
  bull_spread_put: "Bull Spread (Put)",
  bear_spread_call: "Bear Spread (Call)",
  bear_spread_put: "Bear Spread (Put)",
  butterfly: "Butterfly",
  calendar_spread: "Calendar Spread",
};
