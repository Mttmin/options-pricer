use axum::Json;
use options::black_scholes::{black_scholes_approx_american, black_scholes_price};
use options::exotics::{ChooserOption, ConvertibleBond};
use options::optionspreads::{Direction, OptionSpreads};
use options::{BlackScholesPrice, Call, Greeks, Options, Put};
use monte_carlo::binomial::{binomial_price, BinomialParameters, ExerciseStyle};
use monte_carlo::{monte_carlo, price_european};

use crate::error::AppError;
use crate::models::request::*;
use crate::models::response::*;

pub async fn calculate_price(
    Json(request): Json<PriceRequest>,
) -> Result<Json<PriceResponse>, AppError> {
    let response = tokio::task::spawn_blocking(move || match request {
        PriceRequest::Single(req) => price_single(req),
        PriceRequest::Spread(req) => price_spread(req),
        PriceRequest::Exotic(req) => price_exotic(req),
    })
    .await
    .map_err(|e| AppError::Internal(format!("Pricing task failed: {}", e)))??;

    Ok(Json(response))
}

fn price_single(req: SingleOptionRequest) -> Result<PriceResponse, AppError> {
    let direction_sign = match req.direction {
        DirectionInput::Long => 1.0,
        DirectionInput::Short => -1.0,
    };

    let option = match req.option_type {
        OptionTypeInput::Call => Options::new_call(
            req.strike_price, req.spot_price, req.volatility,
            req.risk_free_rate, req.time_to_maturity, req.dividend_yield,
        ),
        OptionTypeInput::Put => Options::new_put(
            req.strike_price, req.spot_price, req.volatility,
            req.risk_free_rate, req.time_to_maturity, req.dividend_yield,
        ),
    };

    // Black-Scholes
    let bs_price = black_scholes_price(option) * direction_sign;

    // Monte Carlo (monomorphized for Call/Put)
    let mc_result = match req.option_type {
        OptionTypeInput::Call => price_european::<Call>(&option, req.mc_simulations),
        OptionTypeInput::Put => price_european::<Put>(&option, req.mc_simulations),
    };
    let (ci_lower, ci_upper) = mc_result.confidence_interval_95();

    // Binomial tree
    let params = BinomialParameters::new(
        req.spot_price, req.binomial_depth, req.time_to_maturity,
        req.volatility, req.risk_free_rate,
    );

    let binomial_euro = match req.option_type {
        OptionTypeInput::Call => {
            let call = Call::new(req.strike_price, req.spot_price, req.volatility,
                req.risk_free_rate, req.time_to_maturity, req.dividend_yield);
            binomial_price(&call, &params, ExerciseStyle::European)
        }
        OptionTypeInput::Put => {
            let put = Put::new(req.strike_price, req.spot_price, req.volatility,
                req.risk_free_rate, req.time_to_maturity, req.dividend_yield);
            binomial_price(&put, &params, ExerciseStyle::European)
        }
    };

    let binomial_amer = match req.option_type {
        OptionTypeInput::Call => {
            let call = Call::new(req.strike_price, req.spot_price, req.volatility,
                req.risk_free_rate, req.time_to_maturity, req.dividend_yield);
            binomial_price(&call, &params, ExerciseStyle::American)
        }
        OptionTypeInput::Put => {
            let put = Put::new(req.strike_price, req.spot_price, req.volatility,
                req.risk_free_rate, req.time_to_maturity, req.dividend_yield);
            binomial_price(&put, &params, ExerciseStyle::American)
        }
    };

    // Black's American approximation
    let bs_american = black_scholes_approx_american(option, None) * direction_sign;

    // Greeks
    let vol = req.volatility;
    let spot = req.spot_price;
    let rfr = req.risk_free_rate;
    let greeks = GreeksResult {
        delta: option.delta(vol, spot) * direction_sign,
        gamma: option.gamma(vol, spot) * direction_sign,
        theta: option.theta(vol, spot) * direction_sign,
        vega: option.vega(spot) * direction_sign,
        rho: option.rho(vol, spot, rfr) * direction_sign,
    };

    Ok(PriceResponse {
        structure_type: "single".to_string(),
        pricing: PricingResult {
            black_scholes: bs_price,
            monte_carlo: Some(MonteCarloResult {
                price: mc_result.price * direction_sign,
                std_error: mc_result.std_error,
                ci_lower: ci_lower * direction_sign,
                ci_upper: ci_upper * direction_sign,
            }),
            binomial_european: Some(binomial_euro * direction_sign),
            binomial_american: Some(binomial_amer * direction_sign),
            bs_american_approx: Some(bs_american),
        },
        greeks: Some(greeks),
    })
}

fn to_direction(d: DirectionInput) -> Direction {
    match d {
        DirectionInput::Long => Direction::LONG,
        DirectionInput::Short => Direction::SHORT,
    }
}

fn price_spread(req: SpreadRequest) -> Result<PriceResponse, AppError> {
    let direction = to_direction(req.direction);
    let s = &req.strikes;
    let spot = req.spot_price;
    let vol = req.volatility;
    let rfr = req.risk_free_rate;
    let ttm = req.time_to_maturity;
    let div = req.dividend_yield;

    let spread = match req.spread_type {
        SpreadTypeInput::Straddle => {
            let strike = s.strike.ok_or(AppError::BadRequest("straddle requires 'strike'".into()))?;
            OptionSpreads::new_straddle(direction, strike, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::Strangle => {
            let sc = s.strike_call.ok_or(AppError::BadRequest("strangle requires 'strike_call'".into()))?;
            let sp = s.strike_put.ok_or(AppError::BadRequest("strangle requires 'strike_put'".into()))?;
            OptionSpreads::new_strangle(direction, sc, sp, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::Strip => {
            let strike = s.strike.ok_or(AppError::BadRequest("strip requires 'strike'".into()))?;
            OptionSpreads::new_strip(strike, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::Strap => {
            let strike = s.strike.ok_or(AppError::BadRequest("strap requires 'strike'".into()))?;
            OptionSpreads::new_strap(strike, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::SyntheticStock => {
            let strike = s.strike.ok_or(AppError::BadRequest("synthetic_stock requires 'strike'".into()))?;
            OptionSpreads::new_synthetic_stock(direction, strike, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::BullSpreadCall => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("bull_spread_call requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("bull_spread_call requires 'strike_high'".into()))?;
            OptionSpreads::new_bull_spread_call(sl, sh, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::BullSpreadPut => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("bull_spread_put requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("bull_spread_put requires 'strike_high'".into()))?;
            OptionSpreads::new_bull_spread_put(sl, sh, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::BearSpreadCall => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("bear_spread_call requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("bear_spread_call requires 'strike_high'".into()))?;
            OptionSpreads::new_bear_spread_call(sl, sh, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::BearSpreadPut => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("bear_spread_put requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("bear_spread_put requires 'strike_high'".into()))?;
            OptionSpreads::new_bear_spread_put(sl, sh, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::Butterfly => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("butterfly requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("butterfly requires 'strike_high'".into()))?;
            let sm = s.strike_medium.ok_or(AppError::BadRequest("butterfly requires 'strike_medium'".into()))?;
            OptionSpreads::new_butterfly(direction, sl, sh, sm, spot, vol, rfr, ttm, div)
        }
        SpreadTypeInput::CalendarSpread => {
            let strike = s.strike.ok_or(AppError::BadRequest("calendar_spread requires 'strike'".into()))?;
            let st = req.short_term_maturity.ok_or(AppError::BadRequest("calendar_spread requires 'short_term_maturity'".into()))?;
            let lt = req.long_term_maturity.ok_or(AppError::BadRequest("calendar_spread requires 'long_term_maturity'".into()))?;
            OptionSpreads::new_calendar_spread(direction, strike, spot, vol, rfr, st, lt, div)
        }
    };

    let bs_price = spread.bs_pricing();

    let mc_result = monte_carlo(&spread, req.mc_simulations);
    let (ci_lower, ci_upper) = mc_result.confidence_interval_95();

    let greeks = GreeksResult {
        delta: spread.delta(vol, spot),
        gamma: spread.gamma(vol, spot),
        theta: spread.theta(vol, spot),
        vega: spread.vega(spot),
        rho: spread.rho(vol, spot, rfr),
    };

    Ok(PriceResponse {
        structure_type: "spread".to_string(),
        pricing: PricingResult {
            black_scholes: bs_price,
            monte_carlo: Some(MonteCarloResult {
                price: mc_result.price,
                std_error: mc_result.std_error,
                ci_lower,
                ci_upper,
            }),
            binomial_european: None,
            binomial_american: None,
            bs_american_approx: None,
        },
        greeks: Some(greeks),
    })
}

fn price_exotic(req: ExoticRequest) -> Result<PriceResponse, AppError> {
    match req.exotic_type {
        ExoticTypeInput::ConvertibleBond => {
            let params: ConvertibleBondRequest = serde_json::from_value(req.params)
                .map_err(|e| AppError::BadRequest(format!("Invalid convertible bond params: {}", e)))?;

            let cb = ConvertibleBond {
                face_value: params.face_value,
                coupon_rate: params.coupon_rate,
                maturity: params.maturity,
                payment_frequency: params.payment_frequency,
                credit_spread: params.credit_spread,
                conversion_price: params.conversion_price,
                stock_price: params.stock_price,
                volatility: params.volatility,
                time_to_maturity: params.time_to_maturity,
                risk_free_rate: params.risk_free_rate,
                dividend_yield: params.dividend_yield,
            };

            Ok(PriceResponse {
                structure_type: "exotic".to_string(),
                pricing: PricingResult {
                    black_scholes: cb.bs_pricing(),
                    monte_carlo: None,
                    binomial_european: None,
                    binomial_american: None,
                    bs_american_approx: None,
                },
                greeks: None,
            })
        }
        ExoticTypeInput::ChooserOption => {
            let params: ChooserOptionRequest = serde_json::from_value(req.params)
                .map_err(|e| AppError::BadRequest(format!("Invalid chooser option params: {}", e)))?;

            let chooser = ChooserOption::new(
                params.spot_price, params.strike_price, params.volatility,
                params.risk_free_rate, params.time_to_maturity, params.dividend_yield,
                params.choose_time,
            );

            Ok(PriceResponse {
                structure_type: "exotic".to_string(),
                pricing: PricingResult {
                    black_scholes: chooser.bs_pricing(),
                    monte_carlo: None,
                    binomial_european: None,
                    binomial_american: None,
                    bs_american_approx: None,
                },
                greeks: None,
            })
        }
    }
}
