use axum::extract::State;
use axum::Json;
use std::sync::Arc;
use std::time::Instant;
use options::black_scholes::{black_scholes_approx_american, black_scholes_price, spread_black_scholes_approx_american};
use options::exotics::{AsianKind, AsianOption, AsianPooling, AutocallableNote, BarrierKind, BarrierOption, BarrierType, ChooserOption, CliquetKind, CliquetOption, ConvertibleBond};
use options::optionspreads::{Direction, OptionSpreads};
use options::{BlackScholesPrice, Call, Greeks, MonteCarloParameters, Options, Payoff, Put};
use numerical_methods::binomial::{binomial_price, BinomialParameters, ExerciseStyle};
use numerical_methods::{monte_carlo, monte_carlo_asian, monte_carlo_autocallable, monte_carlo_barrier, monte_carlo_cliquet, price_european, PenaltySolver, TimestepMode};

use crate::error::AppError;
use crate::models::request::*;
use crate::models::response::*;
use crate::state::AppState;

async fn fetch_dividend_for_request(
    fetcher: &cli::fetcher::DataFetcher,
    symbol: Option<&String>,
    spot_price: f64,
) -> Option<(f64, f64)> {
    if let Some(sym) = symbol {
        match fetcher.get_dividend_info(sym, spot_price).await {
            Ok(Some(div)) => Some((
                div.time_to_dividend,
                div.dividend_amount / spot_price
            )),
            Ok(None) => None,
            Err(e) => {
                eprintln!("Failed to fetch dividend info for {}: {}", sym, e);
                None
            }
        }
    } else {
        None
    }
}

pub async fn calculate_price(
    State(state): State<Arc<AppState>>,
    Json(request): Json<PriceRequest>,
) -> Result<Json<PriceResponse>, AppError> {
    let dividend_info = match &request {
        PriceRequest::Single(req) => fetch_dividend_for_request(
            &state.fetcher,
            req.symbol.as_ref(),
            req.spot_price
        ).await,
        PriceRequest::Spread(req) => fetch_dividend_for_request(
            &state.fetcher,
            req.symbol.as_ref(),
            req.spot_price
        ).await,
        PriceRequest::Exotic(_) => None,
    };

    let response = tokio::task::spawn_blocking(move || match request {
        PriceRequest::Single(req) => price_single(req, dividend_info),
        PriceRequest::Spread(req) => price_spread(req, dividend_info),
        PriceRequest::Exotic(req) => price_exotic(req),
    })
    .await
    .map_err(|e| AppError::Internal(format!("Pricing task failed: {}", e)))??;

    Ok(Json(response))
}

fn price_single(req: SingleOptionRequest, dividend_info: Option<(f64, f64)>) -> Result<PriceResponse, AppError> {
    let total_start = Instant::now();
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
    let bs_start = Instant::now();
    let bs_price = black_scholes_price(option) * direction_sign;
    let bs_ms = bs_start.elapsed().as_secs_f64() * 1000.0;

    // Monte Carlo (monomorphized for Call/Put)
    let mc_start = Instant::now();
    let mc_result = match req.option_type {
        OptionTypeInput::Call => price_european::<Call>(&option, req.mc_simulations),
        OptionTypeInput::Put => price_european::<Put>(&option, req.mc_simulations),
    };
    let mc_ms = mc_start.elapsed().as_secs_f64() * 1000.0;
    let (ci_lower, ci_upper) = mc_result.confidence_interval_95();

    // Binomial tree
    let params = BinomialParameters::new(
        req.spot_price, req.binomial_depth, req.time_to_maturity,
        req.volatility, req.risk_free_rate,
    );

    let binomial_euro_start = Instant::now();
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
    let binomial_euro_ms = binomial_euro_start.elapsed().as_secs_f64() * 1000.0;

    let binomial_amer_start = Instant::now();
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
    let binomial_amer_ms = binomial_amer_start.elapsed().as_secs_f64() * 1000.0;

    // Black's American approximation
    let bs_amer_start = Instant::now();
    let bs_american = black_scholes_approx_american(option, dividend_info) * direction_sign;
    let bs_amer_ms = bs_amer_start.elapsed().as_secs_f64() * 1000.0;

    // PenaltySolver (PDE method for American options)
    let penalty_start = Instant::now();
    let penalty_result = compute_penalty_solver(&req, direction_sign).ok();
    let penalty_ms = if penalty_result.is_some() {
        Some(penalty_start.elapsed().as_secs_f64() * 1000.0)
    } else {
        None
    };

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

    let mut spot_prices = Vec::with_capacity(100);
    let mut payoffs = Vec::with_capacity(100);
    let min_spot = spot * 0.5;
    let max_spot = spot * 1.5;
    let step = (max_spot - min_spot) / 99.0;
    for i in 0..100 {
        let s = min_spot + (i as f64) * step;
        spot_prices.push(s);
        payoffs.push((option.compute(s) - bs_price.abs()) * direction_sign);
    }

    Ok(PriceResponse {
        structure_type: "single".to_string(),
        pricing: PricingResult {
            black_scholes: Some(bs_price),
            monte_carlo: Some(MonteCarloResult {
                price: mc_result.price * direction_sign,
                std_error: mc_result.std_error,
                ci_lower: ci_lower * direction_sign,
                ci_upper: ci_upper * direction_sign,
            }),
            binomial_european: Some(binomial_euro * direction_sign),
            binomial_american: Some(binomial_amer * direction_sign),
            bs_american_approx: Some(bs_american),
            penalty_solver: penalty_result,
            timings: Some(PricingTimings {
                total_ms: total_start.elapsed().as_secs_f64() * 1000.0,
                black_scholes_ms: Some(bs_ms),
                monte_carlo_ms: Some(mc_ms),
                binomial_european_ms: Some(binomial_euro_ms),
                binomial_american_ms: Some(binomial_amer_ms),
                bs_american_approx_ms: Some(bs_amer_ms),
                penalty_solver_ms: penalty_ms,
            }),
        },
        greeks: Some(greeks),
        payoff_curve: Some(PayoffCurve { spot_prices, payoffs }),
    })
}

fn compute_penalty_solver(
    req: &SingleOptionRequest,
    direction_sign: f64,
) -> Result<PenaltySolverResult, AppError> {
    let spatial_nodes = 200;
    let timesteps = 100;
    let penalty_tolerance = 1e-6;
    let timestep_mode = TimestepMode::Adaptive { dnorm: 0.02, scale_d: 1.0 };

    match req.option_type {
        OptionTypeInput::Call => {
            let call = Call::new(
                req.strike_price, req.spot_price, req.volatility,
                req.risk_free_rate, req.time_to_maturity, req.dividend_yield
            );
            solve_and_extract(call, spatial_nodes, timesteps, penalty_tolerance, timestep_mode, direction_sign)
        }
        OptionTypeInput::Put => {
            let put = Put::new(
                req.strike_price, req.spot_price, req.volatility,
                req.risk_free_rate, req.time_to_maturity, req.dividend_yield
            );
            solve_and_extract(put, spatial_nodes, timesteps, penalty_tolerance, timestep_mode, direction_sign)
        }
    }
}

fn solve_and_extract<T>(
    derivative: T,
    spatial_nodes: usize,
    timesteps: usize,
    penalty_tolerance: f64,
    timestep_mode: TimestepMode,
    direction_sign: f64,
) -> Result<PenaltySolverResult, AppError>
where
    T: Payoff + MonteCarloParameters,
{
    let mut solver = PenaltySolver::new(
        derivative,
        spatial_nodes,
        timesteps,
        penalty_tolerance,
        timestep_mode,
    );

    solver.initialize();
    solver.solve();

    let price = solver.option_value() * direction_sign;
    let avg_iters = solver.total_iterations as f64 / solver.time_grid.num_steps as f64;
    let grid_min = solver.spatial_grid.prices[0];
    let grid_max = solver.spatial_grid.prices[solver.spatial_grid.num_nodes - 1];

    Ok(PenaltySolverResult {
        price,
        diagnostics: PenaltySolverDiagnostics {
            total_iterations: solver.total_iterations,
            avg_iterations_per_step: avg_iters,
            max_american_error: solver.max_american_error,
            spatial_nodes: solver.spatial_grid.num_nodes,
            timesteps: solver.time_grid.num_steps,
            grid_min_price: grid_min,
            grid_max_price: grid_max,
        },
    })
}

fn to_direction(d: DirectionInput) -> Direction {
    match d {
        DirectionInput::Long => Direction::LONG,
        DirectionInput::Short => Direction::SHORT,
    }
}

fn price_spread(req: SpreadRequest, dividend_info: Option<(f64, f64)>) -> Result<PriceResponse, AppError> {
    let total_start = Instant::now();
    let direction = to_direction(req.direction);
    let s = &req.strikes;
    let spot = req.spot_price;
    let vol = req.volatility;
    let rfr = req.risk_free_rate;
    let ttm = req.time_to_maturity;
    let div = req.dividend_yield;

    let (spread, allow_pde) = match req.spread_type {
        SpreadTypeInput::Straddle => {
            let strike = s.strike.ok_or(AppError::BadRequest("straddle requires 'strike'".into()))?;
            (OptionSpreads::new_straddle(direction, strike, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::Strangle => {
            let sc = s.strike_call.ok_or(AppError::BadRequest("strangle requires 'strike_call'".into()))?;
            let sp = s.strike_put.ok_or(AppError::BadRequest("strangle requires 'strike_put'".into()))?;
            (OptionSpreads::new_strangle(direction, sc, sp, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::Strip => {
            let strike = s.strike.ok_or(AppError::BadRequest("strip requires 'strike'".into()))?;
            (OptionSpreads::new_strip(strike, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::Strap => {
            let strike = s.strike.ok_or(AppError::BadRequest("strap requires 'strike'".into()))?;
            (OptionSpreads::new_strap(strike, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::SyntheticStock => {
            let strike = s.strike.ok_or(AppError::BadRequest("synthetic_stock requires 'strike'".into()))?;
            (OptionSpreads::new_synthetic_stock(direction, strike, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::BullSpreadCall => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("bull_spread_call requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("bull_spread_call requires 'strike_high'".into()))?;
            (OptionSpreads::new_bull_spread_call(sl, sh, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::BullSpreadPut => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("bull_spread_put requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("bull_spread_put requires 'strike_high'".into()))?;
            (OptionSpreads::new_bull_spread_put(sl, sh, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::BearSpreadCall => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("bear_spread_call requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("bear_spread_call requires 'strike_high'".into()))?;
            (OptionSpreads::new_bear_spread_call(sl, sh, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::BearSpreadPut => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("bear_spread_put requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("bear_spread_put requires 'strike_high'".into()))?;
            (OptionSpreads::new_bear_spread_put(sl, sh, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::Butterfly => {
            let sl = s.strike_low.ok_or(AppError::BadRequest("butterfly requires 'strike_low'".into()))?;
            let sh = s.strike_high.ok_or(AppError::BadRequest("butterfly requires 'strike_high'".into()))?;
            let sm = s.strike_medium.ok_or(AppError::BadRequest("butterfly requires 'strike_medium'".into()))?;
            (OptionSpreads::new_butterfly(direction, sl, sh, sm, spot, vol, rfr, ttm, div), true)
        }
        SpreadTypeInput::CalendarSpread => {
            let strike = s.strike.ok_or(AppError::BadRequest("calendar_spread requires 'strike'".into()))?;
            let st = req.short_term_maturity.ok_or(AppError::BadRequest("calendar_spread requires 'short_term_maturity'".into()))?;
            let lt = req.long_term_maturity.ok_or(AppError::BadRequest("calendar_spread requires 'long_term_maturity'".into()))?;
            (OptionSpreads::new_calendar_spread(direction, strike, spot, vol, rfr, st, lt, div), false)
        }
    };

    let binomial_depth: u16 = 1000;
    let binomial_euro_start = Instant::now();
    let binomial_euro = spread_binomial(&spread, ExerciseStyle::European, binomial_depth);
    let binomial_euro_ms = binomial_euro_start.elapsed().as_secs_f64() * 1000.0;

    let binomial_amer_start = Instant::now();
    let binomial_amer = spread_binomial(&spread, ExerciseStyle::American, binomial_depth);
    let binomial_amer_ms = binomial_amer_start.elapsed().as_secs_f64() * 1000.0;

    let bs_start = Instant::now();
    let bs_price = spread.bs_pricing();
    let bs_ms = bs_start.elapsed().as_secs_f64() * 1000.0;

    let mc_start = Instant::now();
    let mc_result = monte_carlo(&spread, req.mc_simulations);
    let mc_ms = mc_start.elapsed().as_secs_f64() * 1000.0;
    let (ci_lower, ci_upper) = mc_result.confidence_interval_95();

    let greeks = GreeksResult {
        delta: spread.delta(vol, spot),
        gamma: spread.gamma(vol, spot),
        theta: spread.theta(vol, spot),
        vega: spread.vega(spot),
        rho: spread.rho(vol, spot, rfr),
    };

    let bs_amer_spread_start = Instant::now();
    let bs_american_approx_spread = spread_black_scholes_approx_american(
        spread.components(),
        dividend_info
    );
    let bs_amer_spread_ms = bs_amer_spread_start.elapsed().as_secs_f64() * 1000.0;

    let penalty_start = Instant::now();
    let penalty_result = if allow_pde {
        compute_spread_penalty_solver(spread.clone()).ok()
    } else {
        None
    };
    let penalty_ms = if penalty_result.is_some() {
        Some(penalty_start.elapsed().as_secs_f64() * 1000.0)
    } else {
        None
    };

    let mut spot_prices = Vec::with_capacity(100);
    let mut payoffs = Vec::with_capacity(100);
    let min_spot = spot * 0.5;
    let max_spot = spot * 1.5;
    let step = (max_spot - min_spot) / 99.0;

    let is_calendar = matches!(req.spread_type, SpreadTypeInput::CalendarSpread);
    let short_term_maturity = req.short_term_maturity.unwrap_or(0.0);

    for i in 0..100 {
        let s = min_spot + (i as f64) * step;
        spot_prices.push(s);
        
        let mut payoff = 0.0;
        if is_calendar {
            for pos in spread.components() {
                let opt = pos.option;
                if (opt.time_to_maturity() - short_term_maturity).abs() < 1e-6 {
                    // Short-term leg expires, intrinsic value
                    payoff += opt.compute(s) * pos.weight as f64;
                } else {
                    // Long-term leg, price with remaining time
                    let remaining_ttm = opt.time_to_maturity() - short_term_maturity;
                    let new_opt = match opt {
                        Options::Call(_) => Options::new_call(
                            opt.strike_price(), s, opt.volatility(), opt.risk_free_rate(), remaining_ttm, opt.dividend_yield()
                        ),
                        Options::Put(_) => Options::new_put(
                            opt.strike_price(), s, opt.volatility(), opt.risk_free_rate(), remaining_ttm, opt.dividend_yield()
                        ),
                    };
                    payoff += new_opt.bs_pricing() * pos.weight as f64;
                }
            }
        } else {
            payoff = spread.compute(s);
        }
        
        payoffs.push(payoff - bs_price);
    }

    Ok(PriceResponse {
        structure_type: "spread".to_string(),
        pricing: PricingResult {
            black_scholes: Some(bs_price),
            monte_carlo: Some(MonteCarloResult {
                price: mc_result.price,
                std_error: mc_result.std_error,
                ci_lower,
                ci_upper,
            }),
            binomial_european: Some(binomial_euro),
            binomial_american: Some(binomial_amer),
            bs_american_approx: Some(bs_american_approx_spread),
            penalty_solver: penalty_result,
            timings: Some(PricingTimings {
                total_ms: total_start.elapsed().as_secs_f64() * 1000.0,
                black_scholes_ms: Some(bs_ms),
                monte_carlo_ms: Some(mc_ms),
                binomial_european_ms: Some(binomial_euro_ms),
                binomial_american_ms: Some(binomial_amer_ms),
                bs_american_approx_ms: Some(bs_amer_spread_ms),
                penalty_solver_ms: penalty_ms,
            }),
        },
        greeks: Some(greeks),
        payoff_curve: Some(PayoffCurve { spot_prices, payoffs }),
    })
}

fn compute_spread_penalty_solver(
    spread: OptionSpreads,
) -> Result<PenaltySolverResult, AppError> {
    let spatial_nodes = 200;
    let timesteps = 100;
    let penalty_tolerance = 1e-6;
    let timestep_mode = TimestepMode::Adaptive { dnorm: 0.02, scale_d: 1.0 };

    solve_and_extract(
        spread,
        spatial_nodes,
        timesteps,
        penalty_tolerance,
        timestep_mode,
        1.0,
    )
}

fn spread_binomial(spread: &OptionSpreads, style: ExerciseStyle, depth: u16) -> f64 {
    spread
        .components()
        .iter()
        .map(|position| {
            let option = position.option;
            let params = BinomialParameters::new(
                option.spot_price(),
                depth,
                option.time_to_maturity(),
                option.volatility(),
                option.risk_free_rate(),
            );
            binomial_price(&option, &params, style) * position.weight as f64
        })
        .sum()
}

fn price_exotic(req: ExoticRequest) -> Result<PriceResponse, AppError> {
    match req.exotic_type {
        ExoticTypeInput::ConvertibleBond => {
            let total_start = Instant::now();
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

            let bs_start = Instant::now();
            let bs_price = cb.bs_pricing();
            let bs_ms = bs_start.elapsed().as_secs_f64() * 1000.0;

            let mut spot_prices = Vec::with_capacity(100);
            let mut payoffs = Vec::with_capacity(100);
            let min_spot = params.stock_price * 0.5;
            let max_spot = params.stock_price * 1.5;
            let step = (max_spot - min_spot) / 99.0;
            for i in 0..100 {
                let s = min_spot + (i as f64) * step;
                spot_prices.push(s);
                let conversion_value = (params.face_value / params.conversion_price) * s;
                let payoff = params.face_value.max(conversion_value);
                payoffs.push(payoff - bs_price);
            }

            Ok(PriceResponse {
                structure_type: "exotic".to_string(),
                pricing: PricingResult {
                    black_scholes: Some(bs_price),
                    monte_carlo: None,
                    binomial_european: None,
                    binomial_american: None,
                    bs_american_approx: None,
                    penalty_solver: None,
                    timings: Some(PricingTimings {
                        total_ms: total_start.elapsed().as_secs_f64() * 1000.0,
                        black_scholes_ms: Some(bs_ms),
                        monte_carlo_ms: None,
                        binomial_european_ms: None,
                        binomial_american_ms: None,
                        bs_american_approx_ms: None,
                        penalty_solver_ms: None,
                    }),
                },
                greeks: None,
                payoff_curve: Some(PayoffCurve { spot_prices, payoffs }),
            })
        }
        ExoticTypeInput::ChooserOption => {
            let total_start = Instant::now();
            let params: ChooserOptionRequest = serde_json::from_value(req.params)
                .map_err(|e| AppError::BadRequest(format!("Invalid chooser option params: {}", e)))?;

            let chooser = ChooserOption::new(
                params.spot_price, params.strike_price, params.volatility,
                params.risk_free_rate, params.time_to_maturity, params.dividend_yield,
                params.choose_time,
            );

            let bs_start = Instant::now();
            let bs_price = chooser.bs_pricing();
            let bs_ms = bs_start.elapsed().as_secs_f64() * 1000.0;

            let mut spot_prices = Vec::with_capacity(100);
            let mut payoffs = Vec::with_capacity(100);
            let min_spot = params.spot_price * 0.5;
            let max_spot = params.spot_price * 1.5;
            let step = (max_spot - min_spot) / 99.0;
            let remaining_ttm = params.time_to_maturity - params.choose_time;
            
            for i in 0..100 {
                let s = min_spot + (i as f64) * step;
                spot_prices.push(s);
                
                let call = Options::new_call(
                    params.strike_price, s, params.volatility, params.risk_free_rate, remaining_ttm, params.dividend_yield
                );
                let put = Options::new_put(
                    params.strike_price, s, params.volatility, params.risk_free_rate, remaining_ttm, params.dividend_yield
                );
                
                let payoff = call.bs_pricing().max(put.bs_pricing());
                payoffs.push(payoff - bs_price);
            }

            Ok(PriceResponse {
                structure_type: "exotic".to_string(),
                pricing: PricingResult {
                    black_scholes: Some(bs_price),
                    monte_carlo: None,
                    binomial_european: None,
                    binomial_american: None,
                    bs_american_approx: None,
                    penalty_solver: None,
                    timings: Some(PricingTimings {
                        total_ms: total_start.elapsed().as_secs_f64() * 1000.0,
                        black_scholes_ms: Some(bs_ms),
                        monte_carlo_ms: None,
                        binomial_european_ms: None,
                        binomial_american_ms: None,
                        bs_american_approx_ms: None,
                        penalty_solver_ms: None,
                    }),
                },
                greeks: None,
                payoff_curve: Some(PayoffCurve { spot_prices, payoffs }),
            })
        }
        ExoticTypeInput::AsianOption => {
            let total_start = Instant::now();
            let params: AsianOptionRequest = serde_json::from_value(req.params)
                .map_err(|e| AppError::BadRequest(format!("Invalid asian option params: {}", e)))?;

            let kind = match params.option_type {
                OptionTypeInput::Call => AsianKind::Call,
                OptionTypeInput::Put => AsianKind::Put,
            };
            let pooling = match params.pooling {
                AsianPoolingInput::Average => AsianPooling::Average,
                AsianPoolingInput::Max => AsianPooling::Max,
            };
            let asian = AsianOption::new(
                kind, pooling,
                params.strike_price, params.spot_price, params.volatility,
                params.risk_free_rate, params.time_to_maturity, params.dividend_yield,
                params.num_observations,
            );

            let bs_start = Instant::now();
            let cf_price = asian.bs_pricing();
            let bs_ms = bs_start.elapsed().as_secs_f64() * 1000.0;

            let use_cv = pooling == AsianPooling::Average;
            let mc_start = Instant::now();
            let mc = monte_carlo_asian(&asian, params.mc_simulations, params.num_observations, use_cv);
            let mc_ms = mc_start.elapsed().as_secs_f64() * 1000.0;
            let mc_result = MonteCarloResult {
                price: mc.price,
                std_error: mc.std_error,
                ci_lower: mc.price - 1.96 * mc.std_error,
                ci_upper: mc.price + 1.96 * mc.std_error,
            };

            // Payoff curve: closed-form price as a function of S0.
            let mut spot_prices = Vec::with_capacity(100);
            let mut payoffs = Vec::with_capacity(100);
            let min_spot = params.spot_price * 0.5;
            let max_spot = params.spot_price * 1.5;
            let step = (max_spot - min_spot) / 99.0;
            for i in 0..100 {
                let s = min_spot + (i as f64) * step;
                spot_prices.push(s);
                let probe = AsianOption::new(
                    kind, pooling,
                    params.strike_price, s, params.volatility,
                    params.risk_free_rate, params.time_to_maturity, params.dividend_yield,
                    params.num_observations,
                );
                payoffs.push(probe.bs_pricing() - cf_price);
            }

            Ok(PriceResponse {
                structure_type: "exotic".to_string(),
                pricing: PricingResult {
                    black_scholes: Some(cf_price),
                    monte_carlo: Some(mc_result),
                    binomial_european: None,
                    binomial_american: None,
                    bs_american_approx: None,
                    penalty_solver: None,
                    timings: Some(PricingTimings {
                        total_ms: total_start.elapsed().as_secs_f64() * 1000.0,
                        black_scholes_ms: Some(bs_ms),
                        monte_carlo_ms: Some(mc_ms),
                        binomial_european_ms: None,
                        binomial_american_ms: None,
                        bs_american_approx_ms: None,
                        penalty_solver_ms: None,
                    }),
                },
                greeks: None,
                payoff_curve: Some(PayoffCurve { spot_prices, payoffs }),
            })
        }
        ExoticTypeInput::CliquetOption => {
            let total_start = Instant::now();
            let params: CliquetOptionRequest = serde_json::from_value(req.params)
                .map_err(|e| AppError::BadRequest(format!("Invalid cliquet option params: {}", e)))?;

            let kind = match params.option_type {
                OptionTypeInput::Call => CliquetKind::Call,
                OptionTypeInput::Put => CliquetKind::Put,
            };
            let cliquet = CliquetOption::new(
                kind, params.spot_price, params.volatility, params.risk_free_rate,
                params.time_to_maturity, params.dividend_yield, params.num_resets,
                params.local_cap, params.local_floor, params.global_cap, params.global_floor,
            );

            let bs_start = Instant::now();
            let cf_price = cliquet.bs_pricing();
            let bs_ms = bs_start.elapsed().as_secs_f64() * 1000.0;

            let mc_start = Instant::now();
            let mc = monte_carlo_cliquet(&cliquet, params.mc_simulations, params.num_resets);
            let mc_ms = mc_start.elapsed().as_secs_f64() * 1000.0;
            let mc_result = MonteCarloResult {
                price: mc.price,
                std_error: mc.std_error,
                ci_lower: mc.price - 1.96 * mc.std_error,
                ci_upper: mc.price + 1.96 * mc.std_error,
            };

            // Payoff curve: closed-form benchmark vs S0 (S0-independent, flat).
            let mut spot_prices = Vec::with_capacity(100);
            let mut payoffs = Vec::with_capacity(100);
            let min_spot = params.spot_price * 0.5;
            let max_spot = params.spot_price * 1.5;
            let step = (max_spot - min_spot) / 99.0;
            for i in 0..100 {
                let s = min_spot + (i as f64) * step;
                spot_prices.push(s);
                let probe = CliquetOption::new(
                    kind, s, params.volatility, params.risk_free_rate,
                    params.time_to_maturity, params.dividend_yield, params.num_resets,
                    params.local_cap, params.local_floor, params.global_cap, params.global_floor,
                );
                payoffs.push(probe.bs_pricing() - cf_price);
            }

            Ok(PriceResponse {
                structure_type: "exotic".to_string(),
                pricing: PricingResult {
                    black_scholes: Some(cf_price),
                    monte_carlo: Some(mc_result),
                    binomial_european: None,
                    binomial_american: None,
                    bs_american_approx: None,
                    penalty_solver: None,
                    timings: Some(PricingTimings {
                        total_ms: total_start.elapsed().as_secs_f64() * 1000.0,
                        black_scholes_ms: Some(bs_ms),
                        monte_carlo_ms: Some(mc_ms),
                        binomial_european_ms: None,
                        binomial_american_ms: None,
                        bs_american_approx_ms: None,
                        penalty_solver_ms: None,
                    }),
                },
                greeks: None,
                payoff_curve: Some(PayoffCurve { spot_prices, payoffs }),
            })
        }
        ExoticTypeInput::BarrierOption => {
            let total_start = Instant::now();
            let params: BarrierOptionRequest = serde_json::from_value(req.params)
                .map_err(|e| AppError::BadRequest(format!("Invalid barrier option params: {}", e)))?;

            let kind = match params.option_type {
                BarrierKindInput::Call => BarrierKind::Call,
                BarrierKindInput::Put => BarrierKind::Put,
            };
            let barrier_type = match params.barrier_type {
                BarrierTypeInput::UpAndIn => BarrierType::UpAndIn,
                BarrierTypeInput::UpAndOut => BarrierType::UpAndOut,
                BarrierTypeInput::DownAndIn => BarrierType::DownAndIn,
                BarrierTypeInput::DownAndOut => BarrierType::DownAndOut,
            };
            let barrier = BarrierOption {
                kind,
                barrier_type,
                spot_price: params.spot_price,
                strike_price: params.strike_price,
                barrier: params.barrier,
                rebate: params.rebate,
                volatility: params.volatility,
                risk_free_rate: params.risk_free_rate,
                time_to_maturity: params.time_to_maturity,
                dividend_yield: params.dividend_yield,
            };

            let bs_start = Instant::now();
            let bs_price = barrier.bs_pricing();
            let bs_ms = bs_start.elapsed().as_secs_f64() * 1000.0;

            let mc_start = Instant::now();
            let mc = monte_carlo_barrier(&barrier, params.mc_simulations, params.num_steps);
            let mc_ms = mc_start.elapsed().as_secs_f64() * 1000.0;
            let mc_result = MonteCarloResult {
                price: mc.price,
                std_error: mc.std_error,
                ci_lower: mc.price - 1.96 * mc.std_error,
                ci_upper: mc.price + 1.96 * mc.std_error,
            };

            let mut spot_prices = Vec::with_capacity(100);
            let mut payoffs = Vec::with_capacity(100);
            let min_spot = params.spot_price * 0.5;
            let max_spot = params.spot_price * 1.5;
            let step = (max_spot - min_spot) / 99.0;
            for i in 0..100 {
                let s = min_spot + (i as f64) * step;
                spot_prices.push(s);
                let probe = BarrierOption { spot_price: s, ..barrier };
                payoffs.push(probe.bs_pricing() - bs_price);
            }

            Ok(PriceResponse {
                structure_type: "exotic".to_string(),
                pricing: PricingResult {
                    black_scholes: Some(bs_price),
                    monte_carlo: Some(mc_result),
                    binomial_european: None,
                    binomial_american: None,
                    bs_american_approx: None,
                    penalty_solver: None,
                    timings: Some(PricingTimings {
                        total_ms: total_start.elapsed().as_secs_f64() * 1000.0,
                        black_scholes_ms: Some(bs_ms),
                        monte_carlo_ms: Some(mc_ms),
                        binomial_european_ms: None,
                        binomial_american_ms: None,
                        bs_american_approx_ms: None,
                        penalty_solver_ms: None,
                    }),
                },
                greeks: None,
                payoff_curve: Some(PayoffCurve { spot_prices, payoffs }),
            })
        }
        ExoticTypeInput::AutocallableNote => {
            let total_start = Instant::now();
            let params: AutocallableNoteRequest = serde_json::from_value(req.params)
                .map_err(|e| AppError::BadRequest(format!("Invalid autocallable note params: {}", e)))?;

            let note = AutocallableNote {
                spot_price: params.spot_price,
                volatility: params.volatility,
                risk_free_rate: params.risk_free_rate,
                time_to_maturity: params.time_to_maturity,
                dividend_yield: params.dividend_yield,
                notional: params.notional,
                autocall_barrier: params.autocall_barrier,
                coupon_rate: params.coupon_rate,
                protection_barrier: params.protection_barrier,
                num_observations: params.num_observations,
                memory_coupon: params.memory_coupon,
            };

            let mc_start = Instant::now();
            let mc = monte_carlo_autocallable(&note, params.mc_simulations);
            let mc_ms = mc_start.elapsed().as_secs_f64() * 1000.0;
            let mc_result = MonteCarloResult {
                price: mc.price,
                std_error: mc.std_error,
                ci_lower: mc.price - 1.96 * mc.std_error,
                ci_upper: mc.price + 1.96 * mc.std_error,
            };

            // Terminal payoff profile: analytical payout at maturity by spot level.
            let s0 = params.spot_price;
            let ac_level = params.autocall_barrier * s0;
            let ki_level = params.protection_barrier * s0;
            let t = params.time_to_maturity;
            let mut spot_prices = Vec::with_capacity(100);
            let mut payoffs = Vec::with_capacity(100);
            let min_spot = s0 * 0.5;
            let max_spot = s0 * 1.5;
            let step = (max_spot - min_spot) / 99.0;
            for i in 0..100 {
                let s = min_spot + (i as f64) * step;
                spot_prices.push(s);
                let terminal_payout = if s >= ac_level {
                    params.notional * (1.0 + params.coupon_rate * t)
                } else if s < ki_level {
                    params.notional * s / s0
                } else {
                    params.notional
                };
                payoffs.push(terminal_payout - mc.price);
            }

            Ok(PriceResponse {
                structure_type: "exotic".to_string(),
                pricing: PricingResult {
                    black_scholes: None,
                    monte_carlo: Some(mc_result),
                    binomial_european: None,
                    binomial_american: None,
                    bs_american_approx: None,
                    penalty_solver: None,
                    timings: Some(PricingTimings {
                        total_ms: total_start.elapsed().as_secs_f64() * 1000.0,
                        black_scholes_ms: None,
                        monte_carlo_ms: Some(mc_ms),
                        binomial_european_ms: None,
                        binomial_american_ms: None,
                        bs_american_approx_ms: None,
                        penalty_solver_ms: None,
                    }),
                },
                greeks: None,
                payoff_curve: Some(PayoffCurve { spot_prices, payoffs }),
            })
        }
    }
}
