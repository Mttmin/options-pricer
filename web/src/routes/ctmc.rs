use std::sync::Arc;
use std::time::Instant;
use axum::extract::State;
use axum::Json;
use cli::options_chain::{CalibrationDataBuildConfig, OptionsEndpoint};
use numerical_methods::deep_cal::{
    BatesCalibrator, CalibrationConfig, calibrate_heston_from_dataset, default_onnx_path,
};
use numerical_methods::deep_cal_four_two::{
    FourTwoCalibrationConfig, FourTwoCalibrator, calibrate_four_two_from_dataset,
    default_onnx_path as default_four_two_onnx_path,
};
use numerical_methods::ctmc::{price_american_option_four_two, price_american_option_heston};
use options::Options;

use crate::error::AppError;
use crate::models::request::{CtmcPriceRequest, DirectionInput, OptionTypeInput};
use crate::models::response::{CtmcPriceResponse, FourTwoResult, HestonResult, ModelQuote};
use crate::state::AppState;

pub async fn calculate_ctmc_price(
    State(state): State<Arc<AppState>>,
    Json(request): Json<CtmcPriceRequest>,
) -> Result<Json<CtmcPriceResponse>, AppError> {
    let symbol = request.symbol.to_uppercase();

    let market_data = state
        .fetcher
        .update_symbol(&symbol, 90)
        .await
        .map_err(|e| AppError::NotFound(format!("Failed to fetch market data for {}: {}", symbol, e)))?;

    let chain = state
        .fetcher
        .fetch_option_chain_with_underlying(
            &symbol,
            OptionsEndpoint::Historical { date: None },
            market_data.spot_price,
        )
        .await
        .map_err(|e| AppError::NotFound(format!("Failed to fetch option chain for {}: {}", symbol, e)))?;

    let spot_price = market_data.spot_price;
    let rfr = request.risk_free_rate;
    let div = request
        .dividend_yield
        .or(market_data.dividend_yield)
        .unwrap_or(0.0);

    let cal_config = CalibrationDataBuildConfig {
        min_days_to_expiry: 1,
        fallback_risk_free_rate: rfr,
        fallback_dividend_yield: div,
    };
    let (dataset, _coverage) = chain.to_heston_calibration_dataset(&cal_config, None, None);
    let n_instruments = dataset.instruments.len();

    if n_instruments == 0 {
        return Err(AppError::BadRequest(format!(
            "No valid option contracts found for {}",
            symbol
        )));
    }

    let option = match request.option_type {
        OptionTypeInput::Call => Options::new_call(
            request.strike_price,
            spot_price,
            0.2,
            rfr,
            request.time_to_maturity,
            Some(div),
        ),
        OptionTypeInput::Put => Options::new_put(
            request.strike_price,
            spot_price,
            0.2,
            rfr,
            request.time_to_maturity,
            Some(div),
        ),
    };

    let direction_sign: f64 = match request.direction {
        DirectionInput::Long => 1.0,
        DirectionInput::Short => -1.0,
    };
    let n_restarts = request.n_restarts.unwrap_or(3);
    let n_x = request.n_x.unwrap_or(80);
    let m_v = request.m_v.unwrap_or(20);
    let n_time = request.n_time.unwrap_or(50);
    let model = request.model.to_ascii_lowercase();
    let run_four = model == "four_two" || model == "4/2" || model == "both";
    let run_heston = !run_four || model == "both";

    let start = Instant::now();
    let (heston_res, four_res) = if run_heston && run_four {
        let dataset_four = dataset.clone();
        let heston_task = tokio::task::spawn_blocking(move || {
            price_heston_model(option, dataset, n_restarts, n_x, m_v, n_time)
        });
        let four_task = tokio::task::spawn_blocking(move || {
            price_four_two_model(option, dataset_four, n_restarts, n_x, m_v, n_time)
        });
        let (heston_task, four_task) = tokio::join!(heston_task, four_task);
        (Some(flatten_task(heston_task)), Some(flatten_task(four_task)))
    } else if run_heston {
        let heston_task = tokio::task::spawn_blocking(move || {
            price_heston_model(option, dataset, n_restarts, n_x, m_v, n_time)
        })
        .await;
        (Some(flatten_task(heston_task)), None)
    } else {
        let four_task = tokio::task::spawn_blocking(move || {
            price_four_two_model(option, dataset, n_restarts, n_x, m_v, n_time)
        })
        .await;
        (None, Some(flatten_task(four_task)))
    };
    let timing_ms = start.elapsed().as_secs_f64() * 1000.0;

    let heston_quote = heston_res.as_ref().and_then(|r| r.as_ref().ok()).cloned();
    let four_quote = four_res.as_ref().and_then(|r| r.as_ref().ok()).cloned();
    let heston_error = heston_res.and_then(|r| r.err());
    let four_error = four_res.and_then(|r| r.err());
    if heston_quote.is_none() && four_quote.is_none() {
        let msg = heston_error.or(four_error).unwrap_or_else(|| "CTMC pricing failed".to_string());
        return Err(AppError::Internal(msg));
    }

    let primary = heston_quote.as_ref().or(four_quote.as_ref()).expect("one quote");
    Ok(Json(CtmcPriceResponse {
        symbol,
        spot_price,
        price: primary.price * direction_sign,
        heston: primary.heston.clone(),
        four_two: four_quote.as_ref().and_then(|q| q.four_two.clone()),
        model: if run_heston && run_four {
            "both".to_string()
        } else if run_four {
            "four_two".to_string()
        } else {
            "heston".to_string()
        },
        ivrmse: primary.ivrmse,
        n_evals: primary.n_evals,
        n_instruments,
        timing_ms,
        heston_quote: heston_quote.map(|mut q| {
            q.price *= direction_sign;
            q
        }),
        four_two_quote: four_quote.map(|mut q| {
            q.price *= direction_sign;
            q
        }),
        heston_error,
        four_two_error: four_error,
    }))
}

fn flatten_task(
    joined: Result<Result<ModelQuote, AppError>, tokio::task::JoinError>,
) -> Result<ModelQuote, String> {
    match joined {
        Ok(Ok(quote)) => Ok(quote),
        Ok(Err(AppError::BadRequest(msg) | AppError::NotFound(msg) | AppError::Internal(msg))) => Err(msg),
        Err(err) => Err(format!("CTMC task failed: {err}")),
    }
}

fn price_heston_model(
    option: Options,
    dataset: numerical_methods::slv::CalibrationDataset,
    n_restarts: usize,
    n_x: usize,
    m_v: usize,
    n_time: usize,
) -> Result<ModelQuote, AppError> {
    let start = Instant::now();
    let onnx_path = default_onnx_path()
        .ok_or_else(|| AppError::Internal("ONNX model weights not found".to_string()))?;
    let mut calibrator = BatesCalibrator::new(&onnx_path)
        .map_err(|e| AppError::Internal(format!("Failed to load ONNX model: {e}")))?;
    let cal = calibrate_heston_from_dataset(
        &mut calibrator, &dataset, n_restarts, &CalibrationConfig::default(),
    )
    .map_err(|e| AppError::Internal(format!("Calibration failed: {e}")))?;
    let h = &cal.heston;
    let ctmc = price_american_option_heston(
        option, h.v0, h.kappa, h.v_bar, h.sigma, h.rho, n_x, m_v, n_time,
    );
    Ok(ModelQuote {
        price: ctmc.price,
        ivrmse: cal.ivrmse,
        n_evals: cal.n_evals,
        timing_ms: start.elapsed().as_secs_f64() * 1000.0,
        heston: HestonResult {
            kappa: h.kappa,
            theta: h.v_bar,
            sigma: h.sigma,
            rho: h.rho,
            v0: h.v0,
        },
        four_two: None,
    })
}

fn price_four_two_model(
    option: Options,
    dataset: numerical_methods::slv::CalibrationDataset,
    n_restarts: usize,
    n_x: usize,
    m_v: usize,
    n_time: usize,
) -> Result<ModelQuote, AppError> {
    let start = Instant::now();
    let onnx_path = default_four_two_onnx_path()
        .ok_or_else(|| AppError::Internal("4/2 ONNX weights not found".to_string()))?;
    let mut calibrator = FourTwoCalibrator::new(&onnx_path)
        .map_err(|e| AppError::Internal(format!("Failed to load 4/2 ONNX model: {e}")))?;
    let cal = calibrate_four_two_from_dataset(
        &mut calibrator, &dataset, n_restarts, &FourTwoCalibrationConfig::default(),
    )
    .map_err(|e| AppError::Internal(format!("4/2 calibration failed: {e}")))?;
    let p = &cal.params;
    let ctmc = price_american_option_four_two(
        option, p.v0, p.kappa, p.v_bar, p.sigma, p.rho, p.a, p.b, n_x, m_v, n_time,
    );
    let fitted = HestonResult {
        kappa: p.kappa,
        theta: p.v_bar,
        sigma: p.sigma,
        rho: p.rho,
        v0: p.v0,
    };
    Ok(ModelQuote {
        price: ctmc.price,
        ivrmse: cal.ivrmse,
        n_evals: cal.n_evals,
        timing_ms: start.elapsed().as_secs_f64() * 1000.0,
        heston: fitted.clone(),
        four_two: Some(FourTwoResult {
            kappa: fitted.kappa,
            theta: fitted.theta,
            sigma: fitted.sigma,
            rho: fitted.rho,
            v0: fitted.v0,
            a: p.a,
            b: p.b,
        }),
    })
}
