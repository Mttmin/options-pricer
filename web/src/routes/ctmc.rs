use std::sync::Arc;
use std::time::Instant;
use axum::extract::State;
use axum::Json;
use cli::options_chain::{CalibrationDataBuildConfig, OptionsEndpoint};
use numerical_methods::deep_cal::{
    BatesCalibrator, CalibrationConfig, calibrate_heston_from_dataset, default_onnx_path,
};
use numerical_methods::ctmc::price_american_option_heston;
use options::Options;

use crate::error::AppError;
use crate::models::request::{CtmcPriceRequest, DirectionInput, OptionTypeInput};
use crate::models::response::{CtmcPriceResponse, HestonResult};
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

    let response = tokio::task::spawn_blocking(move || -> Result<CtmcPriceResponse, AppError> {
        let start = Instant::now();

        let onnx_path = default_onnx_path()
            .ok_or_else(|| AppError::Internal("ONNX model weights not found".to_string()))?;

        let mut calibrator = BatesCalibrator::new(&onnx_path)
            .map_err(|e| AppError::Internal(format!("Failed to load ONNX model: {}", e)))?;

        let cal = calibrate_heston_from_dataset(
            &mut calibrator,
            &dataset,
            n_restarts,
            &CalibrationConfig::default(),
        )
        .map_err(|e| AppError::Internal(format!("Calibration failed: {}", e)))?;

        let h = &cal.heston;
        let ctmc = price_american_option_heston(
            option, h.v0, h.kappa, h.v_bar, h.sigma, h.rho, n_x, m_v, n_time,
        );

        let timing_ms = start.elapsed().as_secs_f64() * 1000.0;

        Ok(CtmcPriceResponse {
            symbol,
            spot_price,
            price: ctmc.price * direction_sign,
            heston: HestonResult {
                kappa: h.kappa,
                theta: h.v_bar,
                sigma: h.sigma,
                rho: h.rho,
                v0: h.v0,
            },
            ivrmse: cal.ivrmse,
            n_evals: cal.n_evals,
            n_instruments,
            timing_ms,
        })
    })
    .await
    .map_err(|e| AppError::Internal(format!("CTMC task failed: {}", e)))??;

    Ok(Json(response))
}
