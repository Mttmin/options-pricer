use std::sync::Arc;
use axum::extract::{Path, Query, State};
use axum::Json;
use cli::fetcher::DataFetcher;
use cli::options_chain::OptionsEndpoint;
use cli::volatility;
use crate::error::AppError;
use crate::models::request::MarketQuery;
use crate::models::response::VolatilityResponse;
use crate::state::AppState;

pub async fn get_volatility(
    State(state): State<Arc<AppState>>,
    Path(symbol): Path<String>,
    Query(query): Query<MarketQuery>,
) -> Result<Json<VolatilityResponse>, AppError> {
    let symbol = symbol.to_uppercase();

    // Fetch historical data for volatility calculations
    let (prices, _) = state
        .fetcher
        .get_historical_data(&symbol, query.lookback_days)
        .await
        .map_err(|e| AppError::Internal(format!("Failed to fetch historical data: {}", e)))?;

    let historical = DataFetcher::calculate_volatility(&prices);

    // EMA volatility
    let ema = volatility::ema_volatility(&prices, 30, 0.9).ok();

    // VIX-correlated volatility (requires cloning DataFetcher since vix_volatility takes ownership)
    let fetcher_clone = state.fetcher.clone();
    let symbol_clone = symbol.clone();
    let (vix_correlated, vix_fallback_used) = tokio::task::spawn(async move {
        match volatility::vix_volatility(&symbol_clone, fetcher_clone, 120).await {
            Ok(vix) => (Some(vix), Some(false)),
            Err(fallback) => (Some(fallback), Some(true)),
        }
    })
    .await
    .ok()
    .unwrap_or((None, None));

    // Implied volatility from options chain using live spot for ATM selection
    let implied = if let Ok(market) = state
        .fetcher
        .update_symbol(&symbol, query.lookback_days)
        .await
    {
        state
            .fetcher
            .fetch_option_chain_with_underlying(
                &symbol,
                OptionsEndpoint::Historical { date: None },
                market.spot_price,
            )
            .await
            .ok()
            .and_then(|chain| chain.summary_implied_volatility)
    } else {
        None
    };

    Ok(Json(VolatilityResponse {
        symbol,
        historical,
        ema,
        vix_correlated,
        vix_fallback_used,
        implied,
    }))
}
