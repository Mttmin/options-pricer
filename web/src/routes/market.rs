use std::sync::Arc;
use axum::extract::{Path, Query, State};
use axum::Json;
use cli::options_chain::OptionsEndpoint;
use crate::error::AppError;
use crate::models::request::MarketQuery;
use crate::models::response::MarketDataResponse;
use crate::state::AppState;

pub async fn get_market_data(
    State(state): State<Arc<AppState>>,
    Path(symbol): Path<String>,
    Query(query): Query<MarketQuery>,
) -> Result<Json<MarketDataResponse>, AppError> {
    let symbol = symbol.to_uppercase();

    let market_data = state
        .fetcher
        .update_symbol(&symbol, query.lookback_days)
        .await
        .map_err(|e| AppError::NotFound(format!("Failed to fetch data for {}: {}", symbol, e)))?;

    // Fetch IV in background, don't fail if unavailable
    let iv = state
        .fetcher
        .fetch_implied_volatility(&symbol, OptionsEndpoint::Historical { date: None })
        .await
        .ok();

    Ok(Json(MarketDataResponse {
        symbol,
        spot_price: market_data.spot_price,
        historical_volatility: market_data.volatility,
        implied_volatility: iv.or(market_data.implied_volatility),
        corporate_action_detected: market_data.corporate_action_detected,
        dividend_yield: market_data.dividend_yield,
    }))
}
