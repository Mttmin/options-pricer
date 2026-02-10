use std::sync::Arc;
use axum::extract::{Path, State};
use axum::Json;
use cli::options_chain::{OptionsEndpoint, VolatilitySmile};
use crate::error::AppError;
use crate::state::AppState;

pub async fn get_option_chain(
    State(state): State<Arc<AppState>>,
    Path(symbol): Path<String>,
) -> Result<Json<VolatilitySmile>, AppError> {
    let symbol = symbol.to_uppercase();

    let market_data = state
        .fetcher
        .update_symbol(&symbol, 90)
        .await
        .map_err(|e| {
            AppError::NotFound(format!(
                "Failed to fetch spot price for {}: {}",
                symbol, e
            ))
        })?;

    let chain = state
        .fetcher
        .fetch_option_chain_with_underlying(
            &symbol,
            OptionsEndpoint::Historical { date: None },
            market_data.spot_price,
        )
        .await
        .map_err(|e| {
            AppError::NotFound(format!(
                "Failed to fetch option chain for {}: {}",
                symbol, e
            ))
        })?;

    let smile = chain.extract_smile(None, 50);

    Ok(Json(smile))
}
