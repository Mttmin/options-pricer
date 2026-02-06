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

    let chain = state
        .fetcher
        .fetch_option_chain(&symbol, OptionsEndpoint::Historical { date: None })
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
