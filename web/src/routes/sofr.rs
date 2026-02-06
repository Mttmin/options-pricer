use std::sync::Arc;
use axum::extract::State;
use axum::Json;
use serde::Serialize;
use crate::error::AppError;
use crate::state::AppState;

#[derive(Serialize)]
pub struct SofrResponse {
    pub rate: f64,
    pub date: String,
}

pub async fn get_sofr(
    State(state): State<Arc<AppState>>,
) -> Result<Json<SofrResponse>, AppError> {
    let rate = state
        .fetcher
        .fetch_fred_sofr()
        .await
        .map_err(|e| AppError::Internal(format!("Failed to fetch SOFR: {}", e)))?;

    Ok(Json(SofrResponse {
        rate,
        date: chrono::Utc::now().format("%Y-%m-%d").to_string(),
    }))
}
