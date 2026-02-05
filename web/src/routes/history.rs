use std::sync::Arc;
use axum::extract::{Path, Query, State};
use axum::Json;
use serde::{Deserialize, Serialize};
use crate::error::AppError;
use crate::state::AppState;

#[derive(Deserialize)]
pub struct HistoryQuery {
    #[serde(default = "default_days")]
    pub days: usize,
}

fn default_days() -> usize {
    90
}

#[derive(Serialize)]
pub struct PricePoint {
    pub date: String,
    pub price: f64,
}

pub async fn get_price_history(
    State(state): State<Arc<AppState>>,
    Path(symbol): Path<String>,
    Query(query): Query<HistoryQuery>,
) -> Result<Json<Vec<PricePoint>>, AppError> {
    let symbol = symbol.to_uppercase();

    let prices = state
        .fetcher
        .get_historical_prices(&symbol, query.days)
        .await
        .map_err(|e| AppError::NotFound(format!("Failed to fetch history for {}: {}", symbol, e)))?;

    let points: Vec<PricePoint> = prices
        .into_iter()
        .map(|(date, price)| PricePoint { date, price })
        .collect();

    Ok(Json(points))
}
