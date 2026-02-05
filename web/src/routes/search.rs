use std::sync::Arc;
use axum::extract::{Query, State};
use axum::Json;
use serde::{Deserialize, Serialize};
use crate::error::AppError;
use crate::state::AppState;

#[derive(Deserialize)]
pub struct SearchQuery {
    pub q: String,
}

#[derive(Serialize)]
pub struct SearchResult {
    pub symbol: String,
    pub description: String,
    #[serde(rename = "type")]
    pub security_type: String,
}

pub async fn search_tickers(
    State(state): State<Arc<AppState>>,
    Query(query): Query<SearchQuery>,
) -> Result<Json<Vec<SearchResult>>, AppError> {
    if query.q.len() < 1 {
        return Ok(Json(vec![]));
    }

    let results = state
        .fetcher
        .search_symbols(&query.q)
        .await
        .map_err(|e| AppError::Internal(format!("Search failed: {}", e)))?;

    Ok(Json(results.into_iter().map(|(symbol, description, security_type)| {
        SearchResult { symbol, description, security_type }
    }).collect()))
}
