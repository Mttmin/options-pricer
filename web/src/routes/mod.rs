pub mod history;
pub mod market;
pub mod options_chain;
pub mod price;
pub mod search;
pub mod volatility;

use std::sync::Arc;
use axum::Router;
use axum::routing::{get, post};
use crate::state::AppState;

pub fn api_router() -> Router<Arc<AppState>> {
    Router::new()
        .route("/api/health", get(health))
        .route("/api/search", get(search::search_tickers))
        .route("/api/market/{symbol}", get(market::get_market_data))
        .route("/api/history/{symbol}", get(history::get_price_history))
        .route("/api/volatility/{symbol}", get(volatility::get_volatility))
        .route("/api/options-chain/{symbol}", get(options_chain::get_option_chain))
        .route("/api/price", post(price::calculate_price))
}

async fn health() -> axum::Json<serde_json::Value> {
    axum::Json(serde_json::json!({ "status": "ok" }))
}
