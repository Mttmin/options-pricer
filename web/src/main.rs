mod error;
mod models;
mod routes;
mod state;

use std::env;
use std::sync::Arc;
use cli::fetcher::DataFetcher;
use tower_http::cors::{Any, CorsLayer};
use tower_http::services::ServeDir;
use state::AppState;

#[tokio::main]
async fn main() {
    let fetcher = DataFetcher::new();
    let state = Arc::new(AppState { fetcher });

    let cors = CorsLayer::new()
        .allow_origin(Any)
        .allow_methods(Any)
        .allow_headers(Any);

    let app = routes::api_router()
        .fallback_service(ServeDir::new("frontend/dist").fallback(
            tower_http::services::ServeFile::new("frontend/dist/index.html"),
        ))
        .with_state(state)
        .layer(cors);

    let port = env::var("PORT")
        .ok()
        .and_then(|value| value.parse::<u16>().ok())
        .unwrap_or(3000);
    let bind_addr = format!("0.0.0.0:{port}");
    let listener = tokio::net::TcpListener::bind(&bind_addr)
        .await
        .unwrap_or_else(|_| panic!("Failed to bind to {bind_addr}"));

    println!("Server running on http://localhost:{port}");

    axum::serve(listener, app)
        .await
        .expect("Server error");
}
