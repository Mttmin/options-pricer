mod error;
mod models;
mod routes;
mod state;

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

    let listener = tokio::net::TcpListener::bind("0.0.0.0:3000")
        .await
        .expect("Failed to bind to port 3000");

    println!("Server running on http://localhost:3000");

    axum::serve(listener, app)
        .await
        .expect("Server error");
}
