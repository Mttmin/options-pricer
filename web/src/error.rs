use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use serde_json::json;

pub enum AppError {
    BadRequest(String),
    NotFound(String),
    Internal(String),
}

impl IntoResponse for AppError {
    fn into_response(self) -> Response {
        let (status, message) = match self {
            AppError::BadRequest(msg) => (StatusCode::BAD_REQUEST, msg),
            AppError::NotFound(msg) => (StatusCode::NOT_FOUND, msg),
            AppError::Internal(msg) => (StatusCode::INTERNAL_SERVER_ERROR, msg),
        };

        let body = axum::Json(json!({ "error": message }));
        (status, body).into_response()
    }
}

impl From<Box<dyn std::error::Error + Send + Sync>> for AppError {
    fn from(err: Box<dyn std::error::Error + Send + Sync>) -> Self {
        let msg = err.to_string();
        if msg.contains("rate limit") {
            AppError::Internal(msg)
        } else if msg.contains("Invalid") || msg.contains("missing") {
            AppError::NotFound(msg)
        } else {
            AppError::Internal(msg)
        }
    }
}
