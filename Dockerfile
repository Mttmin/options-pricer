FROM rust:1.92-bookworm AS rust-builder
WORKDIR /app
ENV CARGO_REGISTRIES_CRATES_IO_PROTOCOL=sparse
RUN apt-get update \
    && apt-get install -y pkg-config libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Cache dependencies first
COPY Cargo.toml Cargo.lock ./
COPY options/Cargo.toml options/Cargo.toml
COPY numerical-methods/Cargo.toml numerical-methods/Cargo.toml
COPY cli/Cargo.toml cli/Cargo.toml
COPY web/Cargo.toml web/Cargo.toml
COPY options/src options/src
COPY numerical-methods/src numerical-methods/src
COPY cli/src cli/src
COPY web/src web/src
RUN cargo fetch

# Build the web binary
COPY . .
RUN cargo build --release --bin web

FROM node:20-bookworm-slim AS frontend-builder
WORKDIR /app/frontend
COPY frontend/package.json frontend/package-lock.json ./
RUN npm ci
COPY frontend/ ./
RUN npm run build

FROM debian:bookworm-slim AS runtime
WORKDIR /app
RUN apt-get update \
    && apt-get install -y ca-certificates libssl3 \
    && rm -rf /var/lib/apt/lists/*
ENV PORT=3001
COPY --from=rust-builder /app/target/release/web /app/web
COPY --from=frontend-builder /app/frontend/dist /app/frontend/dist
EXPOSE 3001
CMD ["/app/web"]
