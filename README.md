# Options Pricer

Rust-based options pricing system implementing quantitative finance models. Personal learning project for Rust and derivatives pricing.

## Pricing Methods

### European Options

| Method                                      | Model                           | Applies To                                     |
| ------------------------------------------- | ------------------------------- | ---------------------------------------------- |
| Black-Scholes (analytical)                  | Lognormal GBM                   | Vanilla call/put on any equity/ETF/index       |
| Black-Scholes spreads                       | Lognormal GBM                   | Straddle, strangle, butterfly, and 10+         |
| Monte Carlo + antithetic variates           | GBM                             | Vanilla call/put and spreads                   |
| Binomial CRR tree                           | Lognormal GBM                   | Vanilla call/put on any equity/ETF/index       |
| **Heston model (COS Fourier method)** | **Stochastic volatility** | **Vanilla call/put on any equity/index** |

### American Options

| Method                                   | Model                           | Applies To                                                      |
| ---------------------------------------- | ------------------------------- | --------------------------------------------------------------- |
| Black's approximation                    | Lognormal GBM                   | Single call/put on dividend-paying stock                        |
| Binomial CRR tree (backward induction)   | Lognormal GBM                   | Call/put on any equity/ETF/index                                |
| PDE penalty solver (Rannacher smoothing) | Lognormal GBM                   | Call/put and non-calendar spreads                               |
| **CTMC Algorithm**                 | **Stochastic volatility** | **Call/put on any equity; spreads via leg decomposition** |

### Exotics

| Method                                 | Applies To                            |
| -------------------------------------- | ------------------------------------- |
| Convertible bond (NPV + Black-Scholes) | Hybrid debt-equity instruments        |
| Chooser option (Black-Scholes)         | Options with deferred call/put choice |

### Volatility Modeling

| Method                                  | Output                                                |
| --------------------------------------- | ----------------------------------------------------- |
| Historical log-return vol               | Annualized vol from adjusted daily prices             |
| EMA volatility                          | Exponentially weighted rolling 30-day vol             |
| VIX correlation blending                | corr*VIX + (1-corr)*EMA_vol over 120 days             |
| Implied vol extraction                  | OI-weighted IV from near-ATM chain contracts          |
| **Heston calibration (COS + LM)** | **Fit v0, v_bar, kappa, sigma, rho to surface** |

## Heston Stochastic Volatility Pipeline

The full SV pipeline for American option pricing:

1. **Calibrate** — fetch the live option chain, build a `CalibrationDataset`, and run `calibrate_heston` with `HestonEuropeanPricer` (COS Fourier method, microseconds per price) to recover the five Heston parameters from the implied vol surface.
2. **Price** — pass the calibrated parameters to `price_american_option_heston` (or the `_put`/`_call` convenience wrappers), which solves the backward integral equation via CTMC (Algorithm 3.1) using either dense Padé matrix exponential (NM ≤ 600) or Krylov subspace approximation (NM > 600). Dispatch between put and call paths is resolved at compile time via the `Options` enum — no runtime branching in the inner loops.
3. **Spreads** — pass an `OptionSpreads` to `price_american_spread_heston`; each leg is priced independently and combined with its weight.

The calibration uses vega-weighted least squares with a Feller-condition soft penalty and Tikhonov anchor regularization for temporal stability.

## Greeks

Analytical Greeks (Δ, Γ, θ, ρ, vega) are available for all Black-Scholes-based pricers, including dividend adjustments. Supported for single options and option spreads.

## API Connectivity

Live market data is pulled automatically from:

- **Finnhub** — real-time spot price (15-min delayed Alpha Vantage as fallback)
- **Alpha Vantage** — adjusted daily historical prices (premium) and options chain
- **FRED** — VIX time series for correlation-blended vol

Keys are loaded from `api_keys.json` at the project root (not committed). If connectivity is not needed, supply vol and spot directly.

## SLV Models (CTMC Engine)

The CTMC engine supports six stochastic-local volatility models via the `SLVModelVariant` enum:

| Model        | Reference                 |
| ------------ | ------------------------- |
| Heston       | Heston 1993               |
| 4/2          | Grasselli 2017            |
| AlphaHyper   | Da Fonseca & Martini 2016 |
| SABR         | Hagan et al. 2002         |
| HestonSABR   | Hybrid                    |
| QuadraticSLV | Lipton 2002               |

Calibration for non-Heston SLV models requires separate parameter fitting (not yet automated), separate project on my github incoming using deep calibration. 

## Project Status

See [TODO.md](TODO.md) for current tasks and roadmap.

## AI Disclaimer

Tests in the code were written with Claude assistance. Main pricing logic and mathematical implementation were not. Safety checks via AI GitHub agents run regularly. Front-end was vibe coded.
