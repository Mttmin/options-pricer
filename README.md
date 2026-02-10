# Options Pricer

Coding a fully fledged options pricer in rust for options, futures and exotics. Serves as a personal learning project for Rust.

## Current functionalities

Currently supports the following pricing algorithms:

- Black-Scholes (European options, spreads, exotics)
- Monte Carlo (European options and spreads)
- Binomial CRR tree (European and American options; spreads via leg aggregation)
- PDE penalty solver (American options and non-calendar spreads)
- Black's American approximation (single options)

Includes Greeks calculations ($\Delta, \Gamma, \theta, \rho, vega$) for Europeans and Option Spreads

APIs call automatically for most information on the underlying

Front-end offers some functionnalities and exploration for options research but is still work in progress

## Project Status

See [TODO.md](TODO.md) for current tasks and roadmap.

## API Connectivity

As implemented right now, the pricer uses Alpha Vantage API for volatility estimates and Finnhub API for real-time stock price.

To get the vol estimates from historical data, and to avoid any problems from stock splits or dividends, the fetcher uses the adjusted daily endpoint, which requires a premium subscription. In case of a Finnhub problem, the fallback is to use the 15 minutes delayed Alpha Vantage global quote.

In case the online connectivity doesn't interest you, simply input your own vol and spot.

## AI Disclaimer

Tests in the code, were written through claude (Opus 4.5) but main logic was not. Checks on safety are also run by AI Github agents regularly for things I could have missed.

As a back-end/math oriented programmer, front-end dev is not my speciality and as such the frontend was vibe coded
