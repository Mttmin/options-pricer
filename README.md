# Options Pricer

Coding a fully fledged options pricer in rust for options, futures and exotics. Serves as a personal learning project for Rust.

## Current functionalities

Currently supports Black-Scholes and Monte-Carlo pricing of European Calls/Puts and frequent option spreads, along with Binomial Trees for Americans and all options.

Includes Greeks calculations ($\Delta, \Gamma, \theta, \rho, vega$) for Europeans and Option Spreads

## Project Status

See [TODO.md](TODO.md) for current tasks and roadmap.

## API Connectivity

As implemented right now, the pricer uses Alpha Vantage API for volatility estimates and Finnhub API for real-time stock price.

To get the vol estimates from historical data, and to avoid any problems from stock splits or dividends, the fetcher uses the adjusted daily endpoint, which requires a premium subscription. In case of a Finnhub problem, the fallback is to use the 15 minutes delayed Alpha Vantage global quote.

In case the online connectivity doesn't interest you, simply input your own vol and spot.
