use chrono::{DateTime, Duration, NaiveDate, Utc};
use reqwest;
use serde::Deserialize;
use std::collections::HashMap;
use std::sync::Arc;
use tokio::sync::RwLock;
use crate::options_chain::{OptionChain, OptionContract, OptionType, OptionsEndpoint};
use options::black_scholes::calculate_iv;

/// Configuration struct for API keys loaded from api_keys.json
#[derive(Debug, Deserialize)]
pub struct ApiConfig {
    #[serde(rename = "Alpha_vantage")]
    pub alpha_vantage: String,
    #[serde(rename = "Finnhub")]
    pub finnhub: String,
    #[serde(rename = "Finnhub_secret")]
    pub finnhub_secret: String,
    #[serde(rename = "FRED")]
    pub fred: String,
}

/// Load API keys from api_keys.json in the project root
pub fn load_api_keys() -> Result<ApiConfig, Box<dyn std::error::Error + Send + Sync>> {
    let possible_paths = vec![
        "api_keys.json",    // Current directory
        "../api_keys.json", // Parent directory (when running from cli/)
    ];

    let mut config_content: Option<String> = None;
    let mut last_error: String = String::new();

    for path in &possible_paths {
        match std::fs::read_to_string(path) {
            Ok(content) => {
                config_content = Some(content);
                break;
            }
            Err(e) => last_error = format!("{}: {}", path, e),
        }
    }

    let config_content = config_content
        .ok_or_else(|| format!("Failed to read api_keys.json. Last error: {}", last_error))?;

    let config: ApiConfig = serde_json::from_str(&config_content)
        .map_err(|e| format!("Failed to parse api_keys.json: {}", e))?;

    Ok(config)
}

#[derive(Debug, Deserialize)]
struct FinnhubQuote {
    c: f64, // current price
    #[serde(rename = "h")]
    _h: f64, // high
    #[serde(rename = "l")]
    _l: f64, // low
    #[serde(rename = "o")]
    _o: f64, // open
    #[serde(rename = "pc")]
    _pc: f64, // previous close
    #[serde(rename = "t")]
    _t: i64, // timestamp
}

#[derive(Clone, Debug)]
pub struct MarketData {
    pub spot_price: f64,
    pub volatility: f64,
    pub last_price_update: DateTime<Utc>,
    pub vol_calculation_date: NaiveDate,
    pub yesterday_adjusted_close: f64,
    pub corporate_action_detected: bool,
    /// Implied volatility from option chain (None if not fetched)
    pub implied_volatility: Option<f64>,
    /// Dividend yield as decimal (e.g., 0.0044 for 0.44%)
    pub dividend_yield: Option<f64>,
}

#[derive(Clone)]
pub struct DataFetcher {
    alpha_vantage_key: String,
    finnhub_key: String,
    client: reqwest::Client,
    cache: Arc<RwLock<HashMap<String, MarketData>>>,
    fred_key: String,
    options_cache: Arc<RwLock<HashMap<String, OptionChain>>>,
    sofr_cache: Arc<RwLock<Option<(DateTime<Utc>, f64)>>>,
    fundamentals_cache: Arc<RwLock<HashMap<String, (DateTime<Utc>, CompanyFundamentals)>>>,
}

#[derive(Debug, Deserialize)]
struct Observation {
    date: String,
    value: String,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CompanyFundamentals {
    #[serde(rename = "DividendYield")]
    pub dividend_yield: Option<String>,
    #[serde(rename = "DividendPerShare")]
    pub dividend_per_share: Option<String>,
    #[serde(rename = "MarketCapitalization")]
    pub market_cap: Option<String>,
    #[serde(rename = "DividendDate")]
    pub dividend_date: Option<String>,
    #[serde(rename = "ExDividendDate")]
    pub ex_dividend_date: Option<String>,
}

#[derive(Debug, Deserialize)]
struct FredResponse {
    observations: Vec<Observation>,
}

#[derive(Debug, Clone)]
pub struct DividendInfo {
    pub time_to_dividend: f64,
    pub dividend_amount: f64,
    pub dividend_yield: f64,
}

fn calculate_time_to_dividend(
    dividend_date_str: &str
) -> Result<f64, Box<dyn std::error::Error + Send + Sync>> {
    let dividend_date = NaiveDate::parse_from_str(dividend_date_str, "%Y-%m-%d")?;
    let today = Utc::now().date_naive();

    let days_diff = (dividend_date - today).num_days();

    if days_diff < 0 {
        return Err("Dividend date is in the past".into());
    }

    Ok(days_diff as f64 / 365.25)
}

impl Default for DataFetcher {
    fn default() -> Self {
        Self::new()
    }
}

impl DataFetcher {
    pub fn new() -> Self {
        // get the API keys
        let config = load_api_keys().expect("failed to load API keys");
        let alpha_vantage = config.alpha_vantage;
        let finnhub = config.finnhub;
        let fred_key = config.fred;
        Self {
            alpha_vantage_key: alpha_vantage,
            finnhub_key: finnhub,
            fred_key,
            client: reqwest::Client::new(),
            cache: Arc::new(RwLock::new(HashMap::new())),
            options_cache: Arc::new(RwLock::new(HashMap::new())),
            sofr_cache: Arc::new(RwLock::new(None)),
            fundamentals_cache: Arc::new(RwLock::new(HashMap::new())),
        }
    }

    async fn get_current_price_finnhub(
        &self,
        symbol: &str,
    ) -> Result<FinnhubQuote, Box<dyn std::error::Error + Send + Sync>> {
        let url = format!(
            "https://finnhub.io/api/v1/quote?symbol={}",
            symbol
        );

        let response = self
            .client
            .get(&url)
            .header("X-Finnhub-Token", &self.finnhub_key)
            .send()
            .await?
            .json::<FinnhubQuote>()
            .await?;
        Ok(response)
    }

    #[allow(dead_code)]
    /// returns 15 min delayed price from Alpha Vantage
    async fn get_current_price_alpha_vantage(
        &self,
        symbol: &str,
    ) -> Result<f64, Box<dyn std::error::Error + Send + Sync>> {
        let url = format!(
            "https://www.alphavantage.co/query?function=GLOBAL_QUOTE&symbol={}&entitlement=delayed&apikey={}",
            symbol, self.alpha_vantage_key
        );

        let response = self
            .client
            .get(&url)
            .send()
            .await?
            .json::<serde_json::Value>()
            .await?;

        // Find the "Global Quote" key (may have suffix like "- DATA DELAYED BY 15 MINUTES")
        let global_quote_key = response
            .as_object()
            .and_then(|obj| obj.keys().find(|k| k.starts_with("Global Quote")))
            .ok_or("Invalid response format - missing Global Quote key")?;

        let price_str = response[global_quote_key]["05. price"]
            .as_str()
            .ok_or("Invalid response format - missing price")?;
        let price = price_str.parse::<f64>()?;

        Ok(price)
    }

    /// Get current price with fallback mechanism
    /// Tries Finnhub first (real-time), falls back to Alpha Vantage (15-min delayed) if Finnhub fails
    async fn get_current_price_with_fallback(
        &self,
        symbol: &str,
    ) -> Result<f64, Box<dyn std::error::Error + Send + Sync>> {
        // Try Alpha Vantage first (user has premium key = real-time data)
        match self.get_current_price_alpha_vantage(symbol).await {
            Ok(price) => Ok(price),
            Err(av_error) => {
                // Alpha Vantage failed, fall back to Finnhub
                eprintln!(
                    "WARNING: Alpha Vantage quote failed for {} ({}). Falling back to Finnhub real-time.",
                    symbol, av_error
                );

                // Try Finnhub fallback
                self.get_current_price_finnhub(symbol).await.map(|q| q.c)
            }
        }
    }

    pub async fn get_historical_data(
        &self,
        symbol: &str,
        lookback_days: usize,
    ) -> Result<(Vec<(String, f64)>, f64), Box<dyn std::error::Error + Send + Sync>> {
        let url = format!(
            "https://www.alphavantage.co/query?function=TIME_SERIES_DAILY_ADJUSTED&symbol={}&outputsize=full&apikey={}",
            symbol, self.alpha_vantage_key
        );

        let response = self
            .client
            .get(&url)
            .send()
            .await?
            .json::<serde_json::Value>()
            .await?;

        // Check for rate limit or error responses
        if let Some(note) = response.get("Note") {
            return Err(format!(
                "Alpha Vantage rate limit exceeded: {}",
                note.as_str().unwrap_or("unknown")
            )
            .into());
        }
        if let Some(error) = response.get("Error Message") {
            return Err(format!(
                "Alpha Vantage error: {}",
                error.as_str().unwrap_or("unknown")
            )
            .into());
        }
        if let Some(info) = response.get("Information") {
            return Err(format!(
                "Alpha Vantage API info: {}",
                info.as_str().unwrap_or("unknown")
            )
            .into());
        }

        let time_series = response["Time Series (Daily)"].as_object().ok_or_else(|| {
            format!(
                "Invalid response format - missing 'Time Series (Daily)'. Response keys: {:?}",
                response.as_object().map(|o| o.keys().collect::<Vec<_>>())
            )
        })?;

        let mut prices: Vec<(String, f64)> = time_series
            .iter()
            .filter_map(|(date, data)| {
                let price = data["5. adjusted close"].as_str()?.parse::<f64>().ok()?;
                Some((date.clone(), price))
            })
            .collect();

        prices.sort_by(|a, b| b.0.cmp(&a.0)); // Most recent first

        if prices.is_empty() {
            return Err("No historical data available".into());
        }

        let yesterday_close = prices[0].1; // Most recent adjusted close

        prices.truncate(lookback_days);
        prices.reverse(); // Oldest first for volatility calculation

        Ok((prices, yesterday_close))
    }

    pub fn calculate_volatility(prices: &[(String, f64)]) -> f64 {
        if prices.len() < 2 {
            return 0.0;
        }

        let log_returns: Vec<f64> = prices
            .windows(2)
            .map(|pair| (pair[1].1 / pair[0].1).ln())
            .collect();

        // Need at least 2 log returns to calculate variance
        if log_returns.len() < 2 {
            return 0.0;
        }

        let mean = log_returns.iter().sum::<f64>() / log_returns.len() as f64;

        let variance = log_returns.iter().map(|r| (r - mean).powi(2)).sum::<f64>()
            / (log_returns.len() - 1) as f64;

        variance.sqrt() * (252.0_f64).sqrt() // correction for trading days
    }

    /// Fetch the most recent VIX closing value from FRED API
    ///
    /// Uses the VIXCLS series from the Federal Reserve Economic Data API.
    /// Returns VIX value as a percentage
    pub async fn fetch_fred_vix(
        &self,
        num_days: u16,
    ) -> Result<Vec<f64>, Box<dyn std::error::Error + Send + Sync>> {
        let url = format!(
            "https://api.stlouisfed.org/fred/series/observations?series_id=VIXCLS&api_key={}&file_type=json&limit={}&sort_order=desc",
            &self.fred_key, &num_days
        );

        let response = self
            .client
            .get(&url)
            .send()
            .await?
            .json::<FredResponse>()
            .await?;

        // Check if we received any observations
        if response.observations.is_empty() {
            return Err("No VIX observations found in FRED response".into());
        }

        let mut vix_values: Vec<f64> = Vec::new();
        for i in 0..num_days {
            // With &sort_order=desc, latest observation are at the beginning
            let observation = &response.observations[i as usize];

            // Parse the value
            let vix_value: f64 = observation.value.parse().map_err(|e| {
                format!(
                    "Failed to parse VIX value '{}' from date {}: {}",
                    observation.value, observation.date, e
                )
            })?;
            // // Check : VIX typically ranges 10-80, cannot go under 0, but can spike to 150+ in extreme markets
            // if vix_value <= 0.0 || vix_value > 70.0 {
            //     return Err(format!(
            //         "VIX value {:.2} from {} is outside expected range (0-70)",
            //         vix_value, observation.date
            //     )
            //     .into());
            // }
            vix_values.push(vix_value);
        }

        Ok(vix_values)
    }

    /// Fetch the current SOFR (Secured Overnight Financing Rate) from FRED API
    ///
    /// Uses the SOFR30DAYAVG series for stability. Returns rate as a decimal (e.g., 0.0531 for 5.31%)
    /// Caches result for 24 hours to avoid excessive API calls
    pub async fn fetch_fred_sofr(&self) -> Result<f64, Box<dyn std::error::Error + Send + Sync>> {
        // Check cache (24-hour TTL)
        let cache = self.sofr_cache.read().await;
        if let Some((timestamp, rate)) = *cache {
            if Utc::now() - timestamp < Duration::hours(24) {
                return Ok(rate);
            }
        }
        drop(cache);

        // Fetch from FRED using SOFR30DAYAVG series (more stable than daily SOFR)
        let url = format!(
            "https://api.stlouisfed.org/fred/series/observations?series_id=SOFR30DAYAVG&api_key={}&file_type=json&limit=1&sort_order=desc",
            &self.fred_key
        );

        let response = self
            .client
            .get(&url)
            .send()
            .await?
            .json::<FredResponse>()
            .await?;

        // Check if we received any observations
        if response.observations.is_empty() {
            return Err("No SOFR observations found in FRED response".into());
        }

        // Parse rate from latest observation
        let rate_str = &response.observations[0].value;
        let rate_pct: f64 = rate_str.parse().map_err(|e| {
            format!("Failed to parse SOFR value '{}': {}", rate_str, e)
        })?;

        // Convert percentage to decimal (e.g., 5.31 -> 0.0531)
        let rate = rate_pct / 100.0;

        // Validate range (SOFR typically 0-15%)
        if rate < 0.0 || rate > 0.15 {
            return Err(format!("SOFR rate {:.4} outside expected range (0-0.15)", rate).into());
        }

        // Update cache
        let mut cache = self.sofr_cache.write().await;
        *cache = Some((Utc::now(), rate));

        Ok(rate)
    }

    /// Fetch company fundamentals from Alpha Vantage OVERVIEW endpoint
    ///
    /// Returns dividend yield, dividend per share, and market cap
    /// Caches result for 7 days as fundamentals don't change frequently
    pub async fn fetch_company_fundamentals(
        &self,
        symbol: &str,
    ) -> Result<CompanyFundamentals, Box<dyn std::error::Error + Send + Sync>> {
        // Check cache (7-day TTL)
        let cache = self.fundamentals_cache.read().await;
        if let Some((timestamp, fundamentals)) = cache.get(symbol) {
            if Utc::now() - *timestamp < Duration::days(7) {
                return Ok(fundamentals.clone());
            }
        }
        drop(cache);

        // Fetch from Alpha Vantage OVERVIEW endpoint
        let url = format!(
            "https://www.alphavantage.co/query?function=OVERVIEW&symbol={}&apikey={}",
            symbol, &self.alpha_vantage_key
        );

        let response = self.client.get(&url).send().await?;

        // Check for error responses
        let text = response.text().await?;
        let json: serde_json::Value = serde_json::from_str(&text)?;

        if let Some(note) = json.get("Note") {
            return Err(format!("Alpha Vantage rate limit: {}", note).into());
        }

        if let Some(error) = json.get("Error Message") {
            return Err(format!("Alpha Vantage error: {}", error).into());
        }

        // Parse fundamentals
        let fundamentals: CompanyFundamentals = serde_json::from_value(json)?;

        // Update cache
        let mut cache = self.fundamentals_cache.write().await;
        cache.insert(symbol.to_string(), (Utc::now(), fundamentals.clone()));

        Ok(fundamentals)
    }

    /// Fetch dividend information for a given symbol
    ///
    /// Returns processed dividend data including time to next dividend and amount.
    /// Leverages the 7-day cache from fetch_company_fundamentals.
    pub async fn get_dividend_info(
        &self,
        symbol: &str,
        spot_price: f64,
    ) -> Result<Option<DividendInfo>, Box<dyn std::error::Error + Send + Sync>> {
        let fundamentals = self.fetch_company_fundamentals(symbol).await?;

        let date_str = fundamentals.ex_dividend_date
            .as_ref()
            .or(fundamentals.dividend_date.as_ref());

        let div_per_share_str = fundamentals.dividend_per_share.as_ref();

        match (date_str, div_per_share_str) {
            (Some(date), Some(amount_str)) => {
                let dividend_amount = match amount_str.parse::<f64>() {
                    Ok(amt) if amt > 0.0 => amt,
                    _ => return Ok(None),
                };

                let time_to_dividend = match calculate_time_to_dividend(date) {
                    Ok(t) => t,
                    Err(_) => return Ok(None),
                };

                let dividend_yield = fundamentals.dividend_yield
                    .as_ref()
                    .and_then(|s| s.parse::<f64>().ok())
                    .unwrap_or(dividend_amount / spot_price);

                Ok(Some(DividendInfo {
                    time_to_dividend,
                    dividend_amount,
                    dividend_yield,
                }))
            }
            _ => Ok(None),
        }
    }

    fn detect_corporate_action(
        current_price: f64,
        yesterday_adjusted: f64,
        threshold: f64,
    ) -> bool {
        let ratio = current_price / yesterday_adjusted;
        (ratio - 1.0).abs() > threshold
    }

    async fn update_single_symbol(
        &self,
        symbol: &str,
        lookback_days: usize,
    ) -> Result<MarketData, Box<dyn std::error::Error + Send + Sync>> {
        let now = Utc::now();
        let today = now.date_naive();

        // Check what we need to update
        let (needs_vol_update, needs_price_update, existing_data) = {
            let cache = self.cache.read().await;
            if let Some(data) = cache.get(symbol) {
                let needs_vol = data.vol_calculation_date != today;
                let needs_price = now - data.last_price_update > Duration::hours(1);
                (needs_vol, needs_price, Some(data.clone()))
            } else {
                (true, true, None)
            }
        };

        let mut data = existing_data.unwrap_or(MarketData {
            spot_price: 0.0,
            volatility: 0.0,
            last_price_update: DateTime::UNIX_EPOCH,
            vol_calculation_date: NaiveDate::from_ymd_opt(1970, 1, 1).unwrap(),
            yesterday_adjusted_close: 0.0,
            corporate_action_detected: false,
            implied_volatility: None,
            dividend_yield: None,
        });

        // Fetch from both APIs in parallel if we need both
        if needs_vol_update && needs_price_update {
            // println!(
            //     "Fetching both price and volatility for {} in parallel...",
            //     symbol
            // );

            let hist_future = self.get_historical_data(symbol, lookback_days);
            let price_future = self.get_current_price_with_fallback(symbol);
            let fundamentals_future = self.fetch_company_fundamentals(symbol);

            // Execute all requests in parallel
            let (hist_result, price_result, fundamentals_result) =
                tokio::join!(hist_future, price_future, fundamentals_future);

            // Process historical data
            let (prices, yesterday_close) = hist_result?;
            data.volatility = Self::calculate_volatility(&prices);
            data.yesterday_adjusted_close = yesterday_close;
            data.vol_calculation_date = today;

            // Process current price
            let current_price = price_result?;

            let corporate_action =
                Self::detect_corporate_action(current_price, data.yesterday_adjusted_close, 0.05);

            data.corporate_action_detected = corporate_action;

            if corporate_action {
                data.spot_price = data.yesterday_adjusted_close;
                eprintln!(
                    "WARNING: Corporate action detected for {}. Current: {:.2}, Yesterday: {:.2}. Using yesterday's adjusted close.",
                    symbol, current_price, data.yesterday_adjusted_close
                );
            } else {
                data.spot_price = current_price;
            }

            data.last_price_update = now;

            // Process fundamentals (non-critical, don't fail if it errors)
            if let Ok(fundamentals) = fundamentals_result {
                data.dividend_yield = fundamentals.dividend_yield
                    .and_then(|s| s.parse::<f64>().ok())
                    .filter(|&v| v > 0.0);
            }
        } else if needs_vol_update {
            // Only need historical data
            // eprintln!(
            //     "Fetching only volatility for {} (already have recent price)...",
            //     symbol
            // );

            let (prices, yesterday_close) = self.get_historical_data(symbol, lookback_days).await?;
            data.volatility = Self::calculate_volatility(&prices);
            data.yesterday_adjusted_close = yesterday_close;
            data.vol_calculation_date = today;
        } else if needs_price_update {
            // Only need current price
            // eprintln!(
            //     "Fetching only price for {} (already have today's volatility)...",
            //     symbol
            // );

            let current_price = self.get_current_price_with_fallback(symbol).await?;

            let corporate_action =
                Self::detect_corporate_action(current_price, data.yesterday_adjusted_close, 0.05);

            data.corporate_action_detected = corporate_action;

            if corporate_action {
                data.spot_price = data.yesterday_adjusted_close;
                eprintln!(
                    "WARNING: Corporate action detected for {}. Current: {:.2}, Yesterday: {:.2}. Using yesterday's adjusted close.",
                    symbol, current_price, data.yesterday_adjusted_close
                );
            } else {
                data.spot_price = current_price;
            }

            data.last_price_update = now;
        } else {
            // eprintln!("Using cached data for {} (no updates needed)", symbol);
        }

        // Update cache
        {
            let mut cache = self.cache.write().await;
            cache.insert(symbol.to_string(), data.clone());
        }

        Ok(data)
    }

    #[allow(dead_code)]
    async fn update_multiple_symbols(
        &self,
        symbols: &[String],
        lookback_days: usize,
    ) -> HashMap<String, Result<MarketData, String>> {
        // Fetch all symbols in parallel
        let futures: Vec<_> = symbols
            .iter()
            .map(|symbol| {
                let symbol = symbol.clone();
                async move {
                    let result = self.update_single_symbol(&symbol, lookback_days).await;
                    (symbol.clone(), result.map_err(|e| e.to_string()))
                }
            })
            .collect();

        let results = futures::future::join_all(futures).await;
        results.into_iter().collect()
    }

    pub async fn get_market_data(&self, symbol: &str) -> Option<MarketData> {
        let cache = self.cache.read().await;
        cache.get(symbol).cloned()
    }

    pub async fn update_symbol(
        &self,
        symbol: &str,
        lookback_days: usize,
    ) -> Result<MarketData, Box<dyn std::error::Error + Send + Sync>> {
        self.update_single_symbol(symbol, lookback_days).await
    }

    /// Fetch option chain data from Alpha Vantage.
    ///
    /// Supports both HISTORICAL_OPTIONS (all tiers) and REALTIME_OPTIONS (premium).
    /// Returns the full chain with per-contract IV and Greeks.
    /// Automatically computes summary_implied_volatility on the result.
    pub async fn fetch_option_chain(
        &self,
        symbol: &str,
        endpoint: OptionsEndpoint,
    ) -> Result<OptionChain, Box<dyn std::error::Error + Send + Sync>> {
        let cache_key = match endpoint {
            OptionsEndpoint::Historical { date } => match date {
                Some(d) => format!("{}:hist:{}", symbol, d),
                None => format!("{}:hist:latest", symbol),
            },
            OptionsEndpoint::Realtime => format!("{}:realtime", symbol),
        };

        // Check cache (15-minute TTL)
        {
            let cache = self.options_cache.read().await;
            if let Some(cached) = cache.get(&cache_key) {
                let now = Utc::now();
                if (now - cached.fetch_timestamp).num_minutes() < 15 {
                    return Ok(cached.clone());
                }
            }
        }

        let url = match endpoint {
            OptionsEndpoint::Historical { date } => {
                let mut u = format!(
                    "https://www.alphavantage.co/query?function=HISTORICAL_OPTIONS&symbol={}&require_greeks=true&apikey={}",
                    symbol, self.alpha_vantage_key
                );
                if let Some(d) = date {
                    u.push_str(&format!("&date={}", d));
                }
                u
            }
            OptionsEndpoint::Realtime => {
                format!(
                    "https://www.alphavantage.co/query?function=REALTIME_OPTIONS&symbol={}&require_greeks=true&apikey={}",
                    symbol, self.alpha_vantage_key
                )
            }
        };

        let response = self
            .client
            .get(&url)
            .send()
            .await?
            .json::<serde_json::Value>()
            .await?;

        // Standard Alpha Vantage error checking
        if let Some(note) = response.get("Note") {
            return Err(format!(
                "Alpha Vantage rate limit exceeded: {}",
                note.as_str().unwrap_or("unknown")
            )
            .into());
        }
        if let Some(error) = response.get("Error Message") {
            return Err(format!(
                "Alpha Vantage error: {}",
                error.as_str().unwrap_or("unknown")
            )
            .into());
        }
        if let Some(info) = response.get("Information") {
            return Err(format!(
                "Alpha Vantage API info (possible premium-only endpoint): {}",
                info.as_str().unwrap_or("unknown")
            )
            .into());
        }

        let chain = Self::parse_options_response(symbol, &response)?;

        // Cache the result
        {
            let mut cache = self.options_cache.write().await;
            cache.insert(cache_key, chain.clone());
        }

        Ok(chain)
    }

    /// Parse the Alpha Vantage options JSON response into OptionChain.
    fn parse_options_response(
        symbol: &str,
        response: &serde_json::Value,
    ) -> Result<OptionChain, Box<dyn std::error::Error + Send + Sync>> {
        let data_array = response
            .get("data")
            .and_then(|d| d.as_array())
            .ok_or("Invalid response format - missing 'data' array")?;

        let mut contracts: Vec<OptionContract> = Vec::with_capacity(data_array.len());

        for item in data_array {
            match Self::parse_option_contract(item) {
                Ok(contract) => contracts.push(contract),
                Err(e) => {
                    eprintln!("WARNING: Skipping contract due to parse error: {}", e);
                }
            }
        }

        let underlying_price = Self::estimate_underlying_from_chain(&contracts);

        // Try to extract data date from response metadata
        let data_date = response
            .get("meta")
            .and_then(|m| m.get("date"))
            .and_then(|d| d.as_str())
            .and_then(|s| NaiveDate::parse_from_str(s, "%Y-%m-%d").ok())
            .unwrap_or_else(|| Utc::now().date_naive());

        let mut chain = OptionChain {
            symbol: symbol.to_string(),
            underlying_price,
            fetch_timestamp: Utc::now(),
            data_date,
            contracts,
            summary_implied_volatility: None,
        };

        // Fallback IV Calculation
        for c in chain.contracts.iter_mut() {
            if c.implied_volatility == 0.0 || c.implied_volatility.abs() < 1e-6 {
                let t = (c.expiration - chain.data_date).num_days() as f64 / 365.0;
                let r = 0.05;
                let q = 0.0;
                let is_call = match c.option_type {
                    OptionType::Call => true,
                    OptionType::Put => false,
                };

                let price = if c.mark > 0.0 {
                    c.mark
                } else if c.bid > 0.0 && c.ask > 0.0 {
                    (c.bid + c.ask) / 2.0
                } else {
                    c.last
                };

                if t > 0.0 && price > 0.0 {
                    c.implied_volatility = calculate_iv(
                        price,
                        chain.underlying_price,
                        c.strike,
                        t,
                        r,
                        q,
                        is_call,
                    );
                }
            }
        }

        chain.compute_summary_iv();

        Ok(chain)
    }

    /// Parse a single contract JSON object from the Alpha Vantage response.
    fn parse_option_contract(
        item: &serde_json::Value,
    ) -> Result<OptionContract, Box<dyn std::error::Error + Send + Sync>> {
        let parse_f64 = |key: &str| -> f64 {
            item.get(key)
                .and_then(|v| {
                    v.as_str()
                        .and_then(|s| s.parse::<f64>().ok())
                        .or_else(|| v.as_f64())
                })
                .unwrap_or(0.0)
        };

        let parse_optional_f64 = |key: &str| -> Option<f64> {
            item.get(key).and_then(|v| {
                if let Some(s) = v.as_str() {
                    if s == "N/A" || s.is_empty() {
                        None
                    } else {
                        s.parse::<f64>().ok()
                    }
                } else {
                    v.as_f64()
                }
            })
        };

        let parse_u64 = |key: &str| -> u64 {
            item.get(key)
                .and_then(|v| {
                    v.as_str()
                        .and_then(|s| s.parse::<u64>().ok())
                        .or_else(|| v.as_u64())
                })
                .unwrap_or(0)
        };

        let type_str = item
            .get("type")
            .and_then(|v| v.as_str())
            .unwrap_or("");

        let option_type = match type_str.to_lowercase().as_str() {
            "call" => OptionType::Call,
            "put" => OptionType::Put,
            other => return Err(format!("Unknown option type: '{}'", other).into()),
        };

        let expiration_str = item
            .get("expiration")
            .and_then(|v| v.as_str())
            .ok_or("Missing expiration field")?;

        let expiration = NaiveDate::parse_from_str(expiration_str, "%Y-%m-%d")
            .map_err(|e| format!("Failed to parse expiration '{}': {}", expiration_str, e))?;

        let in_the_money_str = item
            .get("in_the_money")
            .and_then(|v| v.as_str())
            .unwrap_or("FALSE");

        Ok(OptionContract {
            contract_id: item
                .get("contractID")
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string(),
            symbol: item
                .get("symbol")
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string(),
            expiration,
            strike: parse_f64("strike"),
            option_type,
            last: parse_f64("last"),
            mark: parse_f64("mark"),
            bid: parse_f64("bid"),
            ask: parse_f64("ask"),
            volume: parse_u64("volume"),
            open_interest: parse_u64("open_interest"),
            implied_volatility: parse_f64("implied_volatility"),
            delta: parse_optional_f64("delta"),
            gamma: parse_optional_f64("gamma"),
            theta: parse_optional_f64("theta"),
            vega: parse_optional_f64("vega"),
            rho: parse_optional_f64("rho"),
            in_the_money: in_the_money_str.eq_ignore_ascii_case("TRUE"),
        })
    }

    /// Estimate underlying price from the chain using the highest open-interest contract.
    fn estimate_underlying_from_chain(contracts: &[OptionContract]) -> f64 {
        if contracts.is_empty() {
            return 0.0;
        }
        contracts
            .iter()
            .max_by_key(|c| c.open_interest)
            .map(|c| c.strike)
            .unwrap_or(0.0)
    }

    /// Fetch option chain implied volatility and optionally update MarketData cache.
    ///
    /// This is a standalone volatility source. The caller decides whether to use
    /// historical vol, EMA vol, VIX-correlated vol, or this market IV.
    pub async fn fetch_implied_volatility(
        &self,
        symbol: &str,
        endpoint: OptionsEndpoint,
    ) -> Result<f64, Box<dyn std::error::Error + Send + Sync>> {
        let chain = self.fetch_option_chain(symbol, endpoint).await?;

        let iv = chain
            .summary_implied_volatility
            .ok_or("Could not compute summary implied volatility from chain")?;

        // Update the cached MarketData with IV if present
        {
            let mut cache = self.cache.write().await;
            if let Some(market_data) = cache.get_mut(symbol) {
                market_data.implied_volatility = Some(iv);
            }
        }

        Ok(iv)
    }

    /// Fetch option chain with a known underlying price (e.g. from a live Finnhub quote).
    /// Recomputes summary IV with the correct underlying.
    pub async fn fetch_option_chain_with_underlying(
        &self,
        symbol: &str,
        endpoint: OptionsEndpoint,
        underlying_price: f64,
    ) -> Result<OptionChain, Box<dyn std::error::Error + Send + Sync>> {
        let mut chain = self.fetch_option_chain(symbol, endpoint).await?;
        chain.underlying_price = underlying_price;
        chain.compute_summary_iv();
        Ok(chain)
    }

    /// Try REALTIME_OPTIONS first, fall back to HISTORICAL_OPTIONS on premium error.
    pub async fn fetch_option_chain_with_fallback(
        &self,
        symbol: &str,
    ) -> Result<OptionChain, Box<dyn std::error::Error + Send + Sync>> {
        match self
            .fetch_option_chain(symbol, OptionsEndpoint::Realtime)
            .await
        {
            Ok(chain) => Ok(chain),
            Err(e) => {
                let err_str = e.to_string();
                if err_str.contains("premium") || err_str.contains("Information") {
                    eprintln!(
                        "WARNING: REALTIME_OPTIONS not available for {} ({}). Falling back to HISTORICAL_OPTIONS.",
                        symbol, err_str
                    );
                    self.fetch_option_chain(symbol, OptionsEndpoint::Historical { date: None })
                        .await
                } else {
                    Err(e)
                }
            }
        }
    }

    /// Search for ticker symbols using Finnhub API
    /// Returns a vector of (symbol, description, type) tuples
    pub async fn search_symbols(
        &self,
        query: &str,
    ) -> Result<Vec<(String, String, String)>, Box<dyn std::error::Error + Send + Sync>> {
        let url = format!(
            "https://finnhub.io/api/v1/search?q={}",
            query
        );

        let response = self
            .client
            .get(&url)
            .header("X-Finnhub-Token", &self.finnhub_key)
            .send()
            .await?
            .json::<serde_json::Value>()
            .await?;

        let results = response["result"]
            .as_array()
            .map(|arr| {
                arr.iter()
                    .filter_map(|item| {
                        let symbol = item["symbol"].as_str()?.to_string();
                        let description = item["description"].as_str().unwrap_or("").to_string();
                        let security_type = item["type"].as_str().unwrap_or("").to_string();
                        // Filter to only US stocks for simplicity
                        if security_type == "Common Stock" || security_type.is_empty() {
                            Some((symbol, description, security_type))
                        } else {
                            None
                        }
                    })
                    .take(10)
                    .collect()
            })
            .unwrap_or_default();

        Ok(results)
    }

    /// Get historical prices for charting
    /// Returns a vector of (date, price) tuples, most recent first
    pub async fn get_historical_prices(
        &self,
        symbol: &str,
        days: usize,
    ) -> Result<Vec<(String, f64)>, Box<dyn std::error::Error + Send + Sync>> {
        let (mut prices, _) = self.get_historical_data(symbol, days.max(30)).await?;
        // get_historical_data returns oldest first after reversal, we want newest first
        prices.reverse();
        prices.truncate(days);
        Ok(prices)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;

    fn pretty_print_json(json: &Value, _label: &str) {
        println!("{}", serde_json::to_string_pretty(json).unwrap());
    }

    /// Helper to inspect JSON structure for debugging
    fn inspect_json_structure(json: &Value, path: String, depth: usize) {
        let indent = "  ".repeat(depth);
        match json {
            Value::Object(map) => {
                println!("{}{}Object with {} keys:", indent, path, map.len());
                for (key, value) in map {
                    let new_path = if path.is_empty() {
                        key.clone()
                    } else {
                        format!("{}.{}", path, key)
                    };
                    inspect_json_structure(value, new_path, depth + 1);
                }
            }
            Value::Array(arr) => {
                println!("{}{}Array with {} elements", indent, path, arr.len());
                if !arr.is_empty() {
                    inspect_json_structure(&arr[0], format!("{}[0]", path), depth + 1);
                }
            }
            Value::String(s) => {
                println!(
                    "{}{}String: {:?}",
                    indent,
                    path,
                    s.chars().take(50).collect::<String>()
                );
            }
            Value::Number(n) => {
                println!("{}{}Number: {}", indent, path, n);
            }
            Value::Bool(b) => {
                println!("{}{}Bool: {}", indent, path, b);
            }
            Value::Null => {
                println!("{}{}Null", indent, path);
            }
        }
    }

    /// Helper to load API keys for tests
    fn get_test_api_keys() -> ApiConfig {
        load_api_keys().expect("Failed to load API keys for testing")
    }

    #[tokio::test]
    #[ignore]
    async fn test_finnhub_raw_response() {
        let config = get_test_api_keys();
        let client = reqwest::Client::new();

        let symbols = vec!["AAPL", "MSFT", "SPY"];

        for symbol in symbols {
            let url = format!(
                "https://finnhub.io/api/v1/quote?symbol={}&token={}",
                symbol, config.finnhub
            );

            println!("\nFetching Finnhub quote for {}...", symbol);
            println!("URL: {}", url.replace(&config.finnhub, "***"));

            let response = client.get(&url).send().await.expect("Request failed");
            let status = response.status();
            println!("Status: {}", status);

            let json: Value = response.json().await.expect("Failed to parse as JSON");
            pretty_print_json(&json, &format!("Finnhub Response for {}", symbol));
            inspect_json_structure(&json, String::new(), 0);

            // Check for expected fields
            println!("\nField presence check:");
            println!("  'c' (current): {}", json.get("c").is_some());
            println!("  'h' (high): {}", json.get("h").is_some());
            println!("  'l' (low): {}", json.get("l").is_some());
            println!("  'o' (open): {}", json.get("o").is_some());
            println!("  'pc' (prev close): {}", json.get("pc").is_some());
            println!("  't' (timestamp): {}", json.get("t").is_some());

            // Small delay to be respectful of API
            tokio::time::sleep(tokio::time::Duration::from_millis(500)).await;
        }
    }

    #[tokio::test]
    #[ignore]
    async fn test_alpha_vantage_quote_raw() {
        let config = get_test_api_keys();
        let client = reqwest::Client::new();

        let symbol = "AAPL";

        let url = format!(
            "https://www.alphavantage.co/query?function=GLOBAL_QUOTE&symbol={}&entitlement=delayed&apikey={}",
            symbol, config.alpha_vantage
        );

        println!("\nFetching Alpha Vantage GLOBAL_QUOTE for {}...", symbol);

        let response = client.get(&url).send().await.expect("Request failed");
        let status = response.status();
        println!("Status: {}", status);

        let json: Value = response.json().await.expect("Failed to parse as JSON");
        pretty_print_json(&json, &format!("Alpha Vantage GLOBAL_QUOTE for {}", symbol));
        inspect_json_structure(&json, String::new(), 0);

        // Check for rate limit
        if json.get("Note").is_some() {
            println!("\nRATE LIMITED - Alpha Vantage has a 75 calls/min limit");
            return;
        }

        // Check for expected fields
        println!("\nField presence check:");
        println!("  'Global Quote': {}", json.get("Global Quote").is_some());
        if let Some(quote) = json.get("Global Quote") {
            println!(
                "  'Global Quote'.'05. price': {}",
                quote.get("05. price").is_some()
            );
            println!(
                "  'Global Quote'.'01. symbol': {}",
                quote.get("01. symbol").is_some()
            );
        }
    }

    #[tokio::test]
    #[ignore]
    async fn test_alpha_vantage_historical_raw() {
        let config = get_test_api_keys();
        let client = reqwest::Client::new();

        let symbol = "AAPL";

        let url = format!(
            "https://www.alphavantage.co/query?function=TIME_SERIES_DAILY_ADJUSTED&symbol={}&outputsize=compact&apikey={}",
            symbol, config.alpha_vantage
        );

        println!(
            "\nFetching Alpha Vantage TIME_SERIES_DAILY_ADJUSTED for {}...",
            symbol
        );
        println!("URL: {}", url.replace(&config.alpha_vantage, "***"));
        println!("Note: Using outputsize=compact (last 100 days) to reduce response size");

        let response = client.get(&url).send().await.expect("Request failed");
        let status = response.status();
        println!("Status: {}", status);

        let json: Value = response.json().await.expect("Failed to parse as JSON");

        // Check for rate limit or errors first
        if let Some(note) = json.get("Note") {
            println!("\nRATE LIMITED: {}", note);
            return;
        }
        if let Some(error) = json.get("Error Message") {
            println!("\nAPI ERROR: {}", error);
            return;
        }
        if let Some(info) = json.get("Information") {
            println!("\nAPI INFO: {}", info);
            return;
        }

        // Print structure overview
        println!("\nTop-level keys:");
        if let Some(obj) = json.as_object() {
            for key in obj.keys() {
                println!("  - {}", key);
            }
        }

        // Print detailed structure (will be large)
        inspect_json_structure(&json, String::new(), 0);

        // Check for expected fields
        println!("\nField presence check:");
        println!(
            "  'Time Series (Daily)': {}",
            json.get("Time Series (Daily)").is_some()
        );

        if let Some(time_series) = json.get("Time Series (Daily)") {
            if let Some(map) = time_series.as_object() {
                println!("  Number of dates: {}", map.len());

                // Show first date's structure
                if let Some((date, data)) = map.iter().next() {
                    println!("\n  Sample date: {}", date);
                    println!("  Fields in date entry:");
                    if let Some(data_obj) = data.as_object() {
                        for key in data_obj.keys() {
                            println!("    - {}", key);
                        }
                    }
                    println!(
                        "  Has '5. adjusted close': {}",
                        data.get("5. adjusted close").is_some()
                    );
                }
            }
        }

        // Print first 3 dates in full for inspection
        if let Some(time_series) = json.get("Time Series (Daily)") {
            if let Some(map) = time_series.as_object() {
                let mut dates: Vec<_> = map.keys().collect();
                dates.sort();
                dates.reverse();

                for (i, date) in dates.iter().take(3).enumerate() {
                    println!("\n--- Date {} Details: {} ---", i + 1, date);
                    pretty_print_json(&map[*date], &format!("Data for {}", date));
                }
            }
        }
    }

    // Deserialization Tests 

    #[tokio::test]
    #[ignore]
    async fn test_finnhub_deserialize() {
        let fetcher = DataFetcher::new();

        let symbols = vec!["AAPL", "MSFT"];

        for symbol in symbols {
            println!("\nTesting Finnhub deserialization for {}...", symbol);

            match fetcher.get_current_price_finnhub(symbol).await {
                Ok(quote) => {
                    println!("Successfully deserialized FinnhubQuote:");
                    // Sanity checks
                    assert!(quote.c > 0.0, "Current price should be positive");
                }
                Err(e) => {
                    println!("Deserialization failed: {}", e);
                    panic!("Finnhub deserialization failed for {}: {}", symbol, e);
                }
            }

            tokio::time::sleep(tokio::time::Duration::from_millis(500)).await;
        }
    }

    #[tokio::test]
    #[ignore]
    async fn test_alpha_vantage_quote_parse() {
        let fetcher = DataFetcher::new();


        let symbol = "AAPL";

        println!(
            "\nTesting Alpha Vantage GLOBAL_QUOTE parsing for {}...",
            symbol
        );

        match fetcher.get_current_price_alpha_vantage(symbol).await {
            Ok(price) => {
                println!("Successfully parsed price: ${:.2}", price);
                assert!(price > 0.0, "Price should be positive");
            }
            Err(e) => {
                println!("Parsing failed: {}", e);
                panic!("Alpha Vantage quote parsing failed: {}", e);
            }
        }
    }

    #[tokio::test]
    #[ignore]
    async fn test_alpha_vantage_historical_parse() {
        let fetcher = DataFetcher::new();


        let symbol = "AAPL";
        let lookback = 30;

        println!(
            "\nTesting Alpha Vantage TIME_SERIES_DAILY_ADJUSTED parsing for {}...",
            symbol
        );
        println!("Requesting {} days of history", lookback);

        match fetcher.get_historical_data(symbol, lookback).await {
            Ok((prices, yesterday_close)) => {
                println!("Successfully parsed historical data:");
                println!("  Number of prices: {}", prices.len());
                println!("  Yesterday's adjusted close: ${:.2}", yesterday_close);

                // Show first 5 prices
                println!("\n  First 5 prices (oldest to newest):");
                for (i, (date, price)) in prices.iter().take(5).enumerate() {
                    println!("    {}: {} -> ${:.2}", i + 1, date, price);
                }

                // Sanity checks
                assert!(!prices.is_empty(), "Should have at least some prices");
                assert!(yesterday_close > 0.0, "Yesterday close should be positive");
            }
            Err(e) => {
                println!("Parsing failed: {}", e);
                panic!("Alpha Vantage historical parsing failed: {}", e);
            }
        }
    }

    // Integration and Smoke Tests

    #[tokio::test]
    #[ignore]
    async fn test_update_symbol_integration() {
        let fetcher = DataFetcher::new();


        let symbol = "AAPL";
        let lookback_days = 30;

        println!("\nTesting full update_symbol integration for {}...", symbol);

        match fetcher.update_symbol(symbol, lookback_days).await {
            Ok(market_data) => {
                println!("Successfully fetched market data:");
                // Sanity checks
                assert!(
                    market_data.spot_price > 0.0,
                    "Spot price should be positive"
                );
                assert!(
                    market_data.volatility > 0.0,
                    "Volatility should be positive"
                );
                assert!(
                    market_data.yesterday_adjusted_close > 0.0,
                    "Yesterday close should be positive"
                );
            }
            Err(e) => {
                println!("Integration test failed: {}", e);
                panic!("update_symbol failed: {}", e);
            }
        }
    }

    #[tokio::test]
    async fn test_api_connectivity_smoke() {
        println!("\nRunning smoke test (basic API connectivity check)");
        println!("   Run full tests with: cargo test --package cli -- --ignored --nocapture\n");

        let config = match load_api_keys() {
            Ok(c) => c,
            Err(e) => {
                println!("Skipping smoke test: {}", e);
                println!("   Make sure api_keys.json exists in project root");
                return;
            }
        };

        let client = reqwest::Client::new();

        // Quick Finnhub check
        let url = format!(
            "https://finnhub.io/api/v1/quote?symbol=AAPL&token={}",
            config.finnhub
        );

        match client.get(&url).send().await {
            Ok(response) => {
                let status = response.status();
                println!("Finnhub API reachable (status: {})", status);
                if !status.is_success() {
                    println!("Non-success status - check API key");
                }
            }
            Err(e) => {
                println!("Finnhub API unreachable: {}", e);
            }
        }
    }
    #[tokio::test]
    async fn test_fetch_fred_vix() {
        let fetcher = DataFetcher::new();


        match fetcher.fetch_fred_vix(1).await {
            Ok(vix_value) => {
                println!("Successfully fetched VIX value from FRED: {:.2}", vix_value[0]);
            }
            Err(e) => {
                println!("Failed to fetch VIX from FRED: {}", e);
                panic!("fetch_fred_vix failed: {}", e);
            }
        }
    }

    // Volatility Tests
    #[test]
    fn test_calculate_volatility() {
        assert_eq!(DataFetcher::calculate_volatility(&vec![("2024-01-01".to_string(), 100.0); 4]), 0.0);
        assert_eq!(DataFetcher::calculate_volatility(&vec![]), 0.0);
        assert_eq!(DataFetcher::calculate_volatility(&vec![("1".to_string(), 100.0), ("2".to_string(), 102.0)]), 0.0);

        let mut trend = vec![("2024-01-01".to_string(), 100.0)];
        for i in 1..31 { trend.push((format!("{}", i), trend.last().unwrap().1 * 1.01)); }
        assert!(DataFetcher::calculate_volatility(&trend) < 0.05);
    }

    // Options Chain Tests

    #[tokio::test]
    #[ignore]
    async fn test_fetch_historical_options_raw() {
        let config = get_test_api_keys();
        let client = reqwest::Client::new();

        let url = "https://www.alphavantage.co/query?function=HISTORICAL_OPTIONS&symbol=IBM";

        println!("\nFetching Alpha Vantage HISTORICAL_OPTIONS for IBM...");

        let response = client
            .get(url)
            .header("X-API-Key", &config.alpha_vantage)
            .send()
            .await
            .expect("Request failed");
        let status = response.status();
        println!("Status: {}", status);

        let json: Value = response.json().await.expect("Failed to parse as JSON");

        // Print top-level structure
        println!("\nTop-level keys:");
        if let Some(obj) = json.as_object() {
            for key in obj.keys() {
                println!("  - {}", key);
            }
        }

        // Inspect first few contracts
        if let Some(data) = json.get("data").and_then(|d| d.as_array()) {
            println!("\nNumber of contracts: {}", data.len());
            for (i, contract) in data.iter().take(3).enumerate() {
                println!("\n--- Contract {} ---", i + 1);
                println!("{}", serde_json::to_string_pretty(contract).unwrap());
            }
        } else {
            println!("\nFull response:");
            println!("{}", serde_json::to_string_pretty(&json).unwrap());
        }
    }

    #[tokio::test]
    #[ignore]
    async fn test_fetch_option_chain_historical() {
        use crate::options_chain::OptionsEndpoint;

        let fetcher = DataFetcher::new();
        let chain = fetcher
            .fetch_option_chain("IBM", OptionsEndpoint::Historical { date: None })
            .await
            .expect("Failed to fetch option chain");

        assert!(!chain.contracts.is_empty(), "Chain should have contracts");
        println!("Fetched {} contracts for {}", chain.contracts.len(), chain.symbol);
        println!("Underlying price estimate: ${:.2}", chain.underlying_price);
        println!("Data date: {}", chain.data_date);

        if let Some(iv) = chain.summary_implied_volatility {
            println!("Summary IV: {:.2}%", iv * 100.0);
        } else {
            println!("Could not compute summary IV");
        }

        // Print first few contracts
        for contract in chain.contracts.iter().take(5) {
            println!(
                "  {:?} {} strike=${:.2} exp={} IV={:.4} OI={}",
                contract.option_type, contract.symbol, contract.strike,
                contract.expiration, contract.implied_volatility, contract.open_interest
            );
        }
    }

    #[tokio::test]
    #[ignore]
    async fn test_fetch_implied_volatility() {
        use crate::options_chain::OptionsEndpoint;

        let fetcher = DataFetcher::new();
        let iv = fetcher
            .fetch_implied_volatility("IBM", OptionsEndpoint::Historical { date: None })
            .await
            .expect("Failed to fetch IV");

        assert!(iv > 0.0 && iv < 2.0, "IV {:.4} outside reasonable range", iv);
        println!("IBM implied volatility: {:.2}%", iv * 100.0);
    }

    #[tokio::test]
    #[ignore]
    async fn test_realtime_options_premium_error() {
        use crate::options_chain::OptionsEndpoint;

        let fetcher = DataFetcher::new();
        let result = fetcher
            .fetch_option_chain("IBM", OptionsEndpoint::Realtime)
            .await;

        match result {
            Ok(chain) => println!("Premium user: got {} contracts", chain.contracts.len()),
            Err(e) => {
                let err = e.to_string();
                assert!(
                    err.contains("premium") || err.contains("Information"),
                    "Error should indicate premium requirement: {}",
                    err
                );
                println!("Expected premium error: {}", err);
            }
        }
    }

    #[tokio::test]
    #[ignore]
    async fn test_volatility_smile_extraction_live() {
        use crate::options_chain::OptionsEndpoint;

        let fetcher = DataFetcher::new();
        let chain = fetcher
            .fetch_option_chain("SPY", OptionsEndpoint::Historical { date: None })
            .await
            .expect("Failed to fetch SPY chain");

        let smile = chain.extract_smile(None, 100);

        println!("SPY Volatility Smile:");
        println!("Underlying: ${:.2}", smile.underlying_price);
        println!("Expirations: {}", smile.smiles_by_expiry.len());

        for (expiry, points) in &smile.smiles_by_expiry {
            println!("\n  Expiration: {}", expiry);
            for p in points.iter().take(5) {
                println!(
                    "    Strike: ${:.2}, IV: {:.2}%, Type: {:?}, Moneyness: {:.3}",
                    p.strike,
                    p.implied_volatility * 100.0,
                    p.option_type,
                    p.moneyness
                );
            }
            if points.len() > 5 {
                println!("    ... ({} more points)", points.len() - 5);
            }
        }
    }

    #[tokio::test]
    #[ignore]
    async fn test_fetch_fred_sofr() {
        let fetcher = DataFetcher::new();
        let rate = fetcher
            .fetch_fred_sofr()
            .await
            .expect("Failed to fetch SOFR");
        assert!(
            rate > 0.0 && rate < 0.15,
            "SOFR rate {} outside reasonable range (0-0.15)",
            rate
        );
        println!("Current SOFR: {:.4}% ({:.6} decimal)", rate * 100.0, rate);
    }

    #[tokio::test]
    #[ignore]
    async fn test_fetch_dividend_yield_aapl() {
        let fetcher = DataFetcher::new();
        let fundamentals = fetcher
            .fetch_company_fundamentals("AAPL")
            .await
            .expect("Failed to fetch AAPL fundamentals");
        assert!(
            fundamentals.dividend_yield.is_some(),
            "AAPL should have a dividend yield"
        );
        let div_str = fundamentals.dividend_yield.as_ref().unwrap();
        let div = div_str.parse::<f64>().expect("Failed to parse dividend yield");
        assert!(div > 0.0, "AAPL dividend yield should be positive");
        println!("AAPL dividend yield: {}% ({} decimal)", div * 100.0, div);
    }

    #[tokio::test]
    #[ignore]
    async fn test_fetch_dividend_yield_no_dividend() {
        let fetcher = DataFetcher::new();
        let fundamentals = fetcher
            .fetch_company_fundamentals("TSLA")
            .await
            .expect("Failed to fetch TSLA fundamentals");
        if let Some(div_str) = fundamentals.dividend_yield.clone() {
            let div = div_str.parse::<f64>().unwrap_or(0.0);
            assert_eq!(
                div, 0.0,
                "TSLA should have no dividend (or zero dividend)"
            );
        }
        println!("TSLA dividend yield: {:?}", fundamentals.dividend_yield);
    }

    #[tokio::test]
    #[ignore]
    async fn test_alpha_vantage_primary() {
        let fetcher = DataFetcher::new();
        let price = fetcher
            .get_current_price_with_fallback("AAPL")
            .await
            .expect("Failed to fetch AAPL price");
        assert!(price > 0.0, "AAPL price should be positive");
        println!("AAPL price (via Alpha Vantage primary): ${:.2}", price);
    }

    #[test]
    fn test_calculate_time_to_dividend_valid() {
        let future_date = (Utc::now().date_naive() + Duration::days(90))
            .format("%Y-%m-%d")
            .to_string();

        let result = calculate_time_to_dividend(&future_date).unwrap();
        assert!((result - 0.246).abs() < 0.01);
    }

    #[test]
    fn test_calculate_time_to_dividend_past_date() {
        let result = calculate_time_to_dividend("2020-01-01");
        assert!(result.is_err());
    }

    #[test]
    fn test_calculate_time_to_dividend_invalid_format() {
        let result = calculate_time_to_dividend("invalid-date");
        assert!(result.is_err());
    }

    #[tokio::test]
    #[ignore]
    async fn test_get_dividend_info_aapl() {
        let fetcher = DataFetcher::new();
        let result = fetcher.get_dividend_info("AAPL", 150.0).await;
        assert!(result.is_ok());
        if let Ok(Some(div_info)) = result {
            assert!(div_info.time_to_dividend >= 0.0);
            assert!(div_info.dividend_amount > 0.0);
            assert!(div_info.dividend_yield > 0.0);
            println!("AAPL Dividend Info: time_to_div={:.3} years, amount=${:.2}, yield={:.2}%",
                div_info.time_to_dividend, div_info.dividend_amount, div_info.dividend_yield * 100.0);
        }
    }
}
