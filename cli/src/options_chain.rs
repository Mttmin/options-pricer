use chrono::{NaiveDate, DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Represents the type of an option contract
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum OptionType {
    Call,
    Put,
}

/// A single option contract as returned by Alpha Vantage
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptionContract {
    pub contract_id: String,
    pub symbol: String,
    pub expiration: NaiveDate,
    pub strike: f64,
    pub option_type: OptionType,
    pub last: f64,
    pub mark: f64,
    pub bid: f64,
    pub ask: f64,
    pub volume: u64,
    pub open_interest: u64,
    pub implied_volatility: f64,
    pub delta: Option<f64>,
    pub gamma: Option<f64>,
    pub theta: Option<f64>,
    pub vega: Option<f64>,
    pub rho: Option<f64>,
    pub in_the_money: bool,
}

/// The full option chain for a symbol
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptionChain {
    pub symbol: String,
    pub underlying_price: f64,
    pub fetch_timestamp: DateTime<Utc>,
    pub data_date: NaiveDate,
    pub contracts: Vec<OptionContract>,
    /// Summary IV computed from ATM contracts
    pub summary_implied_volatility: Option<f64>,
}

/// Distinguishes which Alpha Vantage endpoint to use
#[derive(Debug, Clone, Copy)]
pub enum OptionsEndpoint {
    /// HISTORICAL_OPTIONS - works with all API tiers, previous session data
    Historical { date: Option<NaiveDate> },
    /// REALTIME_OPTIONS - premium tier only
    Realtime,
}

/// A single point on the volatility smile curve
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SmilePoint {
    pub strike: f64,
    pub implied_volatility: f64,
    pub option_type: OptionType,
    pub volume: u64,
    pub open_interest: u64,
    /// Moneyness = strike / underlying_price
    pub moneyness: f64,
}

/// Volatility smile data grouped by expiration, suitable for plotting
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VolatilitySmile {
    pub symbol: String,
    pub underlying_price: f64,
    pub data_date: NaiveDate,
    /// Each key is an expiration date, value is sorted vec of SmilePoints
    pub smiles_by_expiry: HashMap<NaiveDate, Vec<SmilePoint>>,
}

impl OptionChain {
    /// Computes a summary implied volatility from the chain.
    ///
    /// Strategy: Take the nearest expiration, find contracts within ±5%
    /// of ATM (strike near underlying_price), weight by open_interest,
    /// and average. This approximates ATM implied volatility.
    pub fn compute_summary_iv(&mut self) -> Option<f64> {
        if self.contracts.is_empty() || self.underlying_price <= 0.0 {
            self.summary_implied_volatility = None;
            return None;
        }

        // Find nearest expiration
        let min_expiry = self.contracts.iter().map(|c| c.expiration).min()?;

        // Filter to nearest expiry, near-ATM (within 5%), with valid IV
        let atm_tolerance = 0.05;
        let near_atm: Vec<&OptionContract> = self
            .contracts
            .iter()
            .filter(|c| {
                c.expiration == min_expiry
                    && c.implied_volatility > 0.0
                    && ((c.strike / self.underlying_price) - 1.0).abs() <= atm_tolerance
                    && c.open_interest > 0
            })
            .collect();

        if near_atm.is_empty() {
            self.summary_implied_volatility = None;
            return None;
        }

        // Open-interest-weighted average
        let total_oi: f64 = near_atm.iter().map(|c| c.open_interest as f64).sum();

        let weighted_iv: f64 = near_atm
            .iter()
            .map(|c| c.implied_volatility * c.open_interest as f64)
            .sum::<f64>()
            / total_oi;

        self.summary_implied_volatility = Some(weighted_iv);
        Some(weighted_iv)
    }

    /// Extracts volatility smile data suitable for plotting.
    ///
    /// Groups contracts by expiration and sorts by strike.
    /// Filters out contracts with zero IV or insufficient open interest.
    pub fn extract_smile(
        &self,
        expirations: Option<&[NaiveDate]>,
        min_open_interest: u64,
    ) -> VolatilitySmile {
        let mut smiles_by_expiry: HashMap<NaiveDate, Vec<SmilePoint>> = HashMap::new();

        for contract in &self.contracts {
            if let Some(expiries) = expirations
                && !expiries.contains(&contract.expiration) {
                    continue;
            }

            if contract.implied_volatility <= 0.0 || contract.open_interest < min_open_interest {
                continue;
            }

            let point = SmilePoint {
                strike: contract.strike,
                implied_volatility: contract.implied_volatility,
                option_type: contract.option_type,
                volume: contract.volume,
                open_interest: contract.open_interest,
                moneyness: contract.strike / self.underlying_price,
            };

            smiles_by_expiry
                .entry(contract.expiration)
                .or_default()
                .push(point);
        }

        // Sort each expiry's points by strike
        for points in smiles_by_expiry.values_mut() {
            points.sort_by(|a, b| a.strike.partial_cmp(&b.strike).unwrap());
        }

        VolatilitySmile {
            symbol: self.symbol.clone(),
            underlying_price: self.underlying_price,
            data_date: self.data_date,
            smiles_by_expiry,
        }
    }

    /// Returns a list of all available expiration dates, sorted
    pub fn available_expirations(&self) -> Vec<NaiveDate> {
        let mut expiries: Vec<NaiveDate> = self
            .contracts
            .iter()
            .map(|c| c.expiration)
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        expiries.sort();
        expiries
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_contracts() -> Vec<OptionContract> {
        vec![
            OptionContract {
                contract_id: "TEST1".into(),
                symbol: "TEST".into(),
                expiration: NaiveDate::from_ymd_opt(2026, 3, 20).unwrap(),
                strike: 95.0,
                option_type: OptionType::Put,
                last: 1.5,
                mark: 1.5,
                bid: 1.4,
                ask: 1.6,
                volume: 100,
                open_interest: 500,
                implied_volatility: 0.28,
                delta: Some(-0.3),
                gamma: Some(0.02),
                theta: Some(-0.05),
                vega: Some(15.0),
                rho: Some(-2.0),
                in_the_money: false,
            },
            OptionContract {
                contract_id: "TEST2".into(),
                symbol: "TEST".into(),
                expiration: NaiveDate::from_ymd_opt(2026, 3, 20).unwrap(),
                strike: 100.0,
                option_type: OptionType::Call,
                last: 3.0,
                mark: 3.0,
                bid: 2.9,
                ask: 3.1,
                volume: 200,
                open_interest: 1000,
                implied_volatility: 0.22,
                delta: Some(0.5),
                gamma: Some(0.03),
                theta: Some(-0.08),
                vega: Some(20.0),
                rho: Some(5.0),
                in_the_money: false,
            },
            OptionContract {
                contract_id: "TEST3".into(),
                symbol: "TEST".into(),
                expiration: NaiveDate::from_ymd_opt(2026, 3, 20).unwrap(),
                strike: 105.0,
                option_type: OptionType::Call,
                last: 1.2,
                mark: 1.2,
                bid: 1.1,
                ask: 1.3,
                volume: 150,
                open_interest: 800,
                implied_volatility: 0.25,
                delta: Some(0.3),
                gamma: Some(0.02),
                theta: Some(-0.04),
                vega: Some(12.0),
                rho: Some(3.0),
                in_the_money: false,
            },
            OptionContract {
                contract_id: "TEST4".into(),
                symbol: "TEST".into(),
                expiration: NaiveDate::from_ymd_opt(2026, 6, 19).unwrap(),
                strike: 100.0,
                option_type: OptionType::Call,
                last: 8.0,
                mark: 8.0,
                bid: 7.8,
                ask: 8.2,
                volume: 50,
                open_interest: 300,
                implied_volatility: 0.24,
                delta: Some(0.55),
                gamma: Some(0.015),
                theta: Some(-0.03),
                vega: Some(30.0),
                rho: Some(10.0),
                in_the_money: false,
            },
        ]
    }

    fn make_test_chain() -> OptionChain {
        OptionChain {
            symbol: "TEST".into(),
            underlying_price: 100.0,
            fetch_timestamp: Utc::now(),
            data_date: NaiveDate::from_ymd_opt(2026, 2, 3).unwrap(),
            contracts: make_test_contracts(),
            summary_implied_volatility: None,
        }
    }

    #[test]
    fn test_compute_summary_iv() {
        let mut chain = make_test_chain();
        let iv = chain.compute_summary_iv();
        assert!(iv.is_some());
        let iv_val = iv.unwrap();
        // Should be between 0.20 and 0.30 given our test data
        assert!(
            iv_val > 0.20 && iv_val < 0.30,
            "Summary IV {:.4} out of expected range",
            iv_val
        );
    }

    #[test]
    fn test_compute_summary_iv_empty_chain() {
        let mut chain = OptionChain {
            symbol: "TEST".into(),
            underlying_price: 100.0,
            fetch_timestamp: Utc::now(),
            data_date: NaiveDate::from_ymd_opt(2026, 2, 3).unwrap(),
            contracts: vec![],
            summary_implied_volatility: None,
        };
        assert!(chain.compute_summary_iv().is_none());
    }

    #[test]
    fn test_extract_smile() {
        let chain = make_test_chain();
        let smile = chain.extract_smile(None, 0);

        // Should have 2 expirations
        assert_eq!(smile.smiles_by_expiry.len(), 2);

        // Near expiration should have 3 contracts
        let near_expiry = NaiveDate::from_ymd_opt(2026, 3, 20).unwrap();
        assert_eq!(smile.smiles_by_expiry[&near_expiry].len(), 3);

        // Points should be sorted by strike
        let near_points = &smile.smiles_by_expiry[&near_expiry];
        assert!(near_points[0].strike < near_points[1].strike);
        assert!(near_points[1].strike < near_points[2].strike);
    }

    #[test]
    fn test_extract_smile_filtered_expiration() {
        let chain = make_test_chain();
        let target = NaiveDate::from_ymd_opt(2026, 3, 20).unwrap();
        let smile = chain.extract_smile(Some(&[target]), 0);

        assert_eq!(smile.smiles_by_expiry.len(), 1);
        assert!(smile.smiles_by_expiry.contains_key(&target));
    }

    #[test]
    fn test_extract_smile_min_open_interest_filter() {
        let chain = make_test_chain();
        // Filter to OI >= 600 (should exclude the 500-OI put and 300-OI far contract)
        let smile = chain.extract_smile(None, 600);
        let near_expiry = NaiveDate::from_ymd_opt(2026, 3, 20).unwrap();
        assert_eq!(smile.smiles_by_expiry[&near_expiry].len(), 2);
    }

    #[test]
    fn test_available_expirations() {
        let chain = make_test_chain();
        let expiries = chain.available_expirations();
        assert_eq!(expiries.len(), 2);
        assert_eq!(expiries[0], NaiveDate::from_ymd_opt(2026, 3, 20).unwrap());
        assert_eq!(expiries[1], NaiveDate::from_ymd_opt(2026, 6, 19).unwrap());
    }

    #[test]
    fn test_moneyness_calculation() {
        let chain = make_test_chain();
        let smile = chain.extract_smile(None, 0);
        let near_expiry = NaiveDate::from_ymd_opt(2026, 3, 20).unwrap();
        let points = &smile.smiles_by_expiry[&near_expiry];

        assert!((points[0].moneyness - 0.95).abs() < 0.001);
        assert!((points[1].moneyness - 1.00).abs() < 0.001);
        assert!((points[2].moneyness - 1.05).abs() < 0.001);
    }
}
