use crate::fetcher::DataFetcher;
use ggca::correlation::spearman_correlation;

/// Tries and create an ema estimator for volatility,
/// in case of error falls back to a no estimator and pure historical volatility
pub fn ema_volatility(data: &[(String, f64)],n_days_per: u16, lambda: f64) -> Result<f64, f64> {
    // Check for valid parameters before subtraction to avoid overflow
    if data.len() <= n_days_per as usize || n_days_per < 2 {
        return Err(DataFetcher::calculate_volatility(data));
    }

    let n_working = data.len() as u16 - n_days_per;
    if n_working > 1 {
        // Initialize EMA with the first window's volatility
        let first_vol = DataFetcher::calculate_volatility(&data[0..n_days_per as usize]);
        let mut ema = first_vol;

        // loop over the remaining days where the volatility can be calculated
        for i in 1..n_working{
            let i_usize = i as usize;
            let n_days_usize = n_days_per as usize;
            let vol_daily = DataFetcher::calculate_volatility(
                &data[i_usize..i_usize + n_days_usize]
            );
            ema = lambda*ema + (1.0-lambda)*vol_daily;
        }
        Ok(ema)
    }
    else {
        Err(DataFetcher::calculate_volatility(data))
    }
}

fn ema_volatility_all(data: &[(String, f64)],n_days_per: u16, lambda: f64) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    // Check for valid parameters before subtraction to avoid overflow
    if data.len() <= n_days_per as usize || n_days_per < 2 {
        return Err(format!("Wrong input parameters, data length={} vs n_days_per={}, need data.len() > n_days_per and n_days_per >= 2",
            data.len(), n_days_per).into());
    }

    let n_working = data.len() as u16 - n_days_per;
    if n_working > 1 {
        // Initialize EMA with the first window's volatility
        let first_vol = DataFetcher::calculate_volatility(&data[0..n_days_per as usize]);
        let mut ema = first_vol;
        let mut ema_vector: Vec<f64> = Vec::new();
        ema_vector.push(ema);

        // loop over the remaining days where the volatility can be calculated
        for i in 1..n_working{
            let i_usize = i as usize;
            let n_days_usize = n_days_per as usize;
            let vol_daily = DataFetcher::calculate_volatility(
                &data[i_usize..i_usize + n_days_usize]
            );
            ema = lambda*ema + (1.0-lambda)*vol_daily;
            ema_vector.push(ema);
        }
        Ok(ema_vector)
    }
    else {
        Err(format!("Wrong input parameters, n_working={} must be > 1", n_working).into())
    }
}

pub async fn vix_volatility(symbol: &str, data_fetcher:  DataFetcher, num_days: u16) -> Result<f64,f64>{
    // get historical VIX data and close price
    let mut vix_data: Vec<f64> = data_fetcher.fetch_fred_vix(num_days).await.expect("Vix not possible");
    vix_data.iter_mut().for_each(|x| *x /= 100.0);
    let stock_data: Vec<(String, f64)> = data_fetcher.get_historical_data(symbol, num_days.into()).await.expect("Stock data not available").0;
    let mut ema_vol: Vec<f64> = ema_volatility_all(&stock_data, num_days, 0.9).unwrap();
    // extract latest vix and volatility
    let last_vix = vix_data.pop().unwrap_or(-1.0);
    let latest_ema_vol = ema_vol[0];
    ema_vol.remove(0);
    // calculate the correlation between previous vix and price vol for next day
    let correlation = spearman_correlation(&ema_vol, &vix_data);
    // use the ema and correlation to construct and estimate of today's volatility
    if correlation < 0.2{
        eprintln!("Not using vix correlation as it is too low");
        return Err(latest_ema_vol)
    }
    Ok(last_vix*correlation + (1.0-correlation)*latest_ema_vol)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fetcher::DataFetcher;

    fn gen_prices(n: usize, vol: f64) -> Vec<(String, f64)> {
        let mut p = vec![("2024-01-01".to_string(), 100.0)];
        for i in 1..n {
            p.push((format!("2024-{:02}", i), p.last().unwrap().1 * (1.0 + vol * (i as f64 * 0.5).sin())));
        }
        p
    }

    #[test]
    fn test_ema_volatility() {
        let prices = gen_prices(31, 0.02);
        assert!(ema_volatility(&prices, 10, 0.9).is_ok());
        assert!(ema_volatility(&prices, 10, 0.9).unwrap() > 0.0);
        assert!(ema_volatility(&prices, 50, 0.9).is_err()); // insufficient data
    }

    #[test]
    fn test_ema_volatility_all() {
        let prices = gen_prices(31, 0.02);
        let result = ema_volatility_all(&prices, 10, 0.9).unwrap();
        assert_eq!(result.len(), 21);
    }

    #[test]
    fn test_compare_methods() {
        let prices = gen_prices(91, 0.02);
        let std_vol = DataFetcher::calculate_volatility(&prices);
        let ema_vol = ema_volatility(&prices, 30, 0.9).unwrap();

        println!("Standard: {:.2}%, EMA: {:.2}%", std_vol * 100.0, ema_vol * 100.0);
        assert!(std_vol > 0.0 && ema_vol > 0.0);
        assert!(ema_vol / std_vol > 0.3 && ema_vol / std_vol < 3.0);
    }

    #[test]
    fn test_regime_detection() {
        let stable = gen_prices(51, 0.001);
        let volatile = gen_prices(51, 0.03);

        let stable_vol = DataFetcher::calculate_volatility(&stable);
        let volatile_vol = DataFetcher::calculate_volatility(&volatile);

        println!("Stable: {:.2}%, Volatile: {:.2}%", stable_vol * 100.0, volatile_vol * 100.0);
        assert!(volatile_vol > stable_vol * 2.0);
    }

    #[tokio::test]
    #[ignore]
    async fn test_vix_integration() {
        let result = vix_volatility("SPY", DataFetcher::new(), 60).await;
        match result {
            Ok(vol) => println!("VIX: {:.2}%", vol * 100.0),
            Err(vol) => println!("Fallback: {:.2}%", vol * 100.0),
        }
    }

    #[test]
    #[ignore]
    fn test_performance() {
        use std::time::Instant;
        let prices = gen_prices(2521, 0.01);

        let t = Instant::now();
        DataFetcher::calculate_volatility(&prices);
        println!("Standard: {:?}", t.elapsed());

        let t = Instant::now();
        ema_volatility(&prices, 60, 0.9).unwrap();
        println!("EMA: {:?}", t.elapsed());

        let t = Instant::now();
        ema_volatility_all(&prices, 60, 0.9).unwrap();
        println!("EMA_all: {:?}", t.elapsed());
    }
}
