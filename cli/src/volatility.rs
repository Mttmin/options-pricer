use crate::fetcher::DataFetcher;
use ggca::correlation::spearman_correlation;

/// Tries and create an ema estimator for volatility, 
/// in case of error falls back to a no estimator and pure historical volatility
fn ema_volatility(data: &[(String, f64)],n_days_per: u16, lambda: f64) -> Result<f64, f64> {
    let n_working = data.len() as u16 - n_days_per;
    if n_working > 1 || n_days_per < 2{
        let mut ema = 0.0;
        // loop over the days where the volatility can be calculated and calculate it
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
    let n_working = data.len() as u16 - n_days_per;
    if n_working > 1 || n_days_per < 2{
        let mut ema = 0.0;
        let mut ema_vector: Vec<f64> = Vec::new();
        // loop over the days where the volatility can be calculated and calculate it
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
        Err(format!("Wrong input parameters, length of data vs n_days_per={} or n_days_per = {} has to be over 2",n_working, n_days_per).into())
    }
}

async fn vix_volatility(symbol: &str, data_fetcher:  DataFetcher, num_days: u16) -> Result<f64,f64>{
    // get historical VIX data and close price
    let mut vix_data: Vec<f64> = data_fetcher.fetch_fred_vix(num_days).await.expect("Vix not possible");
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
    return Ok(last_vix*correlation + (1.0-correlation)*latest_ema_vol)
}
mod tests {
    use super::*;

    #[test]
    fn compare_vols(){

    }
}
