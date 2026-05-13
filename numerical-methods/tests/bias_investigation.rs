// Bias investigation: CTMC vs BS / binomial reference pricers.
//
// Strategy: degenerate Heston → GBM (sigma_v → 0, v0 = theta, kappa large).
// Under this limit the short variance stays pinned at v0, so spot follows
// GBM with vol = sqrt(v0). Any gap vs Black-Scholes / binomial is CTMC
// discretization or model bias — not legitimate SV effect.

use numerical_methods::binomial::{binomial_price, BinomialParameters, ExerciseStyle};
use numerical_methods::ctmc::{
    price_american_call_heston, price_american_option_heston, price_american_put_heston,
};
use options::black_scholes::black_scholes_price;
use options::Options;

const R: f64 = 0.05;
const T: f64 = 1.0;
const SPOT: f64 = 100.0;
const V0: f64 = 0.04; // vol = 20%

fn bs_call(strike: f64, q: f64) -> f64 {
    black_scholes_price(Options::new_call(
        strike, SPOT, V0.sqrt(), R, T, Some(q),
    ))
}
fn bs_put(strike: f64, q: f64) -> f64 {
    black_scholes_price(Options::new_put(
        strike, SPOT, V0.sqrt(), R, T, Some(q),
    ))
}
fn binomial_euro_put(strike: f64) -> f64 {
    let put = options::Put::new(strike, SPOT, V0.sqrt(), R, T, None);
    let params = BinomialParameters::new(SPOT, 1000, T, V0.sqrt(), R);
    binomial_price(&put, &params, ExerciseStyle::European)
}
fn binomial_am_put(strike: f64) -> f64 {
    let put = options::Put::new(strike, SPOT, V0.sqrt(), R, T, None);
    let params = BinomialParameters::new(SPOT, 1000, T, V0.sqrt(), R);
    binomial_price(&put, &params, ExerciseStyle::American)
}

// Degenerate-Heston CTMC: set sigma_v tiny, v0 = theta → short variance frozen at v0.
// Expectation: CTMC price ≈ BS/binomial (modulo grid error).
fn ctmc_degenerate_put(strike: f64, q: f64, n_x: usize, m_v: usize, n_time: usize) -> f64 {
    price_american_put_heston(
        strike, SPOT, V0, T, R, q,
        /*kappa=*/ 5.0, /*theta=*/ V0, /*sigma_v=*/ 1e-4, /*rho=*/ 0.0,
        n_x, m_v, n_time,
    ).price
}
fn ctmc_degenerate_call(strike: f64, q: f64, n_x: usize, m_v: usize, n_time: usize) -> f64 {
    price_american_call_heston(
        strike, SPOT, V0, T, R, q,
        5.0, V0, 1e-4, 0.0,
        n_x, m_v, n_time,
    ).price
}

// ---- 1. Degenerate CTMC vs BS: European-equivalent American call, q=0 ----
#[test]
#[ignore]
fn bias_1_degenerate_call_q_zero_vs_bs() {
    println!("\n=== TEST 1: Degenerate Heston CTMC call (q=0) vs BS ===");
    println!("American call with q=0 = European call (no EEP) → should match BS.");
    println!("{:<10} {:<12} {:<12} {:<10}", "strike", "ctmc", "bs", "bias");

    let strikes = [80.0, 90.0, 100.0, 110.0, 120.0];
    let mut max_abs_bias = 0.0_f64;
    for &k in &strikes {
        let ctmc = ctmc_degenerate_call(k, 0.0, 120, 15, 80);
        let bs = bs_call(k, 0.0);
        let bias = ctmc - bs;
        println!("{:<10.2} {:<12.4} {:<12.4} {:<+10.4}", k, ctmc, bs, bias);
        max_abs_bias = max_abs_bias.max(bias.abs());
    }
    println!("Max |bias|: {max_abs_bias:.4}");
}

// ---- 2. Degenerate CTMC put vs BS European put (q=0) ----
// Note: American put CAN early-exercise even at q=0 when r > 0 → real EEP.
// Compare to binomial American at matched vol to separate EEP from bias.
#[test]
#[ignore]
fn bias_2_degenerate_put_vs_binomial_american() {
    println!("\n=== TEST 2: Degenerate CTMC American put vs binomial American put ===");
    println!("Both American; matched vol = sqrt(v0) = {:.4}; q=0", V0.sqrt());
    println!("{:<10} {:<12} {:<12} {:<12} {:<10}",
        "strike", "ctmc_am", "binom_am", "bs_euro", "bias_am");

    let strikes = [80.0, 90.0, 100.0, 110.0, 120.0];
    let mut max_abs_bias = 0.0_f64;
    for &k in &strikes {
        let ctmc = ctmc_degenerate_put(k, 0.0, 120, 15, 80);
        let binom = binomial_am_put(k);
        let bs = bs_put(k, 0.0);
        let bias = ctmc - binom;
        println!("{:<10.2} {:<12.4} {:<12.4} {:<12.4} {:<+10.4}",
            k, ctmc, binom, bs, bias);
        max_abs_bias = max_abs_bias.max(bias.abs());
    }
    println!("Max |bias| vs binomial: {max_abs_bias:.4}");
}

// ---- 3. Grid resolution sweep: does bias shrink with finer grids? ----
#[test]
#[ignore]
fn bias_3_grid_resolution_sweep() {
    println!("\n=== TEST 3: Grid resolution — ATM put, degenerate Heston ===");
    let bs_ref = bs_put(100.0, 0.0);
    let binom_am = binomial_am_put(100.0);
    println!("Reference: BS_euro = {bs_ref:.4}, binom_am = {binom_am:.4}, EEP = {:.4}",
        binom_am - bs_ref);
    println!("{:<8} {:<6} {:<8} {:<10} {:<12}",
        "n_x", "m_v", "n_time", "ctmc", "bias_vs_binom");

    let configs = [
        (60, 10, 50),
        (80, 15, 75),
        (100, 20, 100),
        (150, 25, 150),
        (200, 30, 200),
    ];
    for &(nx, mv, nt) in &configs {
        let p = ctmc_degenerate_put(100.0, 0.0, nx, mv, nt);
        let bias = p - binom_am;
        println!("{:<8} {:<6} {:<8} {:<10.4} {:<+12.4}", nx, mv, nt, p, bias);
    }
}

// ---- 4. Vol-of-vol sweep: does bias exist only with sigma_v>0 or always? ----
#[test]
#[ignore]
fn bias_4_vol_of_vol_sweep() {
    println!("\n=== TEST 4: Vol-of-vol sweep, ATM put, q=0 ===");
    println!("kappa=2, theta=v0=0.04, rho=0. As sigma_v → 0, should → BS/binomial.");
    let bs_ref = bs_put(100.0, 0.0);
    let binom_am = binomial_am_put(100.0);
    println!("BS euro = {bs_ref:.4}, binom American = {binom_am:.4}");
    println!("{:<10} {:<10} {:<12}", "sigma_v", "ctmc", "diff_vs_binom");

    for &sv in &[1e-4_f64, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5] {
        let p = price_american_put_heston(
            100.0, SPOT, V0, T, R, 0.0,
            2.0, V0, sv, 0.0,
            120, 20, 80,
        ).price;
        println!("{:<10.4} {:<10.4} {:<+12.4}", sv, p, p - binom_am);
    }
}

// ---- 5. Heston Fourier vs CTMC: same Heston params, European proxy ----
// American call with q=0 equals European call. Compare CTMC(that) to COS pricer.
#[test]
#[ignore]
fn bias_5_ctmc_vs_heston_fourier_european() {
    use numerical_methods::slv::{
        CalibrationInstrument, CalibrationPricer, HestonEuropeanPricer,
        HestonParameters, OptionSide,
    };
    println!("\n=== TEST 5: CTMC (q=0 call = European) vs Heston COS Fourier ===");
    println!("Same Heston params, same strikes. Difference = CTMC grid bias.");
    let params = HestonParameters { v0: V0, v_bar: V0, kappa: 2.0, sigma: 0.3, rho: -0.7 };
    let pricer = HestonEuropeanPricer::default();

    println!("{:<10} {:<12} {:<12} {:<10}", "strike", "ctmc_call", "cos_call", "bias");
    let strikes = [80.0, 90.0, 100.0, 110.0, 120.0];
    for &k in &strikes {
        let ctmc = price_american_call_heston(
            k, SPOT, params.v0, T, R, 0.0,
            params.kappa, params.v_bar, params.sigma, params.rho,
            150, 25, 100,
        ).price;
        let inst = CalibrationInstrument {
            symbol: "T".into(), side: OptionSide::Call, strike: k, maturity_years: T,
            market_price: 0.0, spot: SPOT, risk_free_rate: R, dividend_yield: 0.0,
            implied_vol: None, bid: None, ask: None, open_interest: None, volume: None,
        };
        let cos = pricer.price(&inst, params).unwrap();
        println!("{:<10.2} {:<12.4} {:<12.4} {:<+10.4}", k, ctmc, cos, ctmc - cos);
    }
}

// ---- 6. Maturity sweep: does bias grow with T? ----
#[test]
#[ignore]
fn bias_6_maturity_sweep_put() {
    println!("\n=== TEST 6: Maturity sweep, ATM put, degenerate Heston vs binomial ===");
    println!("{:<8} {:<12} {:<12} {:<10}", "T(yrs)", "ctmc", "binom_am", "bias");
    for &t in &[0.25_f64, 0.5, 1.0, 2.0] {
        let put = options::Put::new(100.0, SPOT, V0.sqrt(), R, t, None);
        let bp = BinomialParameters::new(SPOT, 1000, t, V0.sqrt(), R);
        let binom = binomial_price(&put, &bp, ExerciseStyle::American);
        let ctmc = price_american_put_heston(
            100.0, SPOT, V0, t, R, 0.0,
            5.0, V0, 1e-4, 0.0,
            120, 20, 100,
        ).price;
        println!("{:<8.2} {:<12.4} {:<12.4} {:<+10.4}", t, ctmc, binom, ctmc - binom);
    }
}

// ---- 7. Sign/direction test: is bias always high or always low? ----
#[test]
#[ignore]
fn bias_7_sign_summary() {
    println!("\n=== TEST 7: Bias direction summary (degenerate Heston, grid 150/25/100) ===");
    let mut positive = 0;
    let mut negative = 0;
    let mut sum = 0.0;
    let mut abs_sum = 0.0;
    let mut n = 0;

    for &k in &[80.0_f64, 90.0, 100.0, 110.0, 120.0] {
        let ctmc = ctmc_degenerate_put(k, 0.0, 150, 25, 100);
        let binom = binomial_am_put(k);
        let bias = ctmc - binom;
        if bias > 0.0 { positive += 1 } else { negative += 1 }
        sum += bias;
        abs_sum += bias.abs();
        n += 1;
    }
    println!("Puts: {positive}/{n} biased HIGH, {negative}/{n} biased LOW");
    println!("Mean bias: {:.4}, Mean |bias|: {:.4}", sum / n as f64, abs_sum / n as f64);

    let mut csum = 0.0;
    let mut cabs = 0.0;
    let mut cn = 0;
    for &k in &[80.0_f64, 90.0, 100.0, 110.0, 120.0] {
        let ctmc = ctmc_degenerate_call(k, 0.0, 150, 25, 100);
        let bs = bs_call(k, 0.0);
        csum += ctmc - bs;
        cabs += (ctmc - bs).abs();
        cn += 1;
    }
    println!("Calls (q=0, American=European): mean bias vs BS: {:.4}, mean |bias|: {:.4}",
        csum / cn as f64, cabs / cn as f64);
}

// ---- 14. Boundary trajectory at n_time=30 (bug regime) — suspect collapse ----
#[test]
#[ignore]
fn bias_14_boundary_collapse_hypothesis() {
    println!("\n=== TEST 14: Boundary trajectory at n_time=30 (bug) vs 20 (good) ===");
    println!("Hypothesis: at dt small, J slightly exceeds intrinsic at deep ITM,");
    println!("so f(0) = excess - J < 0 → find_boundary_for_state returns x_grid[0]");
    println!("→ exercise region collapses to single point → no EEP.");
    for &nt in &[20_usize, 30, 50] {
        let r = price_american_put_heston(
            100.0, SPOT, V0, T, R, 0.0, 5.0, V0, 1e-4, 0.0,
            60, 10, nt,
        );
        println!("\n--- n_time={nt}, price={:.4} ---", r.price);
        let mid = r.boundary_s[0].len() / 2;
        print!("k=0..{nt} boundary_s[mid]: ");
        for b in &r.boundary_s { print!("{:.2} ", b[mid]); }
        println!();
    }
}

// ---- 13. Grid refinement for European Heston call (q=0) vs COS ----
#[test]
#[ignore]
fn bias_13_euro_heston_grid_refinement() {
    use numerical_methods::slv::{
        CalibrationInstrument, CalibrationPricer, HestonEuropeanPricer,
        HestonParameters, OptionSide,
    };
    println!("\n=== TEST 13: Grid refinement — ATM Heston call (q=0) vs COS ===");
    let params = HestonParameters { v0: V0, v_bar: V0, kappa: 2.0, sigma: 0.3, rho: -0.7 };
    let pricer = HestonEuropeanPricer::default();
    let inst = CalibrationInstrument {
        symbol: "T".into(), side: OptionSide::Call, strike: 100.0, maturity_years: T,
        market_price: 0.0, spot: SPOT, risk_free_rate: R, dividend_yield: 0.0,
        implied_vol: None, bid: None, ask: None, open_interest: None, volume: None,
    };
    let cos = pricer.price(&inst, params).unwrap();
    println!("COS reference: {cos:.4}");
    println!("{:<6} {:<6} {:<10} {:<10}", "n_x", "m_v", "ctmc", "bias_vs_cos");
    for &(nx, mv) in &[(60_usize,10), (100,15), (150,25), (200,30), (300,40), (400,50)] {
        let p = price_american_call_heston(
            100.0, SPOT, params.v0, T, R, 0.0,
            params.kappa, params.v_bar, params.sigma, params.rho,
            nx, mv, 50,
        ).price;
        println!("{:<6} {:<6} {:<10.4} {:<+10.4}", nx, mv, p, p - cos);
    }
}

// ---- 12. Re-run test 5 (CTMC Euro call vs COS) under n_time sweep ----
#[test]
#[ignore]
fn bias_12_test5_under_n_time_sweep() {
    use numerical_methods::slv::{
        CalibrationInstrument, CalibrationPricer, HestonEuropeanPricer,
        HestonParameters, OptionSide,
    };
    println!("\n=== TEST 12: CTMC ATM call (q=0) vs Heston COS — n_time sweep ===");
    let params = HestonParameters { v0: V0, v_bar: V0, kappa: 2.0, sigma: 0.3, rho: -0.7 };
    let pricer = HestonEuropeanPricer::default();
    let inst = CalibrationInstrument {
        symbol: "T".into(), side: OptionSide::Call, strike: 100.0, maturity_years: T,
        market_price: 0.0, spot: SPOT, risk_free_rate: R, dividend_yield: 0.0,
        implied_vol: None, bid: None, ask: None, open_interest: None, volume: None,
    };
    let cos = pricer.price(&inst, params).unwrap();
    println!("COS Heston call: {cos:.4}");
    println!("{:<8} {:<10} {:<10}", "n_time", "ctmc", "bias");
    for &nt in &[5_usize, 10, 20, 30, 50, 100, 200] {
        let p = price_american_call_heston(
            100.0, SPOT, params.v0, T, R, 0.0,
            params.kappa, params.v_bar, params.sigma, params.rho,
            150, 25, nt,
        ).price;
        println!("{:<8} {:<10.4} {:<+10.4}", nt, p, p - cos);
    }
}

// ---- 11. n_time sweep with DENSE Padé (force via nm<=600) ----
#[test]
#[ignore]
fn bias_11_n_time_sweep_dense_pade() {
    println!("\n=== TEST 11: n_time sweep with DENSE Pade (nm=60*10=600) ===");
    println!("If anomaly persists here, bias is NOT Krylov-specific.");
    println!("Reference binom_am = {:.4}", binomial_am_put(100.0));
    println!("{:<8} {:<10} {:<10}", "n_time", "ctmc", "bias_vs_binom");
    for &nt in &[5_usize, 10, 20, 30, 50, 80, 120, 200, 400] {
        let p = price_american_put_heston(
            100.0, SPOT, V0, T, R, 0.0,
            5.0, V0, 1e-4, 0.0,
            60, 10, nt,
        ).price;
        println!("{:<8} {:<10.4} {:<+10.4}", nt, p, p - binomial_am_put(100.0));
    }
}

// ---- 10. n_time sweep to confirm anomaly: more time steps → less accurate? ----
#[test]
#[ignore]
fn bias_10_n_time_sweep_confirm_anomaly() {
    println!("\n=== TEST 10: n_time sweep, ATM put, degenerate Heston ===");
    println!("Fixed grid (120, 20). Reference binom_am = {:.4}", binomial_am_put(100.0));
    println!("{:<8} {:<10} {:<10}", "n_time", "ctmc", "bias_vs_binom");
    for &nt in &[5_usize, 10, 20, 30, 50, 80, 120, 200, 400] {
        let p = price_american_put_heston(
            100.0, SPOT, V0, T, R, 0.0,
            5.0, V0, 1e-4, 0.0,
            120, 20, nt,
        ).price;
        println!("{:<8} {:<10.4} {:<+10.4}", nt, p, p - binomial_am_put(100.0));
    }
}

// ---- 8. Diagnose where EEP disappears: inspect boundary trajectory ----
#[test]
#[ignore]
fn bias_8_boundary_trajectory_diagnostic() {
    println!("\n=== TEST 8: Boundary trajectory for ATM American put (degenerate Heston) ===");
    println!("If boundary collapses to x_grid[0] or to K throughout → no EEP.");
    println!("Expected: boundary S starts near K at t=T, drifts down toward ~0.8*K at t=0.");

    let result = price_american_put_heston(
        100.0, SPOT, V0, T, R, 0.0,
        5.0, V0, 1e-4, 0.0,
        120, 20, 20, // fewer time steps to print
    );
    println!("Price: {:.4}", result.price);
    let mid = result.boundary_s[0].len() / 2;
    println!("Boundary S at mid variance state (l={mid}) across time:");
    println!("{:<6} {:<10}", "k", "b_s(k,v_mid)");
    for (k, b) in result.boundary_s.iter().enumerate() {
        println!("{:<6} {:<10.4}", k, b[mid]);
    }
}

// ---- 9. Confirm: CTMC behaves as European for puts even at realistic sigma_v ----
#[test]
#[ignore]
fn bias_9_compare_ctmc_put_to_european_heston() {
    use numerical_methods::slv::{
        CalibrationInstrument, CalibrationPricer, HestonEuropeanPricer,
        HestonParameters, OptionSide,
    };
    println!("\n=== TEST 9: CTMC American put vs Heston European put (same params) ===");
    println!("If CTMC ≈ European (not > European by EEP), EEP term is broken.");
    let params = HestonParameters { v0: V0, v_bar: V0, kappa: 2.0, sigma: 0.3, rho: -0.7 };
    let pricer = HestonEuropeanPricer::default();
    println!("{:<8} {:<12} {:<12} {:<10}", "strike", "ctmc_am_put", "heston_euro", "EEP_est");

    for &k in &[80.0_f64, 90.0, 100.0, 110.0, 120.0] {
        let ctmc = price_american_put_heston(
            k, SPOT, params.v0, T, R, 0.0,
            params.kappa, params.v_bar, params.sigma, params.rho,
            150, 25, 100,
        ).price;
        let inst = CalibrationInstrument {
            symbol: "T".into(), side: OptionSide::Put, strike: k, maturity_years: T,
            market_price: 0.0, spot: SPOT, risk_free_rate: R, dividend_yield: 0.0,
            implied_vol: None, bid: None, ask: None, open_interest: None, volume: None,
        };
        let euro = pricer.price(&inst, params).unwrap();
        println!("{:<8.2} {:<12.4} {:<12.4} {:<+10.4}", k, ctmc, euro, ctmc - euro);
    }
}

#[allow(dead_code)]
fn _keep_imports_alive() {
    let _ = (binomial_euro_put, price_american_option_heston);
}
