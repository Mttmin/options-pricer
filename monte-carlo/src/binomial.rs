use options::Payoff;

pub struct BinomialParameters {
    pub spot_price: f64,
    pub depth: u16,
    pub u: f64,
    pub d: f64,
    pub p: f64,
    pub q: f64,
    pub discount: f64,
}

impl BinomialParameters {
    pub fn new(spot_price: f64, depth: u16, time_to_maturity: f64, volatility: f64, risk_free_rate: f64) -> Self {
        let dt: f64 = time_to_maturity / f64::from(depth);
        let u = (volatility * dt.sqrt()).exp();
        let d = 1.0 / u;
        let discount = (-risk_free_rate * dt).exp();
        let p = ((risk_free_rate * dt).exp() - d) / (u - d);
        let q = 1.0 - p;
        
        Self { spot_price, depth, u, d, p, q, discount }
    }
}
 
fn binomial_european<P: Payoff>(product: &P, parameters: &BinomialParameters) -> f64 {
    let depth = parameters.depth;
    let mut values: Vec<f64> = Vec::with_capacity((depth + 1) as usize); // Buffer
    
    // Terminal payoffs: node i has i down moves, (depth-i) up moves
    let s0 = parameters.spot_price;
    let u_power_n = parameters.u.powi(depth as i32);
    let ratio = parameters.d / parameters.u;
    
    let mut price = s0 * u_power_n;
    for _ in 0..=depth {
        values.push(product.compute(price ));
        price *= ratio;
    }
    // Backward induction
    for step in (0..depth).rev() {
        let num_nodes = (step + 1) as usize;
        for j in 0..num_nodes {
            values[j] = parameters.discount * (parameters.p * values[j] + parameters.q * values[j + 1]);
        }
    }
    values[0]
}

fn binomial_american<P: Payoff>(product: &P, parameters: &BinomialParameters) -> f64 {
    let depth = parameters.depth;
    
    // Buffer
    let mut values: Vec<f64> = Vec::with_capacity((depth + 1) as usize);
    
    // Terminal payoffs
    let s0 = parameters.spot_price;
    let u_power_n = parameters.u.powi(depth as i32);
    let ratio = parameters.d / parameters.u;
    
    let mut price = s0 * u_power_n;
    for _ in 0..=depth {
        values.push(product.compute(price));
        price *= ratio;
    }
    
    // Backward induction with early exercise check
    for step in (0..depth).rev() {
        let num_nodes = (step + 1) as usize;
        let mut stock_price = s0 * parameters.u.powi(step as i32);
        
        for j in 0..num_nodes {
            let continuation = parameters.discount * (parameters.p * values[j] + parameters.q * values[j + 1]);
            let early_exercise = product.compute(stock_price);
            values[j] = continuation.max(early_exercise);
            stock_price *= ratio;
        }
    }
    
    values[0]
}

#[derive(Clone, Copy)]
pub enum ExerciseStyle {
    European,
    American,
}

pub fn binomial_price<P: Payoff>(product: &P, parameters: &BinomialParameters, style: ExerciseStyle) -> f64 {
    match style {
        ExerciseStyle::European => binomial_european(product, parameters),
        ExerciseStyle::American => binomial_american(product, parameters),
    }
}