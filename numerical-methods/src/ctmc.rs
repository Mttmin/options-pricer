
fn omega(S: f64, v: f64) -> f64 {
    todo!()
}

fn m(v:f64) -> f64{
    todo!()
}

fn gamma(S:f64) -> f64{
    todo!()
}

fn mu(v:f64) -> f64{
    todo!()
}
fn sigma(v:f64) -> f64{
    todo!()
}

fn stochastic_vol_model(S: f64, v: f64)-> (f64, f64) {
    let delta_t = 0.02;
    let dWt1 = 0.01;
    let dWt2 = 0.02;
    let ds = omega(S,v)*delta_t + m(v)*gamma(S)*dWt1;
    let dv = mu(v)*delta_t + sigma(v)*dWt2;
    (ds,dv)
}