use crate::{Call};

pub enum ExoticOptions {
    ConvertibleBond(ConvertibleBond),
}

#[derive(Debug, Clone, Copy)]
pub struct ConvertibleBond {
    // Bond parameters
    pub face_value: f64,
    pub coupon_rate: f64,
    pub maturity: f64,
    pub payment_frequency: u32,
    pub credit_spread: f64,

    // Conversion option parameters
    pub conversion_price: f64,
    pub stock_price: f64,
    pub volatility: f64,
    pub time_to_maturity: f64,
    pub risk_free_rate: f64,

}

impl ConvertibleBond {
    fn npv(&self) -> f64 {
        // Placeholder for NPV calculation logic
        let periods = (self.maturity * self.payment_frequency as f64) as u32;
        let coupon_payment = self.face_value * self.coupon_rate / self.payment_frequency as f64;
        let mut npv = 0.0;
        for t in 1..=periods {
            let time = t as f64 / self.payment_frequency as f64;
            npv += coupon_payment / (1.0 + self.risk_free_rate + self.credit_spread).powf(time);
        }
        npv += self.face_value / (1.0 + self.risk_free_rate + self.credit_spread).powf(self.maturity);
        npv
    }
    fn conversion_option_price(&self) -> f64 {
        // Placeholder for conversion option pricing logic
        let underlying_call = Call {
            strike_price: self.conversion_price,
            spot_price: self.stock_price,
            volatility: self.volatility,
            risk_free_rate: self.risk_free_rate,
            time_to_maturity: self.time_to_maturity,
        };
        underlying_call.bs_pricing()
    }
    pub fn bs_pricing(&self) -> f64 {
        self.npv() + self.face_value / self.conversion_price * self.conversion_option_price()
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convertible_bond_pricing() {
        let cb = ConvertibleBond {
            face_value: 1000.0,
            coupon_rate: 0.05,
            maturity: 5.0,
            payment_frequency: 2,
            credit_spread: 0.02,
            risk_free_rate: 0.03,
            conversion_price: 50.0,
            stock_price: 55.0,
            volatility: 0.2,
            time_to_maturity: 5.0,
        };
        let price = cb.bs_pricing();
        println!("Convertible Bond Price: {:.4}", price);
        assert!((price - 1320.66).abs() < 1e-2); // expected value
    }
}