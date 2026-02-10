use crate::{BlackScholesPrice, Options, Payoff};

pub struct Position {
    pub option: Options,
    pub weight: f32,
}

impl Payoff for Position {
    #[inline(always)]
    fn compute(&self, spot_t: f64) -> f64 {
        self.option.compute(spot_t) * self.weight as f64
    }

    fn compute_static(_spot_t: f64, _strike: f64) -> f64 {
        panic!(
            "Position is a composite type containing dynamic Options enum - use compute() instance method instead"
        )
    }
}

impl crate::Greeks for Position {
    fn delta(&self, imply_vol: f64, spot_price: f64) -> f64 {
        self.option.delta(imply_vol, spot_price) * self.weight as f64
    }

    fn theta(&self, imply_vol: f64, spot_price: f64) -> f64 {
        self.option.theta(imply_vol, spot_price) * self.weight as f64
    }

    fn gamma(&self, imply_vol: f64, spot_price: f64) -> f64 {
        self.option.gamma(imply_vol, spot_price) * self.weight as f64
    }

    fn vega(&self, spot_price: f64) -> f64 {
        self.option.vega(spot_price) * self.weight as f64
    }

    fn rho(&self, imply_vol: f64, spot_price: f64, interest_rate: f64) -> f64 {
        self.option.rho(imply_vol, spot_price, interest_rate) * self.weight as f64
    }
}

impl BlackScholesPrice for Position {
    fn bs_pricing(&self) -> f64 {
        Options::bs_pricing(&self.option) * self.weight as f64
    }
}
pub enum Direction {
    SHORT,
    LONG,
}
pub struct OptionSpreads {
    components: Vec<Position>,
}

impl OptionSpreads {
    pub fn new(components: Vec<Position>) -> Self {
        OptionSpreads { components }
    }

    pub fn components(&self) -> &[Position] {
        &self.components
    }

    pub fn new_straddle(
        direction: Direction,
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let atm_call: Options = Options::new_call(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let atm_put: Options = Options::new_put(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let mut spreads: Vec<Position> = Vec::new();

        match direction {
            Direction::LONG => {
                spreads.push(Position {
                    option: atm_call,
                    weight: 1.0,
                });
                spreads.push(Position {
                    option: atm_put,
                    weight: 1.0,
                });
            }
            Direction::SHORT => {
                spreads.push(Position {
                    option: atm_call,
                    weight: -1.0,
                });
                spreads.push(Position {
                    option: atm_put,
                    weight: -1.0,
                });
            }
        }

        OptionSpreads {
            components: spreads,
        }
    }

    pub fn new_strangle(
        direction: Direction,
        strike_price_call: f64,
        strike_price_put: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let otm_call: Options = Options::new_call(
            strike_price_call,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let otm_put: Options = Options::new_put(
            strike_price_put,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let mut spreads: Vec<Position> = Vec::new();

        match direction {
            Direction::LONG => {
                spreads.push(Position {
                    option: otm_call,
                    weight: 1.0,
                });
                spreads.push(Position {
                    option: otm_put,
                    weight: 1.0,
                });
            }
            Direction::SHORT => {
                spreads.push(Position {
                    option: otm_call,
                    weight: -1.0,
                });
                spreads.push(Position {
                    option: otm_put,
                    weight: -1.0,
                });
            }
        }

        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_strip(
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let call: Options = Options::new_call(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let put: Options = Options::new_put(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let spreads = vec![
            Position {
                option: call,
                weight: 1.0,
            },
            Position {
                option: put,
                weight: 2.0,
            },
        ];
        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_strap(
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let call: Options = Options::new_call(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let put: Options = Options::new_put(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let spreads = vec![
            Position {
                option: call,
                weight: 2.0,
            },
            Position {
                option: put,
                weight: 1.0,
            },
        ];
        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_synthetic_stock(
        direction: Direction,
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let call: Options = Options::new_call(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let put: Options = Options::new_put(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let mut spreads: Vec<Position> = Vec::new();

        match direction {
            Direction::LONG => {
                spreads.push(Position {
                    option: call,
                    weight: 1.0,
                });
                spreads.push(Position {
                    option: put,
                    weight: -1.0,
                });
            }
            Direction::SHORT => {
                spreads.push(Position {
                    option: call,
                    weight: -1.0,
                });
                spreads.push(Position {
                    option: put,
                    weight: 1.0,
                });
            }
        }

        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_bull_spread_call(
        strike_price_low: f64,
        strike_price_high: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let call_1: Options = Options::new_call(
            strike_price_low,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let call_2: Options = Options::new_call(
            strike_price_high,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let spreads = vec![
            Position {
                option: call_1,
                weight: 1.0,
            },
            Position {
                option: call_2,
                weight: -1.0,
            },
        ];
        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_bull_spread_put(
        strike_price_low: f64,
        strike_price_high: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let put_1: Options = Options::new_put(
            strike_price_low,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let put_2: Options = Options::new_put(
            strike_price_high,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let spreads = vec![
            Position {
                option: put_1,
                weight: 1.0,
            },
            Position {
                option: put_2,
                weight: -1.0,
            },
        ];
        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_bear_spread_call(
        strike_price_low: f64,
        strike_price_high: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let call_1: Options = Options::new_put(
            strike_price_low,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let call_2: Options = Options::new_put(
            strike_price_high,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let spreads = vec![
            Position {
                option: call_1,
                weight: -1.0,
            },
            Position {
                option: call_2,
                weight: 1.0,
            },
        ];
        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_bear_spread_put(
        strike_price_low: f64,
        strike_price_high: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let put_1: Options = Options::new_put(
            strike_price_low,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let put_2: Options = Options::new_put(
            strike_price_high,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let spreads = vec![
            Position {
                option: put_1,
                weight: -1.0,
            },
            Position {
                option: put_2,
                weight: 1.0,
            },
        ];
        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_butterfly(
        direction: Direction,
        strike_price_low: f64,
        strike_price_high: f64,
        strike_price_medium: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        time_to_maturity: f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let call_1: Options = Options::new_call(
            strike_price_low,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let call_2: Options = Options::new_call(
            strike_price_medium,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let call_3: Options = Options::new_call(
            strike_price_high,
            spot_price,
            volatility,
            risk_free_rate,
            time_to_maturity,
            dividend_yield,
        );
        let mut spreads: Vec<Position> = Vec::new();

        match direction {
            Direction::LONG => {
                spreads.push(Position {
                    option: call_1,
                    weight: 1.0,
                });
                spreads.push(Position {
                    option: call_2,
                    weight: -2.0,
                });
                spreads.push(Position {
                    option: call_3,
                    weight: 1.0,
                });
            }
            Direction::SHORT => {
                spreads.push(Position {
                    option: call_1,
                    weight: -1.0,
                });
                spreads.push(Position {
                    option: call_2,
                    weight: 2.0,
                });
                spreads.push(Position {
                    option: call_3,
                    weight: -1.0,
                });
            }
        }
        OptionSpreads {
            components: spreads,
        }
    }
    pub fn new_calendar_spread(
        direction: Direction,
        strike_price: f64,
        spot_price: f64,
        volatility: f64,
        risk_free_rate: f64,
        short_term_maturity: f64,
        long_term_maturity:f64,
        dividend_yield: Option<f64>,
    ) -> OptionSpreads {
        let call_st: Options = Options::new_call(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            short_term_maturity,
            dividend_yield,
        );
        let call_lt: Options = Options::new_put(
            strike_price,
            spot_price,
            volatility,
            risk_free_rate,
            long_term_maturity,
            dividend_yield,
        );
        let mut spreads: Vec<Position> = Vec::new();

        match direction {
            Direction::LONG => {
                spreads.push(Position {
                    option: call_st,
                    weight: -1.0,
                });
                spreads.push(Position {
                    option: call_lt,
                    weight: 1.0,
                });
            }
            Direction::SHORT => {
                spreads.push(Position {
                    option: call_st,
                    weight: 1.0,
                });
                spreads.push(Position {
                    option: call_lt,
                    weight: -1.0,
                });
            }
        }

        OptionSpreads {
            components: spreads,
        }
    }
}

impl BlackScholesPrice for OptionSpreads {
    fn bs_pricing(&self) -> f64 {
        self.components.iter().map(|pos| pos.bs_pricing()).sum()
    }
}

impl crate::Payoff for OptionSpreads {
    #[inline(always)]
    fn compute(&self, spot_t: f64) -> f64 {
        self.components
            .iter()
            .map(|pos| {
                pos.option.compute(spot_t) * pos.weight as f64
            })
            .sum()
    }

    fn compute_static(_spot_t: f64, _strike: f64) -> f64 {
        panic!(
            "OptionSpreads is a composite type with dynamic components - use compute() instance method instead"
        )
    }
}
pub trait MonteCarloParameters {
    fn spot_price(&self) -> f64;
    fn volatility(&self) -> f64;
    fn risk_free_rate(&self) -> f64;
    fn time_to_maturity(&self) -> f64;
}

impl MonteCarloParameters for OptionSpreads {
    fn spot_price(&self) -> f64 {
        // Extract from first component (all components share these parameters)
        self.components[0].option.spot_price()
    }

    fn volatility(&self) -> f64 {
        self.components[0].option.volatility()
    }

    fn risk_free_rate(&self) -> f64 {
        self.components[0].option.risk_free_rate()
    }

    fn time_to_maturity(&self) -> f64 {
        self.components[0].option.time_to_maturity()
    }
}

// Greeks
impl crate::Greeks for OptionSpreads {
    fn delta(&self, imply_vol: f64, spot_price: f64) -> f64 {
        self.components
            .iter()
            .map(|pos| pos.delta(imply_vol, spot_price))
            .sum()
    }

    fn theta(&self, imply_vol: f64, spot_price: f64) -> f64 {
        self.components
            .iter()
            .map(|pos| pos.theta(imply_vol, spot_price))
            .sum()
    }

    fn gamma(&self, imply_vol: f64, spot_price: f64) -> f64 {
        self.components
            .iter()
            .map(|pos| pos.gamma(imply_vol, spot_price))
            .sum()
    }

    fn vega(&self, spot_price: f64) -> f64 {
        self.components.iter().map(|pos| pos.vega(spot_price)).sum()
    }

    fn rho(&self, imply_vol: f64, spot_price: f64, interest_rate: f64) -> f64 {
        self.components
            .iter()
            .map(|pos| pos.rho(imply_vol, spot_price, interest_rate))
            .sum()
    }
}
