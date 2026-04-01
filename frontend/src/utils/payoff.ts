import type { PriceRequest, PriceResponse, SingleOptionRequest } from "../types/index.ts";

export interface PayoffData {
  spotPrices: number[];
  payoffs: number[];
  breakEvenPoints: number[];
  maxProfit: number | null;
  maxLoss: number | null;
}

export function calculatePayoffGrid(
  request: PriceRequest,
  priceResult: PriceResponse
): PayoffData | null {
  if (priceResult.payoff_curve) {
    const { spot_prices, payoffs } = priceResult.payoff_curve;
    const breakEvenPoints = findBreakEvenPoints(spot_prices, payoffs);
    const { maxProfit, maxLoss } = calculateMaxProfitLoss(payoffs);
    return {
      spotPrices: spot_prices,
      payoffs,
      breakEvenPoints,
      maxProfit,
      maxLoss,
    };
  }

  if (request.structure_type === "exotic") {
    return null;
  }

  const bsPrice = priceResult.pricing.black_scholes;
  const premium = Math.abs(bsPrice);

  if (request.structure_type === "single") {
    return calculateSingleOptionPayoff(request, premium);
  }

  if (request.structure_type === "spread") {
    return calculateSpreadPayoff(request, premium);
  }

  return null;
}

function calculateSingleOptionPayoff(
  request: SingleOptionRequest,
  premium: number
): PayoffData {
  const { option_type, strike_price, spot_price, direction } = request;
  const dirMultiplier = direction === "long" ? 1 : -1;

  const minSpot = spot_price * 0.5;
  const maxSpot = spot_price * 1.5;
  const numPoints = 100;
  const step = (maxSpot - minSpot) / (numPoints - 1);

  const spotPrices: number[] = [];
  const payoffs: number[] = [];

  for (let i = 0; i < numPoints; i++) {
    const spot = minSpot + i * step;
    spotPrices.push(spot);

    let intrinsicValue = 0;
    if (option_type === "call") {
      intrinsicValue = Math.max(spot - strike_price, 0);
    } else {
      intrinsicValue = Math.max(strike_price - spot, 0);
    }

    const payoff = dirMultiplier * (intrinsicValue - premium);
    payoffs.push(payoff);
  }

  const breakEvenPoints = findBreakEvenPoints(spotPrices, payoffs);
  const { maxProfit, maxLoss } = calculateMaxProfitLoss(payoffs);

  return { spotPrices, payoffs, breakEvenPoints, maxProfit, maxLoss };
}

function calculateSpreadPayoff(
  request: Extract<PriceRequest, { structure_type: "spread" }>,
  totalPremium: number
): PayoffData {
  const {
    spread_type,
    strikes,
    spot_price,
    direction,
  } = request;

  const dirMultiplier = direction === "long" ? 1 : -1;
  const minSpot = spot_price * 0.5;
  const maxSpot = spot_price * 1.5;
  const numPoints = 100;
  const step = (maxSpot - minSpot) / (numPoints - 1);

  const spotPrices: number[] = [];
  const payoffs: number[] = [];

  for (let i = 0; i < numPoints; i++) {
    const spot = minSpot + i * step;
    spotPrices.push(spot);

    let payoff = 0;

    switch (spread_type) {
      case "straddle": {
        const K = strikes.strike ?? spot_price;
        payoff = Math.max(spot - K, 0) + Math.max(K - spot, 0) - totalPremium;
        break;
      }

      case "strangle": {
        const Kc = strikes.strike_call ?? spot_price * 1.05;
        const Kp = strikes.strike_put ?? spot_price * 0.95;
        payoff =
          Math.max(spot - Kc, 0) + Math.max(Kp - spot, 0) - totalPremium;
        break;
      }

      case "strip": {
        const K = strikes.strike ?? spot_price;
        payoff =
          Math.max(spot - K, 0) + 2 * Math.max(K - spot, 0) - totalPremium;
        break;
      }

      case "strap": {
        const K = strikes.strike ?? spot_price;
        payoff =
          2 * Math.max(spot - K, 0) + Math.max(K - spot, 0) - totalPremium;
        break;
      }

      case "synthetic_stock": {
        const K = strikes.strike ?? spot_price;
        payoff = Math.max(spot - K, 0) - Math.max(K - spot, 0);
        break;
      }

      case "bull_spread_call": {
        const Klow = strikes.strike_low ?? spot_price * 0.95;
        const Khigh = strikes.strike_high ?? spot_price * 1.05;
        payoff =
          Math.max(spot - Klow, 0) - Math.max(spot - Khigh, 0) - totalPremium;
        break;
      }

      case "bull_spread_put": {
        const Klow = strikes.strike_low ?? spot_price * 0.95;
        const Khigh = strikes.strike_high ?? spot_price * 1.05;
        payoff =
          Math.max(Klow - spot, 0) - Math.max(Khigh - spot, 0) - totalPremium;
        break;
      }

      case "bear_spread_call": {
        const Klow = strikes.strike_low ?? spot_price * 0.95;
        const Khigh = strikes.strike_high ?? spot_price * 1.05;
        payoff =
          Math.max(spot - Khigh, 0) - Math.max(spot - Klow, 0) - totalPremium;
        break;
      }

      case "bear_spread_put": {
        const Klow = strikes.strike_low ?? spot_price * 0.95;
        const Khigh = strikes.strike_high ?? spot_price * 1.05;
        payoff =
          Math.max(Khigh - spot, 0) - Math.max(Klow - spot, 0) - totalPremium;
        break;
      }

      case "butterfly": {
        const Klow = strikes.strike_low ?? spot_price * 0.95;
        const Kmid = strikes.strike_medium ?? spot_price;
        const Khigh = strikes.strike_high ?? spot_price * 1.05;
        payoff =
          Math.max(spot - Klow, 0) -
          2 * Math.max(spot - Kmid, 0) +
          Math.max(spot - Khigh, 0) -
          totalPremium;
        break;
      }

      case "calendar_spread": {
        payoff = 0;
        break;
      }
    }

    payoffs.push(dirMultiplier * payoff);
  }

  const breakEvenPoints = findBreakEvenPoints(spotPrices, payoffs);
  const { maxProfit, maxLoss } = calculateMaxProfitLoss(payoffs);

  return { spotPrices, payoffs, breakEvenPoints, maxProfit, maxLoss };
}

function findBreakEvenPoints(spotPrices: number[], payoffs: number[]): number[] {
  const breakEvenPoints: number[] = [];

  for (let i = 1; i < payoffs.length; i++) {
    const prev = payoffs[i - 1];
    const curr = payoffs[i];

    if ((prev <= 0 && curr >= 0) || (prev >= 0 && curr <= 0)) {
      const spotPrev = spotPrices[i - 1];
      const spotCurr = spotPrices[i];
      const interpolated = spotPrev + ((0 - prev) / (curr - prev)) * (spotCurr - spotPrev);
      breakEvenPoints.push(interpolated);
    }
  }

  return breakEvenPoints;
}

function calculateMaxProfitLoss(payoffs: number[]): {
  maxProfit: number | null;
  maxLoss: number | null;
} {
  const maxProfit = Math.max(...payoffs);
  const maxLoss = Math.min(...payoffs);

  return {
    maxProfit: maxProfit > 0 ? maxProfit : null,
    maxLoss: maxLoss < 0 ? maxLoss : null,
  };
}
