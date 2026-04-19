import type { PriceRequest, PriceResponse, SingleOptionRequest, SpreadRequest } from "../types/index.ts";

export type PayoffPoint = {
  spot: number;
  pnl: number;
};

export function calculatePayoffGrid(req: PriceRequest | null, result: PriceResponse | null): PayoffPoint[] {
  if (!req || !result) return [];

  let S = 100;
  let direction: "long" | "short" = "long";

  if (req.structure_type === "single") {
    const singleReq = req as SingleOptionRequest;
    S = singleReq.spot_price;
    direction = singleReq.direction;
  } else if (req.structure_type === "spread") {
    const spreadReq = req as SpreadRequest;
    S = spreadReq.spot_price;
    direction = spreadReq.direction;
  } else {
    return [];
  }

  const points: PayoffPoint[] = [];
  const lo = S * 0.6;
  const hi = S * 1.4;
  const n = 80;

  for (let i = 0; i <= n; i++) {
    const spot = lo + (hi - lo) * i / n;
    let intrinsic = 0;

    if (req.structure_type === "single") {
      const singleReq = req as SingleOptionRequest;
      const K = singleReq.strike_price;
      const type = singleReq.option_type;
      intrinsic = type === "call" ? Math.max(spot - K, 0) : Math.max(K - spot, 0);
    } else if (req.structure_type === "spread") {
      const spreadReq = req as SpreadRequest;
      const strikes = Object.values(spreadReq.strikes).filter(v => v != null) as number[];
      if (strikes.length >= 2) {
        strikes.sort((a, b) => a - b);
        const K1 = strikes[0];
        const K2 = strikes[1];
        intrinsic = Math.max(spot - K1, 0) - Math.max(spot - K2, 0);
      }
    }

    const premium = typeof result.pricing?.black_scholes === "number"
      ? result.pricing.black_scholes
      : 0;

    const pnl = direction === "long" ? intrinsic - premium : premium - intrinsic;

    points.push({ spot, pnl });
  }

  return points;
}
