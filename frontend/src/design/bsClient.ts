// Tiny client-side Black-Scholes — only used for greek sparkline sensitivity curves.
// Real pricing comes from the backend.

export type BSArgs = {
  S: number;
  K: number;
  r: number;
  q: number;
  sigma: number;
  T: number;
  type: "call" | "put";
};

function ndist(x: number): number {
  return Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI);
}

function ncdf(x: number): number {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.sqrt(2);
  const t = 1.0 / (1.0 + p * x);
  const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return 0.5 * (1.0 + sign * y);
}

export type GreeksVec = {
  price: number;
  delta: number;
  gamma: number;
  vega: number;
  theta: number;
  rho: number;
};

export function bs({ S, K, r, q, sigma, T, type }: BSArgs): GreeksVec {
  if (T <= 0 || sigma <= 0) {
    const intrinsic = type === "call" ? Math.max(S - K, 0) : Math.max(K - S, 0);
    return { price: intrinsic, delta: 0, gamma: 0, vega: 0, theta: 0, rho: 0 };
  }
  const d1 = (Math.log(S / K) + (r - q + sigma * sigma / 2) * T) / (sigma * Math.sqrt(T));
  const d2 = d1 - sigma * Math.sqrt(T);
  const Nd1 = ncdf(d1), Nd2 = ncdf(d2);
  const nd1 = ndist(d1);
  let price: number, delta: number, theta: number, rho: number;
  if (type === "call") {
    price = S * Math.exp(-q * T) * Nd1 - K * Math.exp(-r * T) * Nd2;
    delta = Math.exp(-q * T) * Nd1;
    rho = K * T * Math.exp(-r * T) * Nd2;
    theta = -(S * Math.exp(-q * T) * nd1 * sigma) / (2 * Math.sqrt(T))
      - r * K * Math.exp(-r * T) * Nd2 + q * S * Math.exp(-q * T) * Nd1;
  } else {
    price = K * Math.exp(-r * T) * ncdf(-d2) - S * Math.exp(-q * T) * ncdf(-d1);
    delta = -Math.exp(-q * T) * ncdf(-d1);
    rho = -K * T * Math.exp(-r * T) * ncdf(-d2);
    theta = -(S * Math.exp(-q * T) * nd1 * sigma) / (2 * Math.sqrt(T))
      + r * K * Math.exp(-r * T) * ncdf(-d2) - q * S * Math.exp(-q * T) * ncdf(-d1);
  }
  const gamma = Math.exp(-q * T) * nd1 / (S * sigma * Math.sqrt(T));
  const vega = S * Math.exp(-q * T) * nd1 * Math.sqrt(T) / 100;
  return { price, delta, gamma, vega, theta: theta / 365, rho: rho / 100 };
}

export type GreekKey = "delta" | "gamma" | "theta" | "vega" | "rho";

export function greekCurve(key: GreekKey, args: BSArgs): { x: number; y: number }[] {
  const pts: { x: number; y: number }[] = [];
  const lo = args.S * 0.7, hi = args.S * 1.3;
  const n = 40;
  for (let i = 0; i <= n; i++) {
    const s = lo + (hi - lo) * i / n;
    const g = bs({ ...args, S: s });
    pts.push({ x: s, y: g[key] });
  }
  return pts;
}

export function payoffCurve(
  args: { S: number; K: number; type: "call" | "put"; direction: "long" | "short" },
  priceAtEntry: number
): { spot: number; pnl: number }[] {
  const { S, K, type, direction } = args;
  const pts: { spot: number; pnl: number }[] = [];
  const lo = S * 0.6, hi = S * 1.4;
  const n = 80;
  for (let i = 0; i <= n; i++) {
    const spot = lo + (hi - lo) * i / n;
    const intrinsic = type === "call" ? Math.max(spot - K, 0) : Math.max(K - spot, 0);
    const pnl = direction === "long" ? intrinsic - priceAtEntry : priceAtEntry - intrinsic;
    pts.push({ spot, pnl });
  }
  return pts;
}
