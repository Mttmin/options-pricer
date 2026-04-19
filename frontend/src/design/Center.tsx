import { useMemo } from "react";
import type { PricePoint } from "../api/client.ts";
import type {
  CtmcPriceResponse,
  Direction,
  ExerciseStyle,
  GreeksResult,
  OptionType,
  PriceResponse,
  PricingResult,
  PricingTimings,
} from "../types/index.ts";
import { LineChart, PayoffChart, Sparkline } from "./Charts.tsx";
import { greekCurve } from "./bsClient.ts";

export type PriceRow = {
  id: string;
  method: string;
  sub: string;
  price: number | null;
  ms: number | null;
  detail?: string;
  dim?: boolean;
  loading?: boolean;
};

function fmt(v: number | null | undefined, d = 4): string {
  if (v == null || isNaN(v)) return "—";
  return v.toFixed(d);
}

function fmtMs(v: number | null | undefined): string {
  if (v == null) return "—";
  return v < 1 ? `${(v * 1000).toFixed(0)}µs` : `${v.toFixed(1)}ms`;
}

export function buildRows(
  result: PriceResponse | null,
  exerciseStyle: ExerciseStyle,
  ctmc: CtmcPriceResponse | null,
  ctmcLoading: boolean
): PriceRow[] {
  const rows: PriceRow[] = [];
  if (!result) {
    if (ctmc || ctmcLoading) {
      rows.push({
        id: "ctmc",
        method: "CTMC (Heston)",
        sub: "Stochastic vol · calibrated",
        price: ctmc?.price ?? null,
        ms: ctmc?.timing_ms ?? null,
        loading: ctmcLoading,
        detail: ctmc ? `IV RMSE ${(ctmc.ivrmse * 100).toFixed(2)}vp · ${ctmc.n_instruments} instr.` : undefined,
      });
    }
    return rows;
  }

  const p: PricingResult = result.pricing;
  const timings: Partial<PricingTimings> = p.timings ?? {};

  if (exerciseStyle === "european") {
    rows.push({
      id: "bs",
      method: "Black-Scholes",
      sub: "Analytical · GBM",
      price: p.black_scholes,
      ms: timings.black_scholes_ms ?? null,
    });
    if (p.monte_carlo) {
      rows.push({
        id: "mc",
        method: "Monte Carlo",
        sub: "Antithetic variates",
        price: p.monte_carlo.price,
        ms: timings.monte_carlo_ms ?? null,
        detail: `CI [${fmt(p.monte_carlo.ci_lower, 3)}, ${fmt(p.monte_carlo.ci_upper, 3)}] · SE ${p.monte_carlo.std_error.toExponential(1)}`,
      });
    }
    if (p.binomial_european != null) {
      rows.push({
        id: "bin",
        method: "Binomial",
        sub: "CRR tree",
        price: p.binomial_european,
        ms: timings.binomial_european_ms ?? null,
      });
    }
  } else {
    if (p.binomial_american != null) {
      rows.push({
        id: "bin_am",
        method: "Binomial (American)",
        sub: "Backward induction · early exercise",
        price: p.binomial_american,
        ms: timings.binomial_american_ms ?? null,
      });
    }
    if (p.penalty_solver) {
      const d = p.penalty_solver.diagnostics;
      rows.push({
        id: "pde",
        method: "PDE Penalty",
        sub: "Rannacher smoothing",
        price: p.penalty_solver.price,
        ms: timings.penalty_solver_ms ?? null,
        detail: `${d.spatial_nodes} nodes × ${d.timesteps} steps · ${d.avg_iterations_per_step.toFixed(1)} iters/step`,
      });
    }
    if (p.bs_american_approx != null) {
      rows.push({
        id: "bs_approx",
        method: "BS Approx",
        sub: "Black's formula",
        price: p.bs_american_approx,
        ms: timings.bs_american_approx_ms ?? null,
      });
    }
    if (p.binomial_european != null) {
      rows.push({
        id: "bin_eu",
        method: "Binomial (Euro)",
        sub: "For comparison",
        price: p.binomial_european,
        ms: timings.binomial_european_ms ?? null,
        dim: true,
      });
    }
  }

  if (ctmc || ctmcLoading) {
    rows.push({
      id: "ctmc",
      method: "CTMC (Heston)",
      sub: "Stochastic vol · calibrated",
      price: ctmc?.price ?? null,
      ms: ctmc?.timing_ms ?? null,
      loading: ctmcLoading,
      detail: ctmc ? `IV RMSE ${(ctmc.ivrmse * 100).toFixed(2)}vp · ${ctmc.n_instruments} instr.` : undefined,
    });
  }

  return rows;
}

export function HeroPrice({
  label, primary, secondary, ctmc, pinned, onPin, onClearPin,
}: {
  label: string;
  primary: { method: string; price: number | null; ms: number | null };
  secondary: { method: string; price: number | null }[];
  ctmc: CtmcPriceResponse | null;
  pinned: { price: number; method: string } | null;
  onPin: () => void;
  onClearPin: () => void;
}) {
  return (
    <div className="hero-price">
      <div>
        <div className="label">{label}</div>
        <div className="number tnum mono">${fmt(primary.price, 4)}</div>
        <div className="sub">{primary.method} · {fmtMs(primary.ms)}</div>
      </div>
      <div className="altmethods">
        {secondary.map((s, i) => (
          <div key={i} className="alt">
            <div className="mk">{s.method}</div>
            <div>${fmt(s.price, 4)}</div>
          </div>
        ))}
        {ctmc && (
          <div className="alt">
            <div className="mk">CTMC Heston</div>
            <div>${fmt(ctmc.price, 4)}</div>
          </div>
        )}
        {pinned && (
          <div className="alt pin">
            <div className="mk">Pinned</div>
            <div>${fmt(pinned.price, 4)}</div>
          </div>
        )}
        <div style={{ display: "flex", gap: 6 }}>
          {pinned ? (
            <button className="icon-btn" onClick={onClearPin} title="Clear snapshot">Clear pin</button>
          ) : (
            <button className="icon-btn" onClick={onPin} title="Freeze current price" disabled={primary.price == null}>Pin</button>
          )}
        </div>
      </div>
    </div>
  );
}

export function PricingRows({
  rows, selectedId, onSelect, emphasis,
}: {
  rows: PriceRow[];
  selectedId: string | null;
  onSelect: (id: string) => void;
  emphasis: "hero" | "table";
}) {
  return (
    <table className="pricing-table">
      <thead>
        <tr>
          <th>Method</th>
          <th style={{ textAlign: "right" }}>Price</th>
          <th style={{ textAlign: "right" }}>Detail</th>
          <th style={{ textAlign: "right" }}>Time</th>
        </tr>
      </thead>
      <tbody>
        {rows.map(r => (
          <tr
            key={r.id}
            className={emphasis === "hero" && r.id === selectedId ? "primary" : ""}
            style={{
              ...(r.dim ? { opacity: 0.55 } : {}),
              cursor: r.loading ? "default" : "pointer",
            }}
            onClick={() => { if (!r.loading) onSelect(r.id); }}
            title={r.loading ? "" : "Set as primary"}
          >
            <td className="method">
              {r.method}
              <small>{r.sub}</small>
            </td>
            <td className="price tnum">
              {r.loading ? <span className="dim">calibrating…</span> : (r.price != null ? `$${fmt(r.price, 4)}` : "—")}
            </td>
            <td className="detail">{r.detail || ""}</td>
            <td className="detail">{fmtMs(r.ms)}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}

export function GreeksStrip({
  greeks, bsArgs,
}: {
  greeks: GreeksResult;
  bsArgs: { S: number; K: number; r: number; q: number; sigma: number; T: number; type: OptionType };
}) {
  const curves = useMemo(() => {
    return {
      delta: greekCurve("delta", bsArgs),
      gamma: greekCurve("gamma", bsArgs),
      theta: greekCurve("theta", bsArgs),
      vega: greekCurve("vega", bsArgs),
      rho: greekCurve("rho", bsArgs),
    };
  }, [bsArgs]);

  const items: { key: keyof GreeksResult; name: string; glyph: string; desc: string }[] = [
    { key: "delta", name: "Delta", glyph: "Δ", desc: "per $1 spot" },
    { key: "gamma", name: "Gamma", glyph: "Γ", desc: "Δ per $1" },
    { key: "theta", name: "Theta", glyph: "Θ", desc: "per day" },
    { key: "vega", name: "Vega", glyph: "ν", desc: "per 1% vol" },
    { key: "rho", name: "Rho", glyph: "ρ", desc: "per 1% rate" },
  ];

  return (
    <div className="greeks">
      {items.map(it => {
        const v = greeks[it.key];
        const cls = v > 0 ? "pos" : v < 0 ? "neg" : "";
        const color = v > 0 ? "var(--accent)" : "var(--loss)";
        return (
          <div className="greek" key={it.key}>
            <div className="greek-head">
              <div className="greek-name">
                <span className="glyph">{it.glyph}</span>
                {it.name}
              </div>
              <div className={`greek-val tnum ${cls}`}>
                {v != null ? v.toFixed(it.key === "gamma" ? 5 : 4) : "—"}
              </div>
            </div>
            <div className="greek-spark">
              <Sparkline data={curves[it.key as keyof typeof curves]} color={color} />
            </div>
            <div style={{ fontSize: 9, color: "var(--fg-3)", textTransform: "uppercase", letterSpacing: "0.04em" }}>{it.desc}</div>
          </div>
        );
      })}
    </div>
  );
}

export function PriceChartPanel({
  tickerSym, history, range, setRange,
}: {
  tickerSym: string | null;
  history: PricePoint[];
  range: "1M" | "3M" | "6M" | "1Y";
  setRange: (r: "1M" | "3M" | "6M" | "1Y") => void;
}) {
  if (!tickerSym || history.length < 2) {
    return (
      <div className="price-chart-wrap">
        <div className="chart-head">
          <span style={{ fontSize: 11, color: "var(--fg-3)", textTransform: "uppercase", letterSpacing: "0.06em" }}>
            Price · {range}
          </span>
        </div>
        <div className="empty">Select a ticker to load price history.</div>
      </div>
    );
  }

  const days = { "1M": 30, "3M": 90, "6M": 180, "1Y": 365 }[range];
  const slice = history.slice(-days);
  const first = slice[0]?.price, last = slice[slice.length - 1]?.price;
  const chg = last && first ? (last - first) / first : 0;
  const lo = Math.min(...slice.map(p => p.price));
  const hi = Math.max(...slice.map(p => p.price));

  return (
    <div className="price-chart-wrap">
      <div className="chart-head">
        <div style={{ display: "flex", gap: 12, alignItems: "baseline" }}>
          <span style={{ fontSize: 11, color: "var(--fg-3)", textTransform: "uppercase", letterSpacing: "0.06em" }}>
            {tickerSym} · {range}
          </span>
          <span className="mono" style={{ color: chg >= 0 ? "var(--accent)" : "var(--loss)" }}>
            {chg >= 0 ? "+" : ""}{(chg * 100).toFixed(2)}%
          </span>
          <span className="mono dim" style={{ fontSize: 10 }}>
            range ${lo.toFixed(2)}–${hi.toFixed(2)}
          </span>
        </div>
        <div className="tabs">
          {(["1M", "3M", "6M", "1Y"] as const).map(r => (
            <button key={r} aria-pressed={range === r} onClick={() => setRange(r)}>{r}</button>
          ))}
        </div>
      </div>
      <LineChart data={slice} height={140} yKey="price" />
    </div>
  );
}

export function PayoffPanel({
  bsArgs, direction, premium,
}: {
  bsArgs: { S: number; K: number; r: number; q: number; sigma: number; T: number; type: OptionType };
  direction: Direction;
  premium: number;
}) {
  const curve = useMemo(() => {
    const pts: { spot: number; pnl: number }[] = [];
    const lo = bsArgs.S * 0.6, hi = bsArgs.S * 1.4;
    const n = 80;
    for (let i = 0; i <= n; i++) {
      const spot = lo + (hi - lo) * i / n;
      const intrinsic = bsArgs.type === "call" ? Math.max(spot - bsArgs.K, 0) : Math.max(bsArgs.K - spot, 0);
      const pnl = direction === "long" ? intrinsic - premium : premium - intrinsic;
      pts.push({ spot, pnl });
    }
    return pts;
  }, [bsArgs, direction, premium]);

  const maxPnl = Math.max(...curve.map(p => p.pnl));
  const minPnl = Math.min(...curve.map(p => p.pnl));
  const breakeven = curve.find((p, i) => i > 0 && Math.sign(p.pnl) !== Math.sign(curve[i - 1].pnl));

  return (
    <div className="payoff">
      <div className="payoff-head">
        <span style={{ fontSize: 11, color: "var(--fg-3)", textTransform: "uppercase", letterSpacing: "0.06em" }}>
          Payoff at expiry
        </span>
        <div className="mono" style={{ fontSize: 10, color: "var(--fg-3)", display: "flex", gap: 12 }}>
          <span>max <span style={{ color: "var(--accent)" }}>${maxPnl > 1e6 ? "∞" : fmt(maxPnl, 2)}</span></span>
          <span>min <span style={{ color: "var(--loss)" }}>${fmt(minPnl, 2)}</span></span>
          {breakeven && <span>b/e <span style={{ color: "var(--fg-1)" }}>${breakeven.spot.toFixed(2)}</span></span>}
        </div>
      </div>
      <PayoffChart data={curve} height={160} strike={bsArgs.K} spot={bsArgs.S} />
    </div>
  );
}

