import { useMemo } from "react";
import type { PricePoint } from "../api/client.ts";
import type {
  CtmcPriceResponse,
  ModelQuote,
  ExerciseStyle,
  GreeksResult,
  OptionType,
  PayoffCurve,
  PriceResponse,
  PricingResult,
  PricingTimings,
} from "../types/index.ts";
import { LineChart, PayoffChart, Sparkline } from "./Charts.tsx";
import { greekCurve } from "../utils/blackScholes.ts";

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

function fitDetail(quote: ModelQuote, nInstruments: number): string {
  const base = `IV RMSE ${(quote.ivrmse * 100).toFixed(2)}vp · ${nInstruments} instr.`;
  return quote.four_two
    ? `${base} · a ${quote.four_two.a.toFixed(2)} b ${quote.four_two.b.toFixed(3)}`
    : base;
}

function slvRows(ctmc: CtmcPriceResponse | null, loading: boolean): PriceRow[] {
  if (!ctmc && !loading) return [];
  const n = ctmc?.n_instruments ?? 0;
  const specs: {
    id: string;
    method: string;
    quote?: ModelQuote | null;
    error?: string | null;
  }[] = [
    { id: "ctmc_heston", method: "CTMC (Heston)", quote: ctmc?.heston_quote, error: ctmc?.heston_error },
    { id: "ctmc_four_two", method: "CTMC (4/2)", quote: ctmc?.four_two_quote, error: ctmc?.four_two_error },
  ];
  if (ctmc && !ctmc.heston_quote && !ctmc.four_two_quote && !ctmc.heston_error && !ctmc.four_two_error) {
    const legacy: ModelQuote = {
      price: ctmc.price,
      ivrmse: ctmc.ivrmse,
      n_evals: ctmc.n_evals,
      timing_ms: ctmc.timing_ms,
      heston: ctmc.heston,
      four_two: ctmc.four_two,
    };
    if (ctmc.model === "four_two") specs[1].quote = legacy;
    else specs[0].quote = legacy;
  }
  return specs.map(s => ({
    id: s.id,
    method: s.method,
    sub: "Stochastic vol · calibrated",
    price: s.quote?.price ?? null,
    ms: s.quote?.timing_ms ?? null,
    loading: loading && !s.quote && !s.error,
    detail: s.quote ? fitDetail(s.quote, n) : s.error ?? undefined,
  }));
}

export function buildRows(
  result: PriceResponse | null,
  exerciseStyle: ExerciseStyle,
  ctmc: CtmcPriceResponse | null,
  ctmcLoading: boolean
): PriceRow[] {
  const rows: PriceRow[] = [];
  if (!result) return slvRows(ctmc, ctmcLoading);

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

  rows.push(...slvRows(ctmc, ctmcLoading));

  return rows;
}

export function CalibrationStrip({
  title, quote, error, nInstruments,
}: {
  title: string;
  quote?: ModelQuote | null;
  error?: string | null;
  nInstruments: number;
}) {
  const four = quote?.four_two;
  const p = four ?? quote?.heston;
  const cells: { k: string; v: string }[] = p ? [
    { k: "κ", v: p.kappa.toFixed(2) },
    { k: "θ", v: p.theta.toFixed(3) },
    { k: "σ", v: p.sigma.toFixed(3) },
    { k: "ρ", v: p.rho.toFixed(3) },
    { k: "v₀", v: p.v0.toFixed(3) },
  ] : [];
  if (four) cells.push({ k: "a", v: four.a.toFixed(3) }, { k: "b", v: four.b.toFixed(4) });
  return (
    <div className="cal-strip">
      <div className="cal-head">
        <span className="section-label">{title}</span>
        {quote && (
          <span className="mono dim">
            IV RMSE {(quote.ivrmse * 100).toFixed(2)}vp · {nInstruments} instr · {quote.n_evals} evals
          </span>
        )}
      </div>
      {error && !quote ? (
        <div className="empty error">{error}</div>
      ) : (
        <div className={`cal-grid n${cells.length}`}>
          {cells.map(c => (
            <div className="cal-cell" key={c.k}>
              <div className="cal-k">{c.k}</div>
              <div className="cal-v tnum">{c.v}</div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

export function HeroPrice({
  label, primary, secondary, pinned, stale, onPin, onClearPin,
}: {
  label: string;
  primary: { method: string; price: number | null; ms: number | null };
  secondary: { method: string; price: number | null }[];
  pinned: { price: number; method: string } | null;
  stale?: boolean;
  onPin: () => void;
  onClearPin: () => void;
}) {
  return (
    <div className={`hero-price ${stale ? "stale" : ""}`}>
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
        {pinned && (
          <div className="alt pin">
            <div className="mk">Pinned</div>
            <div>${fmt(pinned.price, 4)}</div>
          </div>
        )}
        {pinned ? (
          <button type="button" className="icon-btn" onClick={onClearPin} title="Clear snapshot">Clear pin</button>
        ) : (
          <button type="button" className="icon-btn" onClick={onPin} title="Freeze current price" disabled={primary.price == null}>Pin</button>
        )}
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
  const reference = rows.find(r => r.id === selectedId)?.price ?? null;

  return (
    <table className="pricing-table">
      <thead>
        <tr>
          <th>Method</th>
          <th className="num">Price</th>
          <th className="num">Δ vs primary</th>
          <th className="num">Detail</th>
          <th className="num">Time</th>
        </tr>
      </thead>
      <tbody>
        {rows.map(r => {
          const selected = r.id === selectedId;
          const diff = !selected && reference != null && r.price != null ? r.price - reference : null;
          return (
            <tr
              key={r.id}
              className={[
                emphasis === "hero" && selected ? "primary" : "",
                r.dim ? "dimmed" : "",
                r.loading ? "loading" : "",
              ].filter(Boolean).join(" ")}
              tabIndex={r.loading ? -1 : 0}
              onClick={() => { if (!r.loading) onSelect(r.id); }}
              onKeyDown={(e) => {
                if (r.loading) return;
                if (e.key === "Enter" || e.key === " ") { e.preventDefault(); onSelect(r.id); }
              }}
              title={r.loading ? "" : "Set as primary"}
            >
              <td className="method">
                {r.method}
                <small>{r.sub}</small>
              </td>
              <td className="price tnum">
                {r.loading ? <span className="dim">calibrating…</span> : (r.price != null ? `$${fmt(r.price, 4)}` : "—")}
              </td>
              <td className="delta tnum">
                {diff != null ? (
                  <span className={Math.abs(diff) < 5e-5 ? "dim" : diff > 0 ? "pos" : "neg"}>
                    {diff > 0 ? "+" : ""}{diff.toFixed(4)}
                  </span>
                ) : selected ? <span className="dim">primary</span> : ""}
              </td>
              <td className="detail">{r.detail || ""}</td>
              <td className="detail">{fmtMs(r.ms)}</td>
            </tr>
          );
        })}
      </tbody>
    </table>
  );
}

export function GreeksStrip({
  greeks, bsArgs, showCurves = true,
}: {
  greeks: GreeksResult;
  bsArgs: { S: number; K: number; r: number; q: number; sigma: number; T: number; type: OptionType };
  // Sparklines come from a vanilla BS curve; meaningless for structures without a single strike.
  showCurves?: boolean;
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
            {showCurves && (
              <div className="greek-spark">
                <Sparkline data={curves[it.key as keyof typeof curves]} color={color} />
              </div>
            )}
            <div className="greek-desc">{it.desc}</div>
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
          <span className="section-label">Price · {range}</span>
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
        <div className="chart-head-meta">
          <span className="section-label">{tickerSym} · {range}</span>
          <span className="mono" style={{ color: chg >= 0 ? "var(--accent)" : "var(--loss)" }}>
            {chg >= 0 ? "+" : ""}{(chg * 100).toFixed(2)}%
          </span>
          <span className="mono dim range">range ${lo.toFixed(2)}–${hi.toFixed(2)}</span>
        </div>
        <div className="tabs">
          {(["1M", "3M", "6M", "1Y"] as const).map(r => (
            <button key={r} aria-pressed={range === r} onClick={() => setRange(r)}>{r}</button>
          ))}
        </div>
      </div>
      <LineChart data={slice} height={140} yKey="price" showAxis />
    </div>
  );
}

export function PayoffPanel({
  curve, strike, spot,
}: {
  curve: PayoffCurve;
  strike?: number;
  spot?: number;
}) {
  const points = useMemo(
    () => curve.spot_prices.map((s, i) => ({ spot: s, pnl: curve.payoffs[i] ?? 0 })),
    [curve]
  );

  if (points.length < 2) return null;

  const maxPnl = Math.max(...points.map(p => p.pnl));
  const minPnl = Math.min(...points.map(p => p.pnl));
  const breakeven = points.find((p, i) => i > 0 && Math.sign(p.pnl) !== Math.sign(points[i - 1].pnl));

  return (
    <div className="payoff">
      <div className="payoff-head">
        <span className="section-label">Payoff at expiry</span>
        <div className="mono payoff-stats" title="Evaluated over the plotted spot range (±50% of spot)">
          <span>max <span style={{ color: "var(--accent)" }}>${fmt(maxPnl, 2)}</span></span>
          <span>min <span style={{ color: "var(--loss)" }}>${fmt(minPnl, 2)}</span></span>
          {breakeven && <span>b/e <span style={{ color: "var(--fg-1)" }}>${breakeven.spot.toFixed(2)}</span></span>}
        </div>
      </div>
      <PayoffChart data={points} height={180} strike={strike} spot={spot} />
    </div>
  );
}

