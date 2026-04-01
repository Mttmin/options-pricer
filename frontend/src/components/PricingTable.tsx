import { useState } from "react";
import type { PriceResponse, ExerciseStyle, CtmcPriceResponse } from "../types/index.ts";
import { Spinner } from "./ui/Spinner.tsx";

interface PricingTableProps {
  result: PriceResponse | null;
  loading: boolean;
  exerciseStyle: ExerciseStyle;
  ctmcResult?: CtmcPriceResponse | null;
  ctmcLoading?: boolean;
  ctmcError?: string | null;
}

function formatPrice(val: number | undefined | null): string {
  if (val == null) return "-";
  return val.toFixed(4);
}

function formatGreek(val: number | undefined | null): string {
  if (val == null) return "-";
  return val.toFixed(6);
}

function formatMs(val: number | undefined | null): string {
  if (val == null) return "-";
  return `${val.toFixed(1)} ms`;
}

type GreekKey = "delta" | "gamma" | "theta" | "vega" | "rho";

function greekStrengthClass(key: GreekKey, val: number | undefined | null): string {
  if (val == null) return "text-slate-400";
  const abs = Math.abs(val);

  const high = {
    delta: 0.6,
    gamma: 0.05,
    theta: 50,
    vega: 10,
    rho: 5,
  } satisfies Record<GreekKey, number>;

  const mid = {
    delta: 0.25,
    gamma: 0.015,
    theta: 20,
    vega: 4,
    rho: 2,
  } satisfies Record<GreekKey, number>;

  if (abs >= high[key]) return "text-emerald-400";
  if (abs >= mid[key]) return "text-amber-400";
  return "text-slate-300";
}

export function PricingTable({
  result,
  loading,
  exerciseStyle,
  ctmcResult,
  ctmcLoading,
  ctmcError,
}: PricingTableProps) {
  const [hestonOpen, setHestonOpen] = useState(false);

  if (loading) {
    return (
      <div className="flex items-center justify-center rounded-lg border border-slate-800 bg-slate-900 p-8">
        <Spinner />
        <span className="ml-3 text-slate-400">Calculating prices...</span>
      </div>
    );
  }

  if (!result) return null;

  const { pricing, greeks } = result;
  const timings = pricing.timings;
  const showCtmc = ctmcLoading || ctmcResult != null || ctmcError != null;

  return (
    <div className="overflow-hidden rounded-lg border border-slate-800">
      <table className="w-full text-sm">
        <thead className="bg-slate-900">
          <tr>
            <th className="px-4 py-3 text-left font-medium text-slate-400">
              Method
            </th>
            <th className="px-4 py-3 text-right font-medium text-slate-400">
              Price
            </th>
            <th className="px-4 py-3 text-right font-medium text-slate-400">
              Details
            </th>
          </tr>
        </thead>
        <tbody className="divide-y divide-slate-800">
          {exerciseStyle === "european" && (
            <>
              <tr className="bg-slate-950 hover:bg-slate-900">
                <td className="px-4 py-3 font-medium text-white">Black-Scholes</td>
                <td className="px-4 py-3 text-right font-mono text-green-400">
                  {formatPrice(pricing.black_scholes)}
                </td>
                <td className="px-4 py-3 text-right text-slate-500">
                  Analytical
                  {timings && (
                    <div className="text-xs text-slate-500">
                      Time: {formatMs(timings.black_scholes_ms)}
                    </div>
                  )}
                </td>
              </tr>

              {pricing.monte_carlo && (
                <tr className="bg-slate-950 hover:bg-slate-900">
                  <td className="px-4 py-3 font-medium text-white">Monte Carlo</td>
                  <td className="px-4 py-3 text-right font-mono text-green-400">
                    {formatPrice(pricing.monte_carlo.price)}
                  </td>
                  <td className="px-4 py-3 text-right text-xs text-slate-500">
                    <span className="text-slate-400">
                      95% CI: [{formatPrice(pricing.monte_carlo.ci_lower)},{" "}
                      {formatPrice(pricing.monte_carlo.ci_upper)}]
                    </span>
                    <br />
                    <span>SE: {pricing.monte_carlo.std_error.toFixed(6)}</span>
                    {timings && (
                      <>
                        <br />
                        <span>Time: {formatMs(timings.monte_carlo_ms)}</span>
                      </>
                    )}
                  </td>
                </tr>
              )}

              {pricing.binomial_european != null && (
                <tr className="bg-slate-950 hover:bg-slate-900">
                  <td className="px-4 py-3 font-medium text-white">
                    Binomial (European)
                  </td>
                  <td className="px-4 py-3 text-right font-mono text-green-400">
                    {formatPrice(pricing.binomial_european)}
                  </td>
                  <td className="px-4 py-3 text-right text-slate-500">
                    CRR Tree
                    {timings && (
                      <div className="text-xs text-slate-500">
                        Time: {formatMs(timings.binomial_european_ms)}
                      </div>
                    )}
                  </td>
                </tr>
              )}
            </>
          )}

          {exerciseStyle === "american" && (
            <>
              {pricing.binomial_american != null && (
                <tr className="bg-slate-950 hover:bg-slate-900">
                  <td className="px-4 py-3 font-medium text-white">
                    Binomial (American)
                  </td>
                  <td className="px-4 py-3 text-right font-mono text-green-400">
                    {formatPrice(pricing.binomial_american)}
                  </td>
                  <td className="px-4 py-3 text-right text-slate-500">
                    Early Exercise CRR
                    {timings && (
                      <div className="text-xs text-slate-500">
                        Time: {formatMs(timings.binomial_american_ms)}
                      </div>
                    )}
                  </td>
                </tr>
              )}

              {pricing.bs_american_approx != null && (
                <tr className="bg-slate-950 hover:bg-slate-900">
                  <td className="px-4 py-3 font-medium text-white">
                    BS American Approx
                  </td>
                  <td className="px-4 py-3 text-right font-mono text-green-400">
                    {formatPrice(pricing.bs_american_approx)}
                  </td>
                  <td className="px-4 py-3 text-right text-slate-500">
                    Black's Formula
                    {timings && (
                      <div className="text-xs text-slate-500">
                        Time: {formatMs(timings.bs_american_approx_ms)}
                      </div>
                    )}
                  </td>
                </tr>
              )}

              {pricing.penalty_solver && (
                <tr className="bg-slate-950 hover:bg-slate-900">
                  <td className="px-4 py-3 font-medium text-white">
                    PDE Penalty Method
                  </td>
                  <td className="px-4 py-3 text-right font-mono text-green-400">
                    {formatPrice(pricing.penalty_solver.price)}
                  </td>
                  <td className="px-4 py-3 text-right text-xs text-slate-500">
                    <span className="text-slate-400">
                      Nodes: {pricing.penalty_solver.diagnostics.spatial_nodes},{" "}
                      Timesteps: {pricing.penalty_solver.diagnostics.timesteps}
                    </span>
                    <br />
                    <span>
                      Avg iters/step:{" "}
                      {pricing.penalty_solver.diagnostics.avg_iterations_per_step.toFixed(1)}
                    </span>
                    {timings && (
                      <>
                        <br />
                        <span>Time: {formatMs(timings.penalty_solver_ms)}</span>
                      </>
                    )}
                  </td>
                </tr>
              )}

              {pricing.binomial_european != null && (
                <tr className="bg-slate-950 hover:bg-slate-900 opacity-60">
                  <td className="px-4 py-3 font-medium text-white">
                    Binomial (European)
                  </td>
                  <td className="px-4 py-3 text-right font-mono text-slate-400">
                    {formatPrice(pricing.binomial_european)}
                  </td>
                  <td className="px-4 py-3 text-right text-slate-500 text-xs">
                    For comparison
                    {timings && (
                      <div className="text-xs text-slate-500">
                        Time: {formatMs(timings.binomial_european_ms)}
                      </div>
                    )}
                  </td>
                </tr>
              )}
            </>
          )}
          {/* CTMC row — always shown when available, regardless of exercise style */}
          {showCtmc && (
            <>
              <tr className="bg-slate-950 hover:bg-slate-900">
                <td className="px-4 py-3 font-medium text-white">
                  CTMC (Deep Cal)
                </td>
                <td className="px-4 py-3 text-right font-mono text-green-400">
                  {ctmcLoading ? (
                    <span className="inline-flex items-center justify-end gap-2 text-slate-400">
                      <Spinner size="sm" />
                      <span className="text-xs">calibrating…</span>
                    </span>
                  ) : ctmcError ? (
                    <span className="text-red-400 text-xs">{ctmcError}</span>
                  ) : (
                    formatPrice(ctmcResult?.price)
                  )}
                </td>
                <td className="px-4 py-3 text-right text-xs text-slate-500">
                  {ctmcResult && !ctmcLoading && (
                    <>
                      <span>Heston SV · American</span>
                      <br />
                      <span>
                        Time:{" "}
                        {ctmcResult.timing_ms >= 1000
                          ? `${(ctmcResult.timing_ms / 1000).toFixed(1)}s`
                          : `${ctmcResult.timing_ms.toFixed(0)}ms`}
                      </span>
                      <br />
                      <button
                        onClick={() => setHestonOpen((o) => !o)}
                        className="mt-1 text-slate-400 hover:text-slate-200 underline underline-offset-2"
                      >
                        {hestonOpen ? "Hide params ▲" : "Show params ▼"}
                      </button>
                    </>
                  )}
                </td>
              </tr>

              {hestonOpen && ctmcResult && (
                <tr className="bg-slate-900">
                  <td colSpan={3} className="px-4 py-3">
                    <div className="rounded border border-slate-700 bg-slate-950 px-4 py-3">
                      <p className="text-xs font-medium uppercase tracking-wider text-slate-500 mb-3">
                        Calibrated Heston Parameters
                      </p>
                      <div className="grid grid-cols-5 gap-2 text-center mb-3">
                        {(
                          [
                            ["κ (kappa)", ctmcResult.heston.kappa],
                            ["θ (theta)", ctmcResult.heston.theta],
                            ["σ (sigma)", ctmcResult.heston.sigma],
                            ["ρ (rho)", ctmcResult.heston.rho],
                            ["v₀", ctmcResult.heston.v0],
                          ] as [string, number][]
                        ).map(([label, val]) => (
                          <div
                            key={label}
                            className="rounded border border-slate-700 bg-slate-900 px-2 py-2"
                          >
                            <div className="text-xs text-slate-500 mb-1">{label}</div>
                            <div className="font-mono text-xs text-slate-300">
                              {val.toFixed(4)}
                            </div>
                          </div>
                        ))}
                      </div>
                      <div className="flex flex-wrap gap-4 text-xs text-slate-500">
                        <span>
                          IV RMSE:{" "}
                          <span className="text-slate-400">
                            {(ctmcResult.ivrmse * 100).toFixed(2)} vol pts
                          </span>
                        </span>
                        <span>
                          Instruments:{" "}
                          <span className="text-slate-400">{ctmcResult.n_instruments}</span>
                        </span>
                        <span>
                          Surrogate evals:{" "}
                          <span className="text-slate-400">{ctmcResult.n_evals}</span>
                        </span>
                        <span>
                          Spot:{" "}
                          <span className="text-slate-400">{ctmcResult.spot_price.toFixed(2)}</span>
                        </span>
                      </div>
                    </div>
                  </td>
                </tr>
              )}
            </>
          )}
        </tbody>
      </table>

      {timings && (
        <div className="border-t border-slate-800 bg-slate-900 px-4 py-2 text-xs text-slate-500">
          Total compute time: {formatMs(timings.total_ms)}
        </div>
      )}

      {greeks && (
        <div className="border-t border-slate-800 bg-slate-900 px-4 py-4">
          <div className="mb-3 flex items-center justify-between">
            <p className="text-xs font-medium uppercase tracking-wider text-slate-400">
              Greeks
            </p>
            <span className="text-[11px] text-slate-500">
              Color shows magnitude
            </span>
          </div>
          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-5">
            <div className="group rounded-md border border-slate-800 bg-gradient-to-br from-slate-950 via-slate-900 to-slate-950 px-3 py-2">
              <div className="flex items-center justify-between">
                <span className="text-xs uppercase tracking-wide text-slate-500">Delta</span>
                <span className="text-sm font-semibold text-slate-400">Delta</span>
              </div>
              <div className="mt-2 flex items-baseline justify-between">
                <span className="text-2xl text-slate-300">Δ</span>
                <span className={`font-mono text-base ${greekStrengthClass("delta", greeks.delta)}`}>
                  {formatGreek(greeks.delta)}
                </span>
              </div>
            </div>

            <div className="group rounded-md border border-slate-800 bg-gradient-to-br from-slate-950 via-slate-900 to-slate-950 px-3 py-2">
              <div className="flex items-center justify-between">
                <span className="text-xs uppercase tracking-wide text-slate-500">Gamma</span>
                <span className="text-sm font-semibold text-slate-400">Gamma</span>
              </div>
              <div className="mt-2 flex items-baseline justify-between">
                <span className="text-2xl text-slate-300">Γ</span>
                <span className={`font-mono text-base ${greekStrengthClass("gamma", greeks.gamma)}`}>
                  {formatGreek(greeks.gamma)}
                </span>
              </div>
            </div>

            <div className="group rounded-md border border-slate-800 bg-gradient-to-br from-slate-950 via-slate-900 to-slate-950 px-3 py-2">
              <div className="flex items-center justify-between">
                <span className="text-xs uppercase tracking-wide text-slate-500">Theta</span>
                <span className="text-sm font-semibold text-slate-400">Theta</span>
              </div>
              <div className="mt-2 flex items-baseline justify-between">
                <span className="text-2xl text-slate-300">Θ</span>
                <span className={`font-mono text-base ${greekStrengthClass("theta", greeks.theta)}`}>
                  {formatGreek(greeks.theta)}
                </span>
              </div>
            </div>

            <div className="group rounded-md border border-slate-800 bg-gradient-to-br from-slate-950 via-slate-900 to-slate-950 px-3 py-2">
              <div className="flex items-center justify-between">
                <span className="text-xs uppercase tracking-wide text-slate-500">Vega</span>
                <span className="text-sm font-semibold text-slate-400">Vega</span>
              </div>
              <div className="mt-2 flex items-baseline justify-between">
                <span className="text-2xl text-slate-300">ν</span>
                <span className={`font-mono text-base ${greekStrengthClass("vega", greeks.vega)}`}>
                  {formatGreek(greeks.vega)}
                </span>
              </div>
            </div>

            <div className="group rounded-md border border-slate-800 bg-gradient-to-br from-slate-950 via-slate-900 to-slate-950 px-3 py-2">
              <div className="flex items-center justify-between">
                <span className="text-xs uppercase tracking-wide text-slate-500">Rho</span>
                <span className="text-sm font-semibold text-slate-400">Rho</span>
              </div>
              <div className="mt-2 flex items-baseline justify-between">
                <span className="text-2xl text-slate-300">ρ</span>
                <span className={`font-mono text-base ${greekStrengthClass("rho", greeks.rho)}`}>
                  {formatGreek(greeks.rho)}
                </span>
              </div>
            </div>
          </div>
        </div>
      )}

    </div>
  );
}
