import type { PriceResponse, ExerciseStyle } from "../types/index.ts";
import { Spinner } from "./ui/Spinner.tsx";

interface PricingTableProps {
  result: PriceResponse | null;
  loading: boolean;
  exerciseStyle: ExerciseStyle;
}

function formatPrice(val: number | undefined | null): string {
  if (val == null) return "-";
  return val.toFixed(4);
}

function formatGreek(val: number | undefined | null): string {
  if (val == null) return "-";
  return val.toFixed(6);
}

export function PricingTable({ result, loading, exerciseStyle }: PricingTableProps) {
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
                  </td>
                </tr>
              )}
            </>
          )}
        </tbody>
      </table>

      {greeks && (
        <div className="border-t border-slate-800 bg-slate-900 px-4 py-3">
          <p className="mb-2 text-xs font-medium uppercase tracking-wider text-slate-400">
            Greeks
          </p>
          <div className="flex flex-wrap gap-4 text-sm">
            <span>
              <span className="text-slate-500">Delta:</span>{" "}
              <span className="font-mono text-slate-300">
                {formatGreek(greeks.delta)}
              </span>
            </span>
            <span>
              <span className="text-slate-500">Gamma:</span>{" "}
              <span className="font-mono text-slate-300">
                {formatGreek(greeks.gamma)}
              </span>
            </span>
            <span>
              <span className="text-slate-500">Theta:</span>{" "}
              <span className="font-mono text-slate-300">
                {formatGreek(greeks.theta)}
              </span>
            </span>
            <span>
              <span className="text-slate-500">Vega:</span>{" "}
              <span className="font-mono text-slate-300">
                {formatGreek(greeks.vega)}
              </span>
            </span>
            <span>
              <span className="text-slate-500">Rho:</span>{" "}
              <span className="font-mono text-slate-300">
                {formatGreek(greeks.rho)}
              </span>
            </span>
          </div>
        </div>
      )}
    </div>
  );
}
