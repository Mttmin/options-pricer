import type { MonteCarloResult } from "../types/index.ts";

interface MonteCarloPathsChartProps {
  mcResult: MonteCarloResult | null;
  spotPrice: number;
  timeToMaturity: number;
}

export function MonteCarloPathsChart({
  mcResult,
  spotPrice,
  timeToMaturity,
}: MonteCarloPathsChartProps) {
  if (!mcResult) {
    return (
      <div className="rounded-lg border border-slate-800 bg-slate-900 p-6">
        <h3 className="text-lg font-semibold mb-2">Monte Carlo Simulation Paths</h3>
        <div className="h-48 flex items-center justify-center">
          <div className="text-center">
            <div className="text-slate-500 text-sm mb-2">
              No Monte Carlo simulation data available
            </div>
            <div className="text-slate-600 text-xs">
              Run a calculation to generate simulation paths
            </div>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="rounded-lg border border-slate-800 bg-slate-900 p-6">
      <div className="flex items-center justify-between mb-4">
        <div>
          <h3 className="text-lg font-semibold">Monte Carlo Simulation Paths</h3>
          <p className="text-sm text-slate-400">
            1% sample of paths (visualization placeholder)
          </p>
        </div>
        <div className="text-right text-xs text-slate-500">
          <div>Price: ${mcResult.price.toFixed(4)}</div>
          <div>95% CI: [{mcResult.ci_lower.toFixed(4)}, {mcResult.ci_upper.toFixed(4)}]</div>
        </div>
      </div>

      <div className="relative h-48 rounded border border-slate-800 bg-slate-950 flex items-center justify-center">
        <div className="text-center text-slate-500 text-sm">
          <p className="mb-2">Monte Carlo Path Visualization</p>
          <p className="text-xs text-slate-600">
            Backend needs to return path data for full visualization
          </p>
          <p className="text-xs text-slate-600 mt-1">
            Would show {Math.floor(100000 * 0.01).toLocaleString()} sample paths
          </p>
        </div>
      </div>

      <div className="mt-4 grid grid-cols-3 gap-4 text-sm">
        <div className="text-center p-3 rounded bg-slate-950">
          <div className="text-slate-400 text-xs mb-1">Starting Price</div>
          <div className="font-mono text-white">${spotPrice.toFixed(2)}</div>
        </div>
        <div className="text-center p-3 rounded bg-slate-950">
          <div className="text-slate-400 text-xs mb-1">Time to Maturity</div>
          <div className="font-mono text-white">{timeToMaturity.toFixed(3)}y</div>
        </div>
        <div className="text-center p-3 rounded bg-slate-950">
          <div className="text-slate-400 text-xs mb-1">Std Error</div>
          <div className="font-mono text-white">{mcResult.std_error.toFixed(6)}</div>
        </div>
      </div>
    </div>
  );
}
