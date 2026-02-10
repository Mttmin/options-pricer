import type { MarketDataResponse } from "../types/index.ts";

interface MarketDataDisplayProps {
  data: MarketDataResponse | null;
}

function formatVol(value: number | null): string {
  if (value === null || value === undefined) return "N/A";
  return `${(value * 100).toFixed(2)}%`;
}

export function MarketDataDisplay({ data }: MarketDataDisplayProps) {
  if (!data) return null;

  return (
    <div className="rounded border border-slate-700 bg-slate-800/50 p-4">
      <div className="flex flex-wrap items-center gap-6">
        <div>
          <span className="text-xs font-medium uppercase tracking-wider text-slate-400">
            Symbol
          </span>
          <p className="text-lg font-semibold text-white">{data.symbol}</p>
        </div>
        <div>
          <span className="text-xs font-medium uppercase tracking-wider text-slate-400">
            Spot Price
          </span>
          <p className="text-lg font-semibold text-white">
            ${data.spot_price.toFixed(2)}
          </p>
        </div>
        <div>
          <span className="text-xs font-medium uppercase tracking-wider text-slate-400">
            Historical Vol
          </span>
          <p className="text-lg font-semibold text-white">
            {formatVol(data.historical_volatility)}
          </p>
        </div>
        <div>
          <span className="text-xs font-medium uppercase tracking-wider text-slate-400">
            Implied Vol
          </span>
          <p className="text-lg font-semibold text-white">
            {formatVol(data.implied_volatility)}
          </p>
        </div>
      </div>
      {data.corporate_action_detected && (
        <div className="mt-3 rounded border border-yellow-600/50 bg-yellow-900/20 px-3 py-2 text-sm text-yellow-400">
          Corporate action detected. Volatility figures may be unreliable.
        </div>
      )}
    </div>
  );
}
