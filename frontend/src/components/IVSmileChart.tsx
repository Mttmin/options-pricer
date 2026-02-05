import { useState, useMemo } from "react";
import type { IVSmileData } from "../types/index.ts";
import { Spinner } from "./ui/Spinner.tsx";
import { Select } from "./ui/Select.tsx";

interface IVSmileChartProps {
  data: IVSmileData | null;
  loading: boolean;
}

export function IVSmileChart({ data, loading }: IVSmileChartProps) {
  const expirations = useMemo(() => {
    if (!data) return [];
    return Object.keys(data.smiles_by_expiry).sort();
  }, [data]);

  const [selectedExpiration, setSelectedExpiration] = useState<string>("");

  const currentExpiration = selectedExpiration || expirations[0] || "";
  const smileData = data?.smiles_by_expiry[currentExpiration] || [];

  if (!data) {
    return (
      <div className="relative h-64 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <div className="absolute inset-0 flex items-center justify-center">
          <svg
            className="h-full w-full text-slate-800 opacity-30"
            viewBox="0 0 500 200"
            preserveAspectRatio="none"
          >
            <path
              d="M50,150 Q150,100 250,80 T450,120"
              fill="none"
              stroke="currentColor"
              strokeWidth="2"
            />
          </svg>
        </div>
        <p className="text-slate-500 text-sm z-10 bg-slate-900/80 px-3 py-1 rounded">
          Search for a ticker to view IV smile
        </p>
      </div>
    );
  }

  if (loading) {
    return (
      <div className="h-64 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <Spinner />
        <span className="ml-2 text-slate-400">Loading option chain...</span>
      </div>
    );
  }

  if (expirations.length === 0 || smileData.length === 0) {
    return (
      <div className="h-64 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <p className="text-slate-500 text-sm">
          No option chain data available for this symbol
        </p>
      </div>
    );
  }

  const calls = smileData.filter((p) => p.option_type === "Call");
  const puts = smileData.filter((p) => p.option_type === "Put");

  const allStrikes = smileData.map((p) => p.strike);
  const allIVs = smileData.map((p) => p.implied_volatility * 100);

  const minStrike = Math.min(...allStrikes);
  const maxStrike = Math.max(...allStrikes);
  const minIV = Math.min(...allIVs);
  const maxIV = Math.max(...allIVs);

  const strikeRange = maxStrike - minStrike || 1;
  const ivRange = maxIV - minIV || 1;
  const strikePadding = strikeRange * 0.1;
  const ivPadding = ivRange * 0.1;

  const chartWidth = 500;
  const chartHeight = 200;

  const getX = (strike: number) =>
    ((strike - minStrike + strikePadding) /
      (strikeRange + 2 * strikePadding)) *
    chartWidth;

  const getY = (iv: number) =>
    chartHeight -
    ((iv - minIV + ivPadding) / (ivRange + 2 * ivPadding)) * chartHeight;

  const atmX = getX(data.underlying_price);

  const callsSorted = [...calls].sort((a, b) => a.strike - b.strike);
  const putsSorted = [...puts].sort((a, b) => a.strike - b.strike);

  const callPoints = callsSorted
    .map((p) => `${getX(p.strike)},${getY(p.implied_volatility * 100)}`)
    .join(" ");

  const putPoints = putsSorted
    .map((p) => `${getX(p.strike)},${getY(p.implied_volatility * 100)}`)
    .join(" ");

  return (
    <div className="rounded-lg border border-slate-800 bg-slate-900 p-4">
      <div className="flex items-center justify-between mb-4">
        <div>
          <span className="text-lg font-semibold text-white">
            Implied Volatility Smile
          </span>
          <span className="ml-3 text-sm text-slate-400">
            {data.symbol} @ ${data.underlying_price.toFixed(2)}
          </span>
        </div>
        <div className="w-48">
          <Select
            label="Expiration"
            value={currentExpiration}
            onChange={setSelectedExpiration}
            options={expirations.map((exp) => ({ value: exp, label: exp }))}
          />
        </div>
      </div>

      <div className="relative mb-4">
        <svg
          viewBox={`0 0 ${chartWidth} ${chartHeight}`}
          className="w-full h-64"
          preserveAspectRatio="xMidYMid meet"
        >
          <line
            x1={atmX}
            y1="0"
            x2={atmX}
            y2={chartHeight}
            stroke="#8b5cf6"
            strokeWidth="1"
            strokeDasharray="4 4"
            vectorEffect="non-scaling-stroke"
          />

          {callsSorted.length > 1 && (
            <polyline
              points={callPoints}
              fill="none"
              stroke="#3b82f6"
              strokeWidth="2"
              vectorEffect="non-scaling-stroke"
            />
          )}

          {putsSorted.length > 1 && (
            <polyline
              points={putPoints}
              fill="none"
              stroke="#ef4444"
              strokeWidth="2"
              vectorEffect="non-scaling-stroke"
            />
          )}

          {callsSorted.map((point, idx) => (
            <circle
              key={`call-${idx}`}
              cx={getX(point.strike)}
              cy={getY(point.implied_volatility * 100)}
              r="3"
              fill="#3b82f6"
            />
          ))}

          {putsSorted.map((point, idx) => (
            <circle
              key={`put-${idx}`}
              cx={getX(point.strike)}
              cy={getY(point.implied_volatility * 100)}
              r="3"
              fill="#ef4444"
            />
          ))}
        </svg>

        <div className="absolute left-0 top-0 text-xs text-slate-500">
          {maxIV.toFixed(1)}%
        </div>
        <div className="absolute left-0 bottom-0 text-xs text-slate-500">
          {minIV.toFixed(1)}%
        </div>
        <div className="absolute right-0 bottom-0 text-xs text-slate-500">
          ${maxStrike.toFixed(2)}
        </div>
        <div className="absolute left-1/2 -translate-x-1/2 -bottom-6 text-xs text-slate-400">
          Strike Price
        </div>
        <div
          className="absolute -left-12 top-1/2 -translate-y-1/2 -rotate-90 text-xs text-slate-400"
          style={{ transformOrigin: "center" }}
        >
          Implied Volatility (%)
        </div>
      </div>

      <div className="flex items-center justify-center gap-6 text-sm mt-8">
        <div className="flex items-center gap-2">
          <div className="w-3 h-3 rounded-full bg-blue-500"></div>
          <span className="text-slate-400">Calls</span>
        </div>
        <div className="flex items-center gap-2">
          <div className="w-3 h-3 rounded-full bg-red-500"></div>
          <span className="text-slate-400">Puts</span>
        </div>
        <div className="flex items-center gap-2">
          <div className="w-8 h-0.5 bg-purple-500 opacity-50"></div>
          <span className="text-slate-400">ATM ({data.underlying_price.toFixed(2)})</span>
        </div>
      </div>
    </div>
  );
}
