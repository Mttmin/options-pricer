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

  const chartWidth = 500;
  const chartHeight = 200;
  const padding = { left: 50, right: 30, top: 20, bottom: 40 };

  const getX = (strike: number) =>
    padding.left +
    ((strike - minStrike) / strikeRange) * (chartWidth - padding.left - padding.right);

  const getY = (iv: number) =>
    padding.top +
    (1 - (iv - minIV) / ivRange) * (chartHeight - padding.top - padding.bottom);

  const atmX = getX(data.underlying_price);

  const numGridLinesY = 5;
  const numGridLinesX = 5;
  const yTicks = Array.from({ length: numGridLinesY }, (_, i) =>
    minIV + (i * ivRange) / (numGridLinesY - 1)
  );
  const xTicks = Array.from({ length: numGridLinesX }, (_, i) =>
    minStrike + (i * strikeRange) / (numGridLinesX - 1)
  );

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
          {yTicks.map((tick, i) => {
            const y = getY(tick);
            return (
              <g key={`y-grid-${i}`}>
                <line
                  x1={padding.left}
                  y1={y}
                  x2={chartWidth - padding.right}
                  y2={y}
                  stroke="#334155"
                  strokeWidth="0.5"
                  vectorEffect="non-scaling-stroke"
                />
                <text
                  x={padding.left - 5}
                  y={y}
                  textAnchor="end"
                  dominantBaseline="middle"
                  className="text-[8px] fill-slate-400"
                >
                  {tick.toFixed(1)}%
                </text>
              </g>
            );
          })}

          {xTicks.map((tick, i) => {
            const x = getX(tick);
            return (
              <g key={`x-grid-${i}`}>
                <line
                  x1={x}
                  y1={padding.top}
                  x2={x}
                  y2={chartHeight - padding.bottom}
                  stroke="#334155"
                  strokeWidth="0.5"
                  vectorEffect="non-scaling-stroke"
                />
                <text
                  x={x}
                  y={chartHeight - padding.bottom + 15}
                  textAnchor="middle"
                  className="text-[8px] fill-slate-400"
                >
                  ${tick.toFixed(0)}
                </text>
              </g>
            );
          })}

          <line
            x1={padding.left}
            y1={padding.top}
            x2={padding.left}
            y2={chartHeight - padding.bottom}
            stroke="#475569"
            strokeWidth="1.5"
            vectorEffect="non-scaling-stroke"
          />
          <line
            x1={padding.left}
            y1={chartHeight - padding.bottom}
            x2={chartWidth - padding.right}
            y2={chartHeight - padding.bottom}
            stroke="#475569"
            strokeWidth="1.5"
            vectorEffect="non-scaling-stroke"
          />

          <line
            x1={atmX}
            y1={padding.top}
            x2={atmX}
            y2={chartHeight - padding.bottom}
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

          <text
            x={chartWidth / 2}
            y={chartHeight - 5}
            textAnchor="middle"
            className="text-[10px] fill-slate-400 font-medium"
          >
            Strike Price
          </text>
          <text
            x={15}
            y={chartHeight / 2}
            textAnchor="middle"
            transform={`rotate(-90 15 ${chartHeight / 2})`}
            className="text-[10px] fill-slate-400 font-medium"
          >
            Implied Volatility (%)
          </text>
        </svg>
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
