import type { PricePoint } from "../api/client.ts";
import { Spinner } from "./ui/Spinner.tsx";

interface PriceChartProps {
  data: PricePoint[];
  ticker: string | null;
  currentPrice?: number;
  loading: boolean;
}

export function PriceChart({ data, ticker, currentPrice, loading }: PriceChartProps) {
  if (!ticker) {
    return (
      <div className="relative h-48 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <div className="absolute inset-0 flex items-center justify-center">
          <svg className="h-full w-full text-slate-800 opacity-30" viewBox="0 0 400 150" preserveAspectRatio="none">
            <path
              d="M0,100 Q50,80 100,90 T200,70 T300,85 T400,60"
              fill="none"
              stroke="currentColor"
              strokeWidth="2"
            />
          </svg>
        </div>
        <p className="text-slate-500 text-sm z-10 bg-slate-900/80 px-3 py-1 rounded">
          Search for a ticker to view price chart
        </p>
      </div>
    );
  }

  if (loading) {
    return (
      <div className="h-48 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <Spinner />
        <span className="ml-2 text-slate-400">Loading chart...</span>
      </div>
    );
  }

  if (data.length === 0) {
    return (
      <div className="h-48 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <p className="text-slate-500 text-sm">No price data available</p>
      </div>
    );
  }

  const prices = data.map((p) => p.price);
  const minPrice = Math.min(...prices);
  const maxPrice = Math.max(...prices);
  const priceRange = maxPrice - minPrice || 1;
  const padding = priceRange * 0.1;

  const chartHeight = 150;
  const chartWidth = 400;

  const points = data
    .map((point, i) => {
      const x = (i / (data.length - 1)) * chartWidth;
      const y = chartHeight - ((point.price - minPrice + padding) / (priceRange + 2 * padding)) * chartHeight;
      return `${x},${y}`;
    })
    .join(" ");

  const firstPrice = data[data.length - 1]?.price ?? 0;
  const lastPrice = data[0]?.price ?? 0;
  const priceChange = lastPrice - firstPrice;
  const priceChangePercent = firstPrice ? (priceChange / firstPrice) * 100 : 0;
  const isPositive = priceChange >= 0;

  return (
    <div className="rounded-lg border border-slate-800 bg-slate-900 p-4">
      <div className="flex items-center justify-between mb-3">
        <div>
          <span className="text-lg font-semibold text-white">{ticker}</span>
          {currentPrice && (
            <span className="ml-3 text-2xl font-bold text-white">
              ${currentPrice.toFixed(2)}
            </span>
          )}
        </div>
        <div className={`text-sm font-medium ${isPositive ? "text-green-400" : "text-red-400"}`}>
          {isPositive ? "+" : ""}
          {priceChange.toFixed(2)} ({isPositive ? "+" : ""}
          {priceChangePercent.toFixed(2)}%)
          <span className="text-slate-500 ml-1 text-xs">90d</span>
        </div>
      </div>

      <div className="relative">
        <svg
          viewBox={`0 0 ${chartWidth} ${chartHeight}`}
          className="w-full h-40"
          preserveAspectRatio="none"
        >
          <defs>
            <linearGradient id="chartGradient" x1="0" y1="0" x2="0" y2="1">
              <stop
                offset="0%"
                stopColor={isPositive ? "#22c55e" : "#ef4444"}
                stopOpacity="0.3"
              />
              <stop
                offset="100%"
                stopColor={isPositive ? "#22c55e" : "#ef4444"}
                stopOpacity="0"
              />
            </linearGradient>
          </defs>
          <polygon
            points={`0,${chartHeight} ${points} ${chartWidth},${chartHeight}`}
            fill="url(#chartGradient)"
          />
          <polyline
            points={points}
            fill="none"
            stroke={isPositive ? "#22c55e" : "#ef4444"}
            strokeWidth="2"
            vectorEffect="non-scaling-stroke"
          />
        </svg>

        <div className="absolute left-0 top-0 text-xs text-slate-500">
          ${maxPrice.toFixed(2)}
        </div>
        <div className="absolute left-0 bottom-0 text-xs text-slate-500">
          ${minPrice.toFixed(2)}
        </div>
        <div className="absolute right-0 bottom-0 text-xs text-slate-500">
          {data[0]?.date}
        </div>
        <div className="absolute left-1/2 -translate-x-1/2 bottom-0 text-xs text-slate-500">
          {data[Math.floor(data.length / 2)]?.date}
        </div>
      </div>
    </div>
  );
}
