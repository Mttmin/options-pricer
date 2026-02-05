import type { PriceRequest, PriceResponse } from "../types/index.ts";
import { calculatePayoffGrid } from "../utils/payoff.ts";
import { Spinner } from "./ui/Spinner.tsx";

interface PayoffDiagramProps {
  request: PriceRequest | null;
  priceResult: PriceResponse | null;
  loading: boolean;
}

export function PayoffDiagram({
  request,
  priceResult,
  loading,
}: PayoffDiagramProps) {
  if (!request || !priceResult) {
    return (
      <div className="relative h-48 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <div className="absolute inset-0 flex items-center justify-center">
          <svg
            className="h-full w-full text-slate-800 opacity-30"
            viewBox="0 0 400 150"
            preserveAspectRatio="none"
          >
            <line
              x1="0"
              y1="75"
              x2="400"
              y2="75"
              stroke="currentColor"
              strokeWidth="1"
              strokeDasharray="4"
            />
            <path
              d="M0,120 Q100,80 200,75 T400,40"
              fill="none"
              stroke="currentColor"
              strokeWidth="2"
            />
          </svg>
        </div>
        <p className="text-slate-500 text-sm z-10 bg-slate-900/80 px-3 py-1 rounded">
          Calculate option price to view payoff diagram
        </p>
      </div>
    );
  }

  if (loading) {
    return (
      <div className="h-48 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <Spinner />
        <span className="ml-2 text-slate-400">Calculating payoff...</span>
      </div>
    );
  }

  const payoffData = calculatePayoffGrid(request, priceResult);

  if (!payoffData) {
    return (
      <div className="h-48 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <p className="text-slate-500 text-sm text-center px-4">
          Payoff diagram not available for exotic options.
          <br />
          Use Greeks for sensitivity analysis.
        </p>
      </div>
    );
  }

  const { spotPrices, payoffs, breakEvenPoints, maxProfit, maxLoss } =
    payoffData;

  const minPayoff = Math.min(...payoffs);
  const maxPayoff = Math.max(...payoffs);
  const payoffRange = maxPayoff - minPayoff || 1;
  const padding = payoffRange * 0.15;

  const chartHeight = 150;
  const chartWidth = 400;

  const zeroY =
    chartHeight -
    ((0 - minPayoff + padding) / (payoffRange + 2 * padding)) * chartHeight;

  const points = spotPrices
    .map((_, i) => {
      const x = (i / (spotPrices.length - 1)) * chartWidth;
      const y =
        chartHeight -
        ((payoffs[i] - minPayoff + padding) / (payoffRange + 2 * padding)) *
          chartHeight;
      return `${x},${y}`;
    })
    .join(" ");

  const positiveArea: string[] = [];
  const negativeArea: string[] = [];
  let isFirstPositive = true;
  let isFirstNegative = true;

  for (let i = 0; i < spotPrices.length; i++) {
    const x = (i / (spotPrices.length - 1)) * chartWidth;
    const y =
      chartHeight -
      ((payoffs[i] - minPayoff + padding) / (payoffRange + 2 * padding)) *
        chartHeight;

    if (payoffs[i] >= 0) {
      if (isFirstPositive) {
        positiveArea.push(`M${x},${zeroY}`);
        isFirstPositive = false;
      }
      positiveArea.push(`L${x},${y}`);
    }

    if (payoffs[i] <= 0) {
      if (isFirstNegative) {
        negativeArea.push(`M${x},${zeroY}`);
        isFirstNegative = false;
      }
      negativeArea.push(`L${x},${y}`);
    }
  }

  if (positiveArea.length > 0) {
    const lastX = ((spotPrices.length - 1) / (spotPrices.length - 1)) * chartWidth;
    positiveArea.push(`L${lastX},${zeroY} Z`);
  }

  if (negativeArea.length > 0) {
    const lastX = ((spotPrices.length - 1) / (spotPrices.length - 1)) * chartWidth;
    negativeArea.push(`L${lastX},${zeroY} Z`);
  }

  const currentSpot =
    request.structure_type === "single"
      ? request.spot_price
      : request.structure_type === "spread"
        ? request.spot_price
        : 0;
  const currentSpotX =
    ((currentSpot - spotPrices[0]) / (spotPrices[spotPrices.length - 1] - spotPrices[0])) *
    chartWidth;

  return (
    <div className="rounded-lg border border-slate-800 bg-slate-900 p-4">
      <div className="flex items-center justify-between mb-3">
        <div>
          <span className="text-lg font-semibold text-white">
            Payoff at Expiration
          </span>
          {maxProfit !== null && (
            <span className="ml-3 text-sm text-green-400">
              Max Profit: ${maxProfit.toFixed(2)}
            </span>
          )}
          {maxLoss !== null && (
            <span className="ml-3 text-sm text-red-400">
              Max Loss: ${maxLoss.toFixed(2)}
            </span>
          )}
        </div>
        <div className="text-sm text-slate-400">
          {breakEvenPoints.length > 0 && (
            <span>
              Break-even: ${breakEvenPoints.map((p) => p.toFixed(2)).join(", $")}
            </span>
          )}
        </div>
      </div>

      <div className="relative">
        <svg
          viewBox={`0 0 ${chartWidth} ${chartHeight}`}
          className="w-full h-40"
          preserveAspectRatio="none"
        >
          <defs>
            <linearGradient id="profitGradient" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%" stopColor="#22c55e" stopOpacity="0.3" />
              <stop offset="100%" stopColor="#22c55e" stopOpacity="0" />
            </linearGradient>
            <linearGradient id="lossGradient" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%" stopColor="#ef4444" stopOpacity="0" />
              <stop offset="100%" stopColor="#ef4444" stopOpacity="0.3" />
            </linearGradient>
          </defs>

          {positiveArea.length > 0 && (
            <path d={positiveArea.join(" ")} fill="url(#profitGradient)" />
          )}
          {negativeArea.length > 0 && (
            <path d={negativeArea.join(" ")} fill="url(#lossGradient)" />
          )}

          <line
            x1="0"
            y1={zeroY}
            x2={chartWidth}
            y2={zeroY}
            stroke="#64748b"
            strokeWidth="1"
            strokeDasharray="4 4"
            vectorEffect="non-scaling-stroke"
          />

          <polyline
            points={points}
            fill="none"
            stroke="#60a5fa"
            strokeWidth="2"
            vectorEffect="non-scaling-stroke"
          />

          <line
            x1={currentSpotX}
            y1="0"
            x2={currentSpotX}
            y2={chartHeight}
            stroke="#8b5cf6"
            strokeWidth="1"
            strokeDasharray="2 2"
            vectorEffect="non-scaling-stroke"
          />

          {breakEvenPoints.map((breakEven, idx) => {
            const x =
              ((breakEven - spotPrices[0]) /
                (spotPrices[spotPrices.length - 1] - spotPrices[0])) *
              chartWidth;
            return (
              <circle
                key={idx}
                cx={x}
                cy={zeroY}
                r="3"
                fill="#fbbf24"
                stroke="#92400e"
                strokeWidth="1"
              />
            );
          })}
        </svg>

        <div className="absolute left-0 top-0 text-xs text-slate-500">
          ${maxPayoff.toFixed(2)}
        </div>
        <div className="absolute left-0 bottom-0 text-xs text-slate-500">
          ${minPayoff.toFixed(2)}
        </div>
        <div className="absolute right-0 bottom-0 text-xs text-slate-500">
          ${spotPrices[spotPrices.length - 1].toFixed(2)}
        </div>
        <div className="absolute left-1/2 -translate-x-1/2 bottom-0 text-xs text-slate-500">
          ${spotPrices[Math.floor(spotPrices.length / 2)].toFixed(2)}
        </div>
        <div className="absolute left-0 bottom-0 text-xs text-slate-500">
          ${spotPrices[0].toFixed(2)}
        </div>
      </div>
    </div>
  );
}
