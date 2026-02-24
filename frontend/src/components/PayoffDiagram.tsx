import Plot from "react-plotly.js";
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
      <div className="relative h-64 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
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
      <div className="h-64 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
        <Spinner />
        <span className="ml-2 text-slate-400">Calculating payoff...</span>
      </div>
    );
  }

  const payoffData = calculatePayoffGrid(request, priceResult);

  if (!payoffData) {
    return (
      <div className="h-64 rounded-lg border border-slate-800 bg-slate-900 flex items-center justify-center">
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

  const currentSpot =
    request.structure_type === "single"
      ? request.spot_price
      : request.structure_type === "spread"
        ? request.spot_price
        : 0;

  const positivePayoffs = payoffs.map((p) => (p >= 0 ? p : 0));
  const negativePayoffs = payoffs.map((p) => (p <= 0 ? p : 0));

  const plotData: any[] = [
    {
      x: spotPrices,
      y: positivePayoffs,
      type: "scatter",
      mode: "lines",
      line: { width: 0 },
      fill: "tozeroy",
      fillcolor: "rgba(34, 197, 94, 0.2)",
      hoverinfo: "skip",
      showlegend: false,
    },
    {
      x: spotPrices,
      y: negativePayoffs,
      type: "scatter",
      mode: "lines",
      line: { width: 0 },
      fill: "tozeroy",
      fillcolor: "rgba(239, 68, 68, 0.2)",
      hoverinfo: "skip",
      showlegend: false,
    },
    {
      x: spotPrices,
      y: payoffs,
      type: "scatter",
      mode: "lines",
      line: { color: "#60a5fa", width: 2 },
      name: "Payoff",
      hovertemplate: "Spot: $%{x:.2f}<br>Payoff: $%{y:.2f}<extra></extra>",
    },
  ];

  if (breakEvenPoints.length > 0) {
    plotData.push({
      x: breakEvenPoints,
      y: breakEvenPoints.map(() => 0),
      type: "scatter",
      mode: "markers",
      marker: { color: "#fbbf24", size: 8, line: { color: "#92400e", width: 1 } },
      name: "Break-even",
      hovertemplate: "Break-even: $%{x:.2f}<extra></extra>",
    });
  }

  const layout: any = {
    autosize: true,
    margin: { l: 50, r: 20, t: 20, b: 40 },
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
    xaxis: {
      title: {
        text:
          request.structure_type === "spread" &&
          request.spread_type === "calendar_spread"
            ? "Underlying Asset Price at Short-Term Expiration"
            : "Underlying Asset Price at Expiration",
        font: { color: "#94a3b8", size: 12 },
      },
      tickfont: { color: "#94a3b8" },
      gridcolor: "#334155",
      zerolinecolor: "#475569",
    },
    yaxis: {
      title: {
        text: "Payoff ($)",
        font: { color: "#94a3b8", size: 12 },
      },
      tickfont: { color: "#94a3b8" },
      gridcolor: "#334155",
      zeroline: true,
      zerolinecolor: "#64748b",
      zerolinewidth: 1,
    },
    showlegend: false,
    hovermode: "x unified",
    shapes: [
      {
        type: "line",
        x0: currentSpot,
        x1: currentSpot,
        y0: 0,
        y1: 1,
        yref: "paper",
        line: {
          color: "#8b5cf6",
          width: 1,
          dash: "dash",
        },
      },
    ],
  };

  return (
    <div className="rounded-lg border border-slate-800 bg-slate-900 p-4 pb-4">
      <div className="flex items-center justify-between mb-3">
        <div>
          <span className="text-lg font-semibold text-white">
            {request.structure_type === "spread" &&
            request.spread_type === "calendar_spread"
              ? "Payoff at Short-Term Expiration"
              : "Payoff at Expiration"}
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

      <div className="w-full h-64">
        <Plot
          data={plotData}
          layout={layout}
          useResizeHandler={true}
          style={{ width: "100%", height: "100%" }}
          config={{ displayModeBar: false, responsive: true }}
        />
      </div>
    </div>
  );
}
