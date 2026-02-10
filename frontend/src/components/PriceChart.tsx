import { useState, useRef, useCallback, useEffect } from "react";
import type { PricePoint } from "../api/client.ts";
import { Spinner } from "./ui/Spinner.tsx";

interface PriceChartProps {
  data: PricePoint[];
  ticker: string | null;
  currentPrice?: number;
  loading: boolean;
}

export function PriceChart({ data, ticker, currentPrice, loading }: PriceChartProps) {
  const [hoverPoint, setHoverPoint] = useState<{ price: number; date: string; x: number; y: number } | null>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const [dimensions, setDimensions] = useState({ width: 0, height: 0 });

  useEffect(() => {
    if (!containerRef.current) return;

    const resizeObserver = new ResizeObserver((entries) => {
      if (!entries || entries.length === 0) return;
      const { width, height } = entries[0].contentRect;
      setDimensions({ width, height });
    });

    resizeObserver.observe(containerRef.current);

    return () => resizeObserver.disconnect();
  }, []);

  const chartWidth = dimensions.width || 800;
  const chartHeight = dimensions.height || 400;
  const padding = { left: 50, right: 20, top: 20, bottom: 30 };

  const handleMouseMove = useCallback(
    (e: React.MouseEvent<SVGSVGElement>) => {
      if (!svgRef.current || data.length === 0) return;

      const sortedData = [...data].reverse();
      const prices = sortedData.map((p) => p.price);
      const minPrice = Math.min(...prices);
      const maxPrice = Math.max(...prices);
      const priceRange = maxPrice - minPrice || 1;

      const getX = (index: number) =>
        padding.left + (index / (sortedData.length - 1)) * (chartWidth - padding.left - padding.right);

      const getY = (price: number) =>
        padding.top + (1 - (price - minPrice) / priceRange) * (chartHeight - padding.top - padding.bottom);

      const rect = svgRef.current.getBoundingClientRect();
      const svgX = e.clientX - rect.left;
      const svgWidth = rect.width;
      const normalizedX = (svgX / svgWidth) * chartWidth;

      if (normalizedX < padding.left || normalizedX > chartWidth - padding.right) {
        setHoverPoint(null);
        return;
      }

      const dataX = ((normalizedX - padding.left) / (chartWidth - padding.left - padding.right)) * (sortedData.length - 1);
      const index = Math.round(dataX);
      const clampedIndex = Math.max(0, Math.min(sortedData.length - 1, index));

      const point = sortedData[clampedIndex];
      if (point) {
        setHoverPoint({
          price: point.price,
          date: point.date,
          x: getX(clampedIndex),
          y: getY(point.price),
        });
      }
    },
    [data, chartWidth, chartHeight, padding]
  );

  const handleMouseLeave = useCallback(() => {
    setHoverPoint(null);
  }, []);

  const renderContent = () => {
    if (!ticker) {
      return (
        <div className="absolute inset-0 flex items-center justify-center">
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
        <div className="absolute inset-0 flex items-center justify-center">
          <Spinner />
          <span className="ml-2 text-slate-400">Loading chart...</span>
        </div>
      );
    }

    if (data.length === 0) {
      return (
        <div className="absolute inset-0 flex items-center justify-center">
          <p className="text-slate-500 text-sm">No price data available</p>
        </div>
      );
    }

    const sortedData = [...data].reverse();
    const prices = sortedData.map((p) => p.price);
    const minPrice = Math.min(...prices);
    const maxPrice = Math.max(...prices);
    const priceRange = maxPrice - minPrice || 1;

    const getX = (index: number) =>
      padding.left + (index / (data.length - 1)) * (chartWidth - padding.left - padding.right);

    const getY = (price: number) =>
      padding.top + (1 - (price - minPrice) / priceRange) * (chartHeight - padding.top - padding.bottom);

    const points = sortedData
      .map((point, i) => `${getX(i)},${getY(point.price)}`)
      .join(" ");

    const numGridLinesY = 5;
    const yTicks = Array.from({ length: numGridLinesY }, (_, i) =>
      minPrice + (i * priceRange) / (numGridLinesY - 1)
    );

    const firstPrice = sortedData[0]?.price ?? 0;
    const lastPrice = sortedData[sortedData.length - 1]?.price ?? 0;
    const priceChange = lastPrice - firstPrice;
    const priceChangePercent = firstPrice ? (priceChange / firstPrice) * 100 : 0;
    const isPositive = priceChange >= 0;

    return (
      <>
        <div className="absolute top-4 left-4 z-10">
          <div className="flex items-center gap-3">
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
        </div>

        <svg
          ref={svgRef}
          width={chartWidth}
          height={chartHeight}
          viewBox={`0 0 ${chartWidth} ${chartHeight}`}
          className="w-full h-full absolute inset-0"
          preserveAspectRatio="none"
          onMouseMove={handleMouseMove}
          onMouseLeave={handleMouseLeave}
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
                  className="text-[10px] fill-slate-400"
                >
                  ${tick.toFixed(2)}
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

          <polygon
            points={`${padding.left},${chartHeight - padding.bottom} ${points} ${chartWidth - padding.right},${chartHeight - padding.bottom}`}
            fill="url(#chartGradient)"
          />
          <polyline
            points={points}
            fill="none"
            stroke={isPositive ? "#22c55e" : "#ef4444"}
            strokeWidth="2"
            vectorEffect="non-scaling-stroke"
          />

          {hoverPoint && (
            <>
              <line
                x1={hoverPoint.x}
                y1={padding.top}
                x2={hoverPoint.x}
                y2={chartHeight - padding.bottom}
                stroke="#8b5cf6"
                strokeWidth="1"
                strokeDasharray="4 4"
                vectorEffect="non-scaling-stroke"
              />
              <circle
                cx={hoverPoint.x}
                cy={hoverPoint.y}
                r="4"
                fill={isPositive ? "#22c55e" : "#ef4444"}
                stroke="white"
                strokeWidth="2"
              />
            </>
          )}

          <text
            x={chartWidth / 2}
            y={chartHeight - 5}
            textAnchor="middle"
            className="text-[10px] fill-slate-400"
          >
            Date
          </text>
        </svg>

        {hoverPoint && (
          <div
            className="absolute bg-slate-800 text-white text-xs px-2 py-1 rounded shadow-lg pointer-events-none z-10"
            style={{
              left: hoverPoint.x,
              top: hoverPoint.y,
              transform: "translate(-50%, -120%)",
            }}
          >
            <div className="font-semibold">${hoverPoint.price.toFixed(2)}</div>
            <div className="text-slate-400">{hoverPoint.date}</div>
          </div>
        )}
      </>
    );
  };

  return (
    <div
      ref={containerRef}
      className="relative w-full h-full min-h-[400px] rounded-lg border border-slate-800 bg-slate-900 overflow-hidden"
    >
      {renderContent()}
    </div>
  );
}
