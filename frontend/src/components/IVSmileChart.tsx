import { useState, useMemo, useRef, useEffect, useCallback } from "react";
import type { IVSmileData } from "../types/index.ts";
import { Spinner } from "./ui/Spinner.tsx";
import { Select } from "./ui/Select.tsx";
import { IVSurfaceChart } from "./IVSurfaceChart.tsx";

interface IVSmileChartProps {
  data: IVSmileData | null;
  loading: boolean;
  underlyingPrice?: number | null;
}

interface HoverPoint {
  strike: number;
  iv: number;
  optionType: string;
  moneyness: number;
  x: number;
  y: number;
}

interface ZoomState {
  scale: number;
  translateX: number;
  translateY: number;
}

export function IVSmileChart({ data, loading, underlyingPrice }: IVSmileChartProps) {
  const chartContainerRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const [dimensions, setDimensions] = useState({ width: 0, height: 0 });
  const [hoverPoint, setHoverPoint] = useState<HoverPoint | null>(null);
  const [zoom, setZoom] = useState<ZoomState>({ scale: 1, translateX: 0, translateY: 0 });
  const [isDragging, setIsDragging] = useState(false);
  const [dragStart, setDragStart] = useState<{ x: number; y: number } | null>(null);
  const [viewMode, setViewMode] = useState<"2d" | "3d">("2d");

  useEffect(() => {
    if (!chartContainerRef.current) return;

    // Initial size check
    const { offsetWidth, offsetHeight } = chartContainerRef.current;
    if (offsetWidth > 0 && offsetHeight > 0) {
      setDimensions({ width: offsetWidth, height: offsetHeight });
    }

    const resizeObserver = new ResizeObserver((entries) => {
      if (!entries || entries.length === 0) return;
      const { width, height } = entries[0].contentRect;
      setDimensions({ width, height });
    });

    resizeObserver.observe(chartContainerRef.current);
    return () => resizeObserver.disconnect();
  }, []);

  useEffect(() => {
    const container = chartContainerRef.current;
    if (!container) return;

    const preventScroll = (e: WheelEvent) => {
      e.preventDefault();
    };

    container.addEventListener('wheel', preventScroll, { passive: false });
    return () => container.removeEventListener('wheel', preventScroll);
  }, []);

  const chartWidth = dimensions.width || 500;
  const chartHeight = dimensions.height || 200;
  const padding = { left: 50, right: 30, top: 20, bottom: 40 };

  const expirations = useMemo(() => {
    if (!data) return [];
    return Object.keys(data.smiles_by_expiry).sort();
  }, [data]);

  const [selectedExpiration, setSelectedExpiration] = useState<string>("");

  useEffect(() => {
    if (expirations.length === 0) {
      setSelectedExpiration("");
      return;
    }
    if (!selectedExpiration || !expirations.includes(selectedExpiration)) {
      setSelectedExpiration(expirations[0]);
    }
  }, [expirations, selectedExpiration]);

  const currentExpiration = selectedExpiration || expirations[0] || "";
  const smileData = data?.smiles_by_expiry[currentExpiration] || [];

  const handleMouseMove = useCallback(
    (e: React.MouseEvent<SVGSVGElement>) => {
      if (!svgRef.current || smileData.length === 0 || isDragging) return;

      const allStrikes = smileData.map((p) => p.strike);
      const allIVs = smileData.map((p) => p.implied_volatility * 100);
      const minStrike = Math.min(...allStrikes);
      const maxStrike = Math.max(...allStrikes);
      const minIV = Math.min(...allIVs);
      const maxIV = Math.max(...allIVs);
      const strikeRange = maxStrike === minStrike ? 10 : maxStrike - minStrike;
      const ivRange = maxIV === minIV ? 10 : maxIV - minIV;

      const getX = (strike: number) =>
        padding.left +
        ((strike - minStrike) / strikeRange) * (chartWidth - padding.left - padding.right);

      const getY = (iv: number) =>
        padding.top +
        (1 - (iv - minIV) / ivRange) * (chartHeight - padding.top - padding.bottom);

      const rect = svgRef.current.getBoundingClientRect();
      const svgX = e.clientX - rect.left;
      const svgY = e.clientY - rect.top;
      const svgWidth = rect.width;
      const svgHeight = rect.height;
      let normalizedX = (svgX / svgWidth) * chartWidth;
      let normalizedY = (svgY / svgHeight) * chartHeight;

      normalizedX = (normalizedX - zoom.translateX) / zoom.scale;
      normalizedY = (normalizedY - zoom.translateY) / zoom.scale;

      if (
        normalizedX < padding.left ||
        normalizedX > chartWidth - padding.right ||
        normalizedY < padding.top ||
        normalizedY > chartHeight - padding.bottom
      ) {
        setHoverPoint(null);
        return;
      }

      let closestPoint: typeof smileData[0] | null = null;
      let minDistance = Infinity;

      for (const point of smileData) {
        const px = getX(point.strike);
        const py = getY(point.implied_volatility * 100);
        const distance = Math.sqrt(
          Math.pow(normalizedX - px, 2) + Math.pow(normalizedY - py, 2)
        );

        if (distance < minDistance && distance < 30) {
          minDistance = distance;
          closestPoint = point;
        }
      }

      if (closestPoint) {
        setHoverPoint({
          strike: closestPoint.strike,
          iv: closestPoint.implied_volatility * 100,
          optionType: closestPoint.option_type,
          moneyness: closestPoint.moneyness,
          x: getX(closestPoint.strike),
          y: getY(closestPoint.implied_volatility * 100),
        });
      } else {
        setHoverPoint(null);
      }
    },
    [smileData, chartWidth, chartHeight, padding, zoom, isDragging]
  );

  const handleMouseLeave = useCallback(() => {
    setHoverPoint(null);
  }, []);

  const handleWheel = useCallback(
    (e: React.WheelEvent<SVGSVGElement>) => {
      e.preventDefault();
      if (!svgRef.current) return;

      const rect = svgRef.current.getBoundingClientRect();
      const mouseX = e.clientX - rect.left;
      const mouseY = e.clientY - rect.top;

      const delta = e.deltaY > 0 ? 0.9 : 1.1;
      const newScale = Math.max(1, Math.min(10, zoom.scale * delta));

      if (newScale === 1) {
        setZoom({ scale: 1, translateX: 0, translateY: 0 });
      } else {
        const factor = newScale / zoom.scale;
        const newTranslateX = mouseX - factor * (mouseX - zoom.translateX);
        const newTranslateY = mouseY - factor * (mouseY - zoom.translateY);
        setZoom({ scale: newScale, translateX: newTranslateX, translateY: newTranslateY });
      }
    },
    [zoom]
  );

  const handleMouseDown = useCallback((e: React.MouseEvent<SVGSVGElement>) => {
    if (e.button === 0 && zoom.scale > 1) {
      setIsDragging(true);
      setDragStart({ x: e.clientX - zoom.translateX, y: e.clientY - zoom.translateY });
      e.preventDefault();
    }
  }, [zoom]);

  const handleMouseMoveForPan = useCallback(
    (e: React.MouseEvent<SVGSVGElement>) => {
      if (isDragging && dragStart && zoom.scale > 1) {
        setZoom({
          ...zoom,
          translateX: e.clientX - dragStart.x,
          translateY: e.clientY - dragStart.y,
        });
      }
    },
    [isDragging, dragStart, zoom]
  );

  const handleMouseUp = useCallback(() => {
    setIsDragging(false);
    setDragStart(null);
  }, []);

  const handleDoubleClick = useCallback(() => {
    setZoom({ scale: 1, translateX: 0, translateY: 0 });
  }, []);

  const renderContent = () => {
    if (!data) {
      return (
        <div className="absolute inset-0 flex items-center justify-center">
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
        <div className="absolute inset-0 flex items-center justify-center">
          <Spinner />
          <span className="ml-2 text-slate-400">Loading option chain...</span>
        </div>
      );
    }

    if (expirations.length === 0 || smileData.length === 0) {
      return (
        <div className="absolute inset-0 flex items-center justify-center">
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

    const strikeRange = maxStrike === minStrike ? 10 : maxStrike - minStrike;
    const ivRange = maxIV === minIV ? 10 : maxIV - minIV;

    const getX = (strike: number) =>
      padding.left +
      ((strike - minStrike) / strikeRange) * (chartWidth - padding.left - padding.right);

    const getY = (iv: number) =>
      padding.top +
      (1 - (iv - minIV) / ivRange) * (chartHeight - padding.top - padding.bottom);

    const spotPrice = underlyingPrice ?? data.underlying_price;
    const atmX = getX(spotPrice);

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
      <div className="w-full h-full flex flex-col">
        <div className="flex items-center justify-between mb-4 z-10 relative">
          <div>
            <span className="text-lg font-semibold text-white">
              Implied Volatility {viewMode === "2d" ? "Smile" : "Surface"}
            </span>
            <span className="ml-3 text-sm text-slate-400">
              {data.symbol} @ ${spotPrice.toFixed(2)}
            </span>
            {zoom.scale > 1 && viewMode === "2d" && (
              <span className="ml-3 text-xs text-purple-400">
                {zoom.scale.toFixed(1)}x zoom
              </span>
            )}
          </div>
          <div className="flex items-center gap-3">
            <div className="flex bg-slate-800 rounded-md p-1">
              <button
                onClick={() => setViewMode("2d")}
                className={`px-3 py-1 text-xs rounded-sm transition-colors ${
                  viewMode === "2d"
                    ? "bg-purple-600 text-white"
                    : "text-slate-400 hover:text-slate-200"
                }`}
              >
                2D Smile
              </button>
              <button
                onClick={() => setViewMode("3d")}
                className={`px-3 py-1 text-xs rounded-sm transition-colors ${
                  viewMode === "3d"
                    ? "bg-purple-600 text-white"
                    : "text-slate-400 hover:text-slate-200"
                }`}
              >
                3D Surface
              </button>
            </div>
            {zoom.scale > 1 && viewMode === "2d" && (
              <button
                onClick={() => setZoom({ scale: 1, translateX: 0, translateY: 0 })}
                className="px-3 py-1 text-xs text-purple-400 hover:text-purple-300 border border-purple-500/30 hover:border-purple-500/50 rounded transition-colors"
              >
                Reset Zoom
              </button>
            )}
            {viewMode === "2d" && (
              <div className="w-48">
                <Select
                  label="Expiration"
                  value={currentExpiration}
                  onChange={setSelectedExpiration}
                  options={expirations.map((exp) => ({ value: exp, label: exp }))}
                />
              </div>
            )}
          </div>
        </div>

        {viewMode === "3d" ? (
          <div className="flex-1 min-h-[240px]">
            <IVSurfaceChart data={data} loading={loading} />
          </div>
        ) : (
          <>
            <div
              ref={chartContainerRef}
              className="relative flex-1 min-h-[240px]"
            >
              <svg
                ref={svgRef}
                width={chartWidth}
            height={chartHeight}
            viewBox={`0 0 ${chartWidth} ${chartHeight}`}
            className="w-full h-full absolute inset-0 cursor-crosshair"
            style={{ cursor: isDragging ? "grabbing" : zoom.scale > 1 ? "grab" : "crosshair" }}
            preserveAspectRatio="none"
            onMouseMove={(e) => {
              handleMouseMoveForPan(e);
              handleMouseMove(e);
            }}
            onMouseDown={handleMouseDown}
            onMouseUp={handleMouseUp}
            onMouseLeave={() => {
              handleMouseLeave();
              handleMouseUp();
            }}
            onWheel={handleWheel}
            onDoubleClick={handleDoubleClick}
          >
            <g transform={`translate(${zoom.translateX}, ${zoom.translateY}) scale(${zoom.scale})`}>
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
                    className="text-[10px] fill-slate-400"
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
                <line
                  x1={padding.left}
                  y1={hoverPoint.y}
                  x2={chartWidth - padding.right}
                  y2={hoverPoint.y}
                  stroke="#8b5cf6"
                  strokeWidth="1"
                  strokeDasharray="4 4"
                  vectorEffect="non-scaling-stroke"
                />
                <circle
                  cx={hoverPoint.x}
                  cy={hoverPoint.y}
                  r="5"
                  fill={hoverPoint.optionType === "Call" ? "#3b82f6" : "#ef4444"}
                  stroke="white"
                  strokeWidth="2"
                />
              </>
            )}
            </g>

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

          {hoverPoint && (
            <div
              className="absolute bg-slate-800 text-white text-xs px-3 py-2 rounded shadow-lg pointer-events-none z-20"
              style={{
                left: hoverPoint.x,
                top: hoverPoint.y,
                transform: "translate(-50%, -120%)",
              }}
            >
              <div className="flex items-center gap-2 mb-1">
                <div
                  className={`w-2 h-2 rounded-full ${
                    hoverPoint.optionType === "Call" ? "bg-blue-500" : "bg-red-500"
                  }`}
                />
                <span className="font-semibold">{hoverPoint.optionType}</span>
              </div>
              <div className="text-slate-300">
                <div>Strike: ${hoverPoint.strike.toFixed(2)}</div>
                <div>IV: {hoverPoint.iv.toFixed(2)}%</div>
                <div className="text-slate-400 text-[10px] mt-1">
                  Moneyness: {hoverPoint.moneyness.toFixed(3)}
                </div>
              </div>
            </div>
          )}
        </div>

        <div className="flex items-center justify-center gap-6 text-sm mt-4">
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
            <span className="text-slate-400">
              ATM (Underlying: {spotPrice.toFixed(2)})
            </span>
          </div>
        </div>
        </>
        )}
      </div>
    );
  };

  return (
    <div
      className="relative w-full h-full min-h-[400px] rounded-lg border border-slate-800 bg-slate-900 p-4 overflow-hidden"
    >
      {renderContent()}
    </div>
  );
}
