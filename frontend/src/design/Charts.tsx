import { useMemo } from "react";

type HistoryPoint = { date: string; price?: number; close?: number };

export function LineChart({
  data, width = 640, height = 180, color, yKey = "price",
}: {
  data: HistoryPoint[];
  width?: number;
  height?: number;
  color?: string;
  yKey?: "price" | "close";
}) {
  const { path, area, lastY, fill } = useMemo(() => {
    if (!data || data.length < 2) return { path: "", area: "", lastY: null, fill: true };
    const ys = data.map(d => (d[yKey] ?? d.close ?? d.price ?? 0));
    const minY = Math.min(...ys), maxY = Math.max(...ys);
    const pad = (maxY - minY) * 0.1 || 1;
    const yLo = minY - pad, yHi = maxY + pad;
    const sx = (i: number) => (i / (data.length - 1)) * width;
    const sy = (v: number) => height - ((v - yLo) / (yHi - yLo)) * height;
    const path = data.map((d, i) => {
      const y = d[yKey] ?? d.close ?? d.price ?? 0;
      return `${i ? "L" : "M"}${sx(i).toFixed(2)},${sy(y).toFixed(2)}`;
    }).join(" ");
    const area = path + ` L${width},${height} L0,${height} Z`;
    const lastY = sy(ys[ys.length - 1]);
    const up = ys[ys.length - 1] >= ys[0];
    return { path, area, lastY, fill: up };
  }, [data, width, height, yKey]);

  const strokeColor = color || (fill ? "var(--accent)" : "var(--loss)");
  const gradId = `lc-fill-${Math.random().toString(36).slice(2, 8)}`;

  return (
    <svg width="100%" viewBox={`0 0 ${width} ${height}`} preserveAspectRatio="none" style={{ display: "block" }}>
      <defs>
        <linearGradient id={gradId} x1="0" x2="0" y1="0" y2="1">
          <stop offset="0%" stopColor={strokeColor} stopOpacity="0.18" />
          <stop offset="100%" stopColor={strokeColor} stopOpacity="0" />
        </linearGradient>
      </defs>
      <path d={area} fill={`url(#${gradId})`} />
      <path d={path} fill="none" stroke={strokeColor} strokeWidth="1.5" strokeLinejoin="round" />
      {lastY != null && (
        <circle cx={width - 1} cy={lastY} r="2.5" fill={strokeColor} />
      )}
    </svg>
  );
}

export function Sparkline({
  data, width = 120, height = 22, color = "currentColor", zeroLine = true,
}: {
  data: { x: number; y: number }[];
  width?: number;
  height?: number;
  color?: string;
  zeroLine?: boolean;
}) {
  const { path, zeroY } = useMemo(() => {
    if (!data || data.length < 2) return { path: "", zeroY: null as number | null };
    const ys = data.map(d => d.y);
    let minY = Math.min(...ys, zeroLine ? 0 : Infinity);
    let maxY = Math.max(...ys, zeroLine ? 0 : -Infinity);
    const span = maxY - minY || 1;
    const pad = span * 0.12;
    minY -= pad; maxY += pad;
    const sx = (i: number) => (i / (data.length - 1)) * width;
    const sy = (v: number) => height - ((v - minY) / (maxY - minY)) * height;
    const path = data.map((d, i) => `${i ? "L" : "M"}${sx(i).toFixed(2)},${sy(d.y).toFixed(2)}`).join(" ");
    const zeroY = (maxY >= 0 && minY <= 0) ? sy(0) : null;
    return { path, zeroY };
  }, [data, width, height, zeroLine]);

  return (
    <svg width="100%" height={height} viewBox={`0 0 ${width} ${height}`} preserveAspectRatio="none">
      {zeroY != null && (
        <line x1="0" x2={width} y1={zeroY} y2={zeroY} stroke="var(--border-strong)" strokeDasharray="2 2" strokeWidth="0.5" />
      )}
      <path d={path} fill="none" stroke={color} strokeWidth="1.2" strokeLinejoin="round" />
    </svg>
  );
}

export function PayoffChart({
  data, width = 640, height = 200, strike, spot,
}: {
  data: { spot: number; pnl: number }[];
  width?: number;
  height?: number;
  strike?: number;
  spot?: number;
}) {
  const out = useMemo(() => {
    if (!data || !data.length) return null;
    const xs = data.map(d => d.spot), ys = data.map(d => d.pnl);
    const minX = Math.min(...xs), maxX = Math.max(...xs);
    let minY = Math.min(0, ...ys), maxY = Math.max(0, ...ys);
    const pad = (maxY - minY) * 0.1 || 1;
    minY -= pad; maxY += pad;
    const sx = (v: number) => ((v - minX) / (maxX - minX)) * width;
    const sy = (v: number) => height - ((v - minY) / (maxY - minY)) * height;
    const zeroY = sy(0);
    const strikeX = strike != null ? sx(strike) : null;
    const spotX = spot != null ? sx(spot) : null;

    const posPath = data.map((d, i) => `${i ? "L" : "M"}${sx(d.spot).toFixed(2)},${sy(Math.max(0, d.pnl)).toFixed(2)}`).join(" ") +
      ` L${sx(maxX)},${zeroY} L${sx(minX)},${zeroY} Z`;
    const negPath = data.map((d, i) => `${i ? "L" : "M"}${sx(d.spot).toFixed(2)},${sy(Math.min(0, d.pnl)).toFixed(2)}`).join(" ") +
      ` L${sx(maxX)},${zeroY} L${sx(minX)},${zeroY} Z`;
    const linePath = data.map((d, i) => `${i ? "L" : "M"}${sx(d.spot).toFixed(2)},${sy(d.pnl).toFixed(2)}`).join(" ");
    return { posPath, negPath, linePath, zeroY, strikeX, spotX };
  }, [data, width, height, strike, spot]);

  if (!out) return null;

  return (
    <svg width="100%" viewBox={`0 0 ${width} ${height}`} preserveAspectRatio="none" style={{ display: "block" }}>
      <path d={out.posPath} fill="var(--accent)" opacity="0.15" />
      <path d={out.negPath} fill="var(--loss)" opacity="0.15" />
      <line x1="0" x2={width} y1={out.zeroY} y2={out.zeroY} stroke="var(--border-strong)" strokeWidth="0.5" />
      {out.strikeX != null && (
        <g>
          <line x1={out.strikeX} x2={out.strikeX} y1="0" y2={height} stroke="var(--fg-3)" strokeDasharray="3 3" strokeWidth="0.5" />
          <text x={out.strikeX + 4} y="12" fill="var(--fg-3)" fontSize="9" fontFamily="var(--font-mono)">K</text>
        </g>
      )}
      {out.spotX != null && (
        <g>
          <line x1={out.spotX} x2={out.spotX} y1="0" y2={height} stroke="var(--accent)" strokeDasharray="2 2" strokeWidth="0.5" opacity="0.6" />
          <text x={out.spotX + 4} y="24" fill="var(--accent)" fontSize="9" fontFamily="var(--font-mono)">S</text>
        </g>
      )}
      <path d={out.linePath} fill="none" stroke="var(--fg)" strokeWidth="1.5" strokeLinejoin="round" />
    </svg>
  );
}

export type SmileChartPoint = { strike: number; iv: number; oi?: number };

export function IVSmileChart({
  points, width = 600, height = 220, spot,
}: {
  points: SmileChartPoint[];
  width?: number;
  height?: number;
  spot?: number;
}) {
  const out = useMemo(() => {
    if (!points || !points.length) return null;
    const xs = points.map(p => p.strike), ys = points.map(p => p.iv);
    const minX = Math.min(...xs), maxX = Math.max(...xs);
    const minY = Math.max(0, Math.min(...ys) - 0.02);
    const maxY = Math.max(...ys) + 0.02;
    const sx = (v: number) => 40 + ((v - minX) / (maxX - minX || 1)) * (width - 60);
    const sy = (v: number) => height - 24 - ((v - minY) / (maxY - minY || 1)) * (height - 40);
    const path = points.map((p, i) => `${i ? "L" : "M"}${sx(p.strike).toFixed(2)},${sy(p.iv).toFixed(2)}`).join(" ");
    const yTicks = [minY, (minY + maxY) / 2, maxY].map(v => ({ v, y: sy(v) }));
    return { path, sx, sy, yTicks };
  }, [points, width, height]);

  if (!out) return null;

  return (
    <svg width="100%" viewBox={`0 0 ${width} ${height}`}>
      {out.yTicks.map((t, i) => (
        <g key={i}>
          <line x1="40" x2={width - 20} y1={t.y} y2={t.y} stroke="var(--border)" strokeWidth="0.5" strokeDasharray="2 4" />
          <text x="34" y={t.y + 3} textAnchor="end" fontSize="9" fill="var(--fg-3)" fontFamily="var(--font-mono)">
            {(t.v * 100).toFixed(0)}%
          </text>
        </g>
      ))}
      {spot != null && (
        <g>
          <line x1={out.sx(spot)} x2={out.sx(spot)} y1="10" y2={height - 24} stroke="var(--accent)" strokeDasharray="2 2" strokeWidth="0.5" />
          <text x={out.sx(spot)} y="9" textAnchor="middle" fontSize="9" fill="var(--accent)" fontFamily="var(--font-mono)">spot {spot.toFixed(2)}</text>
        </g>
      )}
      <path d={out.path} fill="none" stroke="var(--accent)" strokeWidth="1.5" />
      {points.map((p, i) => (
        <circle
          key={i}
          cx={out.sx(p.strike)}
          cy={out.sy(p.iv)}
          r={Math.min(4, 1.5 + Math.sqrt(p.oi ?? 0) / 80)}
          fill="var(--fg)"
          opacity="0.4"
        />
      ))}
      {points.filter((_, i) => i % 3 === 0).map((p, i) => (
        <text key={i} x={out.sx(p.strike)} y={height - 8} textAnchor="middle" fontSize="9" fill="var(--fg-3)" fontFamily="var(--font-mono)">
          {p.strike}
        </text>
      ))}
    </svg>
  );
}
