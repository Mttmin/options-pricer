import { useMemo } from "react";
import Plot from "react-plotly.js";
import type { IVSmileData } from "../types/index.ts";

interface IVSurfaceChartProps {
  data: IVSmileData | null;
  loading: boolean;
}

export function IVSurfaceChart({ data, loading }: IVSurfaceChartProps) {
  const plotData = useMemo(() => {
    if (!data) return null;

    const expirations = Object.keys(data.smiles_by_expiry).sort();
    if (expirations.length === 0) return null;

    // Collect all unique strikes across all expirations
    const allStrikesSet = new Set<number>();
    expirations.forEach((exp) => {
      data.smiles_by_expiry[exp].forEach((point) => {
        allStrikesSet.add(point.strike);
      });
    });

    const strikes = Array.from(allStrikesSet).sort((a, b) => a - b);
    
    // Create Z matrix (IVs)
    // Z[y][x] where y is expiration index, x is strike index
    const z: (number | null)[][] = [];

    expirations.forEach((exp) => {
      const row: (number | null)[] = new Array(strikes.length).fill(null);
      const points = data.smiles_by_expiry[exp];
      
      // Average IV if there are both Call and Put for the same strike
      const strikeToIV = new Map<number, { sum: number; count: number }>();
      
      points.forEach((p) => {
        const existing = strikeToIV.get(p.strike) || { sum: 0, count: 0 };
        strikeToIV.set(p.strike, {
          sum: existing.sum + p.implied_volatility * 100,
          count: existing.count + 1,
        });
      });

      strikes.forEach((strike, xIdx) => {
        const ivData = strikeToIV.get(strike);
        if (ivData) {
          row[xIdx] = ivData.sum / ivData.count;
        }
      });
      
      z.push(row);
    });

    return {
      x: strikes,
      y: expirations,
      z: z,
    };
  }, [data]);

  if (loading) {
    return (
      <div className="absolute inset-0 flex items-center justify-center">
        <span className="text-slate-400">Loading surface data...</span>
      </div>
    );
  }

  if (!plotData) {
    return (
      <div className="absolute inset-0 flex items-center justify-center">
        <p className="text-slate-500 text-sm">
          No surface data available
        </p>
      </div>
    );
  }

  return (
    <div className="w-full h-full relative">
      <Plot
        data={[
          {
            x: plotData.x,
            y: plotData.y,
            z: plotData.z,
            type: "surface",
            colorscale: "Viridis",
            showscale: true,
            colorbar: {
              title: { 
                text: "IV (%)",
                side: "right",
                font: { color: "#94a3b8" }
              },
              tickfont: { color: "#94a3b8" },
            },
          },
        ]}
        layout={{
          autosize: true,
          margin: { l: 40, r: 40, b: 40, t: 40 },
          paper_bgcolor: "transparent",
          plot_bgcolor: "transparent",
          scene: {
            xaxis: {
              title: { 
                text: "Strike",
                font: { color: "#cbd5e1" }
              },
              gridcolor: "#334155",
              zerolinecolor: "#475569",
              tickfont: { color: "#94a3b8" },
            },
            yaxis: {
              title: { 
                text: "Expiration",
                font: { color: "#cbd5e1" }
              },
              gridcolor: "#334155",
              zerolinecolor: "#475569",
              tickfont: { color: "#94a3b8" },
            },
            zaxis: {
              title: { 
                text: "Implied Volatility (%)",
                font: { color: "#cbd5e1" }
              },
              gridcolor: "#334155",
              zerolinecolor: "#475569",
              tickfont: { color: "#94a3b8" },
            },
            camera: {
              eye: { x: 1.5, y: 1.5, z: 1.2 },
            },
          },
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
        config={{ responsive: true, displayModeBar: false }}
      />
    </div>
  );
}
