import { useState, useMemo } from "react";
import type { IVSmileData } from "../types/index.ts";

interface OptionChainDisplayProps {
  data: IVSmileData;
  underlyingPrice?: number | null;
}

export function OptionChainDisplay({ data, underlyingPrice }: OptionChainDisplayProps) {
  const expirations = useMemo(() => {
    return Object.keys(data.smiles_by_expiry).sort();
  }, [data]);

  const [selectedExpiration, setSelectedExpiration] = useState<string>("");

  const currentExpiration = selectedExpiration || expirations[0] || "";
  const smileData = data.smiles_by_expiry[currentExpiration] || [];

  const calls = useMemo(
    () => smileData.filter((p) => p.option_type === "Call").sort((a, b) => a.strike - b.strike),
    [smileData]
  );

  const puts = useMemo(
    () => smileData.filter((p) => p.option_type === "Put").sort((a, b) => a.strike - b.strike),
    [smileData]
  );

  if (expirations.length === 0) {
    return (
      <div className="rounded-lg border border-slate-800 bg-slate-900 p-6">
        <p className="text-slate-500 text-sm text-center">
          No option chain data available
        </p>
      </div>
    );
  }

  const spot = underlyingPrice ?? data.underlying_price;

  const isATM = (strike: number) => {
    if (!spot) return false;
    const diff = Math.abs(strike - spot);
    const pct = diff / spot;
    return pct < 0.02;
  };

  return (
    <div className="rounded-lg border border-slate-800 bg-slate-900 p-4">
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-lg font-semibold text-white">Available Options</h3>
        <div className="flex gap-2 overflow-x-auto pb-2 scrollbar-thin scrollbar-thumb-slate-700 scrollbar-track-transparent">
          {expirations.map((exp) => (
            <button
              key={exp}
              onClick={() => setSelectedExpiration(exp)}
              className={`px-3 py-1 text-sm rounded-md transition-colors whitespace-nowrap flex-shrink-0 ${exp === currentExpiration
                  ? "bg-blue-600 text-white"
                  : "bg-slate-800 text-slate-400 hover:bg-slate-700"
                }`}
            >
              {exp}
            </button>
          ))}
        </div>
      </div>

      <div className="grid grid-cols-2 gap-4 overflow-y-auto max-h-[600px]">
        {/* Calls Table */}
        <div>
          <h4 className="text-sm font-semibold text-slate-400 mb-2 uppercase sticky top-0 bg-slate-900 z-10 py-1">
            Calls
          </h4>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead className="sticky top-0 bg-slate-900 z-10">
                <tr className="border-b border-slate-800">
                  <th className="text-left py-2 px-2 text-slate-500 font-medium">
                    Strike
                  </th>
                  <th className="text-right py-2 px-2 text-slate-500 font-medium">
                    IV
                  </th>
                  <th className="text-right py-2 px-2 text-slate-500 font-medium">
                    OI
                  </th>
                  <th className="text-right py-2 px-2 text-slate-500 font-medium">
                    Vol
                  </th>
                </tr>
              </thead>
              <tbody>
                {calls.map((call, idx) => (
                  <tr
                    key={idx}
                    className={`border-b border-slate-800 hover:bg-slate-800 transition-colors ${isATM(call.strike) ? "bg-blue-950/30" : ""
                      }`}
                  >
                    <td className="py-2 px-2 text-white font-mono">
                      ${call.strike.toFixed(2)}
                      {isATM(call.strike) && (
                        <span className="ml-1 text-xs text-blue-400">ATM</span>
                      )}
                    </td>
                    <td className="py-2 px-2 text-right text-white font-mono">
                      {(call.implied_volatility * 100).toFixed(1)}%
                    </td>
                    <td className="py-2 px-2 text-right text-slate-300 font-mono">
                      {call.open_interest.toLocaleString()}
                    </td>
                    <td className="py-2 px-2 text-right text-slate-300 font-mono">
                      {call.volume.toLocaleString()}
                    </td>
                  </tr>
                ))}
                {calls.length === 0 && (
                  <tr>
                    <td colSpan={4} className="py-4 text-center text-slate-500">
                      No call options available
                    </td>
                  </tr>
                )}
              </tbody>
            </table>
          </div>
        </div>

        {/* Puts Table */}
        <div>
          <h4 className="text-sm font-semibold text-slate-400 mb-2 uppercase sticky top-0 bg-slate-900 z-10 py-1">
            Puts
          </h4>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead className="sticky top-0 bg-slate-900 z-10">
                <tr className="border-b border-slate-800">
                  <th className="text-left py-2 px-2 text-slate-500 font-medium">
                    Strike
                  </th>
                  <th className="text-right py-2 px-2 text-slate-500 font-medium">
                    IV
                  </th>
                  <th className="text-right py-2 px-2 text-slate-500 font-medium">
                    OI
                  </th>
                  <th className="text-right py-2 px-2 text-slate-500 font-medium">
                    Vol
                  </th>
                </tr>
              </thead>
              <tbody>
                {puts.map((put, idx) => (
                  <tr
                    key={idx}
                    className={`border-b border-slate-800 hover:bg-slate-800 transition-colors ${isATM(put.strike) ? "bg-red-950/30" : ""
                      }`}
                  >
                    <td className="py-2 px-2 text-white font-mono">
                      ${put.strike.toFixed(2)}
                      {isATM(put.strike) && (
                        <span className="ml-1 text-xs text-red-400">ATM</span>
                      )}
                    </td>
                    <td className="py-2 px-2 text-right text-white font-mono">
                      {(put.implied_volatility * 100).toFixed(1)}%
                    </td>
                    <td className="py-2 px-2 text-right text-slate-300 font-mono">
                      {put.open_interest.toLocaleString()}
                    </td>
                    <td className="py-2 px-2 text-right text-slate-300 font-mono">
                      {put.volume.toLocaleString()}
                    </td>
                  </tr>
                ))}
                {puts.length === 0 && (
                  <tr>
                    <td colSpan={4} className="py-4 text-center text-slate-500">
                      No put options available
                    </td>
                  </tr>
                )}
              </tbody>
            </table>
          </div>
        </div>
      </div>

      <div className="mt-4 text-xs text-slate-500 text-center">
        Data as of {data.data_date} • Underlying: {spot ? `$${spot.toFixed(2)}` : "-"}
      </div>
    </div>
  );
}
