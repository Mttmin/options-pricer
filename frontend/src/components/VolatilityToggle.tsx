import type { VolSource, VolatilityResponse } from "../types/index.ts";

interface VolatilityToggleProps {
  source: VolSource;
  onChange: (source: VolSource) => void;
  volData: VolatilityResponse | null;
}

const VOL_SOURCES: { value: VolSource; label: string; shortLabel: string }[] = [
  { value: "ema", label: "EMA Vol", shortLabel: "EMA" },
  { value: "implied", label: "Implied Vol", shortLabel: "IV" },
  { value: "historical", label: "Historical", shortLabel: "Hist" },
  { value: "vix_correlated", label: "VIX Correlated", shortLabel: "VIX" },
];

function formatVol(val: number | null | undefined): string {
  if (val == null) return "-";
  return (val * 100).toFixed(1) + "%";
}

export function VolatilityToggle({
  source,
  onChange,
  volData,
}: VolatilityToggleProps) {
  return (
    <div className="flex items-center gap-2">
      {VOL_SOURCES.map((opt) => {
        const isActive = source === opt.value;
        const vol =
          opt.value === "implied"
            ? volData?.implied
            : opt.value === "historical"
              ? volData?.historical
              : opt.value === "ema"
                ? volData?.ema
                : volData?.vix_correlated;
        const isAvailable = vol != null;
        const showVixFallback =
          opt.value === "vix_correlated" && volData?.vix_fallback_used;

        return (
          <button
            key={opt.value}
            onClick={() => isAvailable && onChange(opt.value)}
            disabled={!isAvailable}
            className={`rounded px-3 py-2 text-sm font-medium transition-colors flex flex-col items-center min-w-[70px] ${
              isActive
                ? "bg-blue-500 text-white"
                : isAvailable
                  ? "bg-slate-800 text-slate-300 hover:bg-slate-700"
                  : "cursor-not-allowed bg-slate-900 text-slate-600"
            }`}
            title={showVixFallback ? `${opt.label} (fallback: EMA)` : opt.label}
          >
            <span className="text-xs opacity-80">{opt.shortLabel}</span>
            <span className="font-semibold">{formatVol(vol)}</span>
            {showVixFallback && (
              <span className="mt-0.5 text-[10px] uppercase tracking-wide text-amber-300">
                EMA fallback
              </span>
            )}
          </button>
        );
      })}
    </div>
  );
}
