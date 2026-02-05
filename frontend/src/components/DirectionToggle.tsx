import type { Direction } from "../types/index.ts";

interface DirectionToggleProps {
  direction: Direction;
  onChange: (d: Direction) => void;
}

export function DirectionToggle({ direction, onChange }: DirectionToggleProps) {
  return (
    <div className="inline-flex rounded border border-slate-700 bg-slate-800">
      <button
        onClick={() => onChange("long")}
        className={`px-5 py-2 text-sm font-medium rounded-l transition-colors ${
          direction === "long"
            ? "bg-emerald-600 text-white"
            : "text-slate-400 hover:text-white"
        }`}
      >
        Long
      </button>
      <button
        onClick={() => onChange("short")}
        className={`px-5 py-2 text-sm font-medium rounded-r transition-colors ${
          direction === "short"
            ? "bg-red-600 text-white"
            : "text-slate-400 hover:text-white"
        }`}
      >
        Short
      </button>
    </div>
  );
}
