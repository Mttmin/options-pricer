import type { ExerciseStyle } from "../types/index.ts";

interface ExerciseStyleToggleProps {
  style: ExerciseStyle;
  onChange: (style: ExerciseStyle) => void;
}

export function ExerciseStyleToggle({ style, onChange }: ExerciseStyleToggleProps) {
  return (
    <div className="inline-flex rounded border border-slate-700 bg-slate-800">
      <button
        onClick={() => onChange("european")}
        className={`px-5 py-2 text-sm font-medium rounded-l transition-colors ${
          style === "european"
            ? "bg-blue-600 text-white"
            : "text-slate-400 hover:text-white"
        }`}
      >
        European
      </button>
      <button
        onClick={() => onChange("american")}
        className={`px-5 py-2 text-sm font-medium rounded-r transition-colors ${
          style === "american"
            ? "bg-purple-600 text-white"
            : "text-slate-400 hover:text-white"
        }`}
      >
        American
      </button>
    </div>
  );
}
