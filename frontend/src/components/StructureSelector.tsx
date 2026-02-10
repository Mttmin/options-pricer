import type { StructureType } from "../types/index.ts";

interface StructureSelectorProps {
  selected: StructureType;
  onChange: (type: StructureType) => void;
}

const TABS: { value: StructureType; label: string }[] = [
  { value: "single", label: "Single" },
  { value: "spread", label: "Spread" },
  { value: "exotic", label: "Exotic" },
];

export function StructureSelector({ selected, onChange }: StructureSelectorProps) {
  return (
    <div className="inline-flex rounded border border-slate-700 bg-slate-800">
      {TABS.map((tab) => (
        <button
          key={tab.value}
          onClick={() => onChange(tab.value)}
          className={`px-5 py-2 text-sm font-medium transition-colors ${
            selected === tab.value
              ? "bg-blue-500 text-white"
              : "text-slate-400 hover:text-white"
          } ${tab.value === "single" ? "rounded-l" : ""} ${
            tab.value === "exotic" ? "rounded-r" : ""
          }`}
        >
          {tab.label}
        </button>
      ))}
    </div>
  );
}
