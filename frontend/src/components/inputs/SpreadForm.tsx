import type { SpreadType, SpreadStrikes } from "../../types/index.ts";
import { SPREAD_LABELS } from "../../types/index.ts";
import { Input } from "../ui/Input.tsx";
import { Select } from "../ui/Select.tsx";

interface SpreadFormProps {
  spreadType: SpreadType;
  onSpreadTypeChange: (type: SpreadType) => void;
  strikes: SpreadStrikes;
  onStrikesChange: (strikes: SpreadStrikes) => void;
  shortTermMaturity: string;
  longTermMaturity: string;
  onMaturityChange: (field: string, value: string) => void;
}

const SPREAD_OPTIONS = Object.entries(SPREAD_LABELS).map(([value, label]) => ({
  value,
  label,
}));

const STRIKE_FIELDS: Record<SpreadType, (keyof SpreadStrikes)[]> = {
  straddle: ["strike"],
  strangle: ["strike_call", "strike_put"],
  strip: ["strike"],
  strap: ["strike"],
  synthetic_stock: ["strike"],
  bull_spread_call: ["strike_low", "strike_high"],
  bull_spread_put: ["strike_low", "strike_high"],
  bear_spread_call: ["strike_low", "strike_high"],
  bear_spread_put: ["strike_low", "strike_high"],
  butterfly: ["strike_low", "strike_medium", "strike_high"],
  calendar_spread: ["strike"],
};

const STRIKE_LABELS: Record<keyof SpreadStrikes, string> = {
  strike: "Strike",
  strike_call: "Call Strike",
  strike_put: "Put Strike",
  strike_low: "Lower Strike",
  strike_high: "Upper Strike",
  strike_medium: "Middle Strike",
};

export function SpreadForm({
  spreadType,
  onSpreadTypeChange,
  strikes,
  onStrikesChange,
  shortTermMaturity,
  longTermMaturity,
  onMaturityChange,
}: SpreadFormProps) {
  const strikeFields = STRIKE_FIELDS[spreadType];
  const needsMaturity = spreadType === "calendar_spread";

  function handleStrikeChange(field: keyof SpreadStrikes, value: string) {
    const num = value === "" ? undefined : parseFloat(value);
    onStrikesChange({ ...strikes, [field]: num });
  }

  return (
    <div className="space-y-4">
      <Select
        label="Spread Type"
        value={spreadType}
        onChange={(v) => onSpreadTypeChange(v as SpreadType)}
        options={SPREAD_OPTIONS}
      />

      <div className="grid grid-cols-2 gap-4">
        {strikeFields.map((field) => (
          <Input
            key={field}
            label={STRIKE_LABELS[field]}
            type="number"
            value={strikes[field] ?? ""}
            onChange={(v) => handleStrikeChange(field, v)}
            placeholder="100.00"
          />
        ))}
      </div>

      {needsMaturity && (
        <div className="grid grid-cols-2 gap-4">
          <Input
            label="Short-Term TTM"
            type="number"
            value={shortTermMaturity}
            onChange={(v) => onMaturityChange("shortTermMaturity", v)}
            placeholder="0.25"
          />
          <Input
            label="Long-Term TTM"
            type="number"
            value={longTermMaturity}
            onChange={(v) => onMaturityChange("longTermMaturity", v)}
            placeholder="0.5"
          />
        </div>
      )}
    </div>
  );
}
