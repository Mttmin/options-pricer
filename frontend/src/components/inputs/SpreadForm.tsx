import type { SpreadStrikes, SpreadType } from "../../types/index.ts";

export type SpreadFormProps = {
  spreadType: SpreadType;
  onSpreadTypeChange: (type: SpreadType) => void;
  strikes: SpreadStrikes;
  onStrikesChange: (strikes: SpreadStrikes) => void;
  shortTermMaturity?: string;
  longTermMaturity?: string;
  onMaturityChange?: (field: string, value: string) => void;
  showTypeSelector?: boolean;
};

const SPREAD_STRIKE_COUNTS: Record<SpreadType, number> = {
  straddle: 1,
  strangle: 1,
  strip: 1,
  strap: 1,
  synthetic_stock: 1,
  bull_spread_call: 2,
  bull_spread_put: 2,
  bear_spread_call: 2,
  bear_spread_put: 2,
  butterfly: 3,
  calendar_spread: 1,
};

const SPREAD_LABELS: Record<SpreadType, string> = {
  straddle: "Straddle",
  strangle: "Strangle",
  strip: "Strip",
  strap: "Strap",
  synthetic_stock: "Synthetic stock",
  bull_spread_call: "Bull call",
  bull_spread_put: "Bull put",
  bear_spread_call: "Bear call",
  bear_spread_put: "Bear put",
  butterfly: "Butterfly",
  calendar_spread: "Calendar",
};

export function SpreadForm({
  spreadType,
  onSpreadTypeChange,
  strikes,
  onStrikesChange,
  showTypeSelector = true,
}: SpreadFormProps) {
  const count = SPREAD_STRIKE_COUNTS[spreadType] ?? 1;

  const getStrikeKey = (index: number): keyof SpreadStrikes => {
    const keys: (keyof SpreadStrikes)[] = ["strike", "strike_call", "strike_put", "strike_low", "strike_high", "strike_medium"];
    return keys[index] || "strike";
  };

  const handleStrikeChange = (index: number, value: string) => {
    const key = getStrikeKey(index);
    onStrikesChange({ ...strikes, [key]: parseFloat(value) || undefined });
  };

  return (
    <div style={{ fontSize: 12, color: "var(--fg-2)" }}>
      {showTypeSelector && (
        <div style={{ marginBottom: 12 }}>
          <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase", marginBottom: 4 }}>Spread type</label>
          <select
            value={spreadType}
            onChange={(e) => onSpreadTypeChange(e.target.value as SpreadType)}
            style={{ width: "100%", padding: "6px 8px", fontSize: 12, borderRadius: 4 }}
          >
            {(Object.keys(SPREAD_LABELS) as SpreadType[]).map(type => (
              <option key={type} value={type}>{SPREAD_LABELS[type]}</option>
            ))}
          </select>
        </div>
      )}

      {Array.from({ length: count }, (_, i) => (
        <div key={i} style={{ marginBottom: 8 }}>
          <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>
            Strike {i + 1}
          </label>
          <input
            type="number"
            value={strikes[getStrikeKey(i)] ?? ""}
            onChange={(e) => handleStrikeChange(i, e.target.value)}
            style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
          />
        </div>
      ))}
    </div>
  );
}
