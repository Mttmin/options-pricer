import type { ExoticType } from "../../types/index.ts";

export type ExoticFormProps = {
  exoticType: ExoticType;
  onExoticTypeChange: (type: ExoticType) => void;
  values: Record<string, string>;
  onChange: (field: string, value: string) => void;
  manualOverride?: boolean;
};

export function ExoticForm({
  exoticType,
  onExoticTypeChange,
  values,
  onChange,
}: ExoticFormProps) {
  return (
    <div style={{ fontSize: 12, color: "var(--fg-2)" }}>
      <div style={{ marginBottom: 12 }}>
        <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase", marginBottom: 4 }}>Type</label>
        <select
          value={exoticType}
          onChange={(e) => onExoticTypeChange(e.target.value as ExoticType)}
          style={{ width: "100%", padding: "6px 8px", fontSize: 12, borderRadius: 4 }}
        >
          <option value="convertible_bond">Convertible bond</option>
          <option value="chooser_option">Chooser option</option>
        </select>
      </div>

      {exoticType === "convertible_bond" && (
        <div style={{ marginBottom: 8 }}>
          <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Coupon rate (%)</label>
          <input
            type="number"
            value={values.coupon_rate || ""}
            onChange={(e) => onChange("coupon_rate", e.target.value)}
            step="0.01"
            style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
          />
        </div>
      )}

      {exoticType === "chooser_option" && (
        <div style={{ marginBottom: 8 }}>
          <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Choice date</label>
          <input
            type="date"
            value={values.choice_date || ""}
            onChange={(e) => onChange("choice_date", e.target.value)}
            style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
          />
        </div>
      )}

      <div style={{ marginBottom: 8 }}>
        <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Risk-free rate (%)</label>
        <input
          type="number"
          value={values.risk_free_rate || ""}
          onChange={(e) => onChange("risk_free_rate", e.target.value)}
          step="0.01"
          style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
        />
      </div>
    </div>
  );
}
