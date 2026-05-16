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
          <option value="asian_option">Asian option</option>
        </select>
      </div>

      {exoticType === "convertible_bond" && (
        <>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Coupon rate (%)</label>
            <input
              type="number"
              value={values.coupon_rate ?? ""}
              onChange={(e) => onChange("coupon_rate", e.target.value)}
              step="0.01"
              placeholder="5.0"
              style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
            />
          </div>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Face value</label>
            <input
              type="number"
              value={values.face_value ?? ""}
              onChange={(e) => onChange("face_value", e.target.value)}
              step="1"
              placeholder="1000"
              style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
            />
          </div>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Maturity (yrs)</label>
            <input
              type="number"
              value={values.maturity ?? ""}
              onChange={(e) => onChange("maturity", e.target.value)}
              step="0.25"
              placeholder="5"
              style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
            />
          </div>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Payments / yr</label>
            <input
              type="number"
              value={values.payment_frequency ?? ""}
              onChange={(e) => onChange("payment_frequency", e.target.value)}
              step="1"
              min="1"
              placeholder="2"
              style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
            />
          </div>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Credit spread (%)</label>
            <input
              type="number"
              value={values.credit_spread ?? ""}
              onChange={(e) => onChange("credit_spread", e.target.value)}
              step="0.01"
              placeholder="2.0"
              style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
            />
          </div>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Conversion price</label>
            <input
              type="number"
              value={values.conversion_price ?? ""}
              onChange={(e) => onChange("conversion_price", e.target.value)}
              step="0.5"
              placeholder="50"
              style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
            />
          </div>
        </>
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

      {exoticType === "asian_option" && (
        <>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase", marginBottom: 4 }}>Pooling</label>
            <select
              value={values.pooling || "average"}
              onChange={(e) => onChange("pooling", e.target.value)}
              style={{ width: "100%", padding: "6px 8px", fontSize: 12, borderRadius: 4 }}
            >
              <option value="average">Average (Asian)</option>
              <option value="max">Max (Lookback)</option>
            </select>
          </div>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase", marginBottom: 4 }}>Option type</label>
            <select
              value={values.asian_option_type || "call"}
              onChange={(e) => onChange("asian_option_type", e.target.value)}
              style={{ width: "100%", padding: "6px 8px", fontSize: 12, borderRadius: 4 }}
            >
              <option value="call">Call</option>
              <option value="put">Put</option>
            </select>
          </div>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Strike</label>
            <input
              type="number"
              value={values.strike || ""}
              onChange={(e) => onChange("strike", e.target.value)}
              step="0.01"
              placeholder="defaults to spot"
              style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
            />
          </div>
          <div style={{ marginBottom: 8 }}>
            <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Monitoring observations</label>
            <input
              type="number"
              value={values.num_observations || "50"}
              onChange={(e) => onChange("num_observations", e.target.value)}
              step="1"
              min="1"
              style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
            />
          </div>
        </>
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
