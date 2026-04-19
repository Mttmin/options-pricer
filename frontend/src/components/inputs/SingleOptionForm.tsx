export type SingleOptionFormProps = {
  values: {
    optionType: "call" | "put";
    strike: string;
    expirationDate: string;
    spotPrice: string;
    volatility: string;
    riskFreeRate: string;
    dividendYield: string;
  };
  onChange: (field: string, value: string) => void;
  manualOverride: boolean;
  isSpread?: boolean;
};

export function SingleOptionForm({ values, onChange }: SingleOptionFormProps) {
  return (
    <div style={{ fontSize: 12, color: "var(--fg-2)" }}>
      <div style={{ marginBottom: 8 }}>
        <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Strike</label>
        <input
          type="number"
          value={values.strike}
          onChange={(e) => onChange("strike", e.target.value)}
          style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
        />
      </div>
      <div style={{ marginBottom: 8 }}>
        <label style={{ display: "block", fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase" }}>Expiration</label>
        <input
          type="date"
          value={values.expirationDate}
          onChange={(e) => onChange("expirationDate", e.target.value)}
          style={{ width: "100%", padding: "4px 8px", fontSize: 12, borderRadius: 4 }}
        />
      </div>
    </div>
  );
}
