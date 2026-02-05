import type { ExoticType } from "../../types/index.ts";
import { Input } from "../ui/Input.tsx";
import { Select } from "../ui/Select.tsx";

interface ExoticFormProps {
  exoticType: ExoticType;
  onExoticTypeChange: (type: ExoticType) => void;
  values: Record<string, string>;
  onChange: (field: string, value: string) => void;
}

const EXOTIC_OPTIONS = [
  { value: "convertible_bond", label: "Convertible Bond" },
  { value: "chooser_option", label: "Chooser Option" },
];

export function ExoticForm({
  exoticType,
  onExoticTypeChange,
  values,
  onChange,
}: ExoticFormProps) {
  return (
    <div className="space-y-4">
      <Select
        label="Exotic Type"
        value={exoticType}
        onChange={(v) => onExoticTypeChange(v as ExoticType)}
        options={EXOTIC_OPTIONS}
      />

      {exoticType === "convertible_bond" ? (
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-3">
          <Input
            label="Face Value"
            type="number"
            value={values.face_value ?? ""}
            onChange={(v) => onChange("face_value", v)}
            placeholder="1000"
          />
          <Input
            label="Coupon Rate"
            type="number"
            value={values.coupon_rate ?? ""}
            onChange={(v) => onChange("coupon_rate", v)}
            placeholder="0.05"
          />
          <Input
            label="Maturity (years)"
            type="number"
            value={values.maturity ?? ""}
            onChange={(v) => onChange("maturity", v)}
            placeholder="5"
          />
          <Input
            label="Payment Frequency"
            type="number"
            value={values.payment_frequency ?? ""}
            onChange={(v) => onChange("payment_frequency", v)}
            placeholder="2"
          />
          <Input
            label="Credit Spread"
            type="number"
            value={values.credit_spread ?? ""}
            onChange={(v) => onChange("credit_spread", v)}
            placeholder="0.02"
          />
          <Input
            label="Conversion Price"
            type="number"
            value={values.conversion_price ?? ""}
            onChange={(v) => onChange("conversion_price", v)}
            placeholder="50"
          />
          <Input
            label="Stock Price"
            type="number"
            value={values.stock_price ?? ""}
            onChange={(v) => onChange("stock_price", v)}
            placeholder="Auto-filled"
          />
          <Input
            label="Time to Maturity"
            type="number"
            value={values.time_to_maturity ?? ""}
            onChange={(v) => onChange("time_to_maturity", v)}
            placeholder="5"
          />
          <Input
            label="Risk-Free Rate"
            type="number"
            value={values.risk_free_rate ?? ""}
            onChange={(v) => onChange("risk_free_rate", v)}
            placeholder="0.05"
          />
          <Input
            label="Dividend Yield"
            type="number"
            value={values.dividend_yield ?? ""}
            onChange={(v) => onChange("dividend_yield", v)}
            placeholder="0.00"
          />
        </div>
      ) : (
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-3">
          <Input
            label="Spot Price"
            type="number"
            value={values.spot_price ?? ""}
            onChange={(v) => onChange("spot_price", v)}
            placeholder="Auto-filled"
          />
          <Input
            label="Strike Price"
            type="number"
            value={values.strike_price ?? ""}
            onChange={(v) => onChange("strike_price", v)}
            placeholder="100"
          />
          <Input
            label="Time to Maturity"
            type="number"
            value={values.time_to_maturity ?? ""}
            onChange={(v) => onChange("time_to_maturity", v)}
            placeholder="1"
          />
          <Input
            label="Choose Time"
            type="number"
            value={values.choose_time ?? ""}
            onChange={(v) => onChange("choose_time", v)}
            placeholder="0.5"
          />
          <Input
            label="Risk-Free Rate"
            type="number"
            value={values.risk_free_rate ?? ""}
            onChange={(v) => onChange("risk_free_rate", v)}
            placeholder="0.05"
          />
          <Input
            label="Dividend Yield"
            type="number"
            value={values.dividend_yield ?? ""}
            onChange={(v) => onChange("dividend_yield", v)}
            placeholder="0.00"
          />
        </div>
      )}
    </div>
  );
}
