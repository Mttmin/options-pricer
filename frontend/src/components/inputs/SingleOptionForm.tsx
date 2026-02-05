import type { OptionType } from "../../types/index.ts";
import { Input } from "../ui/Input.tsx";
import { Select } from "../ui/Select.tsx";

interface SingleOptionValues {
  optionType: OptionType;
  strike: string;
  expirationDate: string;
  spotPrice: string;
  volatility: string;
  riskFreeRate: string;
  dividendYield: string;
}

interface SingleOptionFormProps {
  values: SingleOptionValues;
  onChange: (field: string, value: string) => void;
  manualOverride: boolean;
}

const OPTION_TYPE_OPTIONS = [
  { value: "call", label: "Call" },
  { value: "put", label: "Put" },
];

export function SingleOptionForm({
  values,
  onChange,
  manualOverride,
}: SingleOptionFormProps) {
  return (
    <div className="grid grid-cols-2 gap-4 sm:grid-cols-3">
      <Select
        label="Option Type"
        value={values.optionType}
        onChange={(v) => onChange("optionType", v)}
        options={OPTION_TYPE_OPTIONS}
      />
      <Input
        label="Strike Price"
        type="number"
        value={values.strike}
        onChange={(v) => onChange("strike", v)}
        placeholder="100.00"
      />
      <Input
        label="Expiration Date"
        type="date"
        value={values.expirationDate}
        onChange={(v) => onChange("expirationDate", v)}
      />
      {manualOverride && (
        <>
          <Input
            label="Spot Price"
            type="number"
            value={values.spotPrice}
            onChange={(v) => onChange("spotPrice", v)}
            placeholder="100.00"
          />
          <Input
            label="Volatility"
            type="number"
            value={values.volatility}
            onChange={(v) => onChange("volatility", v)}
            placeholder="0.20"
          />
        </>
      )}
      <Input
        label="Risk-Free Rate"
        type="number"
        value={values.riskFreeRate}
        onChange={(v) => onChange("riskFreeRate", v)}
        placeholder="0.05"
      />
      <Input
        label="Dividend Yield"
        type="number"
        value={values.dividendYield}
        onChange={(v) => onChange("dividendYield", v)}
        placeholder="0.00"
      />
    </div>
  );
}
