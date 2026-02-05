import { Toggle } from "./ui/Toggle.tsx";

interface ManualOverrideToggleProps {
  enabled: boolean;
  onChange: (enabled: boolean) => void;
}

export function ManualOverrideToggle({ enabled, onChange }: ManualOverrideToggleProps) {
  return <Toggle label="Manual Override" checked={enabled} onChange={onChange} />;
}
