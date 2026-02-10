import type { ReactNode } from "react";
import { Spinner } from "./Spinner";

interface ButtonProps {
  children: ReactNode;
  onClick: () => void;
  loading?: boolean;
  disabled?: boolean;
  variant?: "primary" | "secondary";
}

const variantClasses = {
  primary:
    "bg-blue-500 text-white hover:bg-blue-600 active:bg-blue-700",
  secondary:
    "bg-slate-700 text-slate-200 hover:bg-slate-600 active:bg-slate-500",
} as const;

export function Button({
  children,
  onClick,
  loading = false,
  disabled = false,
  variant = "primary",
}: ButtonProps) {
  const isDisabled = disabled || loading;

  return (
    <button
      onClick={onClick}
      disabled={isDisabled}
      className={`inline-flex items-center justify-center gap-2 rounded px-4 py-2 text-sm font-medium transition-colors ${
        variantClasses[variant]
      } ${isDisabled ? "cursor-not-allowed opacity-50" : ""}`}
    >
      {loading && <Spinner size="sm" />}
      {children}
    </button>
  );
}
