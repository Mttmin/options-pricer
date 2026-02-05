interface SpinnerProps {
  size?: "sm" | "md" | "lg";
}

const sizeClasses = {
  sm: "h-4 w-4 border-2",
  md: "h-6 w-6 border-2",
  lg: "h-8 w-8 border-3",
} as const;

export function Spinner({ size = "md" }: SpinnerProps) {
  return (
    <div
      className={`${sizeClasses[size]} animate-spin rounded-full border-slate-600 border-t-blue-500`}
      role="status"
      aria-label="Loading"
    />
  );
}
