interface InputProps {
  label: string;
  value: string | number;
  onChange: (value: string) => void;
  type?: string;
  disabled?: boolean;
  placeholder?: string;
}

export function Input({
  label,
  value,
  onChange,
  type = "text",
  disabled = false,
  placeholder,
}: InputProps) {
  return (
    <div className="flex flex-col gap-1">
      <label className="text-xs font-medium uppercase tracking-wider text-slate-400">
        {label}
      </label>
      <div className="relative">
        <input
          type={type}
          value={value}
          onChange={(e) => onChange(e.target.value)}
          disabled={disabled}
          placeholder={placeholder}
          className={`w-full rounded border bg-slate-800 px-3 py-2 text-sm text-white outline-none transition-colors ${
            disabled
              ? "cursor-not-allowed border-slate-700 text-slate-500"
              : "border-slate-600 hover:border-slate-500 focus:border-blue-500"
          }`}
        />
        {disabled && (
          <svg
            className="absolute right-2.5 top-1/2 h-3.5 w-3.5 -translate-y-1/2 text-slate-500"
            fill="none"
            viewBox="0 0 24 24"
            stroke="currentColor"
            strokeWidth={2}
          >
            <rect x="3" y="11" width="18" height="11" rx="2" ry="2" />
            <path d="M7 11V7a5 5 0 0 1 10 0v4" />
          </svg>
        )}
      </div>
    </div>
  );
}
