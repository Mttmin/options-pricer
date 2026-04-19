import { useEffect, useRef, useState } from "react";
import { searchTickers, type SearchResult } from "../api/client.ts";
import type {
  Direction,
  ExerciseStyle,
  ExoticType,
  OptionType,
  SpreadStrikes,
  SpreadType,
  StructureType,
  VolSource,
} from "../types/index.ts";

export type TickerInfo = {
  sym: string;
  name: string;
  spot: number;
  iv: number | null;
  hv: number;
  div: number;
};

function Seg<T extends string>({
  value, onChange, options, variant = "",
}: {
  value: T;
  onChange: (v: T) => void;
  options: { value: T; label: string }[];
  variant?: string;
}) {
  return (
    <div className={`seg ${variant}`}>
      {options.map(o => (
        <button
          key={o.value}
          aria-pressed={value === o.value}
          data-val={o.value}
          onClick={() => onChange(o.value)}
        >{o.label}</button>
      ))}
    </div>
  );
}

function Field({
  label, hint, children,
}: {
  label: string;
  hint?: string;
  children: React.ReactNode;
}) {
  return (
    <div className="field">
      <label>{label}{hint && <span className="hint">{hint}</span>}</label>
      {children}
    </div>
  );
}

function calcTTMPreview(dateStr: string): number {
  if (!dateStr) return 0;
  const today = new Date(); today.setHours(0, 0, 0, 0);
  const exp = new Date(dateStr);
  return Math.max(0, (exp.getTime() - today.getTime()) / (365 * 24 * 3600 * 1000));
}

function TickerSearchInput({
  value, onSelect,
}: {
  value: TickerInfo | null;
  onSelect: (symbol: string) => void;
}) {
  const [q, setQ] = useState("");
  const [open, setOpen] = useState(false);
  const [focused, setFocused] = useState(0);
  const [results, setResults] = useState<SearchResult[]>([]);
  const ref = useRef<HTMLDivElement>(null);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | undefined>(undefined);

  useEffect(() => {
    const h = (e: MouseEvent) => {
      if (!ref.current?.contains(e.target as Node)) setOpen(false);
    };
    document.addEventListener("mousedown", h);
    return () => document.removeEventListener("mousedown", h);
  }, []);

  useEffect(() => {
    if (debounceRef.current) clearTimeout(debounceRef.current);
    if (!q.trim()) { setResults([]); return; }
    debounceRef.current = setTimeout(() => {
      searchTickers(q).then(setResults).catch(() => setResults([]));
    }, 250);
    return () => { if (debounceRef.current) clearTimeout(debounceRef.current); };
  }, [q]);

  const submit = (sym: string) => {
    onSelect(sym.toUpperCase());
    setQ("");
    setOpen(false);
  };

  return (
    <div className="ticker-search" ref={ref}>
      <div className="input-group">
        <input
          className="input"
          placeholder={value ? value.sym : "Search ticker (AAPL, MSFT...)"}
          value={q}
          onFocus={() => setOpen(true)}
          onChange={(e) => { setQ(e.target.value); setOpen(true); setFocused(0); }}
          onKeyDown={(e) => {
            if (e.key === "ArrowDown") { setFocused(f => Math.min(f + 1, results.length - 1)); e.preventDefault(); }
            if (e.key === "ArrowUp") { setFocused(f => Math.max(f - 1, 0)); e.preventDefault(); }
            if (e.key === "Enter") {
              if (results[focused]) submit(results[focused].symbol);
              else if (q.trim()) submit(q.trim());
            }
            if (e.key === "Escape") setOpen(false);
          }}
        />
        <div className="adornment">↵</div>
      </div>
      {open && results.length > 0 && (
        <div className="ticker-dropdown">
          {results.map((t, i) => (
            <div
              key={t.symbol}
              className="ticker-result"
              aria-selected={i === focused}
              onMouseEnter={() => setFocused(i)}
              onMouseDown={(e) => { e.preventDefault(); submit(t.symbol); }}
            >
              <div className="sym">{t.symbol}</div>
              <div className="name">{t.description}</div>
              <div className="px">{t.type}</div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

export type SidebarProps = {
  ticker: TickerInfo | null;
  onSelectTicker: (symbol: string) => void;
  structure: StructureType;
  setStructure: (s: StructureType) => void;
  direction: Direction;
  setDirection: (d: Direction) => void;
  exerciseStyle: ExerciseStyle;
  setExerciseStyle: (e: ExerciseStyle) => void;
  optionType: OptionType;
  setOptionType: (o: OptionType) => void;
  strike: string;
  setStrike: (s: string) => void;
  expirationDate: string;
  setExpirationDate: (e: string) => void;
  spot: string;
  setSpot: (s: string) => void;
  volatility: string;
  setVolatility: (v: string) => void;
  riskFreeRate: string;
  setRiskFreeRate: (r: string) => void;
  dividendYield: string;
  setDividendYield: (d: string) => void;
  volSource: VolSource;
  setVolSource: (v: VolSource) => void;
  volData: { implied: number | null; historical: number | null; ema: number | null; vix_correlated: number | null } | null;
  manualOverride: boolean;
  setManualOverride: (m: boolean) => void;
  spreadType: SpreadType;
  setSpreadType: (s: SpreadType) => void;
  strikes: SpreadStrikes;
  setStrikes: (s: SpreadStrikes) => void;
  exoticType: ExoticType;
  setExoticType: (e: ExoticType) => void;
  loadingTicker: boolean;
};

export function Sidebar(props: SidebarProps) {
  const {
    ticker, onSelectTicker, structure, setStructure, direction, setDirection,
    exerciseStyle, setExerciseStyle, optionType, setOptionType,
    strike, setStrike, expirationDate, setExpirationDate,
    spot, setSpot, volatility, setVolatility,
    riskFreeRate, setRiskFreeRate, dividendYield, setDividendYield,
    volSource, setVolSource, volData, manualOverride, setManualOverride,
    spreadType, setSpreadType, strikes, setStrikes,
    exoticType, setExoticType, loadingTicker,
  } = props;

  const volOptions: { key: VolSource; name: string; value: number | null }[] = volData ? [
    { key: "implied", name: "Implied", value: volData.implied },
    { key: "historical", name: "Historical", value: volData.historical },
    { key: "ema", name: "EMA (30d)", value: volData.ema },
    { key: "vix_correlated", name: "VIX corr.", value: volData.vix_correlated },
  ] : [];

  const ttm = calcTTMPreview(expirationDate);

  return (
    <div className="col left">
      <div className="panel">
        <div className="panel-h">
          <span>Underlying</span>
          <span className={`status-dot ${ticker ? "live" : ""}`} title={ticker ? "Live data" : "No ticker"}></span>
        </div>
        <div className="panel-body">
          <TickerSearchInput value={ticker} onSelect={onSelectTicker} />
          {loadingTicker && (
            <div style={{ marginTop: 8, fontSize: 10, color: "var(--fg-3)", fontFamily: "var(--font-mono)" }}>
              Loading market data…
            </div>
          )}
        </div>
      </div>

      <div className="panel">
        <div className="panel-h"><span>Structure</span></div>
        <div className="panel-body" style={{ display: "grid", gap: 10 }}>
          <Seg
            value={structure}
            onChange={setStructure}
            options={[
              { value: "single", label: "Single" },
              { value: "spread", label: "Spread" },
              { value: "exotic", label: "Exotic" },
            ]}
          />
          {structure === "single" && (
            <Seg
              value={optionType}
              onChange={setOptionType}
              options={[
                { value: "call", label: "Call" },
                { value: "put", label: "Put" },
              ]}
              variant="accent"
            />
          )}
          {structure !== "exotic" && (
            <Seg
              value={direction}
              onChange={setDirection}
              options={[
                { value: "long", label: "Long" },
                { value: "short", label: "Short" },
              ]}
              variant="direction"
            />
          )}
          {structure !== "exotic" && (
            <Seg
              value={exerciseStyle}
              onChange={setExerciseStyle}
              options={[
                { value: "european", label: "European" },
                { value: "american", label: "American" },
              ]}
            />
          )}
          {structure === "spread" && (
            <select className="select" value={spreadType} onChange={(e) => setSpreadType(e.target.value as SpreadType)}>
              <option value="straddle">Straddle</option>
              <option value="strangle">Strangle</option>
              <option value="strip">Strip</option>
              <option value="strap">Strap</option>
              <option value="synthetic_stock">Synthetic Stock</option>
              <option value="bull_spread_call">Bull Call Spread</option>
              <option value="bull_spread_put">Bull Put Spread</option>
              <option value="bear_spread_call">Bear Call Spread</option>
              <option value="bear_spread_put">Bear Put Spread</option>
              <option value="butterfly">Butterfly</option>
              <option value="calendar_spread">Calendar Spread</option>
            </select>
          )}
          {structure === "exotic" && (
            <select className="select" value={exoticType} onChange={(e) => setExoticType(e.target.value as ExoticType)}>
              <option value="convertible_bond">Convertible Bond</option>
              <option value="chooser_option">Chooser Option</option>
            </select>
          )}
        </div>
      </div>

      {structure === "single" && (
        <div className="panel">
          <div className="panel-h"><span>Contract</span></div>
          <div className="panel-body">
            <div className="fields-grid">
              <Field label="Strike">
                <div className="input-group">
                  <div className="adornment">$</div>
                  <input className="input mono" type="number" value={strike} onChange={(e) => setStrike(e.target.value)} placeholder="230" />
                </div>
              </Field>
              <Field label="Expiry" hint={ttm ? `${(ttm * 365).toFixed(0)}d · T=${ttm.toFixed(3)}` : ""}>
                <input className="input mono" type="date" value={expirationDate} onChange={(e) => setExpirationDate(e.target.value)} />
              </Field>
            </div>
          </div>
        </div>
      )}

      {structure === "spread" && (
        <div className="panel">
          <div className="panel-h"><span>Strikes</span></div>
          <div className="panel-body">
            <div className="fields-grid">
              {(spreadType === "straddle" || spreadType === "strip" || spreadType === "strap" || spreadType === "synthetic_stock") && (
                <Field label="Strike">
                  <input className="input mono" type="number" value={strikes.strike ?? ""} onChange={(e) => setStrikes({ ...strikes, strike: +e.target.value })} />
                </Field>
              )}
              {spreadType === "strangle" && (
                <>
                  <Field label="Strike Put">
                    <input className="input mono" type="number" value={strikes.strike_put ?? ""} onChange={(e) => setStrikes({ ...strikes, strike_put: +e.target.value })} />
                  </Field>
                  <Field label="Strike Call">
                    <input className="input mono" type="number" value={strikes.strike_call ?? ""} onChange={(e) => setStrikes({ ...strikes, strike_call: +e.target.value })} />
                  </Field>
                </>
              )}
              {(spreadType === "bull_spread_call" || spreadType === "bull_spread_put" || spreadType === "bear_spread_call" || spreadType === "bear_spread_put" || spreadType === "calendar_spread") && (
                <>
                  <Field label="Strike Low">
                    <input className="input mono" type="number" value={strikes.strike_low ?? ""} onChange={(e) => setStrikes({ ...strikes, strike_low: +e.target.value })} />
                  </Field>
                  <Field label="Strike High">
                    <input className="input mono" type="number" value={strikes.strike_high ?? ""} onChange={(e) => setStrikes({ ...strikes, strike_high: +e.target.value })} />
                  </Field>
                </>
              )}
              {spreadType === "butterfly" && (
                <>
                  <Field label="Low"><input className="input mono" type="number" value={strikes.strike_low ?? ""} onChange={(e) => setStrikes({ ...strikes, strike_low: +e.target.value })} /></Field>
                  <Field label="Mid"><input className="input mono" type="number" value={strikes.strike_medium ?? ""} onChange={(e) => setStrikes({ ...strikes, strike_medium: +e.target.value })} /></Field>
                  <Field label="High"><input className="input mono" type="number" value={strikes.strike_high ?? ""} onChange={(e) => setStrikes({ ...strikes, strike_high: +e.target.value })} /></Field>
                </>
              )}
              <Field label="Expiry" hint={ttm ? `${(ttm * 365).toFixed(0)}d` : ""}>
                <input className="input mono" type="date" value={expirationDate} onChange={(e) => setExpirationDate(e.target.value)} />
              </Field>
            </div>
          </div>
        </div>
      )}

      <div className="panel">
        <div className="panel-h">
          <span>Volatility</span>
          <label style={{ fontSize: 10, letterSpacing: 0, color: "var(--fg-3)", display: "flex", gap: 4, alignItems: "center", cursor: "pointer", textTransform: "none" }}>
            <input
              type="checkbox"
              checked={manualOverride}
              onChange={(e) => setManualOverride(e.target.checked)}
              style={{ accentColor: "var(--accent)" }}
            />
            manual
          </label>
        </div>
        <div className="panel-body">
          {!manualOverride && volData && (
            <div className="vol-source-list">
              {volOptions.map(o => (
                <div
                  key={o.key}
                  className="vol-source-row"
                  aria-pressed={volSource === o.key}
                  onClick={() => setVolSource(o.key)}
                >
                  <div className="radio" />
                  <div className="name">{o.name}</div>
                  <div className="value">{o.value != null ? (o.value * 100).toFixed(2) + "%" : "—"}</div>
                </div>
              ))}
            </div>
          )}
          {manualOverride && (
            <Field label="Volatility (σ)" hint="decimal, e.g. 0.24">
              <input className="input mono" type="number" step="0.01" value={volatility} onChange={(e) => setVolatility(e.target.value)} />
            </Field>
          )}
          {!manualOverride && !volData && (
            <div className="empty">Select a ticker to see vol sources.</div>
          )}
        </div>
      </div>

      <div className="panel">
        <div className="panel-h"><span>Market</span></div>
        <div className="panel-body">
          <div className="fields-grid">
            <Field label="Spot" hint={ticker ? "live" : ""}>
              <div className="input-group">
                <div className="adornment">$</div>
                <input className="input mono" type="number" value={spot} onChange={(e) => setSpot(e.target.value)} disabled={!manualOverride && !!ticker} />
              </div>
            </Field>
            <Field label="Risk-Free" hint="SOFR">
              <div className="input-group">
                <input className="input mono" type="number" step="0.001" value={riskFreeRate} onChange={(e) => setRiskFreeRate(e.target.value)} />
                <div className="adornment">%</div>
              </div>
            </Field>
            <Field label="Div Yield">
              <div className="input-group">
                <input className="input mono" type="number" step="0.001" value={dividendYield} onChange={(e) => setDividendYield(e.target.value)} />
                <div className="adornment">%</div>
              </div>
            </Field>
          </div>
        </div>
      </div>
    </div>
  );
}
