import type { PricePoint } from "../api/client.ts";
import type { IVSmileData, SmilePoint } from "../types/index.ts";
import { IVSmileChart, LineChart, type SmileChartPoint } from "./Charts.tsx";
import type { TickerInfo } from "./Sidebar.tsx";

function smileExpiries(data: IVSmileData | null): { expiry: string; points: SmilePoint[] }[] {
  if (!data) return [];
  return Object.entries(data.smiles_by_expiry)
    .map(([expiry, points]) => ({ expiry, points }))
    .sort((a, b) => a.expiry.localeCompare(b.expiry));
}

function nearestExpiry(data: IVSmileData | null): SmilePoint[] {
  const ex = smileExpiries(data);
  return ex[0]?.points ?? [];
}

function toSmileChart(points: SmilePoint[]): SmileChartPoint[] {
  return points.map(p => ({ strike: p.strike, iv: p.implied_volatility, oi: p.open_interest }));
}

function ContextRow({
  title, meta, children, onClick, tall,
}: {
  title: string;
  meta: string;
  children: React.ReactNode;
  onClick: () => void;
  tall?: boolean;
}) {
  return (
    <div className="context-row" onClick={onClick}>
      <div className="context-row-head">
        <div className="context-row-title">{title}</div>
        <div className="context-row-meta">{meta}</div>
      </div>
      <div className={`context-row-preview ${tall ? "tall" : ""}`}>{children}</div>
    </div>
  );
}

export type ModalTab = "history" | "smile" | "chain";

export function RightPanel({
  ticker, history, smile, onOpenModal,
}: {
  ticker: TickerInfo | null;
  history: PricePoint[];
  smile: IVSmileData | null;
  onOpenModal: (tab: ModalTab) => void;
}) {
  if (!ticker) {
    return (
      <div className="col right">
        <div className="panel-h"><span>Market Context</span></div>
        <div className="empty">Select a ticker to see market context</div>
      </div>
    );
  }

  const last = history[history.length - 1]?.price;
  const prev = history[history.length - 2]?.price;
  const chg = last && prev ? (last - prev) / prev : 0;
  const atmPoints = nearestExpiry(smile);
  const atmIV = atmPoints.find(p => Math.abs(p.moneyness) < 0.03)?.implied_volatility;
  const expiriesCount = smileExpiries(smile).length;

  return (
    <div className="col right">
      <div className="ticker-hero">
        <div className="sym">{ticker.sym}</div>
        <div className="px tnum">${last != null ? last.toFixed(2) : ticker.spot.toFixed(2)}</div>
        <div className={`chg tnum ${chg >= 0 ? "up" : "dn"}`}>
          {chg >= 0 ? "+" : ""}{(chg * 100).toFixed(2)}%
        </div>
        <div className="meta">{ticker.name}</div>
      </div>

      <ContextRow
        title="Price history"
        meta={history.length ? `${history.length}d` : "—"}
        onClick={() => onOpenModal("history")}
        tall
      >
        {history.length > 1 ? <LineChart data={history.slice(-60)} height={120} yKey="price" /> : <div className="empty">No history</div>}
      </ContextRow>

      <ContextRow
        title="IV Smile"
        meta={`${expiriesCount} expiries · ${atmIV ? `ATM ${(atmIV * 100).toFixed(1)}%` : "—"}`}
        onClick={() => onOpenModal("smile")}
        tall
      >
        {atmPoints.length ? (
          <IVSmileChart points={toSmileChart(atmPoints)} height={120} spot={ticker.spot} width={300} />
        ) : <div className="empty">No smile data</div>}
      </ContextRow>

      <ContextRow
        title="Option Chain"
        meta={`${atmPoints.length} strikes`}
        onClick={() => onOpenModal("chain")}
      >
        {atmPoints.length ? (
          <div style={{ display: "flex", gap: 4, justifyContent: "center", height: 36, alignItems: "center" }}>
            {atmPoints.slice(0, 13).map((p, i) => (
              <div key={i} style={{
                width: 4,
                height: `${Math.min(32, (p.open_interest / 1000) * 0.4 + 4)}px`,
                background: Math.abs(p.moneyness) < 0.03 ? "var(--accent)" : "var(--fg-3)",
                opacity: Math.abs(p.moneyness) < 0.03 ? 1 : 0.5,
                borderRadius: 1,
              }}/>
            ))}
          </div>
        ) : <div className="empty">No chain data</div>}
      </ContextRow>

      <div className="panel">
        <div className="panel-h"><span>Snapshot</span></div>
        <div className="panel-body dense">
          <div className="iv-stat"><span className="k">Spot</span><span className="v tnum">${ticker.spot.toFixed(2)}</span></div>
          <div className="iv-stat"><span className="k">Implied Vol</span><span className="v tnum">{ticker.iv != null ? (ticker.iv * 100).toFixed(2) + "%" : "—"}</span></div>
          <div className="iv-stat"><span className="k">Historical Vol</span><span className="v tnum">{(ticker.hv * 100).toFixed(2)}%</span></div>
          <div className="iv-stat"><span className="k">Div Yield</span><span className="v tnum">{(ticker.div * 100).toFixed(2)}%</span></div>
          {history.length > 0 && (
            <>
              <div className="iv-stat"><span className="k">Period High</span><span className="v tnum">${Math.max(...history.map(p => p.price)).toFixed(2)}</span></div>
              <div className="iv-stat"><span className="k">Period Low</span><span className="v tnum">${Math.min(...history.map(p => p.price)).toFixed(2)}</span></div>
            </>
          )}
        </div>
      </div>
    </div>
  );
}

export function MarketModal({
  open, tab, setTab, onClose, ticker, history, smile,
}: {
  open: boolean;
  tab: ModalTab;
  setTab: (t: ModalTab) => void;
  onClose: () => void;
  ticker: TickerInfo | null;
  history: PricePoint[];
  smile: IVSmileData | null;
}) {
  if (!open) return null;
  const expiries = smileExpiries(smile);
  const firstExp = expiries[0]?.points ?? [];
  const underlying = smile?.underlying_price ?? ticker?.spot;

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <div className="modal-head">
          <h2>Market Context · {ticker?.sym}</h2>
          <button className="icon-btn" onClick={onClose}>Close · Esc</button>
        </div>
        <div className="tabs-row">
          <button aria-pressed={tab === "history"} onClick={() => setTab("history")}>Price History</button>
          <button aria-pressed={tab === "smile"} onClick={() => setTab("smile")}>IV Smile</button>
          <button aria-pressed={tab === "chain"} onClick={() => setTab("chain")}>Option Chain</button>
        </div>
        <div className="modal-body">
          {tab === "history" && (
            <div>
              <div style={{ marginBottom: 16, display: "flex", gap: 24 }}>
                <div>
                  <div style={{ fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase", letterSpacing: "0.06em" }}>Last close</div>
                  <div className="mono" style={{ fontSize: 22, fontWeight: 500 }}>
                    ${history[history.length - 1]?.price.toFixed(2) ?? "—"}
                  </div>
                </div>
                {history.length > 1 && (
                  <div>
                    <div style={{ fontSize: 10, color: "var(--fg-3)", textTransform: "uppercase", letterSpacing: "0.06em" }}>
                      {history.length}d return
                    </div>
                    <div className="mono" style={{
                      fontSize: 22, fontWeight: 500,
                      color: history[history.length - 1].price >= history[0].price ? "var(--accent)" : "var(--loss)",
                    }}>
                      {(((history[history.length - 1].price - history[0].price) / history[0].price) * 100).toFixed(2)}%
                    </div>
                  </div>
                )}
              </div>
              <LineChart data={history} height={360} yKey="price" />
            </div>
          )}
          {tab === "smile" && (
            <div>
              <div style={{ marginBottom: 12, color: "var(--fg-2)", fontSize: 11 }}>
                Implied vol by strike. Point radius ∝ √OI. Dashed line = spot.
              </div>
              {firstExp.length ? (
                <IVSmileChart points={toSmileChart(firstExp)} height={400} spot={underlying} width={1000} />
              ) : <div className="empty">No smile data available</div>}
            </div>
          )}
          {tab === "chain" && (
            <div>
              <div style={{ marginBottom: 8, color: "var(--fg-2)", fontSize: 11 }}>
                Nearest expiry · {ticker?.sym} · {firstExp.length} contracts
              </div>
              {firstExp.length ? (
                <table className="chain-table">
                  <thead>
                    <tr>
                      <th>Strike</th>
                      <th>Side</th>
                      <th>IV</th>
                      <th>Volume</th>
                      <th>OI</th>
                      <th>Moneyness</th>
                    </tr>
                  </thead>
                  <tbody>
                    {firstExp.map((p, i) => {
                      const atm = Math.abs(p.moneyness) < 0.03;
                      return (
                        <tr key={i} className={atm ? "atm" : ""}>
                          <td className="strike">${p.strike}</td>
                          <td className={atm ? "itm" : "otm"}>{p.option_type}</td>
                          <td className={atm ? "itm" : "otm"}>{(p.implied_volatility * 100).toFixed(1)}%</td>
                          <td className={atm ? "itm" : "otm"}>{p.volume.toLocaleString()}</td>
                          <td className={atm ? "itm" : "otm"}>{p.open_interest.toLocaleString()}</td>
                          <td className={atm ? "itm" : "otm"}>{(p.moneyness * 100).toFixed(2)}%</td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              ) : <div className="empty">No chain data available</div>}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
