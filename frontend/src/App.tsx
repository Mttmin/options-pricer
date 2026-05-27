import { useCallback, useEffect, useMemo, useState } from "react";
import type {
  CtmcPriceResponse,
  Direction,
  ExerciseStyle,
  ExoticType,
  IVSmileData,
  MarketDataResponse,
  OptionType,
  PriceRequest,
  PriceResponse,
  SpreadStrikes,
  SpreadType,
  StructureType,
  VolatilityResponse,
  VolSource,
} from "./types/index.ts";
import {
  fetchMarketData,
  fetchOptionChain,
  fetchPriceHistory,
  fetchSofr,
  fetchVolatility,
  submitCtmcPricing,
  submitPricing,
  type PricePoint,
} from "./api/client.ts";
import { calculateTTM } from "./utils/ttm.ts";
import { Sidebar, type TickerInfo } from "./design/Sidebar.tsx";
import {
  buildRows,
  GreeksStrip,
  HeroPrice,
  PayoffPanel,
  PriceChartPanel,
  PricingRows,
} from "./design/Center.tsx";
import { PayoffChart } from "./design/Charts.tsx";
import { MarketModal, RightPanel, type ModalTab } from "./design/MarketContext.tsx";
import { TweaksPanel, TWEAK_DEFAULTS, type Tweaks } from "./design/TweaksPanel.tsx";
import { SpreadForm } from "./components/inputs/SpreadForm.tsx";
import { calculatePayoffGrid } from "./utils/payoff.ts";

function todayPlus(days: number): string {
  const d = new Date();
  d.setDate(d.getDate() + days);
  return d.toISOString().slice(0, 10);
}

export default function App() {
  const [tweaks, setTweaks] = useState<Tweaks>(() => {
    try { return { ...TWEAK_DEFAULTS, ...JSON.parse(localStorage.getItem("op-tweaks") || "{}") }; }
    catch { return TWEAK_DEFAULTS; }
  });
  useEffect(() => { localStorage.setItem("op-tweaks", JSON.stringify(tweaks)); }, [tweaks]);
  useEffect(() => { document.documentElement.setAttribute("data-theme", tweaks.theme); }, [tweaks.theme]);

  const [selectedSymbol, setSelectedSymbol] = useState<string | null>(null);
  const [marketData, setMarketData] = useState<MarketDataResponse | null>(null);
  const [volData, setVolData] = useState<VolatilityResponse | null>(null);
  const [priceHistory, setPriceHistory] = useState<PricePoint[]>([]);
  const [smile, setSmile] = useState<IVSmileData | null>(null);
  const [marketLoading, setMarketLoading] = useState(false);
  const [marketError, setMarketError] = useState<string | null>(null);

  const [structure, setStructure] = useState<StructureType>("single");
  const [direction, setDirection] = useState<Direction>("long");
  const [exerciseStyle, setExerciseStyle] = useState<ExerciseStyle>("european");
  const [optionType, setOptionType] = useState<OptionType>("call");
  const [strike, setStrike] = useState("100");
  const [expirationDate, setExpirationDate] = useState(todayPlus(60));
  const [spot, setSpot] = useState("100");
  const [volatility, setVolatility] = useState("0.24");
  const [riskFreeRate, setRiskFreeRate] = useState("4.75");
  const [dividendYield, setDividendYield] = useState("0");
  const [volSource, setVolSource] = useState<VolSource>("implied");
  const [manualOverride, setManualOverride] = useState(false);
  const [spreadType, setSpreadType] = useState<SpreadType>("straddle");
  const [strikes, setStrikes] = useState<SpreadStrikes>({});
  const [exoticType, setExoticType] = useState<ExoticType>("convertible_bond");
  const [exoticValues, setExoticValues] = useState<Record<string, string>>({
    risk_free_rate: "0.05",
  });

  const [range, setRange] = useState<"1M" | "3M" | "6M" | "1Y">("3M");
  const [modalOpen, setModalOpen] = useState(false);
  const [modalTab, setModalTab] = useState<ModalTab>("history");
  const [pinned, setPinned] = useState<{ price: number; method: string } | null>(null);
  const [selectedMethod, setSelectedMethod] = useState<string | null>(null);

  const [result, setResult] = useState<PriceResponse | null>(null);
  const [priceError, setPriceError] = useState<string | null>(null);
  const [ctmcResult, setCtmcResult] = useState<CtmcPriceResponse | null>(null);
  const [ctmcLoading, setCtmcLoading] = useState(false);

  useEffect(() => {
    fetchSofr()
      .then(r => setRiskFreeRate((r.rate * 100).toFixed(3)))
      .catch(() => {});
  }, []);

  useEffect(() => {
    const h = (e: KeyboardEvent) => { if (e.key === "Escape") setModalOpen(false); };
    document.addEventListener("keydown", h);
    return () => document.removeEventListener("keydown", h);
  }, []);

  const loadSymbol = useCallback(async (sym: string) => {
    setMarketLoading(true);
    setMarketError(null);
    setSelectedSymbol(sym);
    try {
      const [data, history, chain] = await Promise.all([
        fetchMarketData(sym),
        fetchPriceHistory(sym, 365).catch(() => [] as PricePoint[]),
        fetchOptionChain(sym).catch(() => null),
      ]);
      setMarketData(data);
      setPriceHistory(history);
      setSmile(chain);
      setSpot(data.spot_price.toFixed(2));
      const volToUse = data.implied_volatility ?? data.historical_volatility;
      setVolatility(volToUse.toFixed(4));
      setDividendYield(data.dividend_yield ? (data.dividend_yield * 100).toFixed(2) : "0");
      setStrike(String(Math.round(data.spot_price)));
      setStrikes({ strike: Math.round(data.spot_price) });
      fetchVolatility(sym).then(setVolData).catch(() => setVolData(null));
    } catch (err) {
      setMarketError(err instanceof Error ? err.message : "Failed to fetch market data");
    } finally {
      setMarketLoading(false);
    }
  }, []);

  const ticker: TickerInfo | null = useMemo(() => {
    if (!marketData) return null;
    return {
      sym: marketData.symbol,
      name: marketData.symbol,
      spot: marketData.spot_price,
      iv: marketData.implied_volatility,
      hv: marketData.historical_volatility,
      div: marketData.dividend_yield ?? 0,
    };
  }, [marketData]);

  const volDataForSidebar = useMemo(() => {
    if (!volData && !marketData) return null;
    return {
      implied: volData?.implied ?? marketData?.implied_volatility ?? null,
      historical: volData?.historical ?? marketData?.historical_volatility ?? null,
      ema: volData?.ema ?? null,
      vix_correlated: volData?.vix_correlated ?? null,
    };
  }, [volData, marketData]);

  const effectiveVol = useMemo((): number => {
    if (manualOverride) return parseFloat(volatility) || 0.2;
    if (!volDataForSidebar) return parseFloat(volatility) || 0.2;
    const v = volDataForSidebar[volSource];
    return v ?? volDataForSidebar.historical ?? 0.2;
  }, [manualOverride, volatility, volDataForSidebar, volSource]);

  const bsArgs = useMemo(() => {
    const S = parseFloat(spot) || marketData?.spot_price || 100;
    const K = parseFloat(strike) || S;
    const r = (parseFloat(riskFreeRate) || 5) / 100;
    const q = (parseFloat(dividendYield) || 0) / 100;
    const T = calculateTTM(expirationDate) || 0.25;
    return { S, K, r, q, sigma: effectiveVol, T, type: optionType };
  }, [spot, strike, riskFreeRate, dividendYield, expirationDate, effectiveVol, optionType, marketData]);

  const buildRequest = useCallback((): PriceRequest | null => {
    try {
      const vol = effectiveVol;
      const S = bsArgs.S;
      const r = bsArgs.r;
      const q = bsArgs.q;
      const T = bsArgs.T;

      if (structure === "single") {
        const K = parseFloat(strike);
        if (!S || !K || T <= 0) return null;
        return {
          structure_type: "single",
          option_type: optionType,
          direction,
          strike_price: K,
          spot_price: S,
          volatility: vol,
          risk_free_rate: r,
          time_to_maturity: T,
          dividend_yield: q || null,
        };
      }

      if (structure === "spread") {
        return {
          structure_type: "spread",
          spread_type: spreadType,
          direction,
          spot_price: S,
          volatility: vol,
          risk_free_rate: r,
          time_to_maturity: T,
          dividend_yield: q || null,
          strikes,
        };
      }

      if (exoticType === "convertible_bond") {
        const couponPct = parseFloat(exoticValues.coupon_rate);
        const coupon = Number.isFinite(couponPct) ? couponPct / 100 : 0.05;
        const faceValue = parseFloat(exoticValues.face_value) || 1000;
        const maturity = parseFloat(exoticValues.maturity) || 5;
        const paymentFreq = Math.max(1, parseInt(exoticValues.payment_frequency || "2") || 2);
        const creditSpreadPct = parseFloat(exoticValues.credit_spread);
        const creditSpread = Number.isFinite(creditSpreadPct) ? creditSpreadPct / 100 : 0.02;
        const conversionPrice = parseFloat(exoticValues.conversion_price) || 50;
        return {
          structure_type: "exotic",
          exotic_type: "convertible_bond",
          face_value: faceValue,
          coupon_rate: coupon,
          maturity,
          payment_frequency: paymentFreq,
          credit_spread: creditSpread,
          conversion_price: conversionPrice,
          stock_price: S,
          volatility: vol,
          time_to_maturity: T || maturity,
          risk_free_rate: r,
          dividend_yield: q || null,
        };
      }
      if (exoticType === "asian_option") {
        const pooling = (exoticValues.pooling === "max" ? "max" : "average") as "average" | "max";
        const numObs = Math.max(1, parseInt(exoticValues.num_observations || "50") || 50);
        const asianOptType = (exoticValues.asian_option_type === "put" ? "put" : "call") as "call" | "put";
        const asianStrike = parseFloat(exoticValues.strike || strike) || S;
        return {
          structure_type: "exotic",
          exotic_type: "asian_option",
          option_type: asianOptType,
          pooling,
          spot_price: S,
          strike_price: asianStrike,
          volatility: vol,
          risk_free_rate: r,
          time_to_maturity: T || 1,
          dividend_yield: q || null,
          num_observations: numObs,
        };
      }
      if (exoticType === "cliquet_option") {
        const cqType = (exoticValues.cliquet_option_type === "put" ? "put" : "call") as "call" | "put";
        const numResets = Math.max(1, parseInt(exoticValues.num_resets || "12") || 12);
        const pct = (k: string): number | null => {
          const raw = exoticValues[k];
          if (raw === undefined || raw === "") return null;
          const v = parseFloat(raw);
          return Number.isFinite(v) ? v / 100 : null;
        };
        return {
          structure_type: "exotic",
          exotic_type: "cliquet_option",
          option_type: cqType,
          spot_price: S,
          volatility: vol,
          risk_free_rate: r,
          time_to_maturity: T || 1,
          dividend_yield: q || null,
          num_resets: numResets,
          local_cap: pct("local_cap"),
          local_floor: pct("local_floor"),
          global_cap: pct("global_cap"),
          global_floor: pct("global_floor"),
        };
      }
      if (exoticType === "barrier_option") {
        const bKind = (exoticValues.barrier_kind === "put" ? "put" : "call") as "call" | "put";
        const bType = (exoticValues.barrier_type || "down_and_out") as import("./types/index.ts").BarrierType;
        const barrierStrike = parseFloat(exoticValues.barrier_strike) || S;
        const barrierLevel = parseFloat(exoticValues.barrier_level) || S * 0.9;
        return {
          structure_type: "exotic",
          exotic_type: "barrier_option",
          option_type: bKind,
          barrier_type: bType,
          spot_price: S,
          strike_price: barrierStrike,
          barrier: barrierLevel,
          rebate: parseFloat(exoticValues.rebate) || 0,
          volatility: vol,
          risk_free_rate: r,
          time_to_maturity: T || 1,
          dividend_yield: q || null,
        };
      }
      if (exoticType === "autocallable_note") {
        const notional = parseFloat(exoticValues.notional) || 1000;
        const autocallBarrier = parseFloat(exoticValues.autocall_barrier) || 1.05;
        const couponRate = (parseFloat(exoticValues.coupon_rate) || 8) / 100;
        const protectionBarrier = parseFloat(exoticValues.protection_barrier) || 0.80;
        const numObs = Math.max(1, parseInt(exoticValues.num_observations || "4") || 4);
        return {
          structure_type: "exotic",
          exotic_type: "autocallable_note",
          spot_price: S,
          volatility: vol,
          risk_free_rate: r,
          time_to_maturity: T || 1,
          dividend_yield: q || null,
          notional,
          autocall_barrier: autocallBarrier,
          coupon_rate: couponRate,
          protection_barrier: protectionBarrier,
          num_observations: numObs,
          memory_coupon: exoticValues.memory_coupon === "true",
        };
      }
      return {
        structure_type: "exotic",
        exotic_type: "chooser_option",
        spot_price: S,
        strike_price: parseFloat(strike) || S,
        volatility: vol,
        risk_free_rate: r,
        time_to_maturity: T || 1,
        dividend_yield: q || null,
        choose_time: T * 0.5,
      };
    } catch {
      return null;
    }
  }, [structure, optionType, direction, strike, bsArgs, effectiveVol, spreadType, strikes, exoticType, exoticValues]);

  useEffect(() => {
    const req = buildRequest();
    if (!req) { setResult(null); return; }
    const id = setTimeout(() => {
      submitPricing(req)
        .then(r => { setResult(r); setPriceError(null); })
        .catch(err => setPriceError(err instanceof Error ? err.message : "Pricing failed"));
    }, 150);
    return () => clearTimeout(id);
  }, [buildRequest]);

  useEffect(() => {
    if (structure !== "single" || exerciseStyle !== "american" || !selectedSymbol) {
      setCtmcResult(null);
      return;
    }
    const K = parseFloat(strike);
    const T = bsArgs.T;
    if (!K || T <= 0) { setCtmcResult(null); return; }

    setCtmcLoading(true);
    const id = setTimeout(() => {
      submitCtmcPricing({
        symbol: selectedSymbol,
        option_type: optionType,
        direction,
        strike_price: K,
        time_to_maturity: T,
        risk_free_rate: bsArgs.r,
        dividend_yield: bsArgs.q || null,
      })
        .then(setCtmcResult)
        .catch(() => setCtmcResult(null))
        .finally(() => setCtmcLoading(false));
    }, 500);
    return () => { clearTimeout(id); setCtmcLoading(false); };
  }, [selectedSymbol, structure, exerciseStyle, optionType, direction, strike, bsArgs.T, bsArgs.r, bsArgs.q]);

  const rows = useMemo(
    () => buildRows(result, exerciseStyle, ctmcResult, ctmcLoading),
    [result, exerciseStyle, ctmcResult, ctmcLoading]
  );

  const defaultId = useMemo(() => {
    if (exerciseStyle === "american") return "bin_am";
    return "bs";
  }, [exerciseStyle]);

  const activeId = useMemo(() => {
    if (selectedMethod && rows.some(r => r.id === selectedMethod && !r.loading && r.price != null)) {
      return selectedMethod;
    }
    return defaultId;
  }, [selectedMethod, rows, defaultId]);

  const heroPrimary = useMemo(() => {
    if (!rows.length) return null;
    const r = rows.find(x => x.id === activeId) ?? rows[0];
    return { method: `${r.method} · ${r.sub}`, price: r.price, ms: r.ms };
  }, [rows, activeId]);

  const heroSecondary = useMemo(() => {
    return rows
      .filter(r => r.id !== activeId && !r.loading && r.price != null && !r.dim)
      .slice(0, 3)
      .map(r => ({ method: r.method, price: r.price }));
  }, [rows, activeId]);

  const pinCurrent = () => {
    if (heroPrimary?.price != null) {
      setPinned({ price: heroPrimary.price, method: heroPrimary.method });
    }
  };

  const handleExoticChange = (field: string, value: string) => {
    setExoticValues(prev => ({ ...prev, [field]: value }));
  };

  const payoffData = useMemo(() => {
    if (!result) return null;
    const req = buildRequest();
    if (!req) return null;
    return calculatePayoffGrid(req, result);
  }, [result, buildRequest]);

  const lastClose = priceHistory[priceHistory.length - 1]?.price;
  const prevClose = priceHistory[priceHistory.length - 2]?.price;
  const dayChg = lastClose && prevClose ? (lastClose - prevClose) / prevClose : 0;

  return (
    <div className="app">
      <div className="topbar">
        <div className="brand">
          <div className="brand-mark">∫</div>
          <span>Options Pricer</span>
        </div>
        <div className="topbar-sep"></div>
        {ticker && (
          <div className="topbar-meta">
            <span><strong>{ticker.sym}</strong></span>
            {lastClose != null && <span className="tnum">${lastClose.toFixed(2)}</span>}
            {lastClose != null && prevClose != null && (
              <span className="tnum" style={{ color: dayChg >= 0 ? "var(--accent)" : "var(--loss)" }}>
                {dayChg >= 0 ? "+" : ""}{(dayChg * 100).toFixed(2)}%
              </span>
            )}
            <span>σ={(effectiveVol * 100).toFixed(1)}%</span>
            <span>T={bsArgs.T.toFixed(3)}y</span>
            <span>r={(bsArgs.r * 100).toFixed(2)}%</span>
          </div>
        )}
        <div className="topbar-right">
          <button
            className="icon-btn"
            onClick={() => setTweaks({ ...tweaks, theme: tweaks.theme === "dark" ? "light" : "dark" })}
          >
            {tweaks.theme === "dark" ? "☾ Dark" : "☼ Light"}
          </button>
          <button
            className="icon-btn"
            onClick={() => { setModalTab("smile"); setModalOpen(true); }}
            disabled={!ticker}
          >
            Market context ↗
          </button>
          <button
            className="icon-btn primary"
            onClick={pinCurrent}
            disabled={heroPrimary?.price == null}
          >
            Pin snapshot
          </button>
        </div>
      </div>

      <div className="main" data-density={tweaks.density} data-rightpanel={tweaks.rightPanel}>
        <Sidebar
          ticker={ticker}
          onSelectTicker={loadSymbol}
          loadingTicker={marketLoading}
          structure={structure}
          setStructure={setStructure}
          direction={direction}
          setDirection={setDirection}
          exerciseStyle={exerciseStyle}
          setExerciseStyle={setExerciseStyle}
          optionType={optionType}
          setOptionType={setOptionType}
          strike={strike}
          setStrike={setStrike}
          expirationDate={expirationDate}
          setExpirationDate={setExpirationDate}
          spot={spot}
          setSpot={setSpot}
          volatility={volatility}
          setVolatility={setVolatility}
          riskFreeRate={riskFreeRate}
          setRiskFreeRate={setRiskFreeRate}
          dividendYield={dividendYield}
          setDividendYield={setDividendYield}
          volSource={volSource}
          setVolSource={setVolSource}
          volData={volDataForSidebar}
          manualOverride={manualOverride}
          setManualOverride={setManualOverride}
          spreadType={spreadType}
          setSpreadType={setSpreadType}
          strikes={strikes}
          setStrikes={setStrikes}
          exoticType={exoticType}
          setExoticType={setExoticType}
          exoticValues={exoticValues}
          onExoticChange={handleExoticChange}
        />

        <div className="col center">
          <PriceChartPanel
            tickerSym={ticker?.sym ?? null}
            history={priceHistory}
            range={range}
            setRange={setRange}
          />

          {marketError && (
            <div className="empty" style={{ color: "var(--loss)" }}>{marketError}</div>
          )}

          {structure === "spread" && (
            <div style={{ padding: "14px", borderBottom: "1px solid var(--border)" }}>
              <SpreadForm
                spreadType={spreadType}
                onSpreadTypeChange={setSpreadType}
                strikes={strikes}
                onStrikesChange={setStrikes}
                showTypeSelector={true}
              />
            </div>
          )}

          {structure === "single" && heroPrimary && tweaks.emphasis === "hero" && (
            <HeroPrice
              label={`${direction === "long" ? "Long" : "Short"} · ${optionType.toUpperCase()} · ${exerciseStyle}`}
              primary={heroPrimary}
              secondary={heroSecondary}
              ctmc={ctmcResult}
              pinned={pinned}
              onPin={pinCurrent}
              onClearPin={() => setPinned(null)}
            />
          )}

          {result && rows.length > 0 && (
            <div className="results">
              <div className="results-head">
                <div>
                  <div className="results-title">Pricing methods</div>
                  <div className="results-meta mono">
                    <span>σ {(effectiveVol * 100).toFixed(2)}% ({manualOverride ? "manual" : volSource})</span>
                    <span>S ${bsArgs.S.toFixed(2)}</span>
                    <span>K ${bsArgs.K.toFixed(2)}</span>
                    <span>T {bsArgs.T.toFixed(3)}y</span>
                  </div>
                </div>
              </div>
              <PricingRows
                rows={rows}
                selectedId={activeId}
                onSelect={setSelectedMethod}
                emphasis={tweaks.emphasis}
              />
            </div>
          )}

          {priceError && !result && (
            <div className="empty" style={{ color: "var(--loss)" }}>{priceError}</div>
          )}

          {(structure === "single" || structure === "spread") && result?.greeks && (
            <GreeksStrip greeks={result.greeks} bsArgs={bsArgs} />
          )}

          {structure === "single" && result && heroPrimary?.price != null && (
            <PayoffPanel
              bsArgs={bsArgs}
              direction={direction}
              premium={heroPrimary.price}
            />
          )}

          {(structure === "spread" || structure === "exotic") && payoffData && payoffData.length > 0 && (
            <div className="payoff">
              <div className="payoff-head">
                <span style={{ fontSize: 11, color: "var(--fg-3)", textTransform: "uppercase", letterSpacing: "0.06em" }}>
                  Payoff at expiry
                </span>
              </div>
              <PayoffChart data={payoffData} height={160} strike={bsArgs.K} spot={bsArgs.S} />
            </div>
          )}
        </div>

        {tweaks.rightPanel === "visible" && (
          <RightPanel
            ticker={ticker}
            history={priceHistory}
            smile={smile}
            onOpenModal={(tab) => { setModalTab(tab); setModalOpen(true); }}
          />
        )}
      </div>

      <MarketModal
        open={modalOpen}
        tab={modalTab}
        setTab={setModalTab}
        onClose={() => setModalOpen(false)}
        ticker={ticker}
        history={priceHistory}
        smile={smile}
      />

      <TweaksPanel tweaks={tweaks} setTweaks={setTweaks} />
    </div>
  );
}
