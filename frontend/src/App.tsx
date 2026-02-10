import { useState, useCallback, useEffect } from "react";
import type {
  StructureType,
  Direction,
  OptionType,
  SpreadType,
  ExoticType,
  VolSource,
  ExerciseStyle,
  SpreadStrikes,
  PriceRequest,
  MarketDataResponse,
  VolatilityResponse,
  PriceResponse,
  IVSmileData,
} from "./types/index.ts";
import {
  fetchMarketData,
  fetchVolatility,
  submitPricing,
  fetchPriceHistory,
  fetchOptionChain,
  fetchSofr,
  type PricePoint,
} from "./api/client.ts";
import { calculateTTM } from "./utils/ttm.ts";
import { TickerSearch } from "./components/TickerSearch.tsx";
import { MarketDataDisplay } from "./components/MarketDataDisplay.tsx";
import { PriceChart } from "./components/PriceChart.tsx";
import { StructureSelector } from "./components/StructureSelector.tsx";
import { DirectionToggle } from "./components/DirectionToggle.tsx";
import { ManualOverrideToggle } from "./components/ManualOverrideToggle.tsx";
import { ExerciseStyleToggle } from "./components/ExerciseStyleToggle.tsx";
import { SingleOptionForm } from "./components/inputs/SingleOptionForm.tsx";
import { SpreadForm } from "./components/inputs/SpreadForm.tsx";
import { ExoticForm } from "./components/inputs/ExoticForm.tsx";
import { VolatilityToggle } from "./components/VolatilityToggle.tsx";
import { PricingTable } from "./components/PricingTable.tsx";
import { PayoffDiagram } from "./components/PayoffDiagram.tsx";
import { IVSmileChart } from "./components/IVSmileChart.tsx";
import { OptionChainDisplay } from "./components/OptionChainDisplay.tsx";
import { MonteCarloPathsChart } from "./components/MonteCarloPathsChart.tsx";
import { Button } from "./components/ui/Button.tsx";

export default function App() {
  // Ticker state
  const [selectedTicker, setSelectedTicker] = useState<string | null>(null);

  // Market data state
  const [marketData, setMarketData] = useState<MarketDataResponse | null>(null);
  const [volData, setVolData] = useState<VolatilityResponse | null>(null);
  const [priceHistory, setPriceHistory] = useState<PricePoint[]>([]);
  const [ivSmileData, setIVSmileData] = useState<IVSmileData | null>(null);
  const [marketLoading, setMarketLoading] = useState(false);
  const [marketError, setMarketError] = useState<string | null>(null);

  // Structure state
  const [structureType, setStructureType] = useState<StructureType>("single");
  const [direction, setDirection] = useState<Direction>("long");
  const [manualOverride, setManualOverride] = useState(false);
  const [exerciseStyle, setExerciseStyle] = useState<ExerciseStyle>("european");

  // Single option state
  const [singleValues, setSingleValues] = useState({
    optionType: "call" as OptionType,
    strike: "",
    expirationDate: "",
    spotPrice: "",
    volatility: "",
    riskFreeRate: "0.05",
    dividendYield: "",
  });

  // Spread state
  const [spreadType, setSpreadType] = useState<SpreadType>("straddle");
  const [strikes, setStrikes] = useState<SpreadStrikes>({});
  const [shortTermMaturity, setShortTermMaturity] = useState("");
  const [longTermMaturity, setLongTermMaturity] = useState("");

  // Exotic state
  const [exoticType, setExoticType] = useState<ExoticType>("convertible_bond");
  const [exoticValues, setExoticValues] = useState<Record<string, string>>({
    risk_free_rate: "0.05",
  });

  // Volatility source
  const [volSource, setVolSource] = useState<VolSource>("ema");

  // Pricing state
  const [priceResult, setPriceResult] = useState<PriceResponse | null>(null);
  const [priceLoading, setPriceLoading] = useState(false);
  const [priceError, setPriceError] = useState<string | null>(null);

  // Recalculate when volatility source changes (if we have a previous result)
  useEffect(() => {
    if (priceResult && !manualOverride) {
      handleCalculate();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [volSource]);

  // Fetch SOFR rate on app load
  useEffect(() => {
    const loadSofr = async () => {
      try {
        const sofrData = await fetchSofr();
        console.log("SOFR rate loaded:", sofrData.rate);
        setSingleValues((prev) => ({
          ...prev,
          riskFreeRate: sofrData.rate.toFixed(4),
        }));
      } catch (err) {
        console.warn("SOFR fetch failed, using default 0.05:", err);
        setSingleValues((prev) => ({
          ...prev,
          riskFreeRate: "0.05",
        }));
      }
    };
    loadSofr();
  }, []);

  const handleSearch = useCallback(async (symbol: string) => {
    console.log("Searching for symbol:", symbol);
    setMarketLoading(true);
    setMarketError(null);
    setSelectedTicker(symbol);
    try {
      const [data, history, ivSmile] = await Promise.all([
        fetchMarketData(symbol),
        fetchPriceHistory(symbol, 90),
        fetchOptionChain(symbol).catch((err) => {
          console.warn("Option chain fetch failed:", err);
          return null;
        }),
      ]);
      console.log("Market data:", data);
      console.log("Price history:", history);
      console.log("IV Smile data:", ivSmile);

      setMarketData(data);
      setPriceHistory(history);
      setIVSmileData(ivSmile);
      setSingleValues((prev) => ({
        ...prev,
        spotPrice: data.spot_price.toFixed(2),
        volatility: (
          data.implied_volatility ?? data.historical_volatility
        ).toFixed(4),
        dividendYield: data.dividend_yield
          ? (data.dividend_yield * 100).toFixed(2)
          : "",
      }));
      // Fetch vol data in background
      fetchVolatility(symbol).then(setVolData).catch((err) => {
        console.warn("Volatility fetch failed:", err);
      });
    } catch (err) {
      console.error("Search error:", err);
      setMarketError(err instanceof Error ? err.message : "Failed to fetch data");
    } finally {
      setMarketLoading(false);
    }
  }, []);

  const handleSingleChange = useCallback((field: string, value: string) => {
    setSingleValues((prev) => ({ ...prev, [field]: value }));
  }, []);

  const handleExoticChange = useCallback((field: string, value: string) => {
    setExoticValues((prev) => ({ ...prev, [field]: value }));
  }, []);

  const handleMaturityChange = useCallback((field: string, value: string) => {
    if (field === "shortTermMaturity") setShortTermMaturity(value);
    else setLongTermMaturity(value);
  }, []);

  const getSelectedVol = useCallback((): number => {
    if (!volData && !marketData) return parseFloat(singleValues.volatility) || 0;
    if (volSource === "implied")
      return volData?.implied ?? marketData?.implied_volatility ?? marketData?.historical_volatility ?? 0;
    if (volSource === "historical")
      return volData?.historical ?? marketData?.historical_volatility ?? 0;
    if (volSource === "ema") return volData?.ema ?? volData?.historical ?? 0;
    if (volSource === "vix_correlated")
      return volData?.vix_correlated ?? volData?.historical ?? 0;
    return 0;
  }, [volData, marketData, volSource, singleValues.volatility]);

  const getCurrentRequest = useCallback((): PriceRequest | null => {
    if (!priceResult) return null;

    try {
      const vol = manualOverride
        ? parseFloat(singleValues.volatility) || 0.2
        : getSelectedVol() || 0.2;

      if (structureType === "single") {
        const ttm = singleValues.expirationDate
          ? calculateTTM(singleValues.expirationDate)
          : 0.25;
        const spot = parseFloat(singleValues.spotPrice) || marketData?.spot_price || 100;
        const strike = parseFloat(singleValues.strike) || spot;
        const rfr = parseFloat(singleValues.riskFreeRate) || 0.05;
        const divStr = singleValues.dividendYield.trim();
        const div = divStr ? parseFloat(divStr) / 100 : null;

        return {
          structure_type: "single",
          option_type: singleValues.optionType,
          direction,
          strike_price: strike,
          spot_price: spot,
          volatility: vol,
          risk_free_rate: rfr,
          time_to_maturity: ttm,
          dividend_yield: div,
        };
      } else if (structureType === "spread") {
        const spot = parseFloat(singleValues.spotPrice) || marketData?.spot_price || 100;
        const rfr = parseFloat(singleValues.riskFreeRate) || 0.05;
        const ttm = singleValues.expirationDate
          ? calculateTTM(singleValues.expirationDate)
          : 0.25;
        const divStr = singleValues.dividendYield.trim();
        const div = divStr ? parseFloat(divStr) / 100 : null;

        return {
          structure_type: "spread",
          spread_type: spreadType,
          direction,
          spot_price: spot,
          volatility: vol,
          risk_free_rate: rfr,
          time_to_maturity: ttm,
          dividend_yield: div,
          strikes,
          short_term_maturity: shortTermMaturity
            ? parseFloat(shortTermMaturity)
            : undefined,
          long_term_maturity: longTermMaturity
            ? parseFloat(longTermMaturity)
            : undefined,
        };
      }
    } catch {
      return null;
    }

    return null;
  }, [
    priceResult, manualOverride, singleValues, getSelectedVol, structureType,
    direction, marketData, spreadType, strikes, shortTermMaturity, longTermMaturity,
  ]);

  const handleCalculate = useCallback(async () => {
    setPriceLoading(true);
    setPriceError(null);

    try {
      let request: PriceRequest;
      const vol = manualOverride
        ? parseFloat(singleValues.volatility) || 0.2
        : getSelectedVol() || 0.2;

      if (structureType === "single") {
        const ttm = singleValues.expirationDate
          ? calculateTTM(singleValues.expirationDate)
          : 0.25;
        const spot = parseFloat(singleValues.spotPrice) || marketData?.spot_price;
        const strike = parseFloat(singleValues.strike);
        const rfr = parseFloat(singleValues.riskFreeRate) || 0.05;
        const divStr = singleValues.dividendYield.trim();
        const div = divStr ? parseFloat(divStr) / 100 : null;

        if (!spot || isNaN(spot)) {
          throw new Error("Spot price is required. Search for a ticker or enter manually.");
        }
        if (!strike || isNaN(strike)) {
          throw new Error("Strike price is required.");
        }
        if (ttm <= 0) {
          throw new Error("Expiration date must be in the future.");
        }

        request = {
          structure_type: "single",
          option_type: singleValues.optionType,
          direction,
          strike_price: strike,
          spot_price: spot,
          volatility: vol,
          risk_free_rate: rfr,
          time_to_maturity: ttm,
          dividend_yield: div,
        };
      } else if (structureType === "spread") {
        const spot = parseFloat(singleValues.spotPrice) || marketData?.spot_price || 100;
        const rfr = parseFloat(singleValues.riskFreeRate) || 0.05;
        const ttm = singleValues.expirationDate
          ? calculateTTM(singleValues.expirationDate)
          : 0.25;
        const divStr = singleValues.dividendYield.trim();
        const div = divStr ? parseFloat(divStr) / 100 : null;

        request = {
          structure_type: "spread",
          spread_type: spreadType,
          direction,
          spot_price: spot,
          volatility: vol,
          risk_free_rate: rfr,
          time_to_maturity: ttm,
          dividend_yield: div,
          strikes,
          short_term_maturity: shortTermMaturity
            ? parseFloat(shortTermMaturity)
            : undefined,
          long_term_maturity: longTermMaturity
            ? parseFloat(longTermMaturity)
            : undefined,
        };
      } else {
        if (exoticType === "convertible_bond") {
          request = {
            structure_type: "exotic",
            exotic_type: "convertible_bond",
            face_value: parseFloat(exoticValues.face_value || "1000"),
            coupon_rate: parseFloat(exoticValues.coupon_rate || "0.05"),
            maturity: parseFloat(exoticValues.maturity || "5"),
            payment_frequency: parseInt(exoticValues.payment_frequency || "2"),
            credit_spread: parseFloat(exoticValues.credit_spread || "0.02"),
            conversion_price: parseFloat(exoticValues.conversion_price || "50"),
            stock_price: parseFloat(
              exoticValues.stock_price ||
              singleValues.spotPrice ||
              "100"
            ),
            volatility: vol,
            time_to_maturity: parseFloat(
              exoticValues.time_to_maturity || "5"
            ),
            risk_free_rate: parseFloat(
              exoticValues.risk_free_rate || "0.05"
            ),
            dividend_yield: exoticValues.dividend_yield
              ? parseFloat(exoticValues.dividend_yield)
              : null,
          };
        } else {
          request = {
            structure_type: "exotic",
            exotic_type: "chooser_option",
            spot_price: parseFloat(
              exoticValues.spot_price ||
              singleValues.spotPrice ||
              "100"
            ),
            strike_price: parseFloat(exoticValues.strike_price || "100"),
            volatility: vol,
            risk_free_rate: parseFloat(
              exoticValues.risk_free_rate || "0.05"
            ),
            time_to_maturity: parseFloat(
              exoticValues.time_to_maturity || "1"
            ),
            dividend_yield: exoticValues.dividend_yield
              ? parseFloat(exoticValues.dividend_yield)
              : null,
            choose_time: parseFloat(exoticValues.choose_time || "0.5"),
          };
        }
      }

      const result = await submitPricing(request);
      setPriceResult(result);
    } catch (err) {
      setPriceError(
        err instanceof Error ? err.message : "Pricing failed"
      );
    } finally {
      setPriceLoading(false);
    }
  }, [
    structureType, direction, singleValues, manualOverride, getSelectedVol,
    marketData, spreadType, strikes, shortTermMaturity,
    longTermMaturity, exoticType, exoticValues,
  ]);

  return (
    <div className="min-h-screen bg-slate-950 text-white">
      <div className="max-w-5xl mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold text-center mb-2">Options Pricer</h1>
        <p className="text-slate-400 text-center text-sm mb-8">
          Multi-method pricing engine for options, spreads, and exotics
        </p>

        {/* Ticker Search */}
        <section className="mb-6">
          <TickerSearch
            loading={marketLoading}
            onSearch={handleSearch}
          />
          {marketError && (
            <p className="text-red-400 text-sm mt-2">{marketError}</p>
          )}
        </section>

        {/* Price Chart */}
        <section className="mb-6">
          <PriceChart
            data={priceHistory}
            ticker={selectedTicker}
            currentPrice={marketData?.spot_price}
            loading={marketLoading}
          />
        </section>



        {/* Structure Selection */}
        <section className="mb-6 space-y-4">
          <div className="flex items-center gap-4 flex-wrap">
            <StructureSelector
              selected={structureType}
              onChange={setStructureType}
            />
            {structureType !== "exotic" && (
              <DirectionToggle direction={direction} onChange={setDirection} />
            )}
            {structureType !== "exotic" && (
              <ExerciseStyleToggle
                style={exerciseStyle}
                onChange={setExerciseStyle}
              />
            )}
            <ManualOverrideToggle
              enabled={manualOverride}
              onChange={setManualOverride}
            />
          </div>

          {/* Dynamic Form */}
          <div className="bg-slate-900 rounded-lg p-6 border border-slate-800">
            {structureType === "single" && (
              <SingleOptionForm
                values={singleValues}
                onChange={handleSingleChange}
                manualOverride={manualOverride}
              />
            )}
            {structureType === "spread" && (
              <div className="space-y-4">
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <div>
                    <SingleOptionForm
                      values={singleValues}
                      onChange={handleSingleChange}
                      manualOverride={manualOverride}
                      isSpread={true}
                    />
                  </div>
                  <div>
                    <SpreadForm
                      spreadType={spreadType}
                      onSpreadTypeChange={setSpreadType}
                      strikes={strikes}
                      onStrikesChange={setStrikes}
                      shortTermMaturity={shortTermMaturity}
                      longTermMaturity={longTermMaturity}
                      onMaturityChange={handleMaturityChange}
                    />
                  </div>
                </div>
              </div>
            )}
            {structureType === "exotic" && (
              <ExoticForm
                exoticType={exoticType}
                onExoticTypeChange={setExoticType}
                values={exoticValues}
                onChange={handleExoticChange}
              />
            )}
          </div>
        </section>

        {/* Calculate Button */}
        <section className="mb-8 text-center">
          <Button onClick={handleCalculate} loading={priceLoading}>
            Calculate Prices
          </Button>
          {priceError && (
            <p className="text-red-400 text-sm mt-2">{priceError}</p>
          )}
        </section>

        {/* Pricing Results */}
        {(priceResult || priceLoading) && (
          <section className="mb-8">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-xl font-semibold">Pricing Results</h2>
              <VolatilityToggle
                source={volSource}
                onChange={setVolSource}
                volData={volData}
              />
            </div>
            <PricingTable
              result={priceResult}
              loading={priceLoading}
              exerciseStyle={structureType === "exotic" ? "european" : exerciseStyle}
            />
          </section>
        )}

        {/* Payoff Diagram */}
        {(priceResult || priceLoading) && (
          <section className="mb-8">
            <PayoffDiagram
              request={getCurrentRequest()}
              priceResult={priceResult}
              loading={priceLoading}
            />
          </section>
        )}

        {/* Monte Carlo Paths Visualization */}
        {priceResult?.pricing.monte_carlo && structureType === "single" && (
          <section className="mb-8">
            <MonteCarloPathsChart
              mcResult={priceResult.pricing.monte_carlo}
              spotPrice={parseFloat(singleValues.spotPrice) || marketData?.spot_price || 100}
              timeToMaturity={
                singleValues.expirationDate
                  ? calculateTTM(singleValues.expirationDate)
                  : 0.25
              }
            />
          </section>
        )}

        {/* Market Data Display */}
        {marketData && (
          <section className="mb-6">
            <MarketDataDisplay data={marketData} />
          </section>
        )}

        {/* IV Smile Chart */}
        {ivSmileData && (
          <section className="mb-6">
            <IVSmileChart
              data={ivSmileData}
              loading={marketLoading}
              underlyingPrice={marketData?.spot_price ?? null}
            />
          </section>
        )}

        {/* Option Chain Display */}
        {ivSmileData && (
          <section className="mb-6">
            <OptionChainDisplay
              data={ivSmileData}
              underlyingPrice={marketData?.spot_price ?? null}
            />
          </section>
        )}
      </div>
    </div>
  );
}
