import { useState, useCallback } from "react";
import type { MarketDataResponse, VolatilityResponse } from "../types/index.ts";
import { fetchMarketData, fetchVolatility } from "../api/client.ts";

export function useMarketData() {
  const [marketData, setMarketData] = useState<MarketDataResponse | null>(null);
  const [volData, setVolData] = useState<VolatilityResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const fetchSymbol = useCallback(async (symbol: string) => {
    setLoading(true);
    setError(null);
    try {
      const [market, vol] = await Promise.all([
        fetchMarketData(symbol),
        fetchVolatility(symbol),
      ]);
      setMarketData(market);
      setVolData(vol);
    } catch (err) {
      const message =
        err instanceof Error ? err.message : "Failed to fetch market data";
      setError(message);
      setMarketData(null);
      setVolData(null);
    } finally {
      setLoading(false);
    }
  }, []);

  return { marketData, volData, loading, error, fetchSymbol };
}
