import { useState, useCallback } from "react";
import type { PriceRequest, PriceResponse } from "../types/index.ts";
import { submitPricing } from "../api/client.ts";

export function usePricing() {
  const [result, setResult] = useState<PriceResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const calculate = useCallback(async (request: PriceRequest) => {
    setLoading(true);
    setError(null);
    try {
      const response = await submitPricing(request);
      setResult(response);
    } catch (err) {
      const message =
        err instanceof Error ? err.message : "Pricing request failed";
      setError(message);
      setResult(null);
    } finally {
      setLoading(false);
    }
  }, []);

  return { result, loading, error, calculate };
}
