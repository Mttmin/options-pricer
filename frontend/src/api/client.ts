import type {
  MarketDataResponse,
  VolatilityResponse,
  PriceRequest,
  PriceResponse,
  IVSmileData,
  SofrResponse,
  CtmcPriceRequest,
  CtmcPriceResponse,
} from "../types/index.ts";

const BASE_URL = "/api";

async function handleResponse<T>(response: Response): Promise<T> {
  if (!response.ok) {
    const body = await response.text().catch(() => "");
    throw new Error(
      `API error ${response.status}: ${body || response.statusText}`
    );
  }
  return response.json() as Promise<T>;
}

export async function fetchMarketData(
  symbol: string
): Promise<MarketDataResponse> {
  const response = await fetch(
    `${BASE_URL}/market/${encodeURIComponent(symbol)}`
  );
  return handleResponse<MarketDataResponse>(response);
}

export type SearchResult = {
  symbol: string;
  description: string;
  type: string;
};

export async function searchTickers(query: string): Promise<SearchResult[]> {
  const response = await fetch(
    `${BASE_URL}/search?q=${encodeURIComponent(query)}`
  );
  return handleResponse<SearchResult[]>(response);
}

export type PricePoint = {
  date: string;
  price: number;
};

export async function fetchPriceHistory(
  symbol: string,
  days?: number
): Promise<PricePoint[]> {
  const params = days ? `?days=${days}` : "";
  const response = await fetch(
    `${BASE_URL}/history/${encodeURIComponent(symbol)}${params}`
  );
  return handleResponse<PricePoint[]>(response);
}

export async function fetchVolatility(
  symbol: string
): Promise<VolatilityResponse> {
  const response = await fetch(
    `${BASE_URL}/volatility/${encodeURIComponent(symbol)}`
  );
  return handleResponse<VolatilityResponse>(response);
}

export async function submitPricing(
  request: PriceRequest
): Promise<PriceResponse> {
  const response = await fetch(`${BASE_URL}/price`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(request),
  });
  return handleResponse<PriceResponse>(response);
}

export async function fetchOptionChain(
  symbol: string
): Promise<IVSmileData> {
  const response = await fetch(
    `${BASE_URL}/options-chain/${encodeURIComponent(symbol)}`
  );
  return handleResponse<IVSmileData>(response);
}

export async function fetchSofr(): Promise<SofrResponse> {
  const response = await fetch(`${BASE_URL}/sofr`);
  return handleResponse<SofrResponse>(response);
}

export async function submitCtmcPricing(
  request: CtmcPriceRequest
): Promise<CtmcPriceResponse> {
  const response = await fetch(`${BASE_URL}/price/ctmc`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(request),
  });
  return handleResponse<CtmcPriceResponse>(response);
}
