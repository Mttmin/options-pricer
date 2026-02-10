import { useState, useEffect, useRef, useCallback } from "react";
import { searchTickers, type SearchResult } from "../api/client.ts";
import { Button } from "./ui/Button.tsx";

interface TickerSearchProps {
  loading: boolean;
  onSearch: (symbol: string) => void;
}

export function TickerSearch({ loading, onSearch }: TickerSearchProps) {
  const [query, setQuery] = useState("");
  const [suggestions, setSuggestions] = useState<SearchResult[]>([]);
  const [showSuggestions, setShowSuggestions] = useState(false);
  const [highlightedIndex, setHighlightedIndex] = useState(-1);
  const [searchLoading, setSearchLoading] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | undefined>(undefined);

  const fetchSuggestions = useCallback(async (q: string) => {
    if (q.length < 1) {
      setSuggestions([]);
      return;
    }
    setSearchLoading(true);
    try {
      const results = await searchTickers(q);
      setSuggestions(results);
      setShowSuggestions(results.length > 0);
    } catch {
      setSuggestions([]);
    } finally {
      setSearchLoading(false);
    }
  }, []);

  useEffect(() => {
    if (debounceRef.current) {
      clearTimeout(debounceRef.current);
    }
    debounceRef.current = setTimeout(() => {
      fetchSuggestions(query);
    }, 333);
    return () => {
      if (debounceRef.current) {
        clearTimeout(debounceRef.current);
      }
    };
  }, [query, fetchSuggestions]);

  useEffect(() => {
    function handleClickOutside(e: MouseEvent) {
      if (containerRef.current && !containerRef.current.contains(e.target as Node)) {
        setShowSuggestions(false);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  function handleSelect(symbol: string) {
    setQuery(symbol);
    setShowSuggestions(false);
    setSuggestions([]);
    onSearch(symbol);
  }

  function handleSubmit() {
    const trimmed = query.trim().toUpperCase();
    if (trimmed.length === 0) return;
    setShowSuggestions(false);
    onSearch(trimmed);
  }

  function handleKeyDown(e: React.KeyboardEvent) {
    if (!showSuggestions || suggestions.length === 0) {
      if (e.key === "Enter") {
        handleSubmit();
      }
      return;
    }

    if (e.key === "ArrowDown") {
      e.preventDefault();
      setHighlightedIndex((i) => Math.min(i + 1, suggestions.length - 1));
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      setHighlightedIndex((i) => Math.max(i - 1, 0));
    } else if (e.key === "Enter") {
      e.preventDefault();
      if (highlightedIndex >= 0 && highlightedIndex < suggestions.length) {
        handleSelect(suggestions[highlightedIndex].symbol);
      } else {
        handleSubmit();
      }
    } else if (e.key === "Escape") {
      setShowSuggestions(false);
    }
  }

  return (
    <div ref={containerRef} className="relative">
      <div className="flex items-end gap-3">
        <div className="flex-1">
          <label className="text-xs font-medium uppercase tracking-wider text-slate-400 block mb-1">
            Ticker Symbol
          </label>
          <div className="relative">
            <input
              ref={inputRef}
              type="text"
              value={query}
              onChange={(e) => {
                setQuery(e.target.value.toUpperCase());
                setHighlightedIndex(-1);
              }}
              onFocus={() => suggestions.length > 0 && setShowSuggestions(true)}
              onKeyDown={handleKeyDown}
              placeholder="Search ticker..."
              className="w-full rounded border border-slate-600 bg-slate-800 px-3 py-2 text-sm text-white outline-none transition-colors hover:border-slate-500 focus:border-blue-500"
            />
            {searchLoading && (
              <div className="absolute right-3 top-1/2 -translate-y-1/2">
                <svg className="h-4 w-4 animate-spin text-slate-400" viewBox="0 0 24 24" fill="none">
                  <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                  <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                </svg>
              </div>
            )}
          </div>
        </div>
        <Button onClick={handleSubmit} loading={loading} disabled={query.trim().length === 0}>
          Search
        </Button>
      </div>

      {showSuggestions && suggestions.length > 0 && (
        <div className="absolute z-50 mt-1 w-full max-w-md rounded-lg border border-slate-700 bg-slate-800 shadow-xl overflow-hidden">
          {suggestions.map((result, i) => (
            <button
              key={result.symbol}
              onClick={() => handleSelect(result.symbol)}
              onMouseEnter={() => setHighlightedIndex(i)}
              className={`w-full px-4 py-2 text-left text-sm transition-colors ${
                i === highlightedIndex
                  ? "bg-blue-600 text-white"
                  : "text-slate-300 hover:bg-slate-700"
              }`}
            >
              <span className="font-medium">{result.symbol}</span>
              <span className="ml-2 text-xs opacity-70">{result.description}</span>
            </button>
          ))}
        </div>
      )}
    </div>
  );
}
