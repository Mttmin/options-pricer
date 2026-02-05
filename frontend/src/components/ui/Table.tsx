import type { ReactNode } from "react";

interface TableProps {
  headers: string[];
  rows: (string | number | ReactNode)[][];
}

export function Table({ headers, rows }: TableProps) {
  return (
    <div className="overflow-x-auto rounded border border-slate-700">
      <table className="w-full text-left text-sm">
        <thead>
          <tr className="border-b border-slate-700 bg-slate-800">
            {headers.map((header, i) => (
              <th
                key={i}
                className="px-4 py-2.5 text-xs font-medium uppercase tracking-wider text-slate-400"
              >
                {header}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map((row, rowIdx) => (
            <tr
              key={rowIdx}
              className={`border-b border-slate-700/50 transition-colors hover:bg-slate-800/50 ${
                rowIdx % 2 === 0 ? "bg-slate-900" : "bg-slate-900/60"
              }`}
            >
              {row.map((cell, cellIdx) => (
                <td key={cellIdx} className="px-4 py-2 text-slate-200">
                  {cell}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
