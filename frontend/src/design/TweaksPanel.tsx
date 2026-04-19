import { useState } from "react";

export type Tweaks = {
  density: "compact" | "default" | "comfortable";
  emphasis: "hero" | "table";
  rightPanel: "visible" | "hidden";
  theme: "dark" | "light";
};

export const TWEAK_DEFAULTS: Tweaks = {
  density: "default",
  emphasis: "hero",
  rightPanel: "visible",
  theme: "dark",
};

export function TweaksPanel({
  tweaks, setTweaks,
}: {
  tweaks: Tweaks;
  setTweaks: (t: Tweaks) => void;
}) {
  const [visible, setVisible] = useState(false);

  if (!visible) {
    return (
      <button
        className="icon-btn"
        style={{ position: "fixed", bottom: 16, right: 16, zIndex: 40 }}
        onClick={() => setVisible(true)}
        title="Layout tweaks"
      >
        ⚙ Tweaks
      </button>
    );
  }

  const update = <K extends keyof Tweaks>(k: K, v: Tweaks[K]) => {
    setTweaks({ ...tweaks, [k]: v });
  };

  return (
    <div className="tweaks-panel">
      <div className="tweaks-head">
        <h3>Tweaks</h3>
        <button className="icon-btn" onClick={() => setVisible(false)}>×</button>
      </div>

      <div className="tweak-row">
        <label>Density</label>
        <div className="seg">
          {(["compact", "default", "comfortable"] as const).map(d => (
            <button key={d} aria-pressed={tweaks.density === d} onClick={() => update("density", d)}>{d}</button>
          ))}
        </div>
      </div>

      <div className="tweak-row">
        <label>Price emphasis</label>
        <div className="seg">
          {(["hero", "table"] as const).map(e => (
            <button key={e} aria-pressed={tweaks.emphasis === e} onClick={() => update("emphasis", e)}>{e}</button>
          ))}
        </div>
      </div>

      <div className="tweak-row">
        <label>Market context panel</label>
        <div className="seg">
          {(["visible", "hidden"] as const).map(v => (
            <button key={v} aria-pressed={tweaks.rightPanel === v} onClick={() => update("rightPanel", v)}>{v}</button>
          ))}
        </div>
      </div>

      <div className="tweak-row">
        <label>Theme</label>
        <div className="seg">
          {(["dark", "light"] as const).map(t => (
            <button key={t} aria-pressed={tweaks.theme === t} onClick={() => update("theme", t)}>{t}</button>
          ))}
        </div>
      </div>
    </div>
  );
}
