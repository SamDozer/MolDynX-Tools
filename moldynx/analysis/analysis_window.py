"""
Is there a stationary window for binding-related averages? (Chodera equilibration detection)

For each binding-relevant observable (minimum distance, contacts, interface residues,
buried area from ``interface``; RMSD columns from ``rmsd``): the start t₀ that
maximises the number of effectively independent samples, then -- inside [t₀, end] --
the linear drift and a half-vs-half Welch test. The conservative t₀ is the latest
over the binding observables.

If no window is stationary the result says so; the primary window for averages
(e.g. MM-GBSA) is then a **user decision** (``binding_energy.primary_window`` in the
config), never chosen silently. The proposal offered is the full run plus the
final 20 %.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from moldynx.core.base import BaseAnalysis
from moldynx.statistics import detect_equilibration, drift

BINDING = {"min_interface_dist_nm": "Minimum inter-partner distance (nm)",
           "n_inter_contacts": "Residue–residue contacts",
           "n_interface_residues": "Interface residues",
           "buried_area_nm2": "Buried area (nm²)"}


class AnalysisWindow(BaseAnalysis):
    name = "analysis_window"
    label = "Stationarity and analysis-window selection"
    category = "quality"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    order = 160
    default_params = {"alpha": 0.01}
    outputs = ["results/window_selection.json", "results/window_selection.csv"]

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        series: dict[str, tuple[np.ndarray, np.ndarray, bool]] = {}
        f = ctx.csv_path("interface_timeseries.csv")
        if f.exists():
            df = pd.read_csv(f)
            for col, lab in BINDING.items():
                if col in df and df[col].notna().sum() > 10:
                    ok = df[col].notna()
                    series[lab] = (df.time_ns[ok].to_numpy(), df[col][ok].to_numpy(float), True)
        f = ctx.csv_path("rmsd.csv")
        if f.exists():
            df = pd.read_csv(f)
            tcol = "time_ns" if "time_ns" in df else df.columns[0]
            for col in df.columns:
                if col != tcol and df[col].dtype.kind == "f" and "rmsd" in col.lower():
                    series[col] = (df[tcol].to_numpy(), df[col].to_numpy(float), False)
        if not series:
            return {"status": "skipped", "reason": "no interface or RMSD time series to test"}

        rows = []
        for name, (t, x, binding) in series.items():
            eq = detect_equilibration(x, step=max(1, len(x) // 100))
            i0 = eq["t0"]
            d = drift(t[i0:], x[i0:])
            span = float(t[-1] - t[i0]) or 1.0
            rows.append({"observable": name, "binding_relevant": binding,
                         "t0_ns": float(t[i0]), "stat_ineff": eq["g"], "n_eff": eq["n_eff"],
                         "mean_after_t0": float(x[i0:].mean()),
                         "slope_per_10ns": None if d["slope"] is None else d["slope"] * 10,
                         "p_slope": d["p_slope"], "half_vs_half_p": d["half_p"],
                         "relative_drift": None if d["slope"] is None else
                         d["slope"] * span / (abs(x[i0:].mean()) + 1e-12),
                         "stationary": bool(d["p_slope"] is not None and d["p_slope"] > p["alpha"]
                                            and d["half_p"] > p["alpha"])})
        tab = pd.DataFrame(rows)
        tab.to_csv(ctx.csv_path("window_selection.csv"), index=False)
        t_all = next(iter(series.values()))[0]
        t_end = float(t_all[-1])
        bind = tab[tab.binding_relevant]
        t0 = float(bind.t0_ns.max()) if len(bind) else float(tab.t0_ns.max())
        step = max(t_end * 0.05, 1e-9)
        t0 = float(np.ceil(t0 / step) * step)
        stationary = bool(len(bind) and bind.stationary.all())
        chosen = (ctx.config.binding_energy or {}).get("primary_window", "ask")
        result = {
            "conservative_t0_ns": t0, "end_ns": t_end,
            "binding_observables_stationary": stationary,
            "non_stationary": list(tab[~tab.stationary].observable),
            "proposal": {"full": [float(t_all[0]), t_end],
                         "final_20pct": [round(t_end * 0.8, 6), t_end],
                         "after_t0": [t0, t_end]},
            "primary_window": chosen,
            "decision_required": chosen in (None, "ask"),
            "statement": ("binding-related observables are stationary after t0"
                          if stationary else
                          "no stationary window: report time-resolved values; the primary "
                          "window for averages is a user decision"),
        }
        ctx.csv_path("window_selection.json").write_text(json.dumps(
            {"summary": result, "observables": rows}, indent=2, default=float), encoding="utf-8")
        return result
