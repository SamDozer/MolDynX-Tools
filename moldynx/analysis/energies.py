"""Thermodynamic observables from the GROMACS .edr (panedr): E, T, P, density."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from moldynx.core.base import BaseAnalysis
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE

_WANTED = {
    "Potential": ["Potential"], "Kinetic": ["Kinetic En.", "Kinetic"],
    "Total Energy": ["Total Energy"], "Temperature": ["Temperature"],
    "Pressure": ["Pressure"], "Density": ["Density"],
}
_LABELS = {"Potential": "Potential (kJ/mol)", "Kinetic": "Kinetic (kJ/mol)",
           "Total Energy": "Total energy (kJ/mol)", "Temperature": "Temperature (K)",
           "Pressure": "Pressure (bar)", "Density": "Density (kg m$^{-3}$)"}


class Energies(BaseAnalysis):
    name = "energies"
    label = "Thermodynamic stability observables"
    category = "thermodynamics"
    required_files = {"energy"}
    supported_systems = {"*"}
    outputs = ["results/energies.csv", "tables/energy_summary.csv", "figures/energies.png"]

    def run(self, ctx) -> dict:
        import panedr
        plotting.set_style()
        edr = ctx.fileset.energy or ctx.config.energy
        if edr is None or not Path(edr).exists():
            return {"note": "edr not available"}
        edf = panedr.edr_to_df(str(edr))
        time_ns = edf["Time"].to_numpy() / 1000.0

        def find(terms):
            for t in terms:
                for c in edf.columns:
                    if t.lower() in c.lower():
                        return c
            return None

        series = {"time_ns": time_ns}; found = {}
        for name, terms in _WANTED.items():
            col = find(terms)
            if col is not None:
                series[name] = edf[col].to_numpy(); found[name] = col
        out = pd.DataFrame(series)
        ctx.write_csv(out, "energies.csv")
        st.summary_frame({k: out[k] for k in found}).to_csv(
            ctx.table_path("energy_summary.csv"))

        keys = list(found); ncols = 3; nrows = int(np.ceil(len(keys)/ncols))
        fig, axes = plotting.style.plt.subplots(nrows, ncols, figsize=(4.2*ncols, 3.2*nrows))
        axes = np.atleast_1d(axes).ravel()
        colors = [PALETTE["primary"], PALETTE["secondary"], PALETTE["green"],
                  PALETTE["purple"], PALETTE["helix"], PALETTE["accent"]]
        for ax, k, c in zip(axes, keys, colors):
            y = out[k].to_numpy()
            ax.plot(time_ns, y, color=c, alpha=0.3, lw=0.8)
            ax.plot(time_ns, st.running_mean(y), color=c, lw=2.0)
            ax.set_xlabel("Time (ns)"); ax.set_ylabel(_LABELS.get(k, k), fontsize=11)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        for ax in axes[len(keys):]:
            ax.set_visible(False)
        fig.suptitle("Thermodynamic stability observables", fontweight="bold")
        plotting.save_figure(fig, ctx.fig_path("energies"), dpi=ctx.config.dpi)

        summary = {k: {"mean": float(out[k].mean()), "std": float(out[k].std(ddof=1))}
                   for k in found}
        summary["figure"] = "energies"
        return summary
