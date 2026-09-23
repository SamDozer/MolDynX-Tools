"""Radius of gyration over time (universal solute compactness)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from moldynx.core.base import BaseAnalysis
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class RadiusOfGyration(BaseAnalysis):
    name = "rog"
    label = "Radius of gyration"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    outputs = ["results/rog.csv", "figures/rog.png"]
    default_params = {"window": 20}

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        sel = "protein" if ctx.system.flags.get("has_protein", True) else "all"
        solute = u.select_atoms(sel) if sel != "all" else u.atoms

        rg = np.empty(len(u.trajectory)); times = np.empty_like(rg)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[rog]")):
            rg[i] = solute.radius_of_gyration() / 10.0
            times[i] = ts.time / 1000.0

        ma = st.moving_average(rg, p["window"])
        ctx.write_csv(pd.DataFrame({"time_ns": times, "rg_nm": rg,
                                    "rg_movavg_nm": ma}), "rog.csv")
        mean, std = float(rg.mean()), float(rg.std(ddof=1))

        fig, ax = plotting.new_axes()
        ax.plot(times, rg, color=PALETTE["primary"], alpha=0.35, lw=1.0)
        ax.plot(times, ma, color=PALETTE["primary"], lw=2.2, label="MA")
        ax.axhline(mean, ls="--", lw=1.2, color=PALETTE["accent"],
                   label=f"mean = {mean:.3f} nm")
        ax.fill_between(times, mean - std, mean + std, color=PALETTE["muted"],
                        alpha=0.25, label=f"$\\pm\\sigma$ = {std:.3f} nm")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel(r"$R_g$ (nm)")
        ax.set_title(self.label); ax.legend(loc="best"); ax.margins(x=0.01)
        plotting.save_figure(fig, ctx.fig_path("rog"), dpi=ctx.config.dpi)

        return {"rg": st.describe(rg, "rg_nm"),
                "convergence": st.plateau_detection(times, rg),
                "rg_initial_nm": float(rg[:min(20, len(rg))].mean()),
                "rg_final_nm": float(rg[len(rg)//2:].mean()),
                "figure": "rog"}
