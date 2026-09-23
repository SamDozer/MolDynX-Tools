"""Salt-bridge analysis (acidic O ↔ basic N), count over time + occupancy."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from MDAnalysis.analysis.distances import distance_array

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE

_ACIDIC = "(resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2)"
_BASIC = ("(resname LYS and name NZ) or (resname ARG and name NH1 NH2 NE) or "
          "(resname HIS HSD HSE HSP HISD HISE and name ND1 NE2)")


class SaltBridges(BaseAnalysis):
    name = "salt_bridges"
    label = "Salt-bridge analysis"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/salt_bridges_timeseries.csv", "figures/salt_bridges_timeseries.png"]
    default_params = {"cutoff": 4.0, "window": 20}

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        acidic = u.select_atoms(_ACIDIC); basic = u.select_atoms(_BASIC)
        if acidic.n_atoms == 0 or basic.n_atoms == 0:
            return {"n_salt_bridges_mean": 0.0, "note": "no chargeable groups"}
        a_res = np.array([f"{a.resname}{a.resid}" for a in acidic])
        b_res = np.array([f"{b.resname}{b.resid}" for b in basic])
        n = len(u.trajectory)
        counts = np.empty(n); times = np.empty(n); pair_frames = defaultdict(int)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[salt_bridges]")):
            close = distance_array(acidic.positions, basic.positions) < p["cutoff"]
            ai, bi = np.where(close)
            pairs = set(zip(a_res[ai], b_res[bi]))
            for pr in pairs:
                pair_frames[pr] += 1
            counts[i] = len(pairs); times[i] = ts.time / 1000.0

        ctx.write_csv(pd.DataFrame({"time_ns": times, "n_salt_bridges": counts,
                                    "movavg": st.moving_average(counts, p["window"])}),
                      "salt_bridges_timeseries.csv")
        occ = pd.DataFrame([{"acidic": a, "basic": b, "occupancy_percent": 100*c/n}
                            for (a, b), c in pair_frames.items()])
        if len(occ):
            occ = occ.sort_values("occupancy_percent", ascending=False)
        ctx.write_csv(occ, "salt_bridges_occupancy.csv")

        fig, ax = plotting.new_axes()
        ax.plot(times, counts, color=PALETTE["primary"], alpha=0.35, lw=1.0)
        ax.plot(times, st.moving_average(counts, p["window"]),
                color=PALETTE["primary"], lw=2.2, label="MA")
        ax.axhline(counts.mean(), ls="--", color=PALETTE["accent"], lw=1.1,
                   label=f"mean = {counts.mean():.1f}")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Number of salt bridges")
        ax.set_title(self.label); ax.legend(loc="best")
        plotting.save_figure(fig, ctx.fig_path("salt_bridges_timeseries"), dpi=ctx.config.dpi)

        return {"n_salt_bridges_mean": float(counts.mean()),
                "n_stable_bridges": int((occ["occupancy_percent"] >= 50).sum()) if len(occ) else 0,
                "figure": "salt_bridges_timeseries"}
