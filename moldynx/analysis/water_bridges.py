"""
Water-mediated bridges between two partners (partner A – water – partner B).

MDAnalysis ``WaterBridgeAnalysis`` (order 1: one bridging water) on the full
solvated trajectory (it needs the waters and the charges/bonds of the run input),
strided because the search is expensive. Partners are addressed by atom-index
ranges from the chain records of the full topology.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from moldynx import plotting
from moldynx import statistics as st
from moldynx.core.base import BaseAnalysis
from moldynx.core.system import COMPLEX_SYSTEMS
from moldynx.plotting import PALETTE


class WaterBridges(BaseAnalysis):
    name = "water_bridges"
    label = "Water-mediated bridges between partners"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = COMPLEX_SYSTEMS
    order = 60
    default_params = {"stride": 40, "order": 1, "persistent": 0.5}
    outputs = ["results/water_bridges.csv", "tables/water_bridge_summary.csv",
               "figures/water_bridges.png"]

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        chains = sorted(getattr(ctx.system, "chains", []), key=lambda c: c.n_atoms, reverse=True)[:2]
        water = ctx.system.components.get("water") if hasattr(ctx.system, "components") else None
        if len(chains) < 2 or water is None or not water.n_residues:
            return {"status": "skipped",
                    "reason": "needs two protein chains and explicit water in the run input"}
        from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis
        plotting.set_style()
        chains = sorted(chains, key=lambda c: c.atom_start)
        u = ctx.full_universe()
        sel = [f"index {c.atom_start}:{c.atom_stop - 1}" for c in chains]
        wsel = "resname " + " ".join(water.resnames)
        wb = WaterBridgeAnalysis(u, sel[0], sel[1], water_selection=wsel, order=p["order"],
                                 update_selection=True)
        wb.run(step=p["stride"], verbose=False)
        cbt = wb.count_by_time()
        t = np.array([x[0] for x in cbt], float) / 1000.0
        c = np.array([x[1] for x in cbt], float)
        ctx.write_csv(pd.DataFrame({"time_ns": t, "n_water_bridges": c}), "water_bridges.csv")
        rows = []
        try:
            for entry in wb.count_by_type():
                rows.append({"bridge": " | ".join(str(x) for x in entry[:-1]),
                             "frequency": float(entry[-1])})
        except Exception:  # older MDAnalysis: summary by type unavailable
            pass
        tab = pd.DataFrame(rows, columns=["bridge", "frequency"]).sort_values(
            "frequency", ascending=False)
        tab.to_csv(ctx.table_path("water_bridge_summary.csv"), index=False)
        fig, ax = plotting.new_axes()
        ax.plot(t, c, color=PALETTE["primary"], alpha=0.4, lw=1)
        if len(c) > 2:
            ax.plot(t, st.moving_average(c, max(2, len(c) // 10)), color=PALETTE["primary"], lw=2)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Water bridges")
        ax.set_title(f"Water-mediated bridges (order {p['order']}, every {p['stride']} frames)")
        plotting.save_figure(fig, ctx.fig_path("water_bridges"), dpi=ctx.config.dpi)
        return {"partners": [c.segid for c in chains], "stride": p["stride"],
                "frames": int(len(c)), "mean_bridges": float(c.mean()) if len(c) else 0.0,
                "max_bridges": float(c.max()) if len(c) else 0.0,
                "n_bridge_types": int(len(tab)),
                "n_persistent": int((tab.frequency >= p["persistent"]).sum()),
                "figure": "water_bridges"}
