"""
Example plugin: protein end-to-end distance.

This is a complete, working template showing how to add a new analysis to
MolDynX Tools WITHOUT modifying the core. Copy this file, rename the class, implement
``run()``, and it will be auto-discovered and selectable on the CLI as
``moldynx analyze --analyses end_to_end`` (or run automatically for the systems
it declares support for).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import SystemType
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class EndToEndDistance(BaseAnalysis):
    name = "end_to_end"
    label = "End-to-end distance"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = {SystemType.PROTEIN_ONLY, SystemType.PROTEIN_PEPTIDE}
    outputs = ["results/end_to_end.csv", "figures/end_to_end.png"]
    default_params = {"window": 20}

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        ca = u.select_atoms("protein and name CA")
        d = np.empty(len(u.trajectory)); times = np.empty_like(d)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[end_to_end]")):
            pos = ca.positions / 10.0
            d[i] = float(np.linalg.norm(pos[0] - pos[-1]))
            times[i] = ts.time / 1000.0

        ctx.write_csv(pd.DataFrame({"time_ns": times, "end_to_end_nm": d,
                                    "movavg_nm": st.moving_average(d, p["window"])}),
                      "end_to_end.csv")
        fig, ax = plotting.new_axes()
        ax.plot(times, d, color=PALETTE["primary"], alpha=0.4, lw=1.0)
        ax.plot(times, st.moving_average(d, p["window"]),
                color=PALETTE["primary"], lw=2.2, label="MA")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("End-to-end distance (nm)")
        ax.set_title(self.label); ax.legend(loc="best"); ax.margins(x=0.01)
        plotting.save_figure(fig, ctx.fig_path("end_to_end"), dpi=ctx.config.dpi)
        return {"end_to_end": st.describe(d, "end_to_end_nm"), "figure": "end_to_end"}
