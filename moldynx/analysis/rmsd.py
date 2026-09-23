"""Backbone and Cα RMSD versus the reference frame (universal)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from MDAnalysis.analysis import rms

from moldynx.core.base import BaseAnalysis
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class RMSDAnalysis(BaseAnalysis):
    name = "rmsd"
    label = "Backbone & Cα RMSD"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    outputs = ["results/rmsd.csv", "figures/rmsd.png", "figures/rmsd.pdf"]
    default_params = {"ref_frame": 0, "window": 20}

    def _series(self, u, select, ref_frame):
        R = rms.RMSD(u, select=select, ref_frame=ref_frame).run()
        return R.results.rmsd[:, 2] / 10.0  # Å -> nm

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        # protein if present, else nucleic backbone
        has_protein = ctx.system.flags.get("has_protein", True)
        bb_sel = "backbone" if has_protein else "nucleic and name P O3' O5' C3' C4' C5'"
        ca_sel = ctx.selection("ca", "protein and name CA") if has_protein \
            else "nucleic and name P"

        time_ns = ctx.times_ns(u)
        rmsd_bb = self._series(u, bb_sel, p["ref_frame"])
        rmsd_ca = self._series(u, ca_sel, p["ref_frame"])
        ma_bb = st.moving_average(rmsd_bb, p["window"])
        ma_ca = st.moving_average(rmsd_ca, p["window"])

        df = pd.DataFrame({"time_ns": time_ns,
                           "rmsd_backbone_nm": rmsd_bb, "rmsd_backbone_movavg_nm": ma_bb,
                           "rmsd_ca_nm": rmsd_ca, "rmsd_ca_movavg_nm": ma_ca})
        ctx.write_csv(df, "rmsd.csv")

        fig, ax = plotting.new_axes()
        ax.plot(time_ns, rmsd_bb, color=PALETTE["primary"], alpha=0.35, lw=1.0)
        ax.plot(time_ns, ma_bb, color=PALETTE["primary"], lw=2.2, label="Backbone")
        ax.plot(time_ns, rmsd_ca, color=PALETTE["secondary"], alpha=0.35, lw=1.0)
        ax.plot(time_ns, ma_ca, color=PALETTE["secondary"], lw=2.2, label=r"C$\alpha$")
        plateau = float(np.mean(rmsd_bb[len(rmsd_bb)//2:]))
        plotting.add_reference_line(ax, plateau, label=f"plateau = {plateau:.3f} nm")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("RMSD (nm)")
        ax.set_title(self.label); ax.legend(loc="upper left"); ax.margins(x=0.01)
        plotting.save_figure(fig, ctx.fig_path("rmsd"), dpi=ctx.config.dpi)

        conv = st.plateau_detection(time_ns, rmsd_bb)
        return {"backbone": st.describe(rmsd_bb, "rmsd_backbone_nm"),
                "ca": st.describe(rmsd_ca, "rmsd_ca_nm"),
                "plateau_mean_nm": plateau, "convergence": conv,
                "figure": "rmsd"}
