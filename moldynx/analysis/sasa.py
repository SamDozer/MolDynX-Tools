"""Solvent-accessible surface area (total + per-residue, Shrake-Rupley/mdtraj)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import mdtraj as md

from moldynx.core.base import BaseAnalysis
from moldynx.core.surface import shrake_rupley
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class SASAAnalysis(BaseAnalysis):
    name = "sasa"
    label = "Solvent-accessible surface area"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    outputs = ["results/sasa_total.csv", "results/sasa_per_residue.csv",
               "figures/sasa_timeseries.png", "figures/sasa_per_residue.png"]
    default_params = {"window": 20, "stride": 10}   # SASA is computed frame by frame (see core.surface)

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        # mdtraj reads the cached solute trajectory
        ctx.core_universe()  # ensure core.{pdb,xtc} exist
        traj = md.load(str(ctx.config.data_dir / "core.xtc"),
                       top=str(ctx.config.data_dir / "core.pdb"), stride=p["stride"])
        time_ns = traj.time / 1000.0
        sasa_res = shrake_rupley(traj, mode="residue")      # nm^2 (multi-frame-safe)
        total = sasa_res.sum(axis=1)
        resids = np.array([r.resSeq for r in traj.topology.residues])
        resnames = np.array([r.name for r in traj.topology.residues])

        ma = st.moving_average(total, p["window"])
        ctx.write_csv(pd.DataFrame({"time_ns": time_ns, "sasa_total_nm2": total,
                                    "sasa_total_movavg_nm2": ma}), "sasa_total.csv")
        ctx.write_csv(pd.DataFrame({"resid": resids, "resname": resnames,
                                    "sasa_mean_nm2": sasa_res.mean(axis=0),
                                    "sasa_std_nm2": sasa_res.std(axis=0)}),
                      "sasa_per_residue.csv")
        mean = float(total.mean())

        fig, ax = plotting.new_axes()
        ax.plot(time_ns, total, color=PALETTE["primary"], alpha=0.35, lw=1.0)
        ax.plot(time_ns, ma, color=PALETTE["primary"], lw=2.2, label="MA")
        ax.axhline(mean, ls="--", lw=1.2, color=PALETTE["accent"],
                   label=f"mean = {mean:.1f} nm$^2$")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Total SASA (nm$^2$)")
        ax.set_title(self.label); ax.legend(loc="best"); ax.margins(x=0.01)
        plotting.save_figure(fig, ctx.fig_path("sasa_timeseries"), dpi=ctx.config.dpi)

        fig, ax = plotting.new_axes(figsize=(8.2, 4.6))
        ax.bar(resids, sasa_res.mean(axis=0), color=PALETTE["primary"], width=1.0)
        ax.set_xlabel("Residue number"); ax.set_ylabel("Mean SASA (nm$^2$)")
        ax.set_title("Per-residue solvent accessibility"); ax.margins(x=0.01)
        plotting.save_figure(fig, ctx.fig_path("sasa_per_residue"), dpi=ctx.config.dpi)

        return {"total_sasa": st.describe(total, "sasa_total_nm2"),
                "convergence": st.plateau_detection(time_ns, total),
                "figure": "sasa_timeseries"}
