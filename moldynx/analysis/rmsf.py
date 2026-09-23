"""Per-residue Cα RMSF with flexible-region highlighting (protein systems)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from MDAnalysis.analysis import align, rms

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import SystemType
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE

_PROTEIN_SYSTEMS = {SystemType.PROTEIN_ONLY, SystemType.PROTEIN_PROTEIN,
                    SystemType.PROTEIN_PEPTIDE, SystemType.PROTEIN_LIGAND,
                    SystemType.PROTEIN_DNA, SystemType.PROTEIN_RNA,
                    SystemType.PROTEIN_MEMBRANE, SystemType.MIXED}


class RMSFAnalysis(BaseAnalysis):
    name = "rmsf"
    label = "Per-residue Cα RMSF"
    category = "dynamics"
    required_files = {"trajectory", "topology"}
    supported_systems = _PROTEIN_SYSTEMS
    outputs = ["results/rmsf.csv", "figures/rmsf.png"]
    default_params = {"ref_frame": 0}

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        sel = "protein and name CA"
        avg = align.AverageStructure(u, u, select=sel, ref_frame=p["ref_frame"]).run()
        align.AlignTraj(u, avg.results.universe, select=sel, in_memory=True).run()
        ca = u.select_atoms(sel)
        rmsf = rms.RMSF(ca).run().results.rmsf / 10.0  # nm
        resids, resnames = ca.resids, ca.resnames

        df = pd.DataFrame({"resid": resids, "resname": resnames, "rmsf_nm": rmsf})
        ctx.write_csv(df, "rmsf.csv")
        thr = rmsf.mean() + rmsf.std()
        flex = df[df["rmsf_nm"] > thr]

        fig, ax = plotting.new_axes(figsize=(8.2, 4.6))
        ax.plot(resids, rmsf, color=PALETTE["primary"], lw=1.8)
        ax.fill_between(resids, 0, rmsf, color=PALETTE["primary"], alpha=0.15)
        ax.axhline(thr, ls="--", lw=1.1, color=PALETTE["accent"],
                   label=f"mean + 1$\\sigma$ = {thr:.3f} nm")
        ax.scatter(flex["resid"], flex["rmsf_nm"], color=PALETTE["secondary"],
                   s=28, zorder=5, label="Flexible")
        ax.set_xlabel("Residue number"); ax.set_ylabel("RMSF (nm)")
        ax.set_title(self.label); ax.legend(loc="upper right")
        ax.margins(x=0.01); ax.set_ylim(bottom=0)
        plotting.save_figure(fig, ctx.fig_path("rmsf"), dpi=ctx.config.dpi)

        return {"rmsf": st.describe(rmsf, "rmsf_nm"),
                "n_flexible": int(len(flex)),
                "flexible_resids": flex["resid"].tolist(),
                "most_flexible_resid": int(resids[np.argmax(rmsf)]),
                "figure": "rmsf"}
