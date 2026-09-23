"""Dynamic cross-correlation matrix of Cα fluctuations."""

from __future__ import annotations

import numpy as np
import pandas as pd
from MDAnalysis.analysis import align

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx.plotting import DIV_CMAP


class DCCM(BaseAnalysis):
    name = "dccm"
    label = "Dynamic cross-correlation matrix"
    category = "dynamics"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/dccm.csv", "figures/dccm.png"]

    def run(self, ctx) -> dict:
        plotting.set_style()
        u = ctx.core_universe()
        align.AlignTraj(u, u, select="protein and name CA",
                        ref_frame=ctx.config.start or 0, in_memory=True).run()
        ca = u.select_atoms("protein and name CA")
        resids = ca.resids
        coords = np.array([ca.positions.copy() for _ in u.trajectory]) / 10.0
        disp = coords - coords.mean(axis=0)
        cov = np.einsum("tia,tja->ij", disp, disp) / disp.shape[0]
        norm = np.sqrt(np.diag(cov))
        dccm = np.clip(cov / np.outer(norm, norm), -1, 1)

        pd.DataFrame(dccm, index=resids, columns=resids).to_csv(ctx.csv_path("dccm.csv"))
        fig, ax = plotting.new_axes(figsize=(6.6, 5.6))
        im = ax.imshow(dccm, cmap=DIV_CMAP, vmin=-1, vmax=1, origin="lower",
                       extent=[resids[0], resids[-1], resids[0], resids[-1]])
        ax.set_xlabel("Residue"); ax.set_ylabel("Residue")
        ax.set_title(r"Dynamic cross-correlation (C$\alpha$)")
        fig.colorbar(im, ax=ax, shrink=0.85, label="Correlation")
        plotting.save_figure(fig, ctx.fig_path("dccm"), dpi=ctx.config.dpi)

        off = dccm[~np.eye(len(resids), dtype=bool)]
        return {"mean_correlation": float(off.mean()),
                "frac_strong_positive": float((off > 0.5).mean()),
                "frac_strong_negative": float((off < -0.5).mean()),
                "figure": "dccm"}
