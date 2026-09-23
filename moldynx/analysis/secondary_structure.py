"""DSSP secondary-structure timeline and content (mdtraj, no external binary)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import mdtraj as md
from matplotlib.colors import ListedColormap, BoundaryNorm

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE

_MAP = {"C": 0, "H": 1, "E": 2}
_LABELS = ["Coil", "Helix", "Sheet"]
_COLORS = [PALETTE["coil"], PALETTE["helix"], PALETTE["sheet"]]


class SecondaryStructure(BaseAnalysis):
    name = "secondary_structure"
    label = "DSSP secondary structure"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/ss_fractions.csv", "figures/ss_timeline.png"]

    def run(self, ctx) -> dict:
        plotting.set_style()
        ctx.core_universe()
        traj = md.load(str(ctx.config.data_dir / "core.xtc"),
                       top=str(ctx.config.data_dir / "core.pdb"))
        times = traj.time / 1000.0
        dssp = md.compute_dssp(traj, simplified=True)
        resids = np.array([r.resSeq for r in traj.topology.residues])
        fh = (dssp == "H").mean(axis=1) * 100
        fe = (dssp == "E").mean(axis=1) * 100
        fc = (dssp == "C").mean(axis=1) * 100
        ctx.write_csv(pd.DataFrame({"time_ns": times, "helix_pct": fh,
                                    "sheet_pct": fe, "coil_pct": fc}), "ss_fractions.csv")
        coded = np.vectorize(_MAP.get)(dssp).astype(np.int8)
        pd.DataFrame(coded, index=np.round(times, 3), columns=resids).to_csv(
            ctx.csv_path("ss_timeline.csv"))

        fig, ax = plotting.new_axes()
        for y, lab, c in ((fh, "Helix", PALETTE["helix"]), (fe, "Sheet", PALETTE["sheet"]),
                          (fc, "Coil", PALETTE["coil"])):
            ax.plot(times, st.moving_average(y, 20), color=c, lw=2.0, label=lab)
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Content (%)")
        ax.set_title("Secondary-structure content"); ax.legend(ncol=3, loc="best")
        plotting.save_figure(fig, ctx.fig_path("ss_fractions"), dpi=ctx.config.dpi)

        cmap = ListedColormap(_COLORS); norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
        fig, ax = plotting.new_axes(figsize=(8.4, 4.8))
        im = ax.imshow(coded.T, aspect="auto", origin="lower", cmap=cmap, norm=norm,
                       extent=[times[0], times[-1], resids[0], resids[-1]],
                       interpolation="nearest")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Residue")
        ax.set_title("Secondary-structure timeline")
        cb = fig.colorbar(im, ax=ax, ticks=[0, 1, 2], shrink=0.85)
        cb.ax.set_yticklabels(_LABELS)
        plotting.save_figure(fig, ctx.fig_path("ss_timeline"), dpi=ctx.config.dpi)

        return {"mean_helix_pct": float(fh.mean()), "mean_sheet_pct": float(fe.mean()),
                "mean_coil_pct": float(fc.mean()), "figure": "ss_timeline"}
