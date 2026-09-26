"""
Porcupine plot of the dominant motion (PC1 of the Cα covariance) + a mode PDB.

Uses its own copy of the solute trajectory (alignment in memory would otherwise
rotate the shared trajectory, and its box, under the other analyses). PC1 comes
from an SVD of the aligned, centred coordinates -- no 3N × 3N covariance matrix,
so large proteins stay tractable. Arrows are coloured by chain.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from moldynx import plotting
from moldynx.core.base import BaseAnalysis
from moldynx.plotting import PALETTE

ARROW_SCALE = 2.5


class Porcupine(BaseAnalysis):
    name = "porcupine"
    label = "Porcupine plot of PC1"
    category = "dynamics"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    order = 70
    default_params = {"n_models": 20, "mode_scale": 2.0}
    outputs = ["results/porcupine_pc1.csv", "results/porcupine_pc1.pdb",
               "figures/porcupine_pc1.png"]

    def run(self, ctx) -> dict:
        import MDAnalysis as mda
        from MDAnalysis.analysis import align
        p = self.params(ctx)
        plotting.set_style()
        ctx.core_universe()
        u = mda.Universe(str(ctx.config.data_dir / "core.pdb"), str(ctx.config.data_dir / "core.xtc"))
        ca = u.select_atoms("protein and name CA")
        if ca.n_atoms < 3:
            return {"status": "skipped", "reason": "fewer than three protein Cα atoms"}
        align.AlignTraj(u, u, select="protein and name CA", ref_frame=0, in_memory=True).run()
        X = np.array([ca.positions.copy() for _ in u.trajectory])
        mean = X.mean(0)
        Xc = (X - mean).reshape(len(X), -1)
        _, s, vt = np.linalg.svd(Xc, full_matrices=False)
        var = s ** 2 / max(len(X) - 1, 1)
        pc1 = vt[0].reshape(-1, 3)
        amp = float(np.sqrt(var[0]))
        disp = pc1 * amp
        mag = np.linalg.norm(disp, axis=1)

        labels = np.array(["protein"] * ca.n_atoms, dtype=object)
        for rec in ctx.core_meta.get("chains", []):
            m = (ca.indices >= rec["core_start"]) & (ca.indices < rec["core_stop"])
            labels[m] = rec.get("segid") or f"chain {rec.get('index')}"
        ctx.write_csv(pd.DataFrame({"resid": ca.resids, "chain": labels, "dx_A": disp[:, 0],
                                    "dy_A": disp[:, 1], "dz_A": disp[:, 2], "magnitude_A": mag}),
                      "porcupine_pc1.csv")
        fig = plotting.style.plt.figure(figsize=(7.4, 6.6))
        ax = fig.add_subplot(111, projection="3d")
        ax.plot(*mean.T, color=PALETTE["muted"], lw=0.6, alpha=0.5)
        colours = [PALETTE[c] for c in ("primary", "green", "secondary", "purple")]
        for k, lab in enumerate(dict.fromkeys(labels)):
            m = labels == lab
            d = disp[m] * ARROW_SCALE
            ax.quiver(*mean[m].T, *d.T, color=colours[k % 4], linewidth=1.0, label=lab)
        ax.set_title(f"PC1 ({var[0] / var.sum() * 100:.1f}% of variance)")
        ax.legend(loc="upper left")
        plotting.save_figure(fig, ctx.fig_path("porcupine_pc1"), dpi=ctx.config.dpi)
        with mda.Writer(str(ctx.csv_path("porcupine_pc1.pdb")), n_atoms=ca.n_atoms,
                        multiframe=True) as w:
            for ph in np.linspace(0, 2 * np.pi, p["n_models"], endpoint=False):
                ca.positions = (mean + p["mode_scale"] * amp * np.sin(ph) * pc1).astype(np.float32)
                w.write(ca)
        per_chain = {lab: float(mag[labels == lab].mean()) for lab in dict.fromkeys(labels)}
        return {"pc1_variance_fraction": float(var[0] / var.sum()), "pc1_rms_amplitude_A": amp,
                "max_displacement_resid": int(ca.resids[int(np.argmax(mag))]),
                "max_displacement_A": float(mag.max()), "mean_displacement_by_chain_A": per_chain,
                "figure": "porcupine_pc1"}
