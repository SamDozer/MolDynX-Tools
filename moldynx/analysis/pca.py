"""PCA (essential dynamics) on Cα + PC1/PC2 free-energy landscape."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx.plotting import PALETTE, DIV_CMAP, FEL_CMAP

KB_KJ = 0.00831446261815324


class PCA(BaseAnalysis):
    name = "pca"
    label = "PCA + free-energy landscape"
    category = "dynamics"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/pca_projections.csv", "figures/pca_fel.png"]
    default_params = {"temperature": 310.0, "bins": 80}

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        align.AlignTraj(u, u, select="protein and name CA",
                        ref_frame=ctx.config.start or 0, in_memory=True).run()
        ca = u.select_atoms("protein and name CA")
        X, times = [], []
        for ts in u.trajectory:
            X.append(ca.positions.copy().ravel()); times.append(ts.time / 1000.0)
        X = np.asarray(X) / 10.0; times = np.asarray(times)
        Xc = X - X.mean(axis=0); n = len(X)
        evals, evecs = np.linalg.eigh((Xc.T @ Xc) / (n - 1))
        order = np.argsort(evals)[::-1]
        evals, evecs = np.clip(evals[order], 0, None), evecs[:, order]
        var = evals / evals.sum(); cum = np.cumsum(var)
        proj = Xc @ evecs[:, :5]

        pd.DataFrame({"component": np.arange(1, len(evals)+1), "eigenvalue_nm2": evals,
                      "variance_explained": var, "cumulative_variance": cum}).to_csv(
            ctx.csv_path("pca_eigenvalues.csv"), index=False)
        pd.DataFrame({"time_ns": times, **{f"PC{k+1}": proj[:, k] for k in range(5)}}).to_csv(
            ctx.csv_path("pca_projections.csv"), index=False)

        # variance
        fig, ax = plotting.new_axes()
        kk = min(20, len(evals))
        ax.bar(np.arange(1, kk+1), var[:kk]*100, color=PALETTE["primary"], alpha=0.85)
        ax2 = ax.twinx(); ax2.plot(np.arange(1, kk+1), cum[:kk]*100, "-o",
                                   color=PALETTE["secondary"], ms=4)
        ax2.set_ylabel("Cumulative (%)"); ax2.set_ylim(0, 105); ax2.spines["top"].set_visible(False)
        ax.set_xlabel("Principal component"); ax.set_ylabel("Variance (%)")
        ax.set_title("PCA variance explained")
        plotting.save_figure(fig, ctx.fig_path("pca_variance"), dpi=ctx.config.dpi)

        # FEL
        H, xe, ye = np.histogram2d(proj[:, 0], proj[:, 1], bins=p["bins"], density=True)
        P = (H / H.sum()).T
        with np.errstate(divide="ignore"):
            G = -KB_KJ * p["temperature"] * np.log(P)
        G[~np.isfinite(G)] = np.nan; G -= np.nanmin(G)
        xc = 0.5*(xe[:-1]+xe[1:]); yc = 0.5*(ye[:-1]+ye[1:])
        fig, ax = plotting.new_axes(figsize=(6.8, 5.4))
        im = ax.pcolormesh(xc, yc, G, cmap=FEL_CMAP, shading="auto")
        fig.colorbar(im, ax=ax, label=r"$\Delta G$ (kJ/mol)")
        ax.set_xlabel("PC1 (nm)"); ax.set_ylabel("PC2 (nm)")
        ax.set_title(f"Free-energy landscape (T = {p['temperature']:.0f} K)")
        plotting.save_figure(fig, ctx.fig_path("pca_fel"), dpi=ctx.config.dpi)

        return {"pc1_variance": float(var[0]), "pc2_variance": float(var[1]),
                "n_pc_for_80pct": int(np.searchsorted(cum, 0.80) + 1),
                "figure": "pca_fel"}
