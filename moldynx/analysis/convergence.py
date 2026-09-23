"""Convergence diagnostics: plateau, block averaging, halves, essential-dynamics RMSIP."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as sps
from MDAnalysis.analysis import align, rms

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class Convergence(BaseAnalysis):
    name = "convergence"
    label = "Convergence / equilibration"
    category = "quality"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    order = 200  # runs after rmsd/rog so it can reuse their CSVs
    outputs = ["results/convergence_summary.csv", "figures/convergence_running.png"]
    default_params = {"n_blocks": 10, "n_modes": 10}

    def _series(self, ctx):
        rc, gc = ctx.csv_path("rmsd.csv"), ctx.csv_path("rog.csv")
        if rc.exists() and gc.exists():
            r, g = pd.read_csv(rc), pd.read_csv(gc)
            return (r["time_ns"].to_numpy(), r["rmsd_backbone_nm"].to_numpy(),
                    g["rg_nm"].to_numpy())
        u = ctx.core_universe()
        protein = u.select_atoms("protein")
        R = rms.RMSD(u, select="backbone", ref_frame=ctx.config.start or 0).run()
        rmsd = R.results.rmsd[:, 2] / 10.0
        rg, t = [], []
        for ts in u.trajectory:
            rg.append(protein.radius_of_gyration() / 10.0); t.append(ts.time / 1000.0)
        return np.array(t), rmsd, np.array(rg)

    def _rmsip(self, ctx, k):
        u = ctx.core_universe()
        align.AlignTraj(u, u, select="protein and name CA",
                        ref_frame=ctx.config.start or 0, in_memory=True).run()
        ca = u.select_atoms("protein and name CA")
        X = np.array([ca.positions.copy().ravel() for _ in u.trajectory]) / 10.0
        half = len(X) // 2

        def modes(Xh):
            Xc = Xh - Xh.mean(axis=0)
            ev, evec = np.linalg.eigh((Xc.T @ Xc) / (len(Xh) - 1))
            return evec[:, ::-1][:, :k]
        A, B = modes(X[:half]), modes(X[half:])
        ov = A.T @ B
        return float(np.sqrt((ov ** 2).sum() / k)), ov

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        times, rmsd, rg = self._series(ctx)
        pr, pg = st.plateau_detection(times, rmsd), st.plateau_detection(times, rg)
        half = len(times) // 2
        t_rmsd = sps.ttest_ind(rmsd[:half], rmsd[half:], equal_var=False)
        sem = st.block_average_sem(rmsd, p["n_blocks"])
        rmsip, ov = self._rmsip(ctx, p["n_modes"])

        fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True)
        a1.plot(times, rmsd, color=PALETTE["primary"], alpha=0.3, lw=0.8)
        a1.plot(times, st.running_mean(rmsd), color=PALETTE["primary"], lw=2.2)
        a1.set_ylabel("Backbone RMSD (nm)"); a1.set_title("Cumulative running averages")
        a2.plot(times, rg, color=PALETTE["secondary"], alpha=0.3, lw=0.8)
        a2.plot(times, st.running_mean(rg), color=PALETTE["secondary"], lw=2.2)
        a2.set_ylabel(r"$R_g$ (nm)"); a2.set_xlabel("Time (ns)")
        for a in (a1, a2):
            a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
        plotting.save_figure(fig, ctx.fig_path("convergence_running"), dpi=ctx.config.dpi)

        fig, ax = plotting.new_axes(figsize=(5.6, 4.8))
        im = ax.imshow(np.abs(ov), cmap="viridis", origin="lower", vmin=0, vmax=1,
                       extent=[1, ov.shape[1], 1, ov.shape[0]])
        ax.set_xlabel("PC (2nd half)"); ax.set_ylabel("PC (1st half)")
        ax.set_title(f"Essential-dynamics overlap (RMSIP = {rmsip:.2f})")
        fig.colorbar(im, ax=ax, label="|inner product|")
        plotting.save_figure(fig, ctx.fig_path("essential_dynamics_overlap"), dpi=ctx.config.dpi)

        summ = {"rmsd_converged": pr["converged"],
                "rmsd_tail_slope_nm_per_ns": pr["tail_slope"],
                "rmsd_drift_2nd_minus_1st_nm": pr["drift_second_minus_first"],
                "rg_drift_2nd_minus_1st_nm": pg["drift_second_minus_first"],
                "rmsd_block_sem_nm": sem, "rmsd_halves_ttest_p": float(t_rmsd.pvalue),
                "essential_dynamics_rmsip": rmsip, "figure": "convergence_running"}
        pd.DataFrame([summ]).to_csv(ctx.csv_path("convergence_summary.csv"), index=False)
        return summ
