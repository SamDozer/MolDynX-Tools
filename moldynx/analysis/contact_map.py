"""Residue-residue contact occupancy + time-resolved change map."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE, SEQ_CMAP, DIV_CMAP


class ContactMap(BaseAnalysis):
    name = "contact_map"
    label = "Residue-residue contact map"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/contact_occupancy.csv", "figures/contact_occupancy.png"]
    default_params = {"cutoff": 8.0}

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        ca = u.select_atoms("protein and name CA")
        resids = ca.resids; N = len(resids); cutoff = p["cutoff"]
        seq = np.abs(resids[:, None] - resids[None, :]) > 2

        occ = np.zeros((N, N)); occ1 = np.zeros((N, N)); occ2 = np.zeros((N, N))
        n = len(u.trajectory); half = n // 2
        counts = np.empty(n); times = np.empty(n)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[contact_map]")):
            c = (cdist(ca.positions, ca.positions) < cutoff) & seq
            occ += c
            (occ1 if i < half else occ2)[:] += c
            counts[i] = c.sum() / 2.0; times[i] = ts.time / 1000.0
        occ /= n; occ1 /= max(half, 1); occ2 /= max(n - half, 1)
        diff = occ2 - occ1

        pd.DataFrame(occ, index=resids, columns=resids).to_csv(
            ctx.csv_path("contact_occupancy.csv"))
        ctx.write_csv(pd.DataFrame({"time_ns": times, "n_contacts": counts}),
                      "contacts_vs_time.csv")

        fig, ax = plotting.new_axes(figsize=(6.4, 5.6))
        im = ax.imshow(occ, cmap=SEQ_CMAP, origin="lower", vmin=0, vmax=1,
                       extent=[resids[0], resids[-1], resids[0], resids[-1]])
        ax.set_xlabel("Residue"); ax.set_ylabel("Residue")
        ax.set_title("Contact occupancy")
        fig.colorbar(im, ax=ax, shrink=0.85, label="Occupancy")
        plotting.save_figure(fig, ctx.fig_path("contact_occupancy"), dpi=ctx.config.dpi)

        fig, ax = plotting.new_axes(figsize=(6.4, 5.6))
        im = ax.imshow(diff, cmap=DIV_CMAP, origin="lower", vmin=-1, vmax=1,
                       extent=[resids[0], resids[-1], resids[0], resids[-1]])
        ax.set_xlabel("Residue"); ax.set_ylabel("Residue")
        ax.set_title("Contact change (2nd - 1st half)")
        fig.colorbar(im, ax=ax, shrink=0.85, label=r"$\Delta$ occupancy")
        plotting.save_figure(fig, ctx.fig_path("contact_difference"), dpi=ctx.config.dpi)

        return {"mean_contacts": float(counts.mean()),
                "n_persistent_pairs": int((occ >= 0.5).sum() // 2),
                "figure": "contact_occupancy"}
