"""Ramachandran (phi/psi) analysis: before vs after MD + ensemble density."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis.dihedrals import Ramachandran as MDARama

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx.plotting import PALETTE


def _favoured_fraction(phi, psi) -> float:
    """Rough favoured-region fraction (alpha + beta + PPII + left-handed boxes)."""
    phi = np.asarray(phi); psi = np.asarray(psi)
    alpha = (phi > -160) & (phi < -30) & (psi > -80) & (psi < 60)
    beta = (phi > -180) & (phi < -40) & (psi > 90) & (psi < 180)
    ppII = (phi > -90) & (phi < -40) & (psi > 120) & (psi < 180)
    left = (phi > 30) & (phi < 90) & (psi > 0) & (psi < 90)
    return float((alpha | beta | ppII | left).mean())


class Ramachandran(BaseAnalysis):
    name = "ramachandran"
    label = "Ramachandran (before/after MD)"
    category = "quality"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/ramachandran_first_last.csv",
               "figures/ramachandran_before_after.png",
               "figures/ramachandran_density.png"]

    def run(self, ctx) -> dict:
        plotting.set_style()
        u = ctx.core_universe()
        protein = u.select_atoms("protein")
        rama = MDARama(protein).run()
        angles = rama.results.angles          # (n_frames, n_res, 2), degrees
        first, last = angles[0], angles[-1]
        alla = angles.reshape(-1, 2)

        pd.DataFrame({"phi_first": first[:, 0], "psi_first": first[:, 1],
                      "phi_last": last[:, 0], "psi_last": last[:, 1]}).to_csv(
            ctx.csv_path("ramachandran_first_last.csv"), index=False)

        fav_first = _favoured_fraction(first[:, 0], first[:, 1])
        fav_last = _favoured_fraction(last[:, 0], last[:, 1])
        fav_all = _favoured_fraction(alla[:, 0], alla[:, 1])

        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 4.8))
        for ax, ang, title, frac in ((a1, first, "Before MD (start)", fav_first),
                                     (a2, last, "After MD (end)", fav_last)):
            ax.scatter(ang[:, 0], ang[:, 1], s=18, color=PALETTE["primary"], alpha=0.7)
            ax.axhline(0, color="grey", lw=0.6); ax.axvline(0, color="grey", lw=0.6)
            ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
            ax.set_xticks(range(-180, 181, 90)); ax.set_yticks(range(-180, 181, 90))
            ax.set_xlabel(r"$\phi$ (deg)"); ax.set_ylabel(r"$\psi$ (deg)")
            ax.set_title(f"{title}\nfavoured ~ {frac*100:.0f}%")
        fig.suptitle("Ramachandran: before vs after MD", fontweight="bold")
        plotting.save_figure(fig, ctx.fig_path("ramachandran_before_after"), dpi=ctx.config.dpi)

        fig, ax = plotting.new_axes(figsize=(5.8, 5.2))
        h = ax.hist2d(alla[:, 0], alla[:, 1], bins=120, cmap="viridis",
                      range=[[-180, 180], [-180, 180]], cmin=1)
        ax.set_xlabel(r"$\phi$ (deg)"); ax.set_ylabel(r"$\psi$ (deg)")
        ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
        ax.set_title("Ramachandran density (all frames)")
        fig.colorbar(h[3], ax=ax, label="Count")
        plotting.save_figure(fig, ctx.fig_path("ramachandran_density"), dpi=ctx.config.dpi)

        return {"favoured_fraction_before": fav_first,
                "favoured_fraction_after": fav_last,
                "favoured_fraction_ensemble": fav_all,
                "figure": "ramachandran_before_after"}
