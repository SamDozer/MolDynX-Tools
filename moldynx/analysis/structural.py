"""Global shape descriptors: end-to-end, Dmax, gyration-tensor anisotropy, volume."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


def _gyration(pos):
    r = pos - pos.mean(axis=0)
    S = (r[:, :, None] * r[:, None, :]).mean(axis=0)
    ev = np.sort(np.linalg.eigvalsh(S))
    l1, l2, l3 = ev
    tr = ev.sum()
    rg = float(np.sqrt(tr))
    asph = float(l3 - 0.5 * (l1 + l2))
    kappa2 = float(1 - 3 * (l1*l2 + l2*l3 + l3*l1) / (tr*tr)) if tr > 0 else 0.0
    return rg, asph, kappa2


class StructuralDescriptors(BaseAnalysis):
    name = "structural"
    label = "Global shape descriptors"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/structural.csv", "figures/structural_descriptors.png"]

    def run(self, ctx) -> dict:
        plotting.set_style()
        u = ctx.core_universe()
        protein = u.select_atoms("protein")
        ca = u.select_atoms("protein and name CA")
        n = len(u.trajectory)
        e2e = np.empty(n); dmax = np.empty(n); rg = np.empty(n)
        asph = np.empty(n); kappa2 = np.empty(n); vol = np.empty(n); times = np.empty(n)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[structural]")):
            times[i] = ts.time / 1000.0
            cap = ca.positions / 10.0; ap = protein.positions / 10.0
            e2e[i] = np.linalg.norm(cap[0] - cap[-1])
            dmax[i] = pdist(cap).max()
            rg[i], asph[i], kappa2[i] = _gyration(ap)
            try:
                vol[i] = ConvexHull(ap, qhull_options="QJ").volume
            except Exception:
                vol[i] = np.nan

        df = pd.DataFrame({"time_ns": times, "end_to_end_nm": e2e, "dmax_nm": dmax,
                           "rg_tensor_nm": rg, "asphericity_nm2": asph,
                           "rel_shape_anisotropy": kappa2, "hull_volume_nm3": vol,
                           "elongation": dmax / (2 * rg)})
        ctx.write_csv(df, "structural.csv")

        fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0))
        for ax, (col, lab, c) in zip(axes.ravel(), [
                ("end_to_end_nm", "End-to-end (nm)", PALETTE["primary"]),
                ("dmax_nm", r"$D_{max}$ (nm)", PALETTE["secondary"]),
                ("rel_shape_anisotropy", r"$\kappa^2$", PALETTE["purple"]),
                ("hull_volume_nm3", "Volume (nm$^3$)", PALETTE["green"])]):
            y = df[col].to_numpy()
            ax.plot(times, y, color=c, alpha=0.4, lw=1.0)
            ax.plot(times, st.moving_average(y, 20), color=c, lw=2.0)
            ax.set_xlabel("Time (ns)"); ax.set_ylabel(lab)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        fig.suptitle("Global shape descriptors", fontweight="bold")
        plotting.save_figure(fig, ctx.fig_path("structural_descriptors"), dpi=ctx.config.dpi)

        return {"dmax": st.describe(dmax, "dmax_nm"),
                "rel_shape_anisotropy": st.describe(kappa2, "kappa2"),
                "hull_volume": st.describe(vol[np.isfinite(vol)], "hull_volume_nm3"),
                "figure": "structural_descriptors"}
