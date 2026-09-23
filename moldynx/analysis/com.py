"""Centre-of-mass tracking (universal): displacement, path length, 3D trajectory."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)

from moldynx.core.base import BaseAnalysis
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class CenterOfMass(BaseAnalysis):
    name = "com"
    label = "Centre-of-mass tracking"
    category = "dynamics"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    outputs = ["results/com.csv", "figures/com_xyz.png", "figures/com_3d.png"]

    def run(self, ctx) -> dict:
        plotting.set_style()
        u = ctx.core_universe()
        solute = u.select_atoms("protein" if ctx.system.flags.get("has_protein", True) else "all")
        if solute.n_atoms == 0:
            solute = u.atoms
        n = len(u.trajectory)
        coms = np.empty((n, 3)); boxes = np.empty((n, 3)); times = np.empty(n)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[com]")):
            coms[i] = solute.center_of_mass()
            boxes[i] = ts.dimensions[:3]
            times[i] = ts.time / 1000.0

        coms_nm, box_nm = coms / 10.0, boxes / 10.0
        disp = np.diff(coms_nm, axis=0)
        disp -= box_nm[1:] * np.round(disp / np.where(box_nm[1:] == 0, 1, box_nm[1:]))
        cont = np.vstack([coms_nm[0], coms_nm[0] + np.cumsum(disp, axis=0)])
        rel = cont - cont[0]
        cumdist = np.concatenate([[0.0], np.cumsum(np.linalg.norm(disp, axis=1))])
        net = np.linalg.norm(rel, axis=1)

        ctx.write_csv(pd.DataFrame({"time_ns": times, "com_x_nm": rel[:, 0],
                                    "com_y_nm": rel[:, 1], "com_z_nm": rel[:, 2],
                                    "net_displacement_nm": net,
                                    "cumulative_distance_nm": cumdist}), "com.csv")

        fig, ax = plotting.new_axes()
        for k, c in zip("xyz", (PALETTE["primary"], PALETTE["secondary"], PALETTE["green"])):
            ax.plot(times, rel[:, "xyz".index(k)], color=c, lw=1.6, label=k.upper())
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("COM displacement (nm)")
        ax.set_title("Protein centre-of-mass displacement")
        ax.legend(ncol=3, loc="upper left"); ax.margins(x=0.01)
        plotting.save_figure(fig, ctx.fig_path("com_xyz"), dpi=ctx.config.dpi)

        fig = plt.figure(figsize=(6.4, 5.6))
        ax = fig.add_subplot(111, projection="3d")
        sc = ax.scatter(rel[:, 0], rel[:, 1], rel[:, 2], c=times, cmap="viridis", s=6)
        ax.plot(rel[:, 0], rel[:, 1], rel[:, 2], color=PALETTE["muted"], lw=0.6, alpha=0.6)
        ax.set_xlabel("X (nm)"); ax.set_ylabel("Y (nm)"); ax.set_zlabel("Z (nm)")
        ax.set_title("COM 3D trajectory")
        fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1, label="Time (ns)")
        plotting.save_figure(fig, ctx.fig_path("com_3d"), dpi=ctx.config.dpi)

        return {"total_distance_nm": float(cumdist[-1]),
                "net_displacement_nm": float(net[-1]),
                "displacement": st.describe(net, "net_displacement_nm"),
                "figure": "com_xyz"}
