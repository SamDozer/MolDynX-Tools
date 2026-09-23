"""Conformational clustering on pairwise Cα RMSD (optimal k via silhouette)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import mdtraj as md
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from tqdm import tqdm

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting

_COLORS = ["#1f6f8b", "#e07a5f", "#5f9e6e", "#8367c7", "#d4a017", "#c0392b",
           "#2471a3", "#7f8c8d"]


class Clustering(BaseAnalysis):
    name = "clustering"
    label = "Conformational clustering"
    category = "dynamics"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/cluster_populations.csv", "figures/cluster_populations.png"]
    default_params = {"k_range": (2, 3, 4, 5, 6, 7, 8)}

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        ctx.core_universe()
        traj = md.load(str(ctx.config.data_dir / "core.xtc"),
                       top=str(ctx.config.data_dir / "core.pdb"))
        traj.superpose(traj, 0)
        ca = traj.topology.select("name CA")
        nfr = traj.n_frames
        D = np.empty((nfr, nfr), dtype=np.float32)
        for i in tqdm(range(nfr), desc="[clustering] RMSD matrix", unit="frame"):
            D[i] = md.rmsd(traj, traj, i, atom_indices=ca)
        D = 0.5 * (D + D.T).astype(np.float64)
        np.fill_diagonal(D, 0.0); D = np.clip(D, 0, None)

        scores, labels_by_k = {}, {}
        for k in p["k_range"]:
            if k >= nfr:
                continue
            lbl = AgglomerativeClustering(n_clusters=k, metric="precomputed",
                                          linkage="average").fit_predict(D)
            scores[k] = float(silhouette_score(D, lbl, metric="precomputed"))
            labels_by_k[k] = lbl
        best_k = max(scores, key=scores.get); labels = labels_by_k[best_k]
        order = pd.Series(labels).value_counts().index.tolist()
        remap = {o: n for n, o in enumerate(order)}
        labels = np.array([remap[l] for l in labels])

        struct_dir = ctx.config.results_dir / "cluster_structures"
        struct_dir.mkdir(parents=True, exist_ok=True)
        pops, medoids = [], []
        for c in range(best_k):
            members = np.where(labels == c)[0]
            medoid = members[np.argmin(D[np.ix_(members, members)].mean(axis=1))]
            medoids.append(int(medoid)); pops.append(len(members))
            traj[medoid].save_pdb(str(struct_dir / f"cluster_{c+1:02d}.pdb"))
        pop_df = pd.DataFrame({"cluster": np.arange(1, best_k+1), "population": pops,
                               "percent": 100*np.array(pops)/nfr, "medoid_frame": medoids})
        ctx.write_csv(pop_df, "cluster_populations.csv")
        pd.DataFrame({"k": list(scores), "silhouette": list(scores.values())}).to_csv(
            ctx.csv_path("cluster_silhouette.csv"), index=False)

        fig, ax = plotting.new_axes()
        colors = [_COLORS[i % len(_COLORS)] for i in range(best_k)]
        ax.bar(pop_df["cluster"], pop_df["percent"], color=colors)
        for x, pc in zip(pop_df["cluster"], pop_df["percent"]):
            ax.text(x, pc + 0.5, f"{pc:.0f}%", ha="center", fontsize=10)
        ax.set_xlabel("Cluster"); ax.set_ylabel("Population (%)")
        ax.set_title(f"Conformational clusters (k = {best_k})")
        plotting.save_figure(fig, ctx.fig_path("cluster_populations"), dpi=ctx.config.dpi)

        return {"optimal_k": int(best_k), "silhouette": scores[best_k],
                "dominant_cluster_percent": float(pop_df["percent"].iloc[0]),
                "figure": "cluster_populations"}
