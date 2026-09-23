"""Residue interaction network (RIN) + centrality / hub analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.distance import cdist

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx.plotting import PALETTE

_OCC = 0.5


class ResidueInteractionNetwork(BaseAnalysis):
    name = "rin"
    label = "Residue interaction network"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    order = 150  # after contact_map (reuses its occupancy CSV when present)
    outputs = ["results/rin_node_metrics.csv", "figures/rin_network.png"]
    default_params = {"cutoff": 8.0}

    def _adjacency(self, ctx):
        occ_csv = ctx.csv_path("contact_occupancy.csv")
        if occ_csv.exists():
            df = pd.read_csv(occ_csv, index_col=0)
            return df.columns.astype(int).to_numpy(), df.to_numpy()
        u = ctx.core_universe()
        ca = u.select_atoms("protein and name CA")
        resids = ca.resids
        u.trajectory[ctx.config.start or 0]
        D = cdist(ca.positions, ca.positions)
        seq = np.abs(resids[:, None] - resids[None, :]) > 2
        return resids, ((D < self.default_params["cutoff"]) & seq).astype(float)

    def run(self, ctx) -> dict:
        plotting.set_style()
        resids, occ = self._adjacency(ctx)
        N = len(resids)
        G = nx.Graph(); G.add_nodes_from(resids.tolist())
        for i in range(N):
            for j in range(i + 1, N):
                if occ[i, j] >= _OCC:
                    G.add_edge(int(resids[i]), int(resids[j]), weight=float(occ[i, j]))

        deg = dict(G.degree())
        betw = nx.betweenness_centrality(G, weight="weight")
        clos = nx.closeness_centrality(G)
        node_df = pd.DataFrame({"resid": list(G.nodes()),
                                "degree": [deg[n] for n in G.nodes()],
                                "betweenness": [betw[n] for n in G.nodes()],
                                "closeness": [clos[n] for n in G.nodes()]}).sort_values(
            "betweenness", ascending=False).reset_index(drop=True)
        ctx.write_csv(node_df, "rin_node_metrics.csv")

        fig, ax = plotting.new_axes(figsize=(7.2, 6.4))
        if G.number_of_edges():
            pos = nx.spring_layout(G, seed=42, weight="weight", k=0.3)
            nd = np.array([deg[n] for n in G.nodes()])
            nb = np.array([betw[n] for n in G.nodes()])
            nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.25, width=0.8)
            nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_size=30 + 12*nd,
                                           node_color=nb, cmap="plasma")
            fig.colorbar(nodes, ax=ax, label="Betweenness", shrink=0.8)
        ax.set_axis_off(); ax.set_title("Residue interaction network")
        plotting.save_figure(fig, ctx.fig_path("rin_network"), dpi=ctx.config.dpi)

        return {"n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges(),
                "top_hub_resid": int(node_df.iloc[0]["resid"]) if len(node_df) else None,
                "graph_density": float(nx.density(G)), "figure": "rin_network"}
