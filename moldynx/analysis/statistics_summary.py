"""Cross-metric summary table + Pearson/Spearman correlations (runs last)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from moldynx.core.base import BaseAnalysis
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import DIV_CMAP

# (csv file, time col, value col, label)
_SOURCES = [
    ("rmsd.csv", "time_ns", "rmsd_backbone_nm", "RMSD"),
    ("rog.csv", "time_ns", "rg_nm", "Rg"),
    ("sasa_total.csv", "time_ns", "sasa_total_nm2", "SASA"),
    ("hbonds_timeseries.csv", "time_ns", "hbonds_protein_protein", "HB_intra"),
    ("energies.csv", "time_ns", "Potential", "E_pot"),
    ("com.csv", "time_ns", "net_displacement_nm", "COM_disp"),
    ("native_contacts.csv", "time_ns", "Q", "Q_native"),
    ("ss_fractions.csv", "time_ns", "helix_pct", "Helix"),
]


class StatisticsSummary(BaseAnalysis):
    name = "statistics_summary"
    label = "Cross-metric statistics & correlations"
    category = "quality"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    order = 300  # last: reads other analyses' CSVs
    outputs = ["tables/summary_statistics.csv", "figures/correlation_pearson.png"]

    def run(self, ctx) -> dict:
        plotting.set_style()
        master = None
        for f in ("rog.csv", "rmsd.csv"):
            if ctx.csv_path(f).exists():
                master = pd.read_csv(ctx.csv_path(f))["time_ns"].to_numpy(); break
        if master is None:
            return {"note": "no base time-series found"}

        merged = {"time_ns": master}
        for fname, tcol, vcol, label in _SOURCES:
            pth = ctx.csv_path(fname)
            if not pth.exists():
                continue
            df = pd.read_csv(pth)
            if tcol in df and vcol in df:
                good = np.isfinite(df[tcol]) & np.isfinite(df[vcol])
                if good.sum() >= 2:
                    merged[label] = np.interp(master, df[tcol][good], df[vcol][good])
        mdf = pd.DataFrame(merged)
        cols = [c for c in mdf.columns if c != "time_ns"]
        ctx.write_csv(mdf, "metrics_merged.csv")
        st.summary_frame({c: mdf[c] for c in cols}).to_csv(
            ctx.table_path("summary_statistics.csv"))

        pearson, spearman = st.correlation_matrices(mdf[cols])
        pearson.to_csv(ctx.csv_path("correlation_pearson.csv"))
        spearman.to_csv(ctx.csv_path("correlation_spearman.csv"))
        for corr, nm, title in ((pearson, "correlation_pearson", "Pearson"),
                                (spearman, "correlation_spearman", "Spearman")):
            fig, ax = plotting.new_axes(figsize=(6.6, 5.6))
            im = ax.imshow(corr.to_numpy(), cmap=DIV_CMAP, vmin=-1, vmax=1)
            ax.set_xticks(range(len(cols))); ax.set_yticks(range(len(cols)))
            ax.set_xticklabels(cols, rotation=45, ha="right"); ax.set_yticklabels(cols)
            for i in range(len(cols)):
                for j in range(len(cols)):
                    ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center",
                            fontsize=8)
            ax.set_title(f"{title} correlation")
            fig.colorbar(im, ax=ax, shrink=0.85, label="correlation")
            plotting.save_figure(fig, ctx.fig_path(nm), dpi=ctx.config.dpi)

        notable = [(cols[i], cols[j], round(float(pearson.iloc[i, j]), 2))
                   for i in range(len(cols)) for j in range(i+1, len(cols))
                   if abs(pearson.iloc[i, j]) > 0.6]
        return {"metrics": cols, "notable_correlations": notable,
                "figure": "correlation_pearson"}
