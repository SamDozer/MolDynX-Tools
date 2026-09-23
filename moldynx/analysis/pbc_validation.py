"""
PBC validation: the diagnosis of the raw trajectory and the proof that processing
changed nothing but molecular wholeness / continuity.

The measurements are taken during solute extraction (one pass over the raw data,
see :mod:`moldynx.core.pbc`); this analysis reports them and draws the evidence.
"""

from __future__ import annotations

import json

import pandas as pd

from moldynx import plotting
from moldynx.core.base import BaseAnalysis
from moldynx.plotting import PALETTE


class PBCValidation(BaseAnalysis):
    name = "pbc_validation"
    label = "Periodic-boundary diagnosis and proof"
    category = "quality"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    order = 1   # first: every other analysis relies on the processed trajectory
    outputs = ["results/pbc_per_frame.csv", "results/pbc_summary.json",
               "figures/pbc_validation.png"]

    def run(self, ctx) -> dict:
        plotting.set_style()
        ctx.core_universe()                       # extraction records diagnosis + proof
        summary = json.loads(ctx.csv_path("pbc_summary.json").read_text(encoding="utf-8"))
        df = pd.read_csv(ctx.csv_path("pbc_per_frame.csv"))
        t = df.time_ps.to_numpy() / 1000.0

        n_panels = 2 if summary.get("pairs") else 1
        fig, axes = plotting.style.plt.subplots(n_panels, 1, figsize=(7.6, 2.9 * n_panels + 0.6),
                                                sharex=True, squeeze=False)
        ax = axes[0, 0]
        colours = [PALETTE[c] for c in ("primary", "secondary", "green", "purple", "accent")]
        for k, unit in enumerate(summary["units"][:5]):
            col = f"extent_whole_nm_{k}"
            if col in df:
                ax.plot(t, df[col], color=colours[k % 5], lw=1.2,
                        label=f"{unit['label']} ({unit['frames_split_in_raw']} split frames)")
                split = df.get(f"split_raw_{k}", pd.Series(False, index=df.index)).astype(bool)
                if split.any():
                    ax.plot(t[split], df[col][split], "o", ms=3, color=colours[k % 5])
        if summary.get("half_box_min_nm"):
            ax.axhline(summary["half_box_min_nm"], ls="--", lw=1, color=PALETTE["muted"],
                       label=f"half box ({summary['half_box_min_nm']:.2f} nm)")
        ax.set_ylabel("Molecule extent (nm)")
        ax.set_title(f"PBC treatment: {summary['mode']}")
        ax.legend(fontsize=8, loc="best")
        if n_panels == 2:
            ax2 = axes[1, 0]
            for p in summary["pairs"][:3]:
                i, j = summary_labels_index(summary, p["pair"])
                raw = df.get(f"min_dist_raw_pbc_nm_{i}_{j}")
                proc = df.get(f"min_dist_processed_nm_{i}_{j}")
                if raw is None:
                    continue
                ax2.plot(t, raw, color=PALETTE["primary"], lw=1.4,
                         label=f"{p['pair'][0]}–{p['pair'][1]} raw (PBC-aware)")
                ax2.plot(t, proc, color=PALETTE["secondary"], lw=0.8, ls="--",
                         label="processed")
            ax2.axhline(0.45, ls=":", lw=1, color=PALETTE["muted"], label="0.45 nm contact")
            ax2.set_ylabel("Min. heavy-atom\ndistance (nm)")
            ax2.legend(fontsize=8, loc="best")
        axes[-1, 0].set_xlabel("Time (ns)")
        plotting.save_figure(fig, ctx.fig_path("pbc_validation"), dpi=ctx.config.dpi)

        return {"mode": summary["mode"], "made_whole": summary["made_whole"],
                "checks": summary["checks"],
                "all_checks_pass": bool(all(summary["checks"].values())),
                "units": [{k: u[k] for k in ("label", "frames_split_in_raw",
                                             "whole_box_translations_undone")}
                          for u in summary["units"]],
                "pairs": summary.get("pairs", []),
                "max_dev_from_box_translation_A": summary.get("max_dev_from_box_translation_A"),
                "figure": "pbc_validation"}


def summary_labels_index(summary: dict, pair: list[str]) -> tuple[int, int]:
    labels = [u["label"] for u in summary["units"]]
    return labels.index(pair[0]), labels.index(pair[1])

