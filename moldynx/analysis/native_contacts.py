"""Fraction of native contacts Q(t) relative to the reference frame."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import PROTEIN_SYSTEMS
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE, SEQ_CMAP


class NativeContacts(BaseAnalysis):
    name = "native_contacts"
    label = "Native-contact fraction Q"
    category = "structure"
    required_files = {"trajectory", "topology"}
    supported_systems = PROTEIN_SYSTEMS
    outputs = ["results/native_contacts.csv", "figures/native_contacts_Q.png"]
    default_params = {"cutoff": 8.0, "seq_sep": 4, "window": 20}

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        ca = u.select_atoms("protein and name CA")
        resids = ca.resids
        R0 = p["cutoff"]
        u.trajectory[ctx.config.start or 0]
        D0 = cdist(ca.positions, ca.positions)
        seq = np.abs(resids[:, None] - resids[None, :]) >= p["seq_sep"]
        native = np.triu((D0 < R0) & seq)
        n_native = int(native.sum()) or 1

        n = len(u.trajectory)
        Q = np.empty(n); times = np.empty(n)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[native_contacts]")):
            D = cdist(ca.positions, ca.positions)
            Q[i] = ((D < R0) & native).sum() / n_native
            times[i] = ts.time / 1000.0

        ctx.write_csv(pd.DataFrame({"time_ns": times, "Q": Q,
                                    "Q_movavg": st.moving_average(Q, p["window"])}),
                      "native_contacts.csv")
        fig, ax = plotting.new_axes()
        ax.plot(times, Q, color=PALETTE["primary"], alpha=0.35, lw=1.0)
        ax.plot(times, st.moving_average(Q, p["window"]), color=PALETTE["primary"], lw=2.2)
        plateau = float(Q[len(Q)//2:].mean())
        ax.axhline(plateau, ls="--", color=PALETTE["accent"], lw=1.1,
                   label=f"plateau = {plateau:.2f}")
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Fraction of native contacts, Q")
        ax.set_ylim(0, 1.02); ax.set_title(self.label); ax.legend(loc="best")
        plotting.save_figure(fig, ctx.fig_path("native_contacts_Q"), dpi=ctx.config.dpi)

        return {"n_native_contacts": n_native, "Q_final_plateau": plateau,
                "Q": st.describe(Q, "Q"),
                "convergence": st.plateau_detection(times, Q),
                "figure": "native_contacts_Q"}
