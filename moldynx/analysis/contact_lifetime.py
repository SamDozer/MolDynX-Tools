"""
Inter-partner contact lifetimes and native-contact survival (complexes).

One streaming pass tracks every inter-partner residue–residue contact (heavy atoms
< ``cutoff``): how long each contact survives once formed, how often it re-forms,
and what fraction of the frame-0 contacts are still present over time.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from moldynx import plotting
from moldynx.core.base import BaseAnalysis
from moldynx.core.system import COMPLEX_SYSTEMS
from moldynx.plotting import PALETTE


class ContactLifetime(BaseAnalysis):
    name = "contact_lifetime"
    label = "Contact lifetimes and native-contact survival"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = COMPLEX_SYSTEMS
    order = 45
    default_params = {"cutoff": 4.5, "long_lived_ns": 1.0}
    outputs = ["results/contact_lifetime.csv", "results/native_contact_survival.csv",
               "figures/contact_lifetime_distribution.png", "figures/native_contact_survival.png"]

    def run(self, ctx) -> dict:
        from moldynx.analysis.interface import InterfaceAnalysis
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        pa, pb = InterfaceAnalysis()._partners(ctx, u)
        if pa is None:
            return {"status": "skipped", "reason": "fewer than two partners"}
        (nameA, A), (nameB, B) = pa, pb
        hA, hB = A.select_atoms("not name H*"), B.select_atoms("not name H*")
        locA = np.searchsorted(A.residues.resindices, hA.resindices)
        locB = np.searchsorted(B.residues.resindices, hB.resindices)
        rA, rB = A.residues.resids, B.residues.resids

        n = len(u.trajectory)
        times, survival = np.empty(n), np.zeros(n)
        active, lifetimes = {}, defaultdict(list)
        formed, occ = defaultdict(int), defaultdict(int)
        native = None
        for i, ts in enumerate(ctx.iter_frames(u, desc="[contact_lifetime]")):
            times[i] = ts.time / 1000.0
            pairs = capped_distance(hA.positions, hB.positions, max_cutoff=p["cutoff"],
                                    box=ts.dimensions, return_distances=False)
            cur = set(zip(locA[pairs[:, 0]].tolist(), locB[pairs[:, 1]].tolist())) \
                if len(pairs) else set()
            if native is None:
                native = set(cur)
            survival[i] = len(native & cur) / len(native) if native else 0.0
            for q in cur:
                if q in active:
                    active[q] += 1
                else:
                    active[q] = 1
                    formed[q] += 1
                occ[q] += 1
            for q in [q for q in active if q not in cur]:
                lifetimes[q].append(active.pop(q))
        for q, s in active.items():
            lifetimes[q].append(s)

        dt = float(times[1] - times[0]) if n > 1 else 0.0
        rows, all_lt = [], []
        for (a, b), lts in lifetimes.items():
            ns = np.asarray(lts) * dt
            all_lt.extend(ns.tolist())
            rows.append({"a_resid": int(rA[a]), "b_resid": int(rB[b]), "occupancy": occ[(a, b)] / n,
                         "n_events": formed[(a, b)], "mean_lifetime_ns": float(ns.mean()),
                         "max_lifetime_ns": float(ns.max())})
        df = pd.DataFrame(rows, columns=["a_resid", "b_resid", "occupancy", "n_events",
                                         "mean_lifetime_ns", "max_lifetime_ns"])
        df = df.sort_values("max_lifetime_ns", ascending=False)
        ctx.write_csv(df, "contact_lifetime.csv")
        ctx.write_csv(pd.DataFrame({"time_ns": times, "survival_fraction": survival}),
                      "native_contact_survival.csv")
        all_lt = np.asarray(all_lt)

        fig, ax = plotting.new_axes()
        if all_lt.size:
            ax.hist(all_lt, bins=40, color=PALETTE["primary"], alpha=0.85)
            ax.axvline(np.median(all_lt), ls="--", color=PALETTE["accent"],
                       label=f"median = {np.median(all_lt):.2f} ns")
            ax.legend()
            ax.set_yscale("log")
        ax.set_xlabel("Contact lifetime (ns)")
        ax.set_ylabel("Formation events")
        ax.set_title("Inter-partner contact lifetimes")
        plotting.save_figure(fig, ctx.fig_path("contact_lifetime_distribution"), dpi=ctx.config.dpi)
        fig, ax = plotting.new_axes()
        ax.plot(times, survival * 100, color=PALETTE["secondary"], lw=2)
        ax.set_ylim(0, 105)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Surviving initial contacts (%)")
        ax.set_title(f"Native contact survival ({nameA}–{nameB})")
        plotting.save_figure(fig, ctx.fig_path("native_contact_survival"), dpi=ctx.config.dpi)

        return {"n_native_contacts": len(native or ()),
                "final_survival_fraction": float(survival[-1]) if n else None,
                "n_distinct_contacts": int(len(df)),
                "median_lifetime_ns": float(np.median(all_lt)) if all_lt.size else None,
                "max_lifetime_ns": float(all_lt.max()) if all_lt.size else None,
                "n_long_lived_contacts": int((df.max_lifetime_ns >= p["long_lived_ns"]).sum()),
                "figure": "native_contact_survival"}
