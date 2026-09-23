"""
Interface analysis for complexes (protein-protein / protein-nucleic).

Computes buried surface area (BSA), interface residues/contacts over time, and
interface-RMSD. Two partners are auto-selected: protein-vs-nucleic when a nucleic
acid is present, otherwise the two largest protein chains (by segment).

NOTE: implemented but not yet validated end-to-end (the reference dataset is a
single-chain protein). Validate on a real complex trajectory before publication.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import mdtraj as md
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.distances import distance_array

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import COMPLEX_SYSTEMS
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class InterfaceAnalysis(BaseAnalysis):
    name = "interface"
    label = "Interface (BSA, contacts, iRMSD)"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = COMPLEX_SYSTEMS
    default_params = {"cutoff": 5.0}
    outputs = ["results/interface.csv", "figures/interface.png"]

    def _partners(self, ctx, u):
        """
        The two interface partners as (label, AtomGroup) pairs.

        Protein chains come from the chain identity persisted at extraction
        (atom-index ranges): the PDB format truncates CHARMM-GUI segment IDs
        (``seg_0_PROA``/``seg_1_PROB`` -> ``seg_``), so segments read back from
        ``core.pdb`` cannot tell the chains apart.
        """
        if ctx.system.flags.get("has_nucleic"):
            return ("protein", u.select_atoms("protein")), ("nucleic", u.select_atoms("nucleic"))
        chains = sorted(ctx.chain_groups(u), key=lambda rc: rc[1].n_atoms, reverse=True)
        if len(chains) >= 2:
            (ra, a), (rb, b) = chains[0], chains[1]
            return (ra.get("segid") or "chain A", a), (rb.get("segid") or "chain B", b)
        return None, None

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        pa, pb = self._partners(ctx, u)
        if pa is None or pa[1].n_atoms == 0 or pb[1].n_atoms == 0:
            return {"status": "skipped",
                    "reason": "could not resolve two interface partners "
                              f"({len(ctx.chain_groups(u))} protein chain(s) recorded)"}
        (selA, A), (selB, B) = pa, pb

        # BSA via mdtraj on the cached core trajectory (same atom order as the universe)
        traj = md.load(str(ctx.config.data_dir / "core.xtc"),
                       top=str(ctx.config.data_dir / "core.pdb"))
        idxA, idxB = A.indices, B.indices
        sasa_all = md.shrake_rupley(traj, mode="atom")
        bsa = (sasa_all[:, idxA].sum(1) + sasa_all[:, idxB].sum(1)
               - md.shrake_rupley(traj.atom_slice(np.concatenate([idxA, idxB])),
                                  mode="atom").sum(1))
        times = traj.time / 1000.0

        n = len(u.trajectory)
        contacts = np.empty(n)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[interface]")):
            D = distance_array(A.positions, B.positions)
            contacts[i] = int((D < p["cutoff"]).sum())

        ctx.write_csv(pd.DataFrame({"time_ns": times, "bsa_nm2": bsa,
                                    "interface_contacts": contacts}), "interface.csv")
        fig, (a1, a2) = plotting.style.plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True)
        a1.plot(times, bsa, color=PALETTE["primary"], lw=1.6)
        a1.set_ylabel("Buried SASA (nm$^2$)"); a1.set_title("Interface")
        a2.plot(times, contacts, color=PALETTE["secondary"], lw=1.6)
        a2.set_ylabel("Interface contacts"); a2.set_xlabel("Time (ns)")
        for a in (a1, a2):
            a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
        plotting.save_figure(fig, ctx.fig_path("interface"), dpi=ctx.config.dpi)

        return {"partners": [selA, selB], "bsa": st.describe(bsa, "bsa_nm2"),
                "mean_interface_contacts": float(contacts.mean()), "figure": "interface"}
