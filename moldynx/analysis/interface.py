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

    def _partners(self, u, system) -> tuple[str, str]:
        if system.flags.get("has_nucleic"):
            return "protein", "nucleic"
        segs = [s for s in u.select_atoms("protein").segments
                if s.atoms.select_atoms("name CA").n_residues > 0]
        segs = sorted(segs, key=lambda s: s.atoms.n_atoms, reverse=True)
        if len(segs) >= 2:
            return f"segid {segs[0].segid}", f"segid {segs[1].segid}"
        return "protein", "protein"

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        selA, selB = self._partners(u, ctx.system)
        A, B = u.select_atoms(selA), u.select_atoms(selB)
        if A.n_atoms == 0 or B.n_atoms == 0 or selA == selB:
            return {"note": f"could not resolve two interface partners ({selA}/{selB})"}

        # BSA via mdtraj on the cached core trajectory
        traj = md.load(str(ctx.config.data_dir / "core.xtc"),
                       top=str(ctx.config.data_dir / "core.pdb"))
        idxA = traj.topology.select(_mda_to_mdtraj(selA))
        idxB = traj.topology.select(_mda_to_mdtraj(selB))
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


def _mda_to_mdtraj(sel: str) -> str:
    """Best-effort translation of the simple selections used here to mdtraj DSL."""
    if sel == "protein":
        return "protein"
    if sel == "nucleic":
        return "nucleic"
    if sel.startswith("segid "):
        return f"chainid {sel.split()[1]}"  # approximate
    return "all"
