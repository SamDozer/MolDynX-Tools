"""
Protein-ligand analyses: ligand RMSD, ligand-protein contacts, binding pocket.

Runs for protein-ligand systems. The ligand is taken from the detected ligand
component (or overridden with ``--ligand RESNAME``).

NOTE: implemented but not yet validated end-to-end (the reference dataset has no
ligand). Validate on a real protein-ligand trajectory before publication.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from MDAnalysis.analysis import align
from MDAnalysis.analysis.distances import distance_array

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import SystemType
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class LigandAnalysis(BaseAnalysis):
    name = "ligand"
    label = "Ligand RMSD / contacts / pocket"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = {SystemType.PROTEIN_LIGAND}
    default_params = {"cutoff": 4.0}
    outputs = ["results/ligand_rmsd.csv", "figures/ligand.png"]

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        lig_sel = ctx.config.ligand and f"resname {ctx.config.ligand}" \
            or ctx.selection("ligand", "not protein")
        lig = u.select_atoms(lig_sel)
        if lig.n_atoms == 0:
            return {"note": f"no ligand atoms for selection {lig_sel!r}"}

        # align on protein, then ligand RMSD to reference frame
        align.AlignTraj(u, u, select="protein and name CA",
                        ref_frame=ctx.config.start or 0, in_memory=True).run()
        u.trajectory[ctx.config.start or 0]
        ref = lig.positions.copy()
        protein = u.select_atoms("protein")
        n = len(u.trajectory)
        lrmsd = np.empty(n); ncontacts = np.empty(n); times = np.empty(n)
        pocket = defaultdict(int)
        prot_res_labels = None
        for i, ts in enumerate(ctx.iter_frames(u, desc="[ligand]")):
            d = lig.positions - ref
            lrmsd[i] = np.sqrt((d * d).sum(1).mean()) / 10.0
            D = distance_array(protein.positions, lig.positions)
            close_atoms = (D < p["cutoff"]).any(axis=1)
            ncontacts[i] = int(close_atoms.sum())
            near_res = set(protein[close_atoms].resids.tolist())
            for r in near_res:
                pocket[r] += 1
            times[i] = ts.time / 1000.0

        ctx.write_csv(pd.DataFrame({"time_ns": times, "ligand_rmsd_nm": lrmsd,
                                    "ligand_contacts": ncontacts}), "ligand_rmsd.csv")
        pocket_df = pd.DataFrame([{"resid": r, "occupancy": c / n}
                                  for r, c in pocket.items()]).sort_values(
            "occupancy", ascending=False)
        ctx.write_csv(pocket_df, "ligand_pocket.csv")

        fig, (a1, a2) = plotting.style.plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True)
        a1.plot(times, lrmsd, color=PALETTE["primary"], lw=1.6)
        a1.set_ylabel("Ligand RMSD (nm)"); a1.set_title("Protein-ligand")
        a2.plot(times, ncontacts, color=PALETTE["secondary"], lw=1.6)
        a2.set_ylabel("Ligand contacts"); a2.set_xlabel("Time (ns)")
        for a in (a1, a2):
            a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
        plotting.save_figure(fig, ctx.fig_path("ligand"), dpi=ctx.config.dpi)

        return {"ligand_selection": lig_sel, "ligand_rmsd": st.describe(lrmsd, "ligand_rmsd_nm"),
                "pocket_residues": pocket_df.head(10)["resid"].tolist(), "figure": "ligand"}
