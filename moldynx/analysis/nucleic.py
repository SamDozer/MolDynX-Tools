"""
Protein-nucleic analyses: protein-DNA/RNA contacts and nucleic RMSD.

Runs for protein-DNA and protein-RNA systems. Reports the protein-nucleic contact
count over time and the nucleic-acid backbone RMSD. Detailed groove-width analysis
requires specialised tools (Curves+, do_x3dna) and is intentionally left as an
external step, noted in the report.

NOTE: implemented but not yet validated end-to-end (the reference dataset has no
nucleic acid). Validate on a real protein-DNA/RNA trajectory before publication.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.distances import distance_array

from moldynx.core.base import BaseAnalysis
from moldynx.core.system import SystemType
from moldynx import plotting
from moldynx import statistics as st
from moldynx.plotting import PALETTE


class NucleicAnalysis(BaseAnalysis):
    name = "nucleic"
    label = "Protein-nucleic contacts & nucleic RMSD"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = {SystemType.PROTEIN_DNA, SystemType.PROTEIN_RNA}
    default_params = {"cutoff": 4.0}
    outputs = ["results/nucleic.csv", "figures/nucleic.png"]

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        protein = u.select_atoms("protein")
        nucleic = u.select_atoms("nucleic")
        if protein.n_atoms == 0 or nucleic.n_atoms == 0:
            return {"note": "protein or nucleic selection empty"}

        # nucleic backbone RMSD to reference
        na_bb = "nucleic and name P O3' O5' C3' C4' C5'"
        R = rms.RMSD(u, select=na_bb, ref_frame=ctx.config.start or 0).run()
        na_rmsd = R.results.rmsd[:, 2] / 10.0

        n = len(u.trajectory)
        contacts = np.empty(n); times = np.empty(n)
        for i, ts in enumerate(ctx.iter_frames(u, desc="[nucleic]")):
            D = distance_array(protein.positions, nucleic.positions)
            contacts[i] = int((D < p["cutoff"]).sum())
            times[i] = ts.time / 1000.0

        ctx.write_csv(pd.DataFrame({"time_ns": times, "nucleic_rmsd_nm": na_rmsd,
                                    "protein_nucleic_contacts": contacts}), "nucleic.csv")
        fig, (a1, a2) = plotting.style.plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True)
        a1.plot(times, na_rmsd, color=PALETTE["primary"], lw=1.6)
        a1.set_ylabel("Nucleic RMSD (nm)"); a1.set_title("Protein-nucleic")
        a2.plot(times, contacts, color=PALETTE["secondary"], lw=1.6)
        a2.set_ylabel("Protein-nucleic contacts"); a2.set_xlabel("Time (ns)")
        for a in (a1, a2):
            a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
        plotting.save_figure(fig, ctx.fig_path("nucleic"), dpi=ctx.config.dpi)

        return {"nucleic_rmsd": st.describe(na_rmsd, "nucleic_rmsd_nm"),
                "mean_contacts": float(contacts.mean()),
                "note": "groove-width analysis requires Curves+/do_x3dna (external)",
                "figure": "nucleic"}
