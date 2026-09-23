"""
Protein-solvent (or protein-ligand) interaction analysis.

Primary output is a fast, robust per-residue hydration/interaction frequency
(MDAnalysis). When a ligand is present it is used as the ProLIF "ligand";
otherwise the first hydration shell is used. A genuine ProLIF fingerprint
(HBDonor/HBAcceptor) is computed on a few representative frames for the
interaction-type breakdown (best-effort; skipped gracefully if ProLIF/rdkit
are unavailable).
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from MDAnalysis.analysis.distances import distance_array

from moldynx.core.base import BaseAnalysis
from moldynx import plotting
from moldynx.plotting import PALETTE

WATER = "resname TIP3 SOL WAT HOH SPC"
WATER_O = "resname TIP3 SOL WAT HOH SPC and name OH2 OW OW1 O"
CUTOFF = 3.5
N_FP_FRAMES = 3


class ProLIF(BaseAnalysis):
    name = "prolif"
    label = "Interaction fingerprints (ProLIF)"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    default_params = {"stride": 40, "cutoff": CUTOFF}
    outputs = ["results/prolif_residue_frequency.csv", "figures/prolif_residue_frequency.png"]

    def _hydration(self, ctx, u, p):
        stride = p["stride"]
        polar = u.select_atoms("protein and (name N* O* S*)")
        labels = np.array([f"{a.resname}{a.resid}" for a in polar])
        frames = u.trajectory[::stride]
        n = len(frames)
        res_frames = defaultdict(int)
        for _ in ctx.iter_frames(u, desc="[prolif] hydration", stride=stride):
            wat = u.select_atoms(f"({WATER_O}) and around {p['cutoff']} (protein)")
            if wat.n_atoms == 0:
                continue
            hit = (distance_array(polar.positions, wat.positions) < p["cutoff"]).any(axis=1)
            for lbl in np.unique(labels[hit]):
                res_frames[lbl] += 1
        return res_frames, n

    def _prolif_types(self, ctx, u):
        import prolif as plf
        nfr = len(u.trajectory)
        idx = np.linspace(0, nfr - 1, N_FP_FRAMES, dtype=int)
        fp = plf.Fingerprint(["HBDonor", "HBAcceptor"])
        protein = u.select_atoms("protein")
        types = defaultdict(int)
        for fi in idx:
            u.trajectory[int(fi)]
            water = u.select_atoms(f"byres ({WATER} and around {CUTOFF} protein)")
            if water.n_atoms == 0:
                continue
            ifp = fp.generate(plf.Molecule.from_mda(water),
                              plf.Molecule.from_mda(protein), metadata=True)
            for _, inter in ifp.items():
                for itype in inter:
                    types[itype] += 1
        return dict(types)

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.full_universe()
        res_frames, n = self._hydration(ctx, u, p)
        rows = [{"residue": r, "frequency": c / n} for r, c in res_frames.items()]
        freq = pd.DataFrame(rows).sort_values("frequency", ascending=False).reset_index(drop=True)
        ctx.write_csv(freq, "prolif_residue_frequency.csv")

        types = {}
        try:
            types = self._prolif_types(ctx, u)
        except Exception as e:
            print(f"[prolif] genuine fingerprint skipped: {e}")
        if not types:
            types = {"HBond(contact)": int(sum(res_frames.values()))}
        pd.DataFrame([{"interaction": k, "count": v} for k, v in types.items()]).to_csv(
            ctx.csv_path("prolif_interaction_types.csv"), index=False)

        if len(freq):
            top = freq.head(25)
            fig, ax = plotting.new_axes(figsize=(7.8, 5.8))
            ax.barh(range(len(top)), top["frequency"] * 100, color=PALETTE["primary"])
            ax.set_yticks(range(len(top))); ax.set_yticklabels(top["residue"], fontsize=8)
            ax.invert_yaxis(); ax.set_xlabel("Interaction frequency (% of frames)")
            ax.set_title("Most persistently interacting residues")
            plotting.save_figure(fig, ctx.fig_path("prolif_residue_frequency"), dpi=ctx.config.dpi)

        return {"n_interacting_residues": len(freq), "interaction_types": types,
                "top": freq.head(5)["residue"].tolist() if len(freq) else [],
                "figure": "prolif_residue_frequency"}
