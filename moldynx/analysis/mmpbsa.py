"""
MM/PBSA workflow generation + structural energy-decomposition proxy.

MM-GBSA/MM-PBSA computes a *binding* free energy between two partners. For a
single ligand-free solute it is undefined, so this module GENERATES a ready-to-run
``gmx_MMPBSA`` workflow (input file + commands) for when two groups are defined
(a bound ligand, or two domains/chains), auto-detecting trajectory/topology/index.
It also provides a force-field-independent per-residue **structural contribution**
proxy from persistent contact occupancy (unitless, not kcal/mol).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from moldynx.core.base import BaseAnalysis
from moldynx import plotting
from moldynx.plotting import PALETTE

MMPBSA_IN = """\
&general
  sys_name  = "moldynx_system",
  startframe = 1, endframe = 999999, interval = 10, verbose = 2,
  forcefields = "oldff/leaprc.ff99SB, leaprc.gaff",
/
&gb
  igb = 5, saltcon = 0.150,
/
&decomp
  idecomp = 1, dec_verbose = 0, print_res = "within 6",
/
"""


class MMPBSA(BaseAnalysis):
    name = "mmpbsa"
    label = "MM/PBSA workflow + structural proxy"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    order = 160  # after contact_map for the structural proxy
    outputs = ["tables/mmpbsa.in", "tables/mmpbsa_run.sh", "report/mmpbsa_workflow.md"]

    def _run_script(self, ctx) -> str:
        traj = (ctx.fileset.trajectory or ctx.config.trajectory)
        top = (ctx.fileset.topology or ctx.config.topology)
        gtop = ctx.fileset.gmx_top
        return f"""#!/usr/bin/env bash
# gmx_MMPBSA workflow (run where GROMACS + AmberTools + gmx_MMPBSA are installed).
# MM/PBSA needs TWO groups -- edit the make_ndx selections to define them.
set -euo pipefail
TPR="{top.name if top else 'topol.tpr'}"
XTC="{traj.name if traj else 'traj.xtc'}"
TOP="{gtop.name if gtop else 'topol.top'}"
gmx make_ndx -f "$TPR" -o index.ndx <<'EOF'
r 1-100
name 20 GroupA
r 101-9999
name 21 GroupB
q
EOF
mpirun -np 4 gmx_MMPBSA -O -i mmpbsa.in -cs "$TPR" -ct "$XTC" -ci index.ndx \\
    -cg 20 21 -cp "$TOP" -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv \\
    -do FINAL_DECOMP_MMPBSA.dat -deo FINAL_DECOMP_MMPBSA.csv
"""

    def run(self, ctx) -> dict:
        plotting.set_style()
        ctx.table_path("mmpbsa.in").write_text(MMPBSA_IN, encoding="utf-8")
        ctx.table_path("mmpbsa_run.sh").write_text(self._run_script(ctx), encoding="utf-8")
        has_ligand = ctx.system.flags.get("has_ligand", False)
        meaningful = has_ligand or ctx.system.n_protein_chains >= 2
        (ctx.config.report_dir / "mmpbsa_workflow.md").write_text(
            f"# MM/PBSA workflow\n\nBinding-ΔG is "
            f"{'applicable (define receptor/ligand groups)' if meaningful else 'undefined for this single ligand-free solute'}. "
            f"Generated `tables/mmpbsa.in` and `tables/mmpbsa_run.sh` "
            f"(edit the group selections). Requires GROMACS + AmberTools + gmx_MMPBSA.\n",
            encoding="utf-8")

        summary = {"workflow_generated": True, "binding_dg_meaningful": bool(meaningful)}
        occ_csv = ctx.csv_path("contact_occupancy.csv")
        if occ_csv.exists():
            occ = pd.read_csv(occ_csv, index_col=0)
            resids = occ.columns.astype(int).to_numpy()
            contribution = occ.to_numpy().sum(axis=1)
            ctx.write_csv(pd.DataFrame({"resid": resids,
                                        "structural_contribution": contribution}),
                          "residue_structural_contribution.csv")
            fig, ax = plotting.new_axes(figsize=(8.2, 4.4))
            ax.bar(resids, contribution, color=PALETTE["accent"], width=1.0)
            ax.set_xlabel("Residue"); ax.set_ylabel("Structural contribution")
            ax.set_title("Per-residue structural contribution (MM/PBSA proxy)")
            plotting.save_figure(fig, ctx.fig_path("residue_structural_contribution"),
                                 dpi=ctx.config.dpi)
            summary["figure"] = "residue_structural_contribution"
            summary["top_residues"] = resids[np.argsort(contribution)[::-1][:5]].tolist()
        return summary
