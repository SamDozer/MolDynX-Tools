"""Shared fixtures: miniature copies of two real CHARMM-GUI/GROMACS simulation folders."""

from __future__ import annotations

import csv
import json
import os
import shutil
import time
from pathlib import Path

import pytest

FIX = Path(__file__).parent / "fixtures"
TREES = FIX / "raw_trees"
LOGS = FIX / "gromacs_logs"

A8_SCRIPT = """#!/bin/bash -l
module load GROMACS/2023.3
IN_DIR="../../input/complex_A8HNE1"
init=step3_input
mini_prefix=step4.0_minimization
gmx grompp -f ${IN_DIR}/${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ${IN_DIR}/${init}.gro -r ${IN_DIR}/${init}.gro -p ${IN_DIR}/topol.top -n ${IN_DIR}/index.ndx -maxwarn 1
gmx mdrun -v -deffnm ${mini_prefix} -ntmpi 1 -ntomp 8
"""
A8_APPEND = """#!/bin/bash -l
if [ ! -f step5_production_extended.tpr ]; then
    gmx convert-tpr -s step5_production.tpr -o step5_production_extended.tpr -extend 99000
fi
gmx mdrun -v -deffnm step5_production -s step5_production_extended.tpr -cpi step5_production.cpt -append
"""
# real (trimmed) logs dropped into the miniature tree at their real locations
A8_LOGS = {
    "step4.0_minimization.log": "A8HNE1_step4.0_minimization_machine_precision.log",
    "step4.1_equilibration.log": "A8HNE1_step4.1_equilibration_NVT.log",
    "step4.2_equilibration.log": "A8HNE1_step4.2_equilibration_NPT_303K.log",
    "step5_production.log": "A8HNE1_step5_production_11_sessions.log",
}


def materialise(root: Path, entries, scripts=None, logs=None):
    """Create a tree; each file gets ceil(bytes / 1e6) bytes so relative sizes survive."""
    for i, (rel, nbytes) in enumerate(entries):
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        if scripts and p.name in scripts:
            p.write_text(scripts[p.name], encoding="utf-8")
        elif logs and rel in logs:
            shutil.copyfile(LOGS / logs[rel], p)
        else:
            p.write_bytes(b"\0" * (nbytes // 1_000_000 + 1))
        stamp = time.time() - 10_000 + i
        os.utime(p, (stamp, stamp))


@pytest.fixture
def a8_tree(tmp_path):
    data = json.loads((TREES / "A8HNE1_raw_tree.json").read_text(encoding="utf-8"))
    root = tmp_path / "complex_A8HNE1"
    materialise(root, [(f["path"], f["bytes"]) for f in data["files"]],
                scripts={"run_A8HNE1_gpu.sh": A8_SCRIPT, "run_append_100ns.sh": A8_APPEND},
                logs=A8_LOGS)
    return root, data["expected"]


@pytest.fixture
def q9_tree(tmp_path):
    rows = list(csv.DictReader((TREES / "Q946V6_raw_archive_inventory.csv").open(encoding="utf-8")))
    root = tmp_path / "Complex_Q94v6"
    materialise(root, [(r["path"].replace("\\", "/"), int(r["bytes"])) for r in rows])
    return root, rows
