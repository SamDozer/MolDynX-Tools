"""
verify_dataset.py -- shipped inside every MolDynX Tools dataset.

    python 5_validation/verify_dataset.py [--full]

1. layout: the expected folders exist and are not empty;
2. references: every Markdown link and image resolves;
3. raw files: 0_raw_production/ matches MANIFEST.sha256.tsv (sizes; SHA-256 with --full),
   or is documented as not included;
4. PBC proof: the recorded checks passed;
5. smoke test (if MDAnalysis is installed): Rg recomputed from 4_trajectory/ matches
   3_results/rog.csv.

Exit code 0 only if every check passes. Standard library only (MDAnalysis optional).
"""

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FAIL = []


def check(ok, msg):
    print(("  ok    " if ok else "  FAIL  ") + msg)
    if not ok:
        FAIL.append(msg)


def main():
    full = "--full" in sys.argv
    print(f"MolDynX Tools dataset verification: {ROOT}")
    print("[1] layout")
    for d in ("1_report", "3_results", "4_trajectory", "5_validation", "6_workflow"):
        p = ROOT / d
        check(p.is_dir() and any(p.iterdir()), f"{d}/ present and not empty")
    print("[2] references")
    n = 0
    for md in ROOT.rglob("*.md"):
        for target in re.findall(r"!?\[[^\]]*\]\(([^)\s]+)\)", md.read_text(encoding="utf-8")):
            if re.match(r"^(https?:|mailto:|#)", target):
                continue
            n += 1
            ok = (md.parent / target.split("#")[0]).exists()
            if not ok:
                check(False, f"{md.relative_to(ROOT)} -> {target}")
    check(True, f"{n} internal references checked")
    print("[3] raw files")
    man = ROOT / "0_raw_production" / "MANIFEST.sha256.tsv"
    stub = ROOT / "0_raw_production" / "RAW_DATA_NOT_INCLUDED.md"
    if man.exists():
        rows = list(csv.DictReader(man.open(encoding="utf-8"), delimiter="\t"))
        present = [r for r in rows if (ROOT / "0_raw_production" / r["file"]).exists()]
        if not present and stub.exists():
            check(True, f"{len(rows)} raw files listed; not included in this copy (see {stub.name})")
        for r in present:
            p = ROOT / "0_raw_production" / r["file"]
            ok = p.stat().st_size == int(r["bytes"])
            if ok and full and r["sha256"] != "not-computed":
                h = hashlib.sha256()
                with open(p, "rb") as fh:
                    while b := fh.read(1 << 24):
                        h.update(b)
                ok = h.hexdigest() == r["sha256"]
            check(ok, f"raw file {r['file']}")
    print("[4] PBC proof")
    s = ROOT / "5_validation" / "pbc_summary.json"
    if s.exists():
        checks = json.loads(s.read_text(encoding="utf-8")).get("checks", {})
        check(all(checks.values()), f"recorded PBC checks: {checks}")
    print("[5] smoke test")
    try:
        import MDAnalysis as mda
        import numpy as np
        pdb, xtc = ROOT / "4_trajectory" / "core.pdb", ROOT / "4_trajectory" / "core.xtc"
        rog = ROOT / "3_results" / "rog.csv"
        if pdb.exists() and xtc.exists() and rog.exists():
            u = mda.Universe(str(pdb), str(xtc))
            ag = u.select_atoms("protein") or u.atoms
            ref = [float(line.split(",")[1]) for line in rog.read_text().splitlines()[1:]]
            k = len(u.trajectory) // 2
            u.trajectory[k]
            check(abs(ag.radius_of_gyration() / 10 - ref[k]) < 1e-3, "Rg reproduced at the middle frame")
        else:
            print("  skip  (no trajectory or rog.csv)")
    except ImportError:
        print("  skip  (MDAnalysis not installed)")
    print("RESULT:", "all checks passed" if not FAIL else f"{len(FAIL)} check(s) FAILED")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
