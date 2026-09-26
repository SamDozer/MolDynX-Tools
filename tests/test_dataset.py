"""Dataset assembly, the shipped verifier, and packaging."""

from __future__ import annotations

import json
import subprocess
import sys

from moldynx.dataset import assemble, package


def _fake_run(root):
    (root / "report").mkdir(parents=True)
    (root / "figures").mkdir()
    (root / "results").mkdir()
    (root / "intake").mkdir()
    (root / "figures" / "pbc_validation.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    (root / "report" / "PBC_VALIDATION.md").write_text(
        "# PBC\n\n## Observations\n\n![evidence](../figures/pbc_validation.png)\n", encoding="utf-8")
    (root / "report" / "report.md").write_text("# Report\n\n[PBC](PBC_VALIDATION.md)\n",
                                               encoding="utf-8")
    (root / "results" / "pbc_summary.json").write_text(json.dumps({"checks": {"a": True}}))
    (root / "results" / "rmsd.csv").write_text("time_ns,rmsd_nm\n0,0.1\n")
    (root / "intake" / "INTAKE_REPORT.md").write_text("# Intake\n")
    (root / "data").mkdir()
    (root / "data" / "core_meta.json").write_text(json.dumps({"chains": []}))
    (root / "manifest.json").write_text(json.dumps({"moldynx_version": "0.3.0", "system": {},
                                                    "analyses": {"pbc_validation": {
                                                        "status": "ok",
                                                        "summary": {"all_checks_pass": True}}}}))
    return root


def test_assemble_verify_package(tmp_path):
    run = _fake_run(tmp_path / "run")
    ds = assemble(run, tmp_path / "ds", title="Demo dataset")
    assert (ds / "1_report" / "PBC_VALIDATION.html").exists()
    assert "](figures/pbc_validation.png)" in (ds / "1_report" / "PBC_VALIDATION.md").read_text()
    assert (ds / "5_validation" / "pbc_summary.json").exists()
    assert (ds / "0_raw_production" / "RAW_DATA_NOT_INCLUDED.md").exists()
    assert "pbc_validation | ok | all proof checks pass" in (ds / "README.md").read_text(encoding="utf-8")
    r = subprocess.run([sys.executable, str(ds / "5_validation" / "verify_dataset.py")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout
    z = package(ds, tmp_path / "ds.zip")
    assert z["crc_ok"] and z["verify_exit"] == 0


def test_verifier_catches_a_broken_link(tmp_path):
    run = _fake_run(tmp_path / "run")
    ds = assemble(run, tmp_path / "ds")
    (ds / "1_report" / "report.md").write_text("# R\n\n![x](figures/missing.png)\n")
    r = subprocess.run([sys.executable, str(ds / "5_validation" / "verify_dataset.py")],
                       capture_output=True, text=True)
    assert r.returncode == 1 and "missing.png" in r.stdout
