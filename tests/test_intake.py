"""Intake report on a miniature real folder with real (trimmed) stage logs."""

from __future__ import annotations

import json

from moldynx.cli.main import main
from moldynx.io.intake import render_report, run_intake, write_intake


def test_intake_stages_from_real_logs(a8_tree):
    root, _ = a8_tree
    res = run_intake(root)
    em = res.stages["em"]["log"]["minimization"]
    assert em["reached_emtol"] is False and em["steps"] == 3056
    assert res.stages["npt"]["log"]["ref_t_K"] == [303.15, 303.15]
    assert res.stages["production"]["log"]["sessions"] == 11
    pl = res.fileset.evidence["production_log"]
    assert pl["simulated_ps"] == 100000.0
    assert "extended" in pl and "50,000,000 steps were run" in pl["extended"]


def test_intake_report_sections(a8_tree, tmp_path):
    root, _ = a8_tree
    res = run_intake(root)
    text = render_report(res)
    for heading in ("## Canonical production run", "## Stages", "## Capabilities",
                    "## Inputs referenced by job scripts but not in this folder", "## Ignored"):
        assert heading in text
    assert "requested tolerance not reached" in text
    assert "Extended run" in text
    assert "UNVERIFIED" in text                      # the miniature files have no real headers
    js, md = write_intake(res, tmp_path / "out")
    data = json.loads(js.read_text(encoding="utf-8"))
    assert data["validation"]["ok"] is True
    caps = {c["capability"]: c["status"] for c in data["capabilities"]}
    assert caps["MDP-only settings"] == "unavailable"
    assert caps["binding energy (MM-GBSA/PBSA)"] == "available"


def test_intake_never_writes_into_the_simulation_folder(a8_tree, tmp_path, monkeypatch):
    root, _ = a8_tree
    before = sorted(str(p) for p in root.rglob("*"))
    monkeypatch.chdir(tmp_path)
    rc = main(["intake", "--input", str(root)])
    assert rc == 0
    assert sorted(str(p) for p in root.rglob("*")) == before
    assert (tmp_path / "moldynx_results" / root.name / "intake" / "INTAKE_REPORT.md").exists()


def test_temperature_change_between_stages_is_flagged(a8_tree):
    """A later NPT log at 310 K (as in one of the reference datasets) must be surfaced."""
    root, _ = a8_tree
    from conftest import LOGS
    (root / "step4.2_equilibration.log").write_text(
        (LOGS / "Q946V6_step4.2_equilibration_NPT_310K.log").read_text(encoding="utf-8"),
        encoding="utf-8")
    text = render_report(run_intake(root))
    assert "Temperature changes between stages" in text and "310 K" in text
