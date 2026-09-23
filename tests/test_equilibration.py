"""Equilibration audit on a miniature real folder (real trimmed logs, one real .edr)."""

from __future__ import annotations

import json
import shutil

import numpy as np
import pandas as pd
import pytest

from conftest import FIX
from moldynx.analysis.equilibration import (EquilibrationAudit, energy_statistics,
                                            fraction_time, settle_time)
from moldynx.core.config import RunConfig
from moldynx.core.context import AnalysisContext
from moldynx.io import gromacs
from moldynx.io.discovery import discover_files


@pytest.fixture
def ctx(a8_tree, tmp_path, monkeypatch):
    root, _ = a8_tree
    shutil.copyfile(FIX / "edr" / "A8HNE1_step4.1_equilibration.edr",
                    root / "step4.1_equilibration.edr")
    monkeypatch.setattr(gromacs, "find_gmx", lambda *a, **k: None)
    fs = discover_files(root)
    cfg = RunConfig(input_dir=root, output_dir=tmp_path / "out")
    return AnalysisContext(cfg, None, fs)


def test_audit_on_real_logs(ctx):
    pytest.importorskip("panedr")
    out = EquilibrationAudit().run(ctx)
    summary = json.loads(ctx.csv_path("equilibration_summary.json").read_text(encoding="utf-8"))
    em = summary["stages"]["em"]["minimization"]
    assert em["reached_emtol"] is False and em["steps"] == 3056
    nvt = summary["stages"]["nvt"]
    assert nvt["ref_t_K"] == [303.15, 303.15] and nvt["tcoupl"] == "V-rescale"
    assert nvt["temperature_settled_ps"] == pytest.approx(13.0)          # real NVT energy file
    assert nvt["temperature_after_settling"]["mean"] == pytest.approx(303.2, abs=0.1)
    assert summary["stages"]["npt"]["pcoupl"] == "Berendsen"
    assert summary["crash_dumps"]["pairs"] == 4 and summary["crash_dumps"]["all_accounted_for"]
    assert summary["position_restraints"]["em"]["available"] is False    # no GROMACS: reported
    custody = summary["chain_of_custody"]
    grompp = [c for c in custody if c["tool"] == "grompp"][0]
    assert grompp["args"]["restraint_reference"].endswith("step3_input.gro")
    assert any(c["tool"] == "convert-tpr" and c["args"]["extend_ps"] == "99000" for c in custody)
    gaps = [r.get("gap_after_previous_s") for r in summary["timeline"]]
    assert gaps[1:3] == [13.0, 13.0]                                     # EM->NVT->NPT back to back
    v = " | ".join(out["verdicts"])
    assert "without reaching the requested force tolerance" in v
    assert "all matched to 'Wrote pdb' events" in v
    assert "stage .mdp files are not present" in v
    assert (ctx.config.figures_dir / "equilibration_overview.png").exists()


def test_skipped_without_stage_files(tmp_path):
    (tmp_path / "md.xtc").write_bytes(b"")
    (tmp_path / "md.tpr").write_bytes(b"")
    fs = discover_files(tmp_path, deep=False)
    ctx = AnalysisContext(RunConfig(input_dir=tmp_path, output_dir=tmp_path / "o"), None, fs)
    out = EquilibrationAudit().run(ctx)
    assert out["status"] == "skipped"


def test_statistics_helpers():
    t = np.arange(0, 2001, 1.0)
    x = 1000 + 0.1 * t / 1000 + np.random.default_rng(0).normal(0, 0.02, len(t))  # SE ≈ 0.002/ns
    s = energy_statistics(pd.DataFrame({"Time": t, "Density": x}))
    d = s["terms"]["Density"]
    assert s["tail_ps"] == 500.0
    assert d["drift_per_ns"] == pytest.approx(0.1, abs=0.03) and d["drift_p"] < 1e-3
    temp = np.r_[np.linspace(200, 303, 20), np.full(100, 303.0)]
    assert settle_time(np.arange(len(temp), dtype=float), temp, 303.15, 2.0) == pytest.approx(19.0)
    dens = np.r_[np.linspace(950, 1014, 10), np.full(90, 1014.5)]
    assert fraction_time(np.arange(len(dens), dtype=float), dens) == pytest.approx(9.0)


def test_job_script_parsing(tmp_path):
    from moldynx.io.jobscripts import parse_job_script
    s = tmp_path / "run.sh"
    s.write_text('IN="../in"\nequi=step4.1_equilibration\n'
                 'gmx grompp -f ${IN}/${equi}.mdp -o ${equi}.tpr -c step4.0_minimization.gro '
                 '-r ${IN}/step3_input.gro -p ${IN}/topol.top -n ${IN}/index.ndx -maxwarn 1\n'
                 'gmx mdrun -v -deffnm ${equi} -nb gpu -ntmpi 3 -ntomp 8 -pin on\n'
                 '# gmx mdrun -deffnm commented_out\n', encoding="utf-8")
    calls = parse_job_script(s)
    assert [c.tool for c in calls] == ["grompp", "mdrun"]
    g = calls[0].args
    assert g["mdp"] == "../in/step4.1_equilibration.mdp"
    assert g["coordinates"] == "step4.0_minimization.gro"
    assert g["restraint_reference"] == "../in/step3_input.gro" and g["maxwarn"] == "1"
    assert calls[1].args == {"deffnm": "step4.1_equilibration", "nb": "gpu", "ntmpi": "3",
                             "ntomp": "8"}
