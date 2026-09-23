"""Discovery on miniature copies of two real CHARMM-GUI/GROMACS folders (sizes scaled 1e-6)."""

from __future__ import annotations

from pathlib import Path

import pytest

from moldynx.io.discovery import classify_stage, discover_files


def test_stage_classification():
    assert classify_stage("step4.0_minimization.tpr") == "em"
    assert classify_stage("step4.1_equilibration.log") == "nvt"
    assert classify_stage("step4.2_equilibration-003.xtc") == "npt"
    assert classify_stage("step5_production_ext.tpr") == "production"
    assert classify_stage("md.tpr") == "production"
    assert classify_stage("topol.top") == "setup" and classify_stage("step3_input.gro") == "setup"
    assert classify_stage("mdout.mdp") is None           # grompp record, not a stage input
    assert classify_stage("filenames.txt") is None


def test_a8_production_chosen_by_stage_not_size(a8_tree):
    root, expected = a8_tree
    fs = discover_files(root, deep=True)                 # fake bytes -> falls back, labelled
    rel = lambda p: p.relative_to(root).as_posix() if p else None  # noqa: E731
    assert rel(fs.topology) == expected["canonical_topology"]      # NOT step4.2_equilibration.tpr
    assert rel(fs.trajectory) == expected["canonical_trajectory"]
    assert rel(fs.energy) == expected["canonical_energy"]
    assert rel(fs.log) == expected["canonical_log"]
    assert rel(fs.gmx_top) == expected["gmx_top"]
    assert rel(fs.toppar) == expected["toppar_dir"]
    assert rel(fs.structure) == "Production/step5_production.gro"  # not the step33c.pdb dump
    assert len(fs.crash_dumps) == 8
    assert "UNVERIFIED" in fs.evidence["trajectory_choice"]         # honest about fake files
    assert fs.evidence["verified"] is False


def test_a8_ignores_earlier_analysis_outputs(a8_tree):
    root, _ = a8_tree
    fs = discover_files(root, deep=False)
    ignored = {Path(p).name for p, _ in fs.ignored}
    assert {"analysis", "analysis_nopbc", "mmpbsa_BiP2"} <= ignored
    assert all("analysis" not in str(p) for p in fs.all_trajectories)
    assert any("backup" in why for _, why in fs.ignored)            # #step5_production.gro.1#


def test_a8_stages_and_missing_inputs(a8_tree):
    root, _ = a8_tree
    fs = discover_files(root, deep=False)
    for stage in ("em", "nvt", "npt"):
        assert fs.stage_file(stage, "log") is not None
        assert fs.stage_file(stage, "edr") is not None
        assert fs.stage_file(stage, "tpr") is not None
    assert fs.mdp == {}                                    # the stage MDPs are not in the archive
    missing = fs.evidence["job_script_inputs_not_found"]
    assert "step4.0_minimization.mdp" in missing and "step3_input.gro" in missing
    assert "step5_production_extended.tpr" in missing     # referenced, never copied
    keys = fs.available_keys()
    assert {"em", "nvt", "npt", "stage_tprs", "gmx_top", "toppar", "production_log"} <= keys
    assert "mdp" not in keys


def test_q9_archive(q9_tree):
    root, rows = q9_tree
    fs = discover_files(root, deep=True)
    rel = lambda p: p.relative_to(root).as_posix() if p else None  # noqa: E731
    assert rel(fs.trajectory) == "complex_Q946v6 output/step5_production-001.xtc"
    assert rel(fs.topology) == "complex_Q946v6 output/step5_production_ext.tpr"
    assert any("production run inputs" in w for w in fs.warnings)   # two TPRs, choice explained
    assert rel(fs.gmx_top) == "complex_Q946v6 input/topol.top"
    assert set(fs.mdp) == {"em", "nvt", "npt", "production"}
    ignored = {Path(p).name for p, _ in fs.ignored}
    assert {"complex_Q946v6_analysis", "diagnostics_1ns"} <= ignored   # nojump.xtc, step5_fixed.xtc
    assert len(fs.crash_dumps) == 48
    assert any(p.name == "gromacs_complex_Q946v6.zip" for p in fs.setup_package)
    assert {p.name for p in fs.unclassified} >= {"test_init.tpr", "visual_center_test.gro",
                                                 "index_split.ndx"}


def test_ambiguous_equal_trajectories_are_flagged(tmp_path):
    """Two complete production trajectories covering the same span must not be chosen silently."""
    import MDAnalysis as mda
    import numpy as np
    for name in ("md_run1.xtc", "md_run2.xtc"):
        u = mda.Universe.empty(20, trajectory=True)
        with mda.Writer(str(tmp_path / name), 20) as w:
            for i in range(3):
                u.atoms.positions = np.random.default_rng(i).uniform(0, 9, (20, 3))
                u.trajectory.ts.time = i * 10.0
                w.write(u.atoms)
    fs = discover_files(tmp_path, deep=True)
    assert fs.ambiguities and "same span" in fs.ambiguities[0]
