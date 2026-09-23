"""GROMACS evidence readers, checked against real (trimmed) logs from two production datasets."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from moldynx.io.gromacs import parse_log, parse_mdp, scan_xtc, tpr_header, windows_to_wsl

FIX = Path(__file__).parent / "fixtures"
LOGS = FIX / "gromacs_logs"


# --------------------------------------------------------------------------- #
# logs
# --------------------------------------------------------------------------- #
def test_em_stopped_on_machine_precision():
    log = parse_log(LOGS / "A8HNE1_step4.0_minimization_machine_precision.log")
    assert log.trimmed
    assert log.gromacs_version == "2023.3"
    assert log.is_minimization
    assert log.mdp["integrator"] == "steep"
    assert log.value("emtol") == 100.0 and log.value("nsteps") == 50000.0
    em = log.minimization
    assert em.algorithm == "Steepest Descents"
    assert em.outcome == "converged to machine precision"
    assert em.reached_emtol is False                      # tolerance NOT reached
    assert em.steps == 3056
    assert em.fmax == pytest.approx(7232.5469, rel=1e-6)
    assert em.fmax_atom == 6116
    assert em.potential == pytest.approx(-2.6774416e7)
    assert em.fnorm == pytest.approx(11.151116, rel=1e-6)
    assert log.counts["wrote_pdb"] == 4
    assert "-ntmpi 1 -ntomp 8" in log.commands[0]


def test_em_converged():
    log = parse_log(LOGS / "Q946V6_step4.0_minimization_converged.log")
    em = log.minimization
    assert em.outcome == "converged to Fmax < 1000" and em.reached_emtol is True
    assert em.steps == 2027
    assert em.fmax == pytest.approx(714.97339, rel=1e-6) and em.fmax_atom == 6883
    assert log.value("emtol") == 1000.0 and log.value("nsteps") == 5000.0
    assert log.counts["wrote_pdb"] == 17


def test_nvt_stage():
    log = parse_log(LOGS / "A8HNE1_step4.1_equilibration_NVT.log")
    assert not log.is_minimization and log.minimization is None
    assert log.mdp["tcoupl"] == "V-rescale" and log.mdp["pcoupl"] == "No"
    assert log.value("dt") == 0.001 and log.value("nsteps") == 125000.0
    assert log.mdp["continuation"] == "false"
    assert log.ref_t == [303.15, 303.15]
    assert log.grpopts["tau-t"] == ["1", "1"]
    assert log.grpopts["nrdf"] == ["33024", "3.11729e+06"]
    assert log.n_sessions == 1
    assert log.started == ["2026-06-04T04:11:36"] and log.finished == ["2026-06-04T04:44:50"]
    assert log.performance_ns_day == [5.418]
    assert log.last_statistics_steps == 125001
    assert log.simulated_ps == pytest.approx(125.0)
    assert log.counts["lincs_warnings"] == 0 and log.counts["fatal_errors"] == 0


@pytest.mark.parametrize("name, ref_t", [
    ("A8HNE1_step4.2_equilibration_NPT_303K.log", 303.15),
    ("Q946V6_step4.2_equilibration_NPT_310K.log", 310.0),
])
def test_npt_stage_temperatures_differ(name, ref_t):
    log = parse_log(LOGS / name)
    assert log.mdp["pcoupl"] == "Berendsen" and log.mdp["pcoupltype"] == "Isotropic"
    assert log.value("tau-p") == 1.0
    assert float(log.mdp["ref-p_xx"]) == 1.0
    assert float(log.mdp["compressibility_xx"]) == pytest.approx(4.5e-5)
    assert log.mdp["refcoord-scaling"] == "COM" and log.mdp["continuation"] == "true"
    assert log.ref_t == [ref_t, ref_t]
    assert log.simulated_ps == pytest.approx(2000.0)


@pytest.mark.parametrize("name, sessions", [
    ("A8HNE1_step5_production_11_sessions.log", 11),
    ("Q946V6_step5_production_12_sessions.log", 12),
])
def test_production_sessions_and_extension(name, sessions):
    log = parse_log(LOGS / name)
    assert log.n_sessions == sessions
    assert log.value("dt") == 0.002
    assert log.value("nsteps") == 500000.0             # the 1-ns TPR that was later extended
    assert log.mdp["pcoupl"] == "C-rescale" and log.value("tau-p") == 5.0
    assert log.last_statistics_steps == 50000001
    assert log.simulated_ps == pytest.approx(100000.0)  # 100 ns proven from the log
    assert log.mdp_last.get("nsteps") == "500000"       # appended sessions do not reprint params


# --------------------------------------------------------------------------- #
# mdp
# --------------------------------------------------------------------------- #
def test_parse_mdp(tmp_path):
    p = tmp_path / "x.mdp"
    p.write_text("; comment\ndefine = -DPOSRES  ; restraints\nnsteps=125000\ngen_vel = yes\n")
    mdp = parse_mdp(p)
    assert mdp == {"define": "-DPOSRES", "nsteps": "125000", "gen-vel": "yes"}


# --------------------------------------------------------------------------- #
# xtc
# --------------------------------------------------------------------------- #
def _write_xtc(path: Path, n_atoms: int, n_frames: int, dt: float = 5.0):
    import MDAnalysis as mda
    u = mda.Universe.empty(n_atoms, trajectory=True)
    u.dimensions = [30, 30, 30, 90, 90, 90]
    rng = np.random.default_rng(0)
    with mda.Writer(str(path), n_atoms) as w:
        for i in range(n_frames):
            u.atoms.positions = rng.uniform(0, 30, (n_atoms, 3))
            u.trajectory.ts.time = i * dt
            u.trajectory.ts.data["step"] = i * 5000
            w.write(u.atoms)


@pytest.mark.parametrize("n_atoms", [5, 250])          # uncompressed (<= 9) and compressed
def test_scan_xtc_matches_mdanalysis(tmp_path, n_atoms):
    import MDAnalysis as mda
    p = tmp_path / "t.xtc"
    _write_xtc(p, n_atoms, 26)
    info = scan_xtc(p)
    assert info.readable and not info.truncated and info.error is None
    assert info.natoms == n_atoms and info.n_frames == 26
    assert info.t_first_ps == 0.0 and info.t_last_ps == pytest.approx(125.0)
    assert info.dt_ps == pytest.approx(5.0) and info.uniform_spacing
    u = mda.Universe.empty(n_atoms, trajectory=False)
    u.load_new(str(p))
    assert len(u.trajectory) == info.n_frames


def test_scan_xtc_does_not_write_cache_files(tmp_path):
    p = tmp_path / "t.xtc"
    _write_xtc(p, 100, 5)
    before = sorted(x.name for x in tmp_path.iterdir())
    scan_xtc(p)
    assert sorted(x.name for x in tmp_path.iterdir()) == before


def test_scan_xtc_truncated_and_fake(tmp_path):
    p = tmp_path / "t.xtc"
    _write_xtc(p, 100, 5)
    data = p.read_bytes()
    (tmp_path / "cut.xtc").write_bytes(data[: len(data) - 50])
    cut = scan_xtc(tmp_path / "cut.xtc")
    assert cut.n_frames == 4 and cut.truncated
    (tmp_path / "fake.xtc").write_bytes(b"")
    fake = scan_xtc(tmp_path / "fake.xtc")
    assert not fake.readable and fake.n_frames == 0


# --------------------------------------------------------------------------- #
# tpr / paths
# --------------------------------------------------------------------------- #
def test_tpr_header_unreadable_is_reported_not_raised(tmp_path):
    p = tmp_path / "fake.tpr"
    p.write_bytes(b"\x00" * 16)
    h = tpr_header(p)
    assert not h.readable and h.error


def test_position_restraint_parser(monkeypatch):
    import subprocess

    from moldynx.io import gromacs as g
    # gmx dump lines in the real format (filtered by grep at the source)
    dump = "\n".join([
        "           functype[1240]=POSRES, pos0A=( 0.0e+00, 0.0e+00, 0.0e+00), "
        "fcA=( 4.00000000e+02, 4.00000000e+02, 4.00000000e+02), pos0B=( 0, 0, 0), fcB=( 400, 400, 400)",
        "           functype[1242]=POSRES, pos0A=( 0.0e+00, 0.0e+00, 0.0e+00), "
        "fcA=( 4.00000000e+01, 4.00000000e+01, 4.00000000e+01), pos0B=( 0, 0, 0), fcB=( 40, 40, 40)",
        "      Position Rest.:",
        "            0 type=1240 (POSRES)   0",
        "            1 type=1240 (POSRES)   4",
        "            2 type=1242 (POSRES)   6",
        "      Position Rest.:",
        "            0 type=1242 (POSRES)   0",
        "      Position Rest.:",
    ])
    monkeypatch.setattr(g, "run_gmx_pipeline",
                        lambda *a, **k: subprocess.CompletedProcess([], 0, dump, ""))
    r = g.tpr_position_restraints(g.Gmx("native", ["gmx"]), "x.tpr")
    assert r.available and r.n_restrained == 4
    assert r.by_force_constant == {"400": 2, "40": 2}
    assert r.by_molecule_block == [3, 1]
    none = g.tpr_position_restraints(None, "x.tpr")
    assert not none.available and "not available" in none.error


def test_windows_to_wsl():
    assert windows_to_wsl(r"E:\a b\c.tpr") == "/mnt/e/a b/c.tpr"
    assert windows_to_wsl("/already/posix") == "/already/posix"
