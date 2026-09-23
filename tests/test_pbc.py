"""PBC diagnose -> treat -> prove on synthetic two-chain systems with known answers."""

from __future__ import annotations

import numpy as np
import pytest

from moldynx.core.pbc import PBCProcessor

BOX = 30.0
N = 10           # atoms per chain, bonded in a line 1 Å apart


def _universe(frames_true: list[np.ndarray]):
    """Two bonded 10-atom chains (segments A, B); coordinates wrapped into a 30 Å box."""
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader
    u = mda.Universe.empty(2 * N, n_residues=2 * N, n_segments=2,
                           atom_resindex=np.arange(2 * N),
                           residue_segindex=[0] * N + [1] * N, trajectory=True)
    u.add_TopologyAttr("name", ["CA"] * (2 * N))
    u.add_TopologyAttr("resname", ["ALA"] * (2 * N))
    u.add_TopologyAttr("segid", ["A", "B"])
    u.add_TopologyAttr("masses", np.full(2 * N, 12.0))
    bonds = [(i, i + 1) for i in range(N - 1)] + [(N + i, N + i + 1) for i in range(N - 1)]
    u.add_TopologyAttr("bonds", bonds)
    coords = np.array([np.mod(f, BOX) for f in frames_true], dtype=np.float32)
    dims = np.tile([BOX, BOX, BOX, 90, 90, 90], (len(frames_true), 1)).astype(np.float32)
    u.load_new(coords, format=MemoryReader, dimensions=dims)
    return u


def _chain(x0: float) -> np.ndarray:
    return np.column_stack([x0 + np.arange(N), np.full(N, 15.0), np.full(N, 15.0)])


def _run(u, mode="nojump"):
    chains = [u.segments[0].atoms, u.segments[1].atoms]
    proc = PBCProcessor(u.atoms, chains, ["A", "B"], mode=mode, pairs=[(0, 1)])
    out = []
    for ts in u.trajectory:
        out.append(proc.process(ts).copy())
    return proc, proc.summary(), np.array(out)


def test_split_chain_is_made_whole_and_proven():
    # chain B sits across the x boundary (true x 25..34); chain A at 13..22 (3 Å away)
    frames = [np.vstack([_chain(13.0), _chain(25.0)]) for _ in range(5)]
    proc, s, out = _run(_universe(frames))
    b = s["units"][1]
    assert b["frames_split_in_raw"] == 5 and s["units"][0]["frames_split_in_raw"] == 0
    assert b["frames_with_multiple_translations"] == 5          # wrapped atoms moved, rest not
    assert all(s["checks"].values()), s["checks"]
    assert s["max_bond_A"] == pytest.approx(1.0, abs=1e-3)
    pair = s["pairs"][0]
    assert pair["min_dist_raw_pbc_nm"][0] == pytest.approx(0.3, abs=1e-4)
    assert pair["frames_separated_in_processed"] == 0
    np.testing.assert_allclose(np.ptp(out[0, N:, 0]), N - 1, atol=1e-4)   # B whole


def test_drift_through_the_boundary_counts_one_jump_event_and_flags_separation():
    # B moves +1 Å per frame from x 21 to 32 and passes through the boundary;
    # A stays at 10..19, so through the boundary B approaches A's *other* side
    frames = [np.vstack([_chain(10.0), _chain(21.0 + k)]) for k in range(12)]
    proc, s, out = _run(_universe(frames))
    b = s["units"][1]
    assert b["whole_box_translations_undone"] == 1               # one crossing event
    assert b["frames_translated"] >= 1
    # the processed trajectory is continuous: B keeps moving +1 Å per frame
    com_x = out[:, N:, 0].mean(axis=1)
    np.testing.assert_allclose(np.diff(com_x), 1.0, atol=1e-3)
    assert s["checks"]["only_whole_box_translations"] and s["checks"]["molecules_whole"]
    # ...but the nearest image of B is now on A's other side: reported, not hidden
    assert s["checks"]["interchain_distance_preserved"] is False
    assert s["pairs"][0]["frames_separated_in_processed"] > 0


def test_mode_none_records_but_does_not_touch():
    frames = [np.vstack([_chain(13.0), _chain(25.0)]) for _ in range(3)]
    u = _universe(frames)
    raw = u.trajectory.coordinate_array.copy()
    proc, s, out = _run(u, mode="none")
    np.testing.assert_array_equal(out, raw)
    assert s["mode"] == "none" and s["made_whole"] is False
    assert s["units"][1]["frames_split_in_raw"] == 0   # cannot be detected without making whole


def test_chain_identity_survives_pdb_segid_truncation(tmp_path):
    """CHARMM-GUI segids collapse to 'seg_' in PDB files; chain groups must not depend on them."""
    import json

    import MDAnalysis as mda
    from moldynx.core.config import RunConfig
    from moldynx.core.context import AnalysisContext

    frames = [np.vstack([_chain(13.0), _chain(25.0)])]
    u = _universe(frames)
    u.segments.segids = ["seg_0_PROA", "seg_1_PROB"]
    pdb = tmp_path / "core.pdb"
    u.atoms.write(str(pdb))
    back = mda.Universe(str(pdb))
    assert len(set(back.segments.segids)) == 1          # the trap: both chains read back as one

    cfg = RunConfig(output_dir=tmp_path / "out")
    ctx = AnalysisContext.__new__(AnalysisContext)
    ctx.config, ctx._core = cfg, back
    cfg.ensure_dirs()
    (cfg.data_dir / "core_meta.json").write_text(json.dumps({"chains": [
        {"segid": "seg_0_PROA", "core_start": 0, "core_stop": N},
        {"segid": "seg_1_PROB", "core_start": N, "core_stop": 2 * N}]}))
    groups = ctx.chain_groups()
    assert [r["segid"] for r, _ in groups] == ["seg_0_PROA", "seg_1_PROB"]
    assert [g.n_atoms for _, g in groups] == [N, N]


def test_auto_mode_choice():
    frames = [np.vstack([_chain(13.0), _chain(25.0)])]
    u = _universe(frames)
    two = PBCProcessor(u.atoms, [u.segments[0].atoms, u.segments[1].atoms], ["A", "B"])
    one = PBCProcessor(u.atoms, [u.atoms], ["all"])
    assert two.mode == "nojump" and one.mode == "whole"
