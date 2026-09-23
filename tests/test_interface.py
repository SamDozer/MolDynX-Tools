"""Interface analysis on the real two-chain geometry of a reference complex (Cα-only fixture)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from conftest import FIX
from moldynx.analysis.interface import InterfaceAnalysis
from moldynx.core.config import RunConfig
from moldynx.core.context import AnalysisContext


class _System:
    flags = {"has_nucleic": False}
    chains = []


@pytest.fixture
def ctx(tmp_path):
    import MDAnalysis as mda
    src = mda.Universe(str(FIX / "structures" / "A8HNE1_ZmBiP2_CA_only.pdb"))
    assert len(set(src.segments.segids)) == 1                 # truncated 'seg_' in the fixture
    cfg = RunConfig(output_dir=tmp_path / "out",
                    params={"interface": {"contact_cutoff": 8.0, "interface_cutoff": 8.0,
                                          "ever_cutoff": 10.0, "window": 2, "bsa_stride": 1}})
    cfg.ensure_dirs()
    # a 5-frame trajectory: the real complex, with the second chain drifting away by 0.5 Å/frame
    core_pdb, core_xtc = cfg.data_dir / "core.pdb", cfg.data_dir / "core.xtc"
    src.atoms.write(str(core_pdb))
    with mda.Writer(str(core_xtc), src.atoms.n_atoms) as w:
        for k in range(5):
            src.atoms[187:].translate([0.5 * k, 0, 0]) if k else None
            src.trajectory.ts.time = k * 100.0
            w.write(src.atoms)
    (cfg.data_dir / "core_meta.json").write_text(json.dumps({"chains": [
        {"segid": "seg_0_PROA", "index": 0, "core_start": 0, "core_stop": 187},
        {"segid": "seg_1_PROB", "index": 1, "core_start": 187, "core_stop": 850}]}))
    c = AnalysisContext.__new__(AnalysisContext)
    c.config, c.system, c._core = cfg, _System(), mda.Universe(str(core_pdb), str(core_xtc))
    return c


def test_interface_on_real_complex_geometry(ctx):
    out = InterfaceAnalysis().run(ctx)
    assert out["partners"] == ["seg_0_PROA", "seg_1_PROB"]      # resolved despite 'seg_' PDB
    assert out["partner_residue_ranges"] == [[1, 187], [188, 850]]
    ts = pd.read_csv(ctx.csv_path("interface_timeseries.csv"))
    assert len(ts) == 5
    # the partner drifts away -> the minimum distance grows, contacts do not increase
    assert np.all(np.diff(ts.min_interface_dist_nm) > 0)
    assert ts.n_inter_contacts.iloc[-1] <= ts.n_inter_contacts.iloc[0]
    core_a, core_b = out["core_interface_residues"]
    res = pd.read_csv(ctx.csv_path("interface_residues.csv"))
    assert len(core_a) == int(((res.partner == "seg_0_PROA") &
                               (res.interface_occupancy >= 0.5)).sum())   # full list, no [:10]
    ever_a, ever_b = out["residues_within_10A_ever"]
    assert set(core_a) <= set(ever_a) and set(core_b) <= set(ever_b)
    bsa = ts.buried_area_nm2.to_numpy()
    assert np.isfinite(bsa).all() and (bsa > 0).all()
    assert np.all(np.diff(bsa) <= 1e-3)                       # partner leaves -> less buried
    assert out["trend"]["contacts"]["slope_per_10ns"] is not None
    for f in ("interface", "interface_residues", "interface_contact_map"):
        assert (ctx.config.figures_dir / f"{f}.png").exists()


def test_interface_skips_with_reason_for_single_chain(ctx):
    (ctx.config.data_dir / "core_meta.json").write_text(json.dumps({"chains": [
        {"segid": "only", "index": 0, "core_start": 0, "core_stop": 850}]}))
    out = InterfaceAnalysis().run(ctx)
    assert out["status"] == "skipped" and "1 protein chain" in out["reason"]
