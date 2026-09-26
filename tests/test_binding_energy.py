"""MM-GBSA/PBSA: analysis of real gmx_MMPBSA outputs (reference values from a manual analysis)."""

from __future__ import annotations

import pytest

from conftest import FIX
from moldynx.binding import analyse as be
from moldynx.binding import prepare as prep

MM = FIX / "gmx_mmpbsa"


@pytest.fixture(scope="module")
def result():
    return be.analyse(MM / "gb" / "FINAL_RESULTS_MMGBSA.csv", MM / "pb" / "FINAL_RESULTS_MMPBSA.csv",
                      frame_dt_ns=0.1, windows={"0-100": (0, 100), "80-100": (80, 100)},
                      decomp_csv=MM / "gb" / "FINAL_DECOMP_MMGBSA_frames1-3.csv",
                      gb_dat=MM / "gb" / "FINAL_RESULTS_MMGBSA.dat")


def _hl(res, method, window):
    return next(h for h in res["headline"] if h["method"] == method and h["window"] == window)


@pytest.mark.parametrize("method, window, mean, sem", [
    ("GB", "0-100", -33.4, 2.1), ("GB", "80-100", -37.8, 2.4),
    ("PB", "0-100", -40.7, 2.5), ("PB", "80-100", -47.1, 3.3),
])
def test_reproduces_manual_analysis(result, method, window, mean, sem):
    h = _hl(result, method, window)
    assert h["mean"] == pytest.approx(mean, abs=0.05)
    assert h["sem"] == pytest.approx(sem, abs=0.05)          # autocorrelation-corrected


def test_checks_gb_pb_and_entropy_gate(result):
    c = result["checks"]
    assert c["frames"] == 101 and c["single_trajectory_consistent"] and c["same_frames_gb_pb"]
    assert c["max_abs_ggas_gb_minus_pb"] < 0.02
    assert result["gb_vs_pb"]["pearson_r"] == pytest.approx(0.73, abs=0.01)
    assert result["gb_vs_pb"]["mean_offset_pb_minus_gb"] == pytest.approx(-7.4, abs=0.05)
    ie, c2 = result["entropy"]["IE"], result["entropy"]["C2"]
    assert ie["sigma_int_kcal"] == pytest.approx(86.72) and ie["valid"] is False
    assert c2["valid"] is False and c2["limit_kcal"] == 6.0


def test_decomposition_parsing(result):
    assert result["residues_in_decomposition"] > 50
    first = result["decomposition"]["0-100"]
    assert set(first.side) == {"R", "L"} and first.TOTAL.iloc[0] < 0


def test_prepare_helpers():
    ions = prep.derive_ionic_strength({"POT": 1498, "CLA": 1480}, 15605.2)
    assert ions["ionic_strength_M"] == pytest.approx(0.1594, abs=1e-3)
    assert ions["salt_1to1_M"] == pytest.approx(0.1575, abs=1e-3)
    top = ('#include "toppar/forcefield.itp"\n[ system ]\nx\n\n[ molecules ]\n; Compound #mols\n'
           "PROA 1\nPROB 1\nTIP3 518060\nPOT 1498\nCLA 1480\n")
    assert prep.solute_molecules(top) == ["PROA", "PROB"]
    out = prep.complex_topology(top, ["PROA", "PROB"])
    assert "TIP3" not in out and out.rstrip().endswith("PROB\t1") and "#include" in out
    assert prep.ranges([5, 1, 2, 3, 9, 10]) == "1-3,5,9-10"
    gb, pb = prep.input_files(303.15, 0.1594, 1, 1001, 10, "A/1-5 B/10-12")
    assert 'print_res           = "A/1-5 B/10-12"' in gb and "forcefields" not in gb
    assert "PBRadii             = 7" in pb and "istrng              = 0.1594" in pb
