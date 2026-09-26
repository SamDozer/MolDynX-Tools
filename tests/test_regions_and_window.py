"""Phase 5b/6: contact lifetimes, porcupine, annotations, analysis window, autocorrelation stats."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from test_interface import ctx  # noqa: F401  (shared real-geometry complex fixture)
from moldynx.analysis.analysis_window import AnalysisWindow
from moldynx.analysis.annotation import Annotation
from moldynx.analysis.contact_lifetime import ContactLifetime
from moldynx.analysis.interface import InterfaceAnalysis
from moldynx.analysis.porcupine import Porcupine
from moldynx.core import annotations as ann
from moldynx.statistics import describe_correlated, detect_equilibration, statistical_inefficiency


# --------------------------------------------------------------------------- #
def test_statistical_inefficiency():
    rng = np.random.default_rng(1)
    assert statistical_inefficiency(rng.normal(size=5000)) < 1.2          # white noise
    x = np.zeros(20000)
    for i in range(1, len(x)):                                              # AR(1), phi = 0.9
        x[i] = 0.9 * x[i - 1] + rng.normal()
    g = statistical_inefficiency(x)
    assert 14 < g < 25                                                     # theory: 19
    d = describe_correlated(x)
    assert d["sem"] == pytest.approx(d["sem_naive"] * np.sqrt(g), rel=1e-6)


def test_detect_equilibration_skips_transient():
    rng = np.random.default_rng(2)
    x = np.r_[np.linspace(10, 0, 200), rng.normal(0, 0.3, 800)]
    assert 150 <= detect_equilibration(x, step=10)["t0"] <= 260


# --------------------------------------------------------------------------- #
def test_region_transfer_and_motifs():
    pytest.importorskip("Bio")
    ref = "MKVLAAGIVGLLLAQPAVSAQEKEVGTVIGIDLGTTYSCVGVFKNGRVEIIANDQGNRITPSYVAFTDGERLIGD"
    query = ref[4:40] + "GG" + ref[40:]            # N-terminal truncation + an insertion
    mapped = ann.transfer_regions(query, ref, {"core": [26, 60]})[0]
    assert mapped["start"] == 26 - 4 and not mapped["uncertain"]
    assert mapped["end"] == 60 - 4 + 2                     # shifted by the insertion
    assert ann.locate_motif(query, ref[30:40])["match"] == "exact"
    # a motif hanging off the N-terminus (only its last part is in the construct)
    part = ann.locate_motif(query, "XXXXXX" + query[:6])
    assert part["match"] in ("partial", "similar") and part["start"] == 1


def test_numbering_offset():
    c = ann.chain_annotations([{"segid": "B", "resid_first": 188, "resid_last": 850,
                                "sequence": "A"}],
                              {"chains": [{"segid": "B", "display": "Partner", "role": "receptor"}]})
    assert c["B"].bio(188) == 1 and c["B"].bio_range == (1, 663) and c["B"].md(205) == 392
    c = ann.chain_annotations([{"segid": "B", "resid_first": 214, "resid_last": 876,
                                "sequence": "A"}], {"chains": [{"segid": "B",
                                                                "numbering_offset": 213}]})
    assert c["B"].bio(418) == 205


# --------------------------------------------------------------------------- #
def test_lifetime_porcupine_annotation_window(ctx):  # noqa: F811
    InterfaceAnalysis().run(ctx)
    life = ContactLifetime().run(ctx)
    assert life["n_native_contacts"] > 0 and 0 <= life["final_survival_fraction"] <= 1
    por = Porcupine().run(ctx)
    assert 0 < por["pc1_variance_fraction"] <= 1
    assert set(por["mean_displacement_by_chain_A"]) == {"seg_0_PROA", "seg_1_PROB"}

    seq_b = json.loads((ctx.config.data_dir / "core_meta.json").read_text())["chains"][1]["sequence"]
    ctx.config.annotations = {
        "chains": [{"segid": "seg_1_PROB", "display": "Partner B", "role": "receptor"}],
        "domains": {"seg_1_PROB": {"reference_name": "self", "reference_sequence": seq_b,
                                   "regions": {"first half": [1, 300]}}},
        "motifs": [{"name": "m1", "sequence": seq_b[99:111], "chain": "seg_1_PROB"}],
        "docking_site": {"seg_1_PROB": [100, 101, 102]}}
    out = Annotation().run(ctx)
    assert out["chains"]["seg_1_PROB"]["biological_range"] == [1, 663]
    assert out["motifs"][0]["start"] == 100 and out["motifs"][0]["match"] == "exact"
    reg = pd.read_csv(ctx.csv_path("region_summary.csv"))
    assert set(reg.kind) == {"domain", "motif", "docking site"}

    assert AnalysisWindow().run(ctx)["status"] == "skipped"     # 5 frames: too short to judge
    t = np.linspace(0, 100, 101)
    rng = np.random.default_rng(3)
    pd.DataFrame({"time_ns": t,
                  "n_inter_contacts": 60 + 0.4 * t + rng.normal(0, 2, 101),        # still growing
                  "min_interface_dist_nm": 0.27 + rng.normal(0, 0.01, 101)}).to_csv(
        ctx.csv_path("interface_timeseries.csv"), index=False)
    win = AnalysisWindow().run(ctx)
    assert win["binding_observables_stationary"] is False
    assert "Residue–residue contacts" in win["non_stationary"]
    assert win["decision_required"] is True and win["proposal"]["final_20pct"] == [80.0, 100.0]


def test_annotation_without_config_says_so(ctx):  # noqa: F811
    out = Annotation().run(ctx)
    assert out["annotated"] is False and "no annotations" in out["note"]
