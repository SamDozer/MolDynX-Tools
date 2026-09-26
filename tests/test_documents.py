"""Generated documents: content from results, honest wording, self-contained HTML."""

from __future__ import annotations

import json
import re

import numpy as np
import pytest

from conftest import FIX
from moldynx.report import documents as docs
from test_equilibration import ctx as equil_ctx  # noqa: F401
from test_pbc import _chain, _run, _universe


def _html_ok(html: str, md: str):
    text = re.sub(r"<[^>]+>", " ", re.sub(r"<style>.*?</style>", "", html, flags=re.S))
    assert "**" not in text and "`" not in text                      # markdown fully rendered
    assert html.count("<table>") == len(re.findall(r"(?m)^\|[-| ]+\|$", md))
    assert not re.search(r'src="https?://', html)                    # self-contained
    assert html.count("<h2") == len(re.findall(r"(?m)^## ", md))


def test_pbc_document():
    frames = [np.vstack([_chain(13.0), _chain(25.0)]) for _ in range(3)]
    _, s, _ = _run(_universe(frames))
    md = docs.pbc_document(s, fig_rel="missing.png")
    assert "## Proof that processing changed only wholeness" in md
    assert "| B | 10 | 3 | 0 |" in md and "**FAIL**" not in md
    _html_ok(docs.to_html(md, FIX, "PBC"), md)


def test_binding_energy_document():
    from moldynx.binding import CITATIONS, GMX_MMPBSA_URL
    from moldynx.binding.analyse import analyse
    mm = FIX / "gmx_mmpbsa"
    r = analyse(mm / "gb" / "FINAL_RESULTS_MMGBSA.csv", mm / "pb" / "FINAL_RESULTS_MMPBSA.csv",
                frame_dt_ns=0.1, windows={"0-100 ns": (0, 100), "80-100 ns": (80, 100)},
                decomp_csv=mm / "gb" / "FINAL_DECOMP_MMGBSA_frames1-3.csv",
                gb_dat=mm / "gb" / "FINAL_RESULTS_MMGBSA.dat")
    s = {k: r[k] for k in ("checks", "headline", "trend_per_ns", "gb_vs_pb", "entropy", "closure",
                           "residues_in_decomposition")}
    s["hotspots"] = {w: d.head(5).to_dict("records") for w, d in r["decomposition"].items()}
    s["method"] = {"tool": GMX_MMPBSA_URL, "citations": CITATIONS}
    md = docs.binding_energy_document(json.loads(json.dumps(s, default=float)),
                                      {"summary": {"decision_required": True}})
    assert "not experimental affinities" in md and GMX_MMPBSA_URL in md
    assert "| GB | 80-100 ns | -37.8 ± 2.4 |" in md
    assert "invalid — not reported" in md                           # entropy gate in the text
    assert "Please cite" in md and "10.1021/acs.jctc.1c00645" in md
    _html_ok(docs.to_html(md, FIX, "BE"), md)


def test_equilibration_document(equil_ctx):  # noqa: F811
    pytest.importorskip("panedr")
    from moldynx.analysis.equilibration import EquilibrationAudit
    EquilibrationAudit().run(equil_ctx)
    s = json.loads(equil_ctx.csv_path("equilibration_summary.json").read_text(encoding="utf-8"))
    md = docs.equilibration_document(s, fig_rel="../figures/equilibration_overview.png")
    assert "requested tolerance not reached" in md
    assert "## What could not be recovered" in md and "## Interpretation" in md
    out = docs.write_all(equil_ctx)
    assert {p.name for p in out} >= {"EQUILIBRATION.md", "EQUILIBRATION.html"}
    html = (equil_ctx.config.report_dir / "EQUILIBRATION.html").read_text(encoding="utf-8")
    assert "data:image/png;base64," in html                          # figure embedded
    _html_ok(html, md)
