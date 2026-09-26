"""
Data-driven documents: PBC_VALIDATION, EQUILIBRATION, BINDING_ENERGY (+ HTML).

Every number comes from a results file of the same run; no system names or
conclusions live in code. Each document separates Observations, Interpretation
and Limitations, and the wording rules are fixed here: no "stable" from an RMSD
plateau, binding energies are end-point estimates (not affinities), entropy is
quoted only when valid, unknown metadata is stated as not recorded.
"""

from __future__ import annotations

import base64
import html as _html
import json
import re
from pathlib import Path

CSS = ("body{font-family:Arial,Helvetica,sans-serif;max-width:980px;margin:2rem auto;padding:0 1rem;"
       "color:#1d1d1f;line-height:1.55}h1{border-bottom:3px solid #1f6f8b;padding-bottom:.3rem}"
       "h2{color:#1f6f8b;margin-top:2rem}table{border-collapse:collapse;margin:1rem 0;font-size:.92em}"
       "th,td{border:1px solid #ccc;padding:5px 9px;text-align:left;vertical-align:top}"
       "th{background:#f0f4f6}code{background:#f4f4f4;padding:0 3px;border-radius:3px}"
       "img{max-width:100%;border:1px solid #eee}blockquote{border-left:4px solid #e07a5f;"
       "margin:1rem 0;padding:.4rem 1rem;background:#fdf3ef}nav#toc{background:#f7f9fa;"
       "padding:.6rem 1.2rem;border:1px solid #e3e8eb}nav#toc a{text-decoration:none}")


# --------------------------------------------------------------------------- #
# Markdown -> self-contained HTML
# --------------------------------------------------------------------------- #
def to_html(md_text: str, base_dir: Path, title: str) -> str:
    """Render Markdown (tables, fenced code, sub) with a contents list; embed local images."""
    try:
        import markdown
        body = markdown.markdown(md_text, extensions=["tables", "fenced_code", "attr_list",
                                                      "sane_lists", "md_in_html"])
    except ImportError:   # plain fallback keeps the document readable
        body = "<pre>" + _html.escape(md_text) + "</pre>"
    toc, n = [], 0

    def anchor(m):
        nonlocal n
        n += 1
        level, text = m.group(1), m.group(2)
        slug = f"s{n}"
        if level == "2":
            toc.append(f'<li><a href="#{slug}">{re.sub("<[^>]+>", "", text)}</a></li>')
        return f'<h{level} id="{slug}">{text}</h{level}>'
    body = re.sub(r"<h([23])>(.*?)</h\1>", anchor, body)

    def embed(m):
        src = m.group(1)
        p = (base_dir / src).resolve()
        if p.exists() and p.suffix.lower() in (".png", ".jpg", ".jpeg", ".svg"):
            mime = "image/svg+xml" if p.suffix.lower() == ".svg" else f"image/{p.suffix[1:].lower()}"
            return f'src="data:{mime};base64,{base64.b64encode(p.read_bytes()).decode()}"'
        return m.group(0)
    body = re.sub(r'src="([^"]+)"', embed, body)
    nav = f'<nav id="toc"><strong>Contents</strong><ul>{"".join(toc)}</ul></nav>' if toc else ""
    return (f"<!doctype html><html><head><meta charset='utf-8'><title>{_html.escape(title)}"
            f"</title><style>{CSS}</style></head><body>{nav}{body}</body></html>")


def write_doc(path: Path, md_text: str, formats=("md", "html")) -> list[Path]:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = []
    if "md" in formats:
        path.write_text(md_text, encoding="utf-8")
        out.append(path)
    if "html" in formats:
        title = (re.search(r"^# (.+)$", md_text, re.M) or [None, path.stem])[1]
        h = path.with_suffix(".html")
        h.write_text(to_html(md_text, path.parent, title), encoding="utf-8")
        out.append(h)
    return out


def _load(p: Path) -> dict | None:
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _f(x, nd=2, unit=""):
    try:
        return f"{float(x):.{nd}f}{unit}"
    except (TypeError, ValueError):
        return "—"


def _pm(d: dict, key_mean="mean", key_err="sem", nd=1):
    return f"{_f(d.get(key_mean), nd)} ± {_f(d.get(key_err), nd)}"


# --------------------------------------------------------------------------- #
# PBC
# --------------------------------------------------------------------------- #
def pbc_document(s: dict, fig_rel: str = "../figures/pbc_validation.png") -> str:
    c = s["checks"]
    ok = all(c.values())
    L = ["# Periodic-boundary validation", "",
         f"Treatment: **{s['mode']}** · molecules made whole: {'yes' if s['made_whole'] else '**no**'}"
         f" · frames: {s['n_frames']}", "",
         "## Observations", "",
         "| Molecule | Atoms | Frames split in the raw data | Whole-box translations undone | "
         "Extent mean / max (nm) |", "|---|---|---|---|---|"]
    for u in s["units"]:
        L.append(f"| {u['label']} | {u['n_atoms']:,} | {u['frames_split_in_raw']} | "
                 f"{u['whole_box_translations_undone']} | {_f(u['extent_whole_nm_mean'])} / "
                 f"{_f(u['extent_whole_nm_max'])} |")
    L.append("")
    for p in s.get("pairs", []):
        r = p["min_dist_raw_pbc_nm"]
        L.append(f"- Minimum heavy-atom distance {p['pair'][0]}–{p['pair'][1]} (PBC-aware, raw): "
                 f"{_f(r[0], 3)}–{_f(r[2], 3)} nm (mean {_f(r[1], 3)}); frames without a "
                 f"heavy-atom contact: {p['frames_without_heavy_atom_contact']}.")
    if s.get("half_box_min_nm"):
        L.append(f"- Half of the smallest box edge (minimum-image limit): {_f(s['half_box_min_nm'])} nm.")
    L += ["", "## Proof that processing changed only wholeness / continuity", "",
          "| Check | Result |", "|---|---|"]
    names = {"only_whole_box_translations": "every atom moved only by whole box vectors",
             "molecules_whole": "no bond longer than 2.5 Å after processing",
             "multiple_translations_only_in_split_frames":
                 "two translation vectors only in frames where the molecule was split",
             "interchain_distance_preserved": "partner distances equal their PBC-aware raw values"}
    L += [f"| {names.get(k, k)} | {'pass' if v else '**FAIL**'} |" for k, v in c.items()]
    L += ["", f"Largest deviation from a pure box translation: "
              f"{_f(s.get('max_dev_from_box_translation_A'), 5)} Å.", "",
          f"![PBC evidence]({fig_rel})", "", "## Interpretation", "",
          ("All checks pass: the processed solute trajectory differs from the raw one only by "
           "whole-box translations that make molecules whole and continuous; distance-based "
           "analyses are unaffected by periodic images." if ok else
           "At least one check failed. Analyses that depend on distances between molecules "
           "must be read with the failed check in mind (see the table above)."), "",
          "## Limitations", "",
          "- The diagnosis and proof cover the analysed frames only (frame slice of this run).",
          "- A centre-of-mass distance is not a proximity measure for extended molecules; "
          "contact is judged from minimum heavy-atom distances."]
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------- #
# Equilibration
# --------------------------------------------------------------------------- #
LABEL = {"em": "Energy minimisation", "nvt": "NVT", "npt": "NPT", "production": "Production"}


def equilibration_document(s: dict, fig_rel: str = "../figures/equilibration_overview.png") -> str:
    st = s["stages"]
    L = ["# Preparation and equilibration", "",
         "Reconstructed from the stage logs, energy files, run inputs and job scripts found at "
         "intake. Missing inputs are listed, never invented.", "",
         "## Stages", "", "| Stage | Integrator / dt | Length | Thermostat | Barostat | Outcome |",
         "|---|---|---|---|---|---|"]
    for k in ("em", "nvt", "npt", "production"):
        r = st.get(k)
        if not r:
            continue
        em = r.get("minimization")
        if em:
            outcome = (f"{em.get('outcome')} in {em.get('steps')} steps; F<sub>max</sub> "
                       f"{_f(em.get('fmax'), 1)} kJ mol⁻¹ nm⁻¹"
                       + (" — **requested tolerance not reached**"
                          if em.get("reached_emtol") is False else ""))
            length, thermo, baro = f"{r.get('nsteps')} steps max", "—", "—"
        else:
            sim = r.get("simulated_ps")
            length = (f"{sim / 1000:g} ns" if sim and sim >= 1000 else f"{sim:g} ps") if sim else "—"
            ref = r.get("ref_t_K") or []
            thermo = f"{r.get('tcoupl')} {ref[0]:g} K" if ref else str(r.get("tcoupl"))
            baro = str(r.get("pcoupl")) + (f", τ<sub>p</sub> {r.get('tau_p_ps')} ps"
                                           if r.get("pcoupl") not in (None, "No") else "")
            c = r.get("counts", {})
            outcome = (f"{r.get('sessions')} session(s); LINCS warnings "
                       f"{c.get('lincs_warnings', 0)}, fatal errors {c.get('fatal_errors', 0)}")
        L.append(f"| {LABEL[k]} | {r.get('integrator')} / {r.get('dt_ps') or '—'} ps | {length} | "
                 f"{thermo} | {baro} | {outcome} |")
    L += ["", "## Observations", ""]
    loc = st.get("em", {}).get("fmax_location")
    if loc:
        L.append(f"- Largest residual force after minimisation: atom {loc['atom_number']} "
                 f"({loc['atom']} of {loc['resname']}, chain {loc.get('chain', loc.get('segid'))}, "
                 f"chain residue {loc.get('chain_residue', loc.get('resid'))}).")
    nvt = st.get("nvt", {})
    if nvt.get("temperature_settled_ps") is not None:
        a = nvt.get("temperature_after_settling", {})
        L.append(f"- NVT temperature within ±2 K of target from {_f(nvt['temperature_settled_ps'], 0)} ps; "
                 f"{_f(a.get('mean'), 2)} ± {_f(a.get('sd'), 2)} K afterwards.")
    npt = st.get("npt", {})
    e = npt.get("energy", {}).get("terms", {})
    if e:
        parts = []
        for term, unit, nd in (("Temperature", "K", 3), ("Pressure", "bar", 2),
                               ("Density", "kg m⁻³", 2)):
            if term in e:
                parts.append(f"{term.lower()} {_f(e[term]['tail_mean'], nd)} ± "
                             f"{_f(e[term]['tail_sd'], nd)} {unit}")
        tail = npt["energy"].get("tail_ps")
        L.append(f"- NPT, last {_f(tail, 0)} ps: " + "; ".join(parts) + ".")
    if npt.get("box_x_nm"):
        b = npt["box_x_nm"]
        v = npt.get("volume_nm3") or [None, None, None]
        L.append(f"- Box {_f(b[0], 3)} → {_f(b[1], 3)} nm (volume change {_f(v[2], 1)} %); density "
                 f"within 0.1 % of its final value after {_f(npt.get('density_reached_99_9pct_ps'), 0)} ps.")
    pr = s.get("position_restraints", {})
    if pr:
        L += ["", "| Run input | Restrained atoms | By force constant (kJ mol⁻¹ nm⁻²) | Per molecule type |",
              "|---|---|---|---|"]
        for k, r in pr.items():
            if r.get("available"):
                fc = ", ".join(f"{n} @ {f}" for f, n in r["by_force_constant"].items()) or "none"
                L.append(f"| {LABEL.get(k, k)} | {r['n_restrained']:,} | {fc} | "
                         f"{', '.join(str(x) for x in r['by_molecule_block']) or '—'} |")
            else:
                L.append(f"| {LABEL.get(k, k)} | not audited | {r.get('error', '')} | — |")
    prot = s.get("protonation", {}).get("non_default", [])
    if prot:
        L += ["", "Non-default protonation states: " + "; ".join(
            f"{p['resname']} {p.get('chain_residue', p['resid'])} ({p.get('chain', '')}, {p['state']})"
            for p in prot) + "."]
    tl = [r for r in s.get("timeline", []) if r.get("gap_after_previous_s") is not None]
    if tl:
        L.append("Gaps between consecutive stages: " + ", ".join(
            f"{LABEL[r['stage']]} started {_f(r['gap_after_previous_s'], 0)} s after the previous "
            "stage ended" for r in tl) + ".")
    L += ["", f"![Equilibration overview]({fig_rel})", "", "## Interpretation", ""]
    L += [f"- {v[0].upper() + v[1:]}." for v in s.get("verdicts", [])]
    L += ["", "## What could not be recovered", ""]
    L += [f"- {m}" for m in s.get("missing", [])] or ["- Nothing: every expected stage file was found."]
    L += ["", "## Limitations", "",
          "- A Berendsen barostat relaxes the box but does not sample the NPT ensemble; "
          "only production statistics should be interpreted thermodynamically.",
          "- A residual drift at the end of equilibration means production started from a system "
          "that was still relaxing; production convergence is assessed separately."]
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------- #
# Binding energy
# --------------------------------------------------------------------------- #
def binding_energy_document(s: dict, window_sel: dict | None = None) -> str:
    L = ["# Binding free energy (MM-GBSA / MM-PBSA)", "",
         "> **End-point estimates, not experimental affinities.** Values are mean ± "
         "autocorrelation-corrected standard error (kcal/mol) with the number of effectively "
         "independent frames.", "",
         f"Computed with [gmx_MMPBSA]({s['method']['tool']}) (AmberTools MMPBSA.py), prepared "
         "and analysed by MolDynX Tools.", "", "## Results", "",
         "| Method | Window | ΔG<sub>bind</sub> | N (effective) |", "|---|---|---|---|"]
    for h in s["headline"]:
        L.append(f"| {h['method']} | {h['window']} | {_pm(h)} | {h['n']} ({_f(h['n_eff'], 0)}) |")
    tr = s.get("trend_per_ns", {})
    drifting = [m for m, d in tr.items() if d.get("p") is not None and d["p"] < 0.01]
    if drifting:
        L += ["", "> **Not stationary:** ΔG<sub>bind</sub> drifts over the run ("
              + "; ".join(f"{m} {_f(tr[m]['slope'] * 10, 2)} kcal/mol per 10 ns, p = "
                          f"{tr[m]['p']:.1g}" for m in drifting)
              + "). Each window describes a different state; no single value is an equilibrium "
                "average."]
    if window_sel and window_sel.get("summary", {}).get("decision_required"):
        L += ["", "> The averaging window was **not** chosen by the user "
                  "(`binding_energy.primary_window`); the full run and its final 20 % are shown."]
    g = s.get("gb_vs_pb")
    if g:
        L += ["", f"GB and PB agree frame by frame with Pearson r = {_f(g['pearson_r'], 2)}; "
                  f"PB − GB = {_f(g['mean_offset_pb_minus_gb'], 1)} kcal/mol on average."]
    ent = s.get("entropy") or {}
    if ent:
        L += ["", "## Entropy", "", "| Method | σ(interaction energy) | Validity limit | Status |",
              "|---|---|---|---|"]
        for k, e in ent.items():
            L.append(f"| {k} | {_f(e['sigma_int_kcal'], 1)} | < {e['limit_kcal']} | "
                     f"{'valid' if e['valid'] else '**invalid — not reported**'} |")
    clo = s.get("closure") or {}
    if clo:
        L += ["", "## Per-residue decomposition", "",
              f"{s.get('residues_in_decomposition')} residues decomposed. Closure (Σ residues vs "
              "ΔG<sub>bind</sub>, GB):", "", "| Window | Σ residues | Total | Unattributed |",
              "|---|---|---|---|"]
        L += [f"| {w} | {_f(c['sum_residues'], 2)} | {_f(c['total'], 2)} | "
              f"{_f(c['unattributed'], 2)}{'' if c['closes'] else ' **(> 1 kcal/mol)**'} |"
              for w, c in clo.items()]
        for w, rows in (s.get("hotspots") or {}).items():
            L += ["", f"Top residues, {w}:", "", "| Partner | Residue | ΔG (kcal/mol) |",
                  "|---|---|---|"]
            L += [f"| {r.get('partner') or '—'} | {r['resname']} {r['residue']} | {_f(r['TOTAL'], 2)} |"
                  for r in rows[:15]]
    c = s["checks"]
    L += ["", "## Validation", "",
          f"- Frames: {c['frames']}; GB and PB on identical frames: {c['same_frames_gb_pb']}.",
          f"- Largest bonded Δ term: {_f(c['max_abs_bonded_delta'], 6)} kcal/mol "
          f"({'single-trajectory consistent' if c['single_trajectory_consistent'] else '**inconsistent**'}).",
          "", "## Limitations", "",
          "- Implicit solvent, single-trajectory protocol, no conformational entropy unless its "
          "validity limit is met.",
          "- Ionic strength and settings are recorded in `prepare_meta.json`.",
          "", "## Please cite", ""] + [f"- {c}" for c in s["method"]["citations"]]
    return "\n".join(L) + "\n"


def write_all(ctx, formats=("md", "html")) -> list[Path]:
    """Write every document whose results exist; returns written paths."""
    rd, out = ctx.config.report_dir, []
    s = _load(ctx.csv_path("pbc_summary.json"))
    if s:
        out += write_doc(rd / "PBC_VALIDATION.md", pbc_document(s), formats)
    s = _load(ctx.csv_path("equilibration_summary.json"))
    if s:
        out += write_doc(rd / "EQUILIBRATION.md", equilibration_document(s), formats)
    s = _load(ctx.csv_path("binding_energy_summary.json"))
    if s:
        out += write_doc(rd / "BINDING_ENERGY.md",
                         binding_energy_document(s, _load(ctx.csv_path("window_selection.json"))),
                         formats)
    return out
