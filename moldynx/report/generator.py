"""
Structured report builder.

Content is assembled as a list of typed blocks (heading/paragraph/table/image),
then rendered to Markdown and self-contained HTML from the same source, so the
two formats never drift.  PDF is attempted via weasyprint if installed.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path


# --------------------------------------------------------------------------- #
# Document model + renderers
# --------------------------------------------------------------------------- #
def _md(blocks) -> str:
    out = []
    for kind, payload in blocks:
        if kind == "h1":
            out.append(f"# {payload}\n")
        elif kind == "h2":
            out.append(f"## {payload}\n")
        elif kind == "p":
            out.append(f"{payload}\n")
        elif kind == "table":
            headers, rows = payload
            out.append("| " + " | ".join(headers) + " |")
            out.append("|" + "|".join(["---"] * len(headers)) + "|")
            for r in rows:
                out.append("| " + " | ".join(str(c) for c in r) + " |")
            out.append("")
        elif kind == "img":
            path, caption = payload
            rel = f"../figures/{Path(path).name}"
            out.append(f"\n![{caption}]({rel})\n\n*Figure. {caption}.*\n")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
def _fmt(x, nd=3):
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return "n/a"


def generate_report(ctx, results: dict, manifest, formats=("md", "html")) -> list[Path]:
    system = ctx.system
    blocks = []
    blocks.append(("h1", f"MD Analysis Report — {system.system_type.value}"))
    blocks.append(("p", f"*Generated {date.today().isoformat()} by MolDynX Tools "
                        f"v{manifest.data.get('moldynx_version')}.*"))

    # -- system ----------------------------------------------------------- #
    blocks.append(("h2", "1. System"))
    comp_rows = [[c.ctype.value, c.n_residues, ", ".join(c.resnames[:8])]
                 for c in system.components.values() if c.n_residues]
    blocks.append(("table", (["Component", "Residues", "Resnames"], comp_rows)))
    blocks.append(("p", f"Detected system type: **{system.system_type.value}** · "
                        f"{system.n_protein_chains} protein chain(s) · "
                        f"{system.n_atoms:,} atoms."))

    # -- methods ---------------------------------------------------------- #
    blocks.append(("h2", "2. Methods & reproducibility"))
    libs = manifest.data.get("libraries", {})
    lib_str = ", ".join(f"{k} {v}" for k, v in libs.items() if v)
    git = manifest.data.get("git", {})
    blocks.append(("p", f"Analyses were run with MolDynX Tools (Python "
                        f"{manifest.data.get('python')}). Key libraries: {lib_str}. "
                        f"Git commit: `{git.get('commit')}`"
                        f"{' (dirty)' if git.get('dirty') else ''}. "
                        f"Full provenance in `manifest.json`."))

    # -- results ---------------------------------------------------------- #
    blocks.append(("h2", "3. Results"))
    key_rows = _key_result_rows(results)
    if key_rows:
        blocks.append(("table", (["Observable", "Value"], key_rows)))

    # -- figures ---------------------------------------------------------- #
    blocks.append(("h2", "4. Figures"))
    for name, summ in results.items():
        if isinstance(summ, dict) and summ.get("figure"):
            fig = ctx.config.figures_dir / f"{summ['figure']}.png"
            if fig.exists():
                blocks.append(("img", (fig, name)))

    # -- interpretation --------------------------------------------------- #
    blocks.append(("h2", "5. Interpretation"))
    for line in _interpretation(results):
        blocks.append(("p", line))

    # -- limitations ------------------------------------------------------ #
    blocks.append(("h2", "6. Limitations & provenance"))
    blocks.append(("p", "Convergence diagnostics are necessary but not sufficient "
                        "evidence of equilibration; consider replicate/longer runs. "
                        "All parameters, seeds, library versions, input fingerprints "
                        "and per-analysis runtimes are recorded in `manifest.json` / "
                        "`manifest.yaml` for exact reproduction."))

    # -- companion documents (written from the results files) ------------- #
    from moldynx.report import documents
    docs = documents.write_all(ctx, formats=[f for f in formats if f in ("md", "html")])
    doc_names = sorted({p.stem for p in docs})
    if doc_names:
        blocks.insert(2, ("h2", "Companion documents"))
        blocks.insert(3, ("p", " · ".join(f"[{n}]({n}.md)" for n in doc_names)))

    # -- render ----------------------------------------------------------- #
    ctx.config.report_dir.mkdir(parents=True, exist_ok=True)
    written = list(docs)
    md_text = _md(blocks)
    if "md" in formats:
        p = ctx.config.report_dir / "report.md"
        p.write_text(md_text, encoding="utf-8")
        written.append(p)
    html_text = documents.to_html(md_text.replace(".md)", ".html)"), ctx.config.report_dir,
                                  "MD analysis report")
    if "html" in formats:
        p = ctx.config.report_dir / "report.html"
        p.write_text(html_text, encoding="utf-8")
        written.append(p)
    if "pdf" in formats:
        try:
            from weasyprint import HTML
            pdf = ctx.config.report_dir / "report.pdf"
            HTML(string=html_text).write_pdf(str(pdf))
            written.append(pdf)
        except Exception:
            pass  # weasyprint not installed; MD/HTML still produced
    return written


def _key_result_rows(results: dict) -> list[list]:
    rows = []
    g = lambda d, *k: _nested(d, k)
    if "rmsd" in results:
        rows.append(["Backbone RMSD (mean)", f"{_fmt(g(results,'rmsd','backbone','mean'))} nm"])
        rows.append(["RMSD converged?", g(results, "rmsd", "convergence", "converged")])
    if "rog" in results:
        rows.append(["Radius of gyration (mean)", f"{_fmt(g(results,'rog','rg','mean'))} nm"])
    if "rmsf" in results:
        rows.append(["Most flexible residue", g(results, "rmsf", "most_flexible_resid")])
    if "sasa" in results:
        rows.append(["Total SASA (mean)", f"{_fmt(g(results,'sasa','total_sasa','mean'),1)} nm²"])
    if "pbc_validation" in results:
        rows.append(["PBC proof (all checks)",
                     "pass" if g(results, "pbc_validation", "all_checks_pass") else "FAIL"])
    if "interface" in results and g(results, "interface", "partners"):
        rows.append(["Interface: min. distance (mean)",
                     f"{_fmt(g(results,'interface','min_interface_dist_nm','mean'))} nm"])
        rows.append(["Interface: core residues",
                     " / ".join(str(x) for x in (g(results, "interface",
                                                     "n_core_interface_residues") or []))])
        rows.append(["Interface: persistent contacts", g(results, "interface", "n_persistent_contacts")])
    if "equilibration" in results:
        for v in (g(results, "equilibration", "verdicts") or [])[:3]:
            rows.append(["Preparation", v])
    return [r for r in rows if r[1] not in (None, "n/a")]


def _nested(d, keys):
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def _interpretation(results: dict) -> list[str]:
    lines = []
    conv = _nested(results, ("rmsd", "convergence", "converged"))
    if conv is not None:
        # an RMSD plateau alone is not evidence of stability or equilibrium
        lines.append(f"**RMSD time course.** The backbone RMSD "
                     f"{'levels off' if conv else 'had not levelled off'} over the analysed "
                     f"window; this alone does not establish stability — see the convergence, "
                     f"interface and stationarity results.")
    win = _nested(results, ("analysis_window", "binding_observables_stationary"))
    if win is not None:
        lines.append("**Stationarity.** " + (
            "Binding-related observables are stationary after "
            f"{_fmt(_nested(results, ('analysis_window', 'conservative_t0_ns')), 1)} ns."
            if win else "Binding-related observables are not stationary over the run: "
                        "averages describe a changing state (see the window selection)."))
    be = _nested(results, ("mmpbsa", "headline"))
    if be:
        lines.append("**Binding energy.** " + "; ".join(
            f"{h['method']} {h['window']}: {_fmt(h['mean'], 1)} ± {_fmt(h['sem'], 1)} kcal/mol"
            for h in be) + " — end-point estimates, not experimental affinities.")
    rg = _nested(results, ("rog",))
    if rg:
        init, fin = rg.get("rg_initial_nm"), rg.get("rg_final_nm")
        if init and fin:
            trend = "expanded" if fin > init + 0.15 else "compacted" if fin < init - 0.15 else "kept a stable size"
            lines.append(f"**Compactness.** The solute {trend} "
                         f"(Rg {_fmt(init)} → {_fmt(fin)} nm).")
    if not lines:
        lines.append("See the per-analysis CSVs and figures for detailed results.")
    return lines
