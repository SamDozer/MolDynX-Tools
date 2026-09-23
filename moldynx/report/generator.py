"""
Structured report builder.

Content is assembled as a list of typed blocks (heading/paragraph/table/image),
then rendered to Markdown and self-contained HTML from the same source, so the
two formats never drift.  PDF is attempted via weasyprint if installed.
"""

from __future__ import annotations

import base64
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


def _html(blocks) -> str:
    css = ("body{font-family:Arial,Helvetica,sans-serif;max-width:960px;margin:2rem auto;"
           "padding:0 1rem;color:#222;line-height:1.5}h1{border-bottom:3px solid #1f6f8b}"
           "h2{color:#1f6f8b;margin-top:2rem}table{border-collapse:collapse;margin:1rem 0}"
           "th,td{border:1px solid #ccc;padding:6px 10px;text-align:left}"
           "th{background:#f0f4f6}img{max-width:100%;border:1px solid #eee}"
           "figure{margin:1.5rem 0}figcaption{color:#666;font-size:.9em}")
    out = [f"<!doctype html><html><head><meta charset='utf-8'>"
           f"<style>{css}</style></head><body>"]
    for kind, payload in blocks:
        if kind == "h1":
            out.append(f"<h1>{payload}</h1>")
        elif kind == "h2":
            out.append(f"<h2>{payload}</h2>")
        elif kind == "p":
            out.append(f"<p>{_inline_html(payload)}</p>")
        elif kind == "table":
            headers, rows = payload
            out.append("<table><tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr>")
            for r in rows:
                out.append("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>")
            out.append("</table>")
        elif kind == "img":
            path, caption = payload
            p = Path(path)
            if p.exists():
                b64 = base64.b64encode(p.read_bytes()).decode()
                out.append(f"<figure><img src='data:image/png;base64,{b64}'/>"
                           f"<figcaption>{caption}</figcaption></figure>")
    out.append("</body></html>")
    return "\n".join(out)


def _inline_html(text: str) -> str:
    return text.replace("**", "")  # markdown bold markers -> plain (kept simple)


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

    # -- render ----------------------------------------------------------- #
    ctx.config.report_dir.mkdir(parents=True, exist_ok=True)
    written = []
    if "md" in formats:
        p = ctx.config.report_dir / "report.md"
        p.write_text(_md(blocks), encoding="utf-8")
        written.append(p)
    if "html" in formats:
        p = ctx.config.report_dir / "report.html"
        p.write_text(_html(blocks), encoding="utf-8")
        written.append(p)
    if "pdf" in formats:
        try:
            from weasyprint import HTML
            pdf = ctx.config.report_dir / "report.pdf"
            HTML(string=_html(blocks)).write_pdf(str(pdf))
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
        lines.append(f"**Stability.** Backbone RMSD "
                     f"{'reached a plateau' if conv else 'had not fully plateaued'} "
                     f"over the analysed window.")
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
