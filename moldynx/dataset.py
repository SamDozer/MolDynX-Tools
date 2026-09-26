"""
Assemble a shareable, self-verifying dataset from a MolDynX run, verify it, and zip it.

Layout::

    <dataset>/
      README.md / .html          generated overview
      0_raw_production/          canonical raw files (only with include_raw) + MANIFEST.sha256.tsv
      1_report/                  report + companion documents (.md/.html) + figures/
      2_tables/  3_results/      summary tables, one CSV/JSON per analysis
      4_trajectory/              solute trajectory (core.pdb/.xtc) + core_meta.json
      5_validation/              intake report, PBC proof, window selection, verify_dataset.py
      6_workflow/                run manifest, configuration
      7_binding_energy/          gmx_MMPBSA inputs, run script and raw outputs (if prepared/run)

``verify_dataset.py`` (standard library; MDAnalysis optional) checks the layout, that
every Markdown link/image resolves, the SHA-256 manifest, the PBC proof record and --
if MDAnalysis is available -- recomputes Rg from the shipped trajectory.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
import zipfile
from pathlib import Path

LAYOUT = ["1_report", "1_report/figures", "2_tables", "3_results", "4_trajectory",
          "5_validation", "6_workflow"]
STUB = "RAW_DATA_NOT_INCLUDED.md"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(1 << 24):
            h.update(b)
    return h.hexdigest()


def _copytree(src: Path, dst: Path) -> None:
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns("*_offsets.npz", "__pycache__", "*.lock"))


def _relink(md: str) -> str:
    """Run layout (report/ next to figures/) -> dataset layout (1_report/figures)."""
    return md.replace("](../figures/", "](figures/")


def assemble(run_dir: str | Path, out_dir: str | Path, include_raw: bool = False,
             title: str | None = None) -> Path:
    from moldynx.report.documents import to_html
    run, out = Path(run_dir), Path(out_dir)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"{out} exists and is not empty")
    for d in LAYOUT:
        (out / d).mkdir(parents=True, exist_ok=True)
    rep = run / "report"
    for f in sorted(rep.glob("*.md")):
        text = _relink(f.read_text(encoding="utf-8"))
        (out / "1_report" / f.name).write_text(text, encoding="utf-8")
    _copytree(run / "figures", out / "1_report" / "figures")
    for f in sorted((out / "1_report").glob("*.md")):
        text = f.read_text(encoding="utf-8").replace(".md)", ".html)")
        f.with_suffix(".html").write_text(to_html(text, f.parent, f.stem), encoding="utf-8")
    _copytree(run / "tables", out / "2_tables")
    for f in sorted((run / "results").glob("*")):
        if f.is_file():
            dst = "5_validation" if f.name.startswith(("pbc_", "window_selection")) else "3_results"
            shutil.copy2(f, out / dst / f.name)
        elif f.is_dir():
            _copytree(f, out / "3_results" / f.name)
    for f in ("core.pdb", "core.xtc", "core_meta.json"):
        if (run / "data" / f).exists():
            shutil.copy2(run / "data" / f, out / "4_trajectory" / f)
    _copytree(run / "intake", out / "5_validation")
    for f in ("manifest.json", "manifest.yaml"):
        if (run / f).exists():
            shutil.copy2(run / f, out / "6_workflow" / f)
    if (run / "binding_energy").is_dir():
        _copytree(run / "binding_energy", out / "7_binding_energy")

    manifest = _load_json(run / "manifest.json") or {}
    intake = _load_json(run / "intake" / "intake_manifest.json") or {}
    raw_rows = _raw_files(intake)
    (out / "0_raw_production").mkdir(exist_ok=True)
    lines = ["file\tbytes\tsha256\trole\toriginal_path"]
    for role, p in raw_rows:
        digest = sha256(p) if include_raw else "not-computed"
        if include_raw:
            shutil.copy2(p, out / "0_raw_production" / p.name)
        lines.append(f"{p.name}\t{p.stat().st_size}\t{digest}\t{role}\t{p}")
    (out / "0_raw_production" / "MANIFEST.sha256.tsv").write_text("\n".join(lines) + "\n",
                                                                  encoding="utf-8")
    if not include_raw:
        (out / "0_raw_production" / STUB).write_text(
            "# Raw simulation files are not included\n\nThe manifest lists the canonical files "
            "this dataset was built from (sizes and original paths). Copy them into this folder "
            "and run `python 5_validation/verify_dataset.py --full` to hash and check them.\n",
            encoding="utf-8")
    shutil.copy2(Path(__file__).parent / "verify_template.py",
                 out / "5_validation" / "verify_dataset.py")
    readme = dataset_readme(title or run.name, manifest, intake, out)
    (out / "README.md").write_text(readme, encoding="utf-8")
    (out / "README.html").write_text(to_html(readme.replace(".md)", ".html)"), out, "README"),
                                     encoding="utf-8")
    return out


def _load_json(p: Path) -> dict | None:
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _raw_files(intake: dict) -> list[tuple[str, Path]]:
    fsd = intake.get("fileset", {})
    rows = []
    for role in ("trajectory", "topology", "energy", "log", "structure", "gmx_top", "index"):
        if fsd.get(role):
            rows.append((role, Path(fsd[role])))
    for stage, kinds in (fsd.get("stages") or {}).items():
        if stage in ("em", "nvt", "npt"):
            for kind in ("log", "edr", "tpr"):
                for p in kinds.get(kind, [])[:1]:
                    rows.append((f"{stage}_{kind}", Path(p)))
    return [(r, p) for r, p in rows if p.exists()]


def dataset_readme(title: str, manifest: dict, intake: dict, out: Path) -> str:
    sysd = manifest.get("system", {})
    an = manifest.get("analyses", {})
    ev = intake.get("evidence", {})
    L = [f"# {title}", "",
         f"Molecular dynamics analysis dataset produced by **MolDynX Tools** "
         f"v{manifest.get('moldynx_version', '?')} (formerly mdforge). Every number in the "
         "documents below is computed from files in this folder.", "",
         "## System at a glance", "", "| Property | Value |", "|---|---|",
         f"| System type | {sysd.get('system_type', '—')} |",
         f"| Atoms (run input) | {sysd.get('n_atoms', '—')} |",
         f"| Protein chains | {sysd.get('n_protein_chains', '—')} "
         f"{sysd.get('protein_chain_lengths', '')} |",
         f"| Ions | {', '.join(f'{k} {v}' for k, v in (sysd.get('ion_counts') or {}).items()) or '—'} |",
         f"| Production run | {ev.get('trajectory_choice', '—')} |", ""]
    pl = ev.get("production_log") or {}
    if pl:
        L += [f"Production log: {pl.get('sessions')} mdrun session(s), "
              f"{(pl.get('simulated_ps') or 0) / 1000:g} ns simulated.", ""]
    L += ["## Key results", "", "| Analysis | Status | Headline |", "|---|---|---|"]
    for name, rec in sorted(an.items()):
        s = rec.get("summary") or {}
        head = ""
        if name == "pbc_validation":
            head = "all proof checks pass" if s.get("all_checks_pass") else "see PBC_VALIDATION"
        elif name == "mmpbsa" and s.get("headline"):
            head = "; ".join(f"{h['method']} {h['window']} {h['mean']:.1f} ± {h['sem']:.1f}"
                             for h in s["headline"])
        elif name == "interface" and s.get("partners"):
            head = f"{s.get('n_persistent_contacts')} persistent contacts"
        elif name == "analysis_window":
            head = s.get("statement", "")
        status = s.get("status", rec.get("status"))
        L.append(f"| {name} | {status} | {head} |")
    docs = sorted(p.name for p in (out / "1_report").glob("*.md"))
    L += ["", "## Documents", ""] + [f"- [{d[:-3]}](1_report/{d})" for d in docs]
    L += ["", "## Folder map", "", "| Folder | Contents |", "|---|---|",
          "| `0_raw_production/` | manifest of the canonical raw files (and the files, if included) |",
          "| `1_report/` | report and companion documents, figures |",
          "| `2_tables/`, `3_results/` | summary tables; one CSV/JSON per analysis |",
          "| `4_trajectory/` | processed solute trajectory and its chain/PBC record |",
          "| `5_validation/` | intake report, PBC proof, window selection, `verify_dataset.py` |",
          "| `6_workflow/` | run manifest: versions, parameters, input fingerprints |",
          "| `7_binding_energy/` | gmx_MMPBSA inputs, run script, outputs (if prepared) |",
          "", "## Verifying this dataset", "",
          "```bash", "python 5_validation/verify_dataset.py          # --full also hashes raw files",
          "```", "", "## Caveats", "",
          "- Binding energies are end-point estimates, not experimental affinities.",
          "- An RMSD plateau alone is not evidence of stability; see the stationarity analysis.",
          "- Metadata not recorded in the simulation files (e.g. force-field variant) is stated "
          "as not recorded, never guessed."]
    return "\n".join(L) + "\n"


def package(dataset_dir: str | Path, zip_path: str | Path) -> dict:
    """Deterministic zip (one top folder), then extract and run the shipped verifier."""
    import subprocess
    import sys
    src, zp = Path(dataset_dir), Path(zip_path)
    files = sorted(p for p in src.rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for p in files:
            z.write(p, f"{src.name}/{p.relative_to(src).as_posix()}")
    with zipfile.ZipFile(zp) as z:
        bad = z.testzip()
        with tempfile.TemporaryDirectory() as tmp:
            z.extractall(tmp)
            r = subprocess.run([sys.executable, str(Path(tmp) / src.name / "5_validation" /
                                                    "verify_dataset.py")],
                               capture_output=True, text=True)
    return {"zip": str(zp), "files": len(files), "bytes": zp.stat().st_size,
            "crc_ok": bad is None, "verify_exit": r.returncode, "verify_output": r.stdout[-2000:]}


def verify(dataset_dir: str | Path, full: bool = False) -> int:
    import subprocess
    import sys
    args = [sys.executable, str(Path(dataset_dir) / "5_validation" / "verify_dataset.py")]
    return subprocess.call(args + (["--full"] if full else []))

