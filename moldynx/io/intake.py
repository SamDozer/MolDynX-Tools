"""
Intake: the first stage of every run.

Discovers the files (by evidence), validates them, summarises every simulation
stage from its log, measures the cluster-clock offset, reports what each missing
file disables, and writes

* ``intake_manifest.json`` -- machine-readable, consumed by later stages;
* ``INTAKE_REPORT.md``     -- what a researcher reads before trusting the run.

Nothing here modifies the simulation folder.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from moldynx.io import gromacs
from moldynx.io.discovery import FileSet, discover_files
from moldynx.io.validation import ValidationResult, validate_fileset

STAGE_LABEL = {"setup": "Setup", "em": "Energy minimisation", "nvt": "NVT equilibration",
               "npt": "NPT equilibration", "production": "Production"}


@dataclass
class IntakeResult:
    fileset: FileSet
    validation: ValidationResult
    stages: dict[str, dict] = field(default_factory=dict)
    clock: dict = field(default_factory=dict)
    system: dict | None = None

    @property
    def ok(self) -> bool:
        return self.validation.ok

    def to_dict(self) -> dict:
        return {"fileset": self.fileset.to_dict(), "evidence": self.fileset.evidence,
                "validation": {"ok": self.validation.ok, "errors": self.validation.errors,
                               "warnings": self.validation.warnings},
                "capabilities": self.validation.capabilities,
                "stages": self.stages, "clock": self.clock, "system": self.system}


# --------------------------------------------------------------------------- #
def _stage_summary(fs: FileSet) -> dict[str, dict]:
    out = {}
    for stage in ("em", "nvt", "npt", "production"):
        log = fs.log if stage == "production" else fs.stage_file(stage, "log")
        files = {k: [str(p) for p in v] for k, v in fs.stages.get(stage, {}).items()}
        entry: dict = {"files": files, "mdp": str(fs.mdp[stage]) if stage in fs.mdp else None}
        if log is not None and fs.deep:
            try:
                info = gromacs.parse_log(log)
                m = info.mdp
                entry["log"] = {
                    "path": str(log), "gromacs": info.gromacs_version,
                    "sessions": info.n_sessions, "started": info.started[:1],
                    "finished": info.finished[-1:], "integrator": m.get("integrator"),
                    "dt_ps": m.get("dt"), "nsteps_first_block": m.get("nsteps"),
                    "nsteps_last_block": info.mdp_last.get("nsteps"),
                    "simulated_ps": info.simulated_ps, "tcoupl": m.get("tcoupl"),
                    "ref_t_K": info.ref_t, "pcoupl": m.get("pcoupl"),
                    "tau_p_ps": m.get("tau-p"), "continuation": m.get("continuation"),
                    "counts": info.counts,
                    "minimization": info.minimization.__dict__ if info.minimization else None,
                }
            except (OSError, ValueError) as exc:
                entry["log"] = {"path": str(log), "error": str(exc)}
        out[stage] = entry
    return out


def _clock_offset(fs: FileSet) -> dict:
    """Cluster clock (log 'Finished mdrun') minus local file mtime, per log, in hours."""
    logs = [fs.log] + [fs.stage_file(s, "log") for s in ("em", "nvt", "npt")]
    rows = []
    for log in [p for p in logs if p is not None]:
        try:
            info = gromacs.parse_log(log)
        except (OSError, ValueError):
            continue
        if not info.finished:
            continue
        cluster = datetime.fromisoformat(info.finished[-1])
        local = datetime.fromtimestamp(log.stat().st_mtime)
        rows.append({"log": log.name, "cluster_finished": info.finished[-1],
                     "file_mtime": local.isoformat(timespec="seconds"),
                     "offset_h": round((cluster - local).total_seconds() / 3600, 3)})
    out: dict = {"per_log": rows}
    if rows:
        offs = [r["offset_h"] for r in rows]
        out["spread_h"] = round(max(offs) - min(offs), 3)
        out["constant"] = out["spread_h"] <= 0.1      # within 6 minutes
        out["offset_h"] = round(sum(offs) / len(offs), 2)
    return out


def run_intake(input_dir: str | Path, deep: bool = True, detect: bool = False,
               allow_ambiguous: bool = False,
               include_dirs: list[str | Path] | None = None) -> IntakeResult:
    fs = discover_files(input_dir, deep=deep, include_dirs=include_dirs)
    val = validate_fileset(fs, allow_ambiguous=allow_ambiguous)
    res = IntakeResult(fileset=fs, validation=val)
    res.stages = _stage_summary(fs)
    if deep:
        res.clock = _clock_offset(fs)
        if res.clock.get("constant") is False:
            val.warnings.append(
                f"the cluster-clock offset varies by {res.clock['spread_h']} h between logs — "
                "files may come from different runs or were copied at different times")
    if detect and val.ok and fs.topology is not None:
        from moldynx.core.system import detect_system
        res.system = detect_system(fs.topology).to_dict()
    return res


# --------------------------------------------------------------------------- #
# writers
# --------------------------------------------------------------------------- #
def _rel(p, root: Path) -> str:
    if not p:
        return "—"
    try:
        return Path(p).relative_to(root).as_posix()
    except ValueError:
        return str(p)


def write_intake(res: IntakeResult, out_dir: str | Path) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    js = out_dir / "intake_manifest.json"
    js.write_text(json.dumps(res.to_dict(), indent=2, default=str), encoding="utf-8")
    md = out_dir / "INTAKE_REPORT.md"
    md.write_text(render_report(res), encoding="utf-8")
    return [js, md]


def render_report(res: IntakeResult) -> str:
    fs, val = res.fileset, res.validation
    root = fs.directory
    ev = fs.evidence
    L = [f"# Intake report — `{root.name}`", ""]
    L.append(f"**Status:** {'ready' if val.ok else 'BLOCKED — resolve the errors below'}"
             f" · evidence {'verified from file headers' if ev.get('verified') else 'NOT verified (names/stages only)'}")
    L.append("")
    if val.errors:
        L += ["## Errors", ""] + [f"- {e}" for e in val.errors] + [""]
    L += ["## Canonical production run", "",
          "| Role | File | Why |", "|---|---|---|",
          f"| trajectory | `{_rel(fs.trajectory, root)}` | {ev.get('trajectory_choice', '—')} |",
          f"| run input | `{_rel(fs.topology, root)}` | {ev.get('topology_choice', '—')} |",
          f"| energies | `{_rel(fs.energy, root)}` | production-stage energy file |",
          f"| log | `{_rel(fs.log, root)}` | production-stage mdrun log |",
          f"| structure | `{_rel(fs.structure, root)}` | production/setup structure (never a crash dump) |",
          f"| topology (#include) | `{_rel(fs.gmx_top, root)}` | nearest `topol.top` to the run input |",
          f"| force-field folder | `{_rel(fs.toppar, root)}` | next to `topol.top` |", ""]
    pl = ev.get("production_log")
    if pl:
        L += [f"Production log: {pl['sessions']} mdrun session(s); final statistics over "
              f"{pl['last_statistics_steps']} steps"
              + (f" = {pl['simulated_ps'] / 1000:g} ns simulated" if pl.get("simulated_ps") else "")
              + (f"; run input declared nsteps {pl['first_block_nsteps']} in the first session "
                 f"and {pl['last_block_nsteps']} in the last (extended)"
                 if pl.get("first_block_nsteps") != pl.get("last_block_nsteps") else "") + ".", ""]
        if pl.get("extended"):
            L += [f"> **Extended run:** {pl['extended']}.", ""]
    L += ["## Stages", "", "| Stage | Log | Ensemble / settings | Outcome |", "|---|---|---|---|"]
    for stage in ("em", "nvt", "npt", "production"):
        s = res.stages.get(stage, {})
        lg = s.get("log")
        if not lg or "error" in lg:
            present = ", ".join(sorted(s.get("files", {}))) or "none"
            L.append(f"| {STAGE_LABEL[stage]} | — | files present: {present} | not audited |")
            continue
        if lg.get("minimization"):
            em = lg["minimization"]
            outcome = (f"{em['outcome']} in {em['steps']} steps; Fmax {em['fmax']:.4g}"
                       if em.get("outcome") else "outcome not found")
            if em.get("reached_emtol") is False:
                outcome += " — **requested tolerance not reached**"
            settings = f"{lg['integrator']}"
        else:
            settings = (f"{lg['integrator']}, dt {lg['dt_ps']} ps, {lg['tcoupl']} "
                        f"{'/'.join(str(t) for t in lg['ref_t_K'][:1])} K, pcoupl {lg['pcoupl']}")
            sim = lg.get("simulated_ps")
            outcome = (f"{sim / 1000:g} ns" if sim and sim >= 1000 else f"{sim:g} ps" if sim else "—")
            outcome += f", {lg['sessions']} session(s)"
        c = lg["counts"]
        if c["lincs_warnings"] or c["fatal_errors"]:
            outcome += f"; LINCS warnings {c['lincs_warnings']}, fatal errors {c['fatal_errors']}"
        L.append(f"| {STAGE_LABEL[stage]} | `{_rel(lg['path'], root)}` | {settings} | {outcome} |")
    L.append("")
    temps = {st: res.stages[st]["log"]["ref_t_K"][:1] for st in ("nvt", "npt", "production")
             if res.stages.get(st, {}).get("log", {}).get("ref_t_K")}
    if len({tuple(v) for v in temps.values()}) > 1:
        L += ["> **Temperature changes between stages:** "
              + ", ".join(f"{STAGE_LABEL[k]} {v[0]:g} K" for k, v in temps.items())
              + ". Production starts from a state equilibrated at a different temperature.", ""]
    L += ["## Capabilities", "", "| Capability | Status | Consequence |", "|---|---|---|"]
    for c in val.capabilities:
        cons = c["enables"] if c["status"] == "available" else \
            f"{c['without']} (missing: {', '.join(c['missing'])})"
        L.append(f"| {c['capability']} | {c['status']} | {cons} |")
    L.append("")
    if res.clock.get("per_log"):
        L += ["## Clock", "",
              f"Cluster clock minus local file time: **{res.clock['offset_h']:+g} h**, "
              f"{'constant' if res.clock['constant'] else 'NOT constant'} across "
              f"{len(res.clock['per_log'])} log(s) (spread {res.clock['spread_h']} h).", ""]
    if val.warnings:
        L += ["## Warnings", ""] + [f"- {w}" for w in val.warnings] + [""]
    miss = ev.get("job_script_inputs_not_found")
    if miss:
        L += ["## Inputs referenced by job scripts but not in this folder", "",
              "| File | Referenced by |", "|---|---|"]
        L += [f"| `{k}` | {', '.join(v)} |" for k, v in sorted(miss.items())] + [""]
    L += ["## Ignored", ""]
    reasons: dict[str, int] = {}
    for _p, why in fs.ignored:
        reasons[why] = reasons.get(why, 0) + 1
    L += [f"- {n} × {why}" for why, n in sorted(reasons.items())] or ["- nothing"]
    if fs.crash_dumps:
        L.append(f"- {len(fs.crash_dumps)} minimisation crash-dump structures (`stepN{{b,c}}.pdb`)"
                 " — accounted for in the equilibration audit, never used as structures")
    if fs.unclassified:
        L.append(f"- {len(fs.unclassified)} unclassified simulation file(s): "
                 + ", ".join(f"`{p.name}`" for p in fs.unclassified[:8])
                 + (" …" if len(fs.unclassified) > 8 else ""))
    if res.system:
        s = res.system
        L += ["", "## System", "",
              f"{s['system_type']} · {s['n_atoms']:,} atoms · "
              f"{s['n_protein_chains']} protein chain(s)", ""]
        if s.get("chains"):
            L += ["| # | segid | residues | atoms | resid range |", "|---|---|---|---|---|"]
            L += [f"| {c['index']} | {c['segid']} | {c['n_residues']} | {c['n_atoms']:,} | "
                  f"{c['resid_first']}–{c['resid_last']} |" for c in s["chains"]]
        if s.get("ion_counts"):
            L += ["", "Ions: " + ", ".join(f"{k} {v}" for k, v in s["ion_counts"].items())]
    return "\n".join(L) + "\n"
