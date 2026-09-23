"""
Recursive discovery of GROMACS simulation files -- decided by evidence, not size.

The user provides a directory. This module

1. classifies every file by **stage** (setup / minimisation / NVT / NPT /
   production) and **kind** (tpr, xtc, edr, log, mdp, ...), ignoring backups,
   caches and folders that hold *earlier analysis outputs*;
2. picks the canonical production run from **evidence**: the trajectory's atom
   count must equal the run input's, the longest complete, uniformly spaced
   trajectory wins, and the production log's final step is recorded -- file
   size and names are only a last-resort tie-breaker, and are labelled so;
3. records every piece of evidence, warning and ambiguity on the returned
   :class:`FileSet`, which intake/validation turn into a report.

Why not "largest file per role"? On a real CHARMM-GUI folder the equilibration
``.tpr`` files are larger than the production one and minimisation crash dumps
(``stepNc.pdb``) are the largest ``.pdb`` files, so size picks the wrong run.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from moldynx.io import gromacs

# --------------------------------------------------------------------------- #
# Classification tables (data, so other naming schemes can be added)
# --------------------------------------------------------------------------- #
_TRAJ_EXT = {".xtc", ".trr", ".dcd", ".nc", ".netcdf", ".h5", ".xyz"}
_TOPOLOGY_EXT = {".tpr", ".psf", ".prmtop", ".parm7"}
_STRUCTURE_EXT = {".gro", ".pdb", ".g96"}

# stage -> regex on the lower-case file name
STAGE_PATTERNS: list[tuple[str, str]] = [
    ("em", r"step4\.0|(^|[_.-])(em|min|minim|minimi[sz]ation)([_.-]|$)"),
    ("nvt", r"step4\.1|(^|[_.-])nvt([_.-]|$)"),
    ("npt", r"step4\.2|(^|[_.-])npt([_.-]|$)"),
    ("production", r"step5|(^|[_.-])(prod|production|md)([_.-]|$)"),
    ("setup", r"^step[123]_|^topol\.top$|^index\.ndx$|charmm-?gui|^readme$"),
]
STAGES = ("setup", "em", "nvt", "npt", "production")

# folders that contain derived / earlier-analysis outputs, never inputs
_DERIVED_DIR = re.compile(r"analysis|diagnostic|^results?$|_results$|^mmpbsa|^__pycache__$|^\.git$",
                          re.I)
_CRASH_DUMP = re.compile(r"^step-?\d+[bc]\.pdb$", re.I)
_BACKUP = re.compile(r"^#.*#$")
_CACHE = re.compile(r"_offsets\.npz$|\.lock$", re.I)
_EXTEND_HINT = re.compile(r"(ext|extend|extended|\d+ns)(?=\.tpr$)|_ext[_.]", re.I)


def classify_stage(name: str) -> str | None:
    low = name.lower()
    for stage, pat in STAGE_PATTERNS:
        if re.search(pat, low):
            return stage
    return None


@dataclass
class FileSet:
    """Resolved semantic file roles for a simulation, with the evidence behind them."""
    directory: Path
    trajectory: Path | None = None
    topology: Path | None = None            # the run input that produced the trajectory
    structure: Path | None = None           # single-frame .gro/.pdb (never a crash dump)
    energy: Path | None = None              # production .edr
    index: Path | None = None               # .ndx
    gmx_top: Path | None = None             # topol.top (#includes)
    toppar: Path | None = None              # folder the topol.top includes from
    log: Path | None = None                 # production md .log
    all_trajectories: list[Path] = field(default_factory=list)
    all_files: dict[str, list[Path]] = field(default_factory=dict)
    # -- new in 0.3 ------------------------------------------------------- #
    stages: dict[str, dict[str, list[Path]]] = field(default_factory=dict)
    mdp: dict[str, Path] = field(default_factory=dict)          # stage -> .mdp
    job_scripts: list[Path] = field(default_factory=list)
    crash_dumps: list[Path] = field(default_factory=list)
    setup_package: list[Path] = field(default_factory=list)     # CHARMM-GUI zip/tgz
    ignored: list[tuple[str, str]] = field(default_factory=list)
    unclassified: list[Path] = field(default_factory=list)
    evidence: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    ambiguities: list[str] = field(default_factory=list)
    deep: bool = False

    def available_keys(self) -> set[str]:
        """Semantic keys that are present (used to gate analyses)."""
        keys = set()
        for role in ("trajectory", "topology", "structure", "energy", "index", "gmx_top",
                     "toppar", "log"):
            if getattr(self, role) is not None:
                keys.add(role)
        if self.topology or self.structure:
            keys.add("topology")
        for stage in ("em", "nvt", "npt"):
            files = self.stages.get(stage, {})
            if files.get("log") or files.get("edr"):
                keys.add(stage)
        if any(self.stages.get(s, {}).get("tpr") for s in ("em", "nvt", "npt")):
            keys.add("stage_tprs")
        if self.job_scripts:
            keys.add("job_scripts")
        if self.mdp:
            keys.add("mdp")
        if self.log is not None:
            keys.add("production_log")
        return keys

    def stage_file(self, stage: str, kind: str) -> Path | None:
        files = self.stages.get(stage, {}).get(kind, [])
        return files[0] if files else None

    def to_dict(self) -> dict:
        out = {role: (str(getattr(self, role)) if getattr(self, role) else None)
               for role in ("directory", "trajectory", "topology", "structure",
                            "energy", "index", "gmx_top", "toppar", "log")}
        out["stages"] = {s: {k: [str(p) for p in v] for k, v in d.items()}
                         for s, d in self.stages.items()}
        out["mdp"] = {k: str(v) for k, v in self.mdp.items()}
        out["job_scripts"] = [str(p) for p in self.job_scripts]
        out["crash_dumps"] = len(self.crash_dumps)
        out["ignored"] = [list(x) for x in self.ignored]
        out["unclassified"] = [str(p) for p in self.unclassified]
        out["warnings"] = list(self.warnings)
        out["ambiguities"] = list(self.ambiguities)
        return out


# --------------------------------------------------------------------------- #
# scanning
# --------------------------------------------------------------------------- #
def _is_derived_dir(d: Path) -> str | None:
    if _DERIVED_DIR.search(d.name):
        return f"folder name '{d.name}' marks derived/analysis output"
    man = d / "manifest.json"
    if man.exists():
        try:
            data = json.loads(man.read_text(encoding="utf-8"))
            if "moldynx_version" in data or "mdforge_version" in data:
                return "contains a MolDynX/mdforge run manifest"
        except (OSError, ValueError):
            pass
    return None


def _walk(root: Path, include_dirs: set[Path], fs: FileSet) -> list[Path]:
    files: list[Path] = []
    stack = [root]
    while stack:
        d = stack.pop()
        try:
            entries = sorted(d.iterdir())
        except OSError as exc:
            fs.warnings.append(f"cannot read folder {d}: {exc}")
            continue
        for p in entries:
            if p.is_dir():
                reason = None if p.resolve() in include_dirs else _is_derived_dir(p)
                if reason:
                    fs.ignored.append((str(p), reason))
                else:
                    stack.append(p)
            elif p.is_file():
                files.append(p)
    return files


def _script_mentions_gromacs(p: Path) -> bool:
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    return bool(re.search(r"\b(grompp|mdrun|convert-tpr|trjconv)\b", text))


def discover_files(directory: str | Path, deep: bool = True,
                   include_dirs: list[str | Path] | None = None) -> FileSet:
    """
    Recursively scan ``directory`` and resolve the canonical file for each role.

    ``deep`` reads trajectory frame headers, run-input headers and logs (read-only;
    a 6 GB XTC takes seconds). With ``deep=False`` only names/stages are used and
    every choice is marked unverified.
    """
    directory = Path(directory)
    if not directory.exists():
        raise FileNotFoundError(f"Input directory does not exist: {directory}")
    fs = FileSet(directory=directory, deep=deep)
    incl = {Path(p).resolve() for p in (include_dirs or [])}

    files = _walk(directory, incl, fs)
    by_ext: dict[str, list[Path]] = {}
    stages: dict[str, dict[str, list[Path]]] = {s: {} for s in STAGES}
    for p in files:
        name = p.name
        if _BACKUP.match(name):
            fs.ignored.append((str(p), "GROMACS backup (#…#)"))
            continue
        if _CACHE.search(name):
            fs.ignored.append((str(p), "cache/lock file"))
            continue
        ext = p.suffix.lower()
        by_ext.setdefault(ext, []).append(p)
        if _CRASH_DUMP.match(name):
            fs.crash_dumps.append(p)
            continue
        if ext in (".zip", ".tgz", ".gz") and re.search(r"charmm|gromacs", name, re.I):
            fs.setup_package.append(p)
            continue
        if ext in (".sh", ".slurm", ".pbs", ".job") and _script_mentions_gromacs(p):
            fs.job_scripts.append(p)
            continue
        if p.parent.name.lower() == "toppar" and ext == ".itp":
            stages["setup"].setdefault("itp", []).append(p)
            continue
        stage = classify_stage(name)
        kind = ext.lstrip(".") or name.lower()
        if stage is None:
            if ext in _TRAJ_EXT | _TOPOLOGY_EXT | {".edr", ".log", ".gro", ".mdp", ".top", ".ndx"}:
                fs.unclassified.append(p)
            continue
        stages[stage].setdefault(kind, []).append(p)

    fs.stages = {s: d for s, d in stages.items() if d}
    fs.all_files = by_ext
    fs.all_trajectories = sorted(p for e in _TRAJ_EXT for p in by_ext.get(e, []))
    for stage in STAGES:
        mdps = fs.stages.get(stage, {}).get("mdp", [])
        # mdout.mdp is grompp's record; the CHARMM-GUI stage file is the input
        inputs = [m for m in mdps if m.name.lower() != "mdout.mdp"] or mdps
        if inputs:
            fs.mdp[stage] = sorted(inputs, key=lambda m: (m.name.lower() == "mdout.mdp", str(m)))[0]
    for p in by_ext.get(".mdp", []):
        if p.name.lower() == "mdout.mdp":
            fs.evidence.setdefault("mdout_mdp", []).append(str(p))

    _resolve_production(fs)
    _resolve_setup(fs)
    _check_job_script_references(fs, files)
    return fs


# --------------------------------------------------------------------------- #
# canonical production run
# --------------------------------------------------------------------------- #
def _candidates(fs: FileSet, kinds: set[str]) -> list[Path]:
    prod = [p for k in kinds for p in fs.stages.get("production", {}).get(k, [])]
    if prod:
        return sorted(prod)
    return sorted(p for p in fs.unclassified if p.suffix.lower().lstrip(".") in kinds)


def _resolve_production(fs: FileSet) -> None:
    ev = fs.evidence
    trajs = _candidates(fs, {e.lstrip(".") for e in _TRAJ_EXT})
    tprs = _candidates(fs, {"tpr"})
    other_tops = _candidates(fs, {"psf", "prmtop", "parm7"})
    # --- evidence ------------------------------------------------------------ #
    traj_info = {str(p): gromacs.scan_xtc(p) if (fs.deep and p.suffix.lower() == ".xtc") else None
                 for p in trajs}
    tpr_info = {str(p): gromacs.tpr_header(p) if fs.deep else None for p in tprs}
    ev["trajectory_candidates"] = {k: (v.to_dict() if v else None) for k, v in traj_info.items()}
    ev["tpr_candidates"] = {k: (v.to_dict() if v else None) for k, v in tpr_info.items()}

    readable = {k: v for k, v in traj_info.items() if v and v.readable and not v.truncated}
    verified = bool(readable)
    if readable:
        tpr_atoms = {k: v.natoms for k, v in tpr_info.items() if v and v.readable}
        matched = {k: v for k, v in readable.items()
                   if not tpr_atoms or v.natoms in tpr_atoms.values()}
        if tpr_atoms and not matched:
            fs.ambiguities.append(
                "no production trajectory has the same atom count as any production run input "
                f"(trajectories: { {Path(k).name: v.natoms for k, v in readable.items()} }; "
                f"run inputs: { {Path(k).name: n for k, n in tpr_atoms.items()} })")
            matched = readable
        span = {k: (v.t_last_ps or 0) - (v.t_first_ps or 0) for k, v in matched.items()}
        best = max(span.values())
        winners = [k for k, s in span.items() if abs(s - best) < 1e-6]
        if len(winners) > 1:
            sizes = {Path(k).name: Path(k).stat().st_size for k in winners}
            fs.ambiguities.append(f"{len(winners)} production trajectories cover the same span "
                                  f"({best:.0f} ps): {sizes}")
        chosen = sorted(winners, key=lambda k: Path(k).stat().st_size, reverse=True)[0]
        fs.trajectory = Path(chosen)
        t = matched[chosen]
        ev["trajectory_choice"] = (f"longest complete trajectory with a matching atom count: "
                                   f"{t.n_frames} frames, {t.t_first_ps:g}-{t.t_last_ps:g} ps, "
                                   f"dt {t.dt_ps:g} ps, uniform={t.uniform_spacing}, "
                                   f"{t.natoms} atoms")
        for k, v in traj_info.items():
            if v and v.truncated:
                fs.warnings.append(f"trajectory {Path(k).name} ends in an incomplete frame")
        # --- run input: same atom count; extension hint only breaks ties -------- #
        same = [k for k, n in tpr_atoms.items() if n == t.natoms]
        if same:
            fs.topology = Path(_pick_tpr(same, fs))
            ev["topology_choice"] = (f"run input with the trajectory's atom count ({t.natoms})"
                                     + ("" if len(same) == 1 else
                                        f"; {len(same)} candidates with identical atom counts "
                                        f"({', '.join(Path(s).name for s in same)}) -- "
                                        "the extended/latest one was taken, its inputrec is "
                                        "compared in the provenance audit"))
            if len(same) > 1:
                fs.warnings.append(ev["topology_choice"])
    elif trajs:
        # unreadable (or not XTC) -> fall back to stage + size, labelled unverified
        fs.trajectory = max(trajs, key=lambda p: p.stat().st_size)
        ev["trajectory_choice"] = "UNVERIFIED: largest production-stage trajectory (headers not read)"
    if fs.topology is None and tprs:
        fs.topology = Path(_pick_tpr([str(p) for p in tprs], fs))
        ev["topology_choice"] = ("UNVERIFIED: production-stage run input chosen by stage/name "
                                 "(atom counts not compared)") if not verified else \
            ev.get("topology_choice", "production-stage run input")
        if len(tprs) > 1:
            fs.warnings.append(f"{len(tprs)} production run inputs "
                               f"({', '.join(p.name for p in tprs)}); took {fs.topology.name}")
    if fs.topology is None and other_tops:
        fs.topology = other_tops[0]
    ev["verified"] = verified

    # --- energy / log / structure / index from the production stage --------- #
    prod = fs.stages.get("production", {})
    edrs = prod.get("edr", [])
    fs.energy = max(edrs, key=lambda p: p.stat().st_size) if edrs else None
    logs = prod.get("log", [])
    if logs:
        fs.log = max(logs, key=lambda p: p.stat().st_size)
        if fs.deep:
            try:
                info = gromacs.parse_log(fs.log)
                ev["production_log"] = {
                    "path": str(fs.log), "sessions": info.n_sessions,
                    "last_statistics_steps": info.last_statistics_steps,
                    "simulated_ps": info.simulated_ps,
                    "first_block_nsteps": info.mdp.get("nsteps"),
                    "last_block_nsteps": info.mdp_last.get("nsteps"),
                    "finished": info.finished[-1:] or None,
                }
                try:
                    declared = int(float(info.mdp.get("nsteps", "-1")))
                except ValueError:
                    declared = -1
                ran = (info.last_statistics_steps or 1) - 1
                if declared >= 0 and ran > declared:
                    note = (f"the run input printed in the log declares nsteps = {declared:,} but "
                            f"{ran:,} steps were run: the run was extended (e.g. convert-tpr "
                            f"-extend) and continued with -cpi/-append, which does not reprint "
                            f"the parameters")
                    ev["production_log"]["extended"] = note
                    fs.warnings.append(note)
                t = fs.trajectory and traj_info.get(str(fs.trajectory))
                if t and t.readable and info.simulated_ps is not None:
                    span = (t.t_last_ps or 0) - (t.t_first_ps or 0)
                    if abs(span - info.simulated_ps) > max(1.0, t.dt_ps or 0):
                        fs.warnings.append(
                            f"trajectory span {span:g} ps differs from the production log "
                            f"({info.simulated_ps:g} ps simulated)")
            except (OSError, ValueError) as exc:
                fs.warnings.append(f"could not parse production log {fs.log.name}: {exc}")
    gros = [p for p in prod.get("gro", []) + prod.get("pdb", [])]
    setup_structs = [p for k in ("gro", "pdb") for p in fs.stages.get("setup", {}).get(k, [])]
    fs.structure = (sorted(gros)[0] if gros else (sorted(setup_structs)[0] if setup_structs
                                                  else None))
    ndx = prod.get("ndx", []) + fs.stages.get("setup", {}).get("ndx", [])
    fs.index = ndx[0] if ndx else None


def _pick_tpr(paths: list[str], fs: FileSet) -> str:
    """Among run inputs of the same run: prefer an extension (ext/extend/Nns), then newest."""
    hinted = [p for p in paths if _EXTEND_HINT.search(Path(p).name)]
    pool = hinted or paths
    return sorted(pool, key=lambda p: Path(p).stat().st_mtime, reverse=True)[0]


# --------------------------------------------------------------------------- #
# setup files (topol.top + toppar) and job-script references
# --------------------------------------------------------------------------- #
def _resolve_setup(fs: FileSet) -> None:
    tops = [p for s in STAGES for p in fs.stages.get(s, {}).get("top", [])] + \
        [p for p in fs.unclassified if p.suffix.lower() == ".top"]
    if tops:
        ref = fs.topology.parent if fs.topology else fs.directory

        def distance(p: Path) -> tuple[int, int]:
            try:
                p.relative_to(ref)
                inside = 0
            except ValueError:
                inside = 1
            return inside, len(p.parts)
        fs.gmx_top = sorted(tops, key=lambda p: (p.name.lower() != "topol.top", *distance(p)))[0]
        tp = fs.gmx_top.parent / "toppar"
        fs.toppar = tp if tp.is_dir() else None
        includes = []
        try:
            includes = re.findall(r'#include\s+"([^"]+)"',
                                  fs.gmx_top.read_text(encoding="utf-8", errors="replace"))
        except OSError:
            pass
        missing = [i for i in includes if not (fs.gmx_top.parent / i).exists()]
        fs.evidence["gmx_top_includes"] = {"n": len(includes), "missing": missing}
        if missing:
            fs.warnings.append(f"{fs.gmx_top.name} includes files that are not present: "
                               f"{', '.join(missing[:6])}{' …' if len(missing) > 6 else ''}")


_ASSIGN = re.compile(r'^\s*(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)=["\']?([^"\'\s#]*)["\']?')
_FILE_TOKEN = re.compile(r"[^\s\"'<>|;]+\.(?:tpr|mdp|gro|top|ndx|cpt|pdb)\b")


def _check_job_script_references(fs: FileSet, files: list[Path]) -> None:
    """Report inputs that job scripts reference but that are not in this folder."""
    names = {p.name for p in files}
    missing: dict[str, list[str]] = {}
    for script in fs.job_scripts:
        try:
            lines = script.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        env: dict[str, str] = {}
        for line in lines:
            m = _ASSIGN.match(line)
            if m:
                env[m.group(1)] = m.group(2)
                continue
            if not re.search(r"\b(grompp|mdrun|convert-tpr)\b", line):
                continue
            expanded = re.sub(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?",
                              lambda mm: env.get(mm.group(1), mm.group(0)), line)
            outputs = set(re.findall(r"-o\s+(\S+)", expanded))
            for tok in _FILE_TOKEN.findall(expanded):
                if tok in outputs or "$" in tok:
                    continue
                base = Path(tok).name
                if base not in names:
                    missing.setdefault(base, [])
                    if script.name not in missing[base]:
                        missing[base].append(script.name)
    if missing:
        fs.evidence["job_script_inputs_not_found"] = missing


def pick_trajectory(fs: FileSet, interactive: bool = False, ask=None) -> Path | None:
    """
    Choose the trajectory to analyse. When discovery found the choice ambiguous and
    ``interactive`` is set, ``ask(question, options) -> index`` decides.
    """
    if interactive and ask is not None and fs.ambiguities and len(fs.all_trajectories) > 1:
        options = [f"{p.name}  ({p.stat().st_size / 1e6:.1f} MB)" for p in fs.all_trajectories]
        idx = ask("Several trajectories are plausible — which one is the production run?",
                  options)
        return fs.all_trajectories[idx]
    return fs.trajectory
