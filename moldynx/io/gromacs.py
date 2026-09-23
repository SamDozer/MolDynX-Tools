"""
Read-only readers for GROMACS files, used as *evidence* by intake and audits.

Everything here only reads. In particular the XTC scanner walks frame headers
without building MDAnalysis offset caches, which would write hidden ``.npz``
files into the (raw, provenance-critical) simulation folder.

* :func:`parse_log`      -- mdrun ``.log``: sessions, the ``Input Parameters``
                            dump, coupling groups, minimisation outcome,
                            warnings, timings. Tolerates trimmed logs.
* :func:`parse_mdp`      -- ``.mdp`` key/value pairs.
* :func:`scan_xtc`       -- atoms, frames, times, steps, spacing of an ``.xtc``.
* :func:`tpr_header`     -- version, atom count, velocities flag of a ``.tpr``.
* :func:`find_gmx` / :func:`run_gmx` -- a GROMACS binary, native or through WSL,
                            for evidence that needs ``gmx dump`` / ``gmx check``.
"""

from __future__ import annotations

import os
import re
import shutil
import struct
import subprocess
import tempfile
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

# --------------------------------------------------------------------------- #
# mdrun log
# --------------------------------------------------------------------------- #
_MINIMIZERS = r"(Steepest Descents|Polak-Ribiere Conjugate Gradients|Low-Memory BFGS Minimizer)"
_TIMESTAMP = "%a %b %d %H:%M:%S %Y"


@dataclass
class MinimizationResult:
    algorithm: str | None = None
    outcome: str | None = None           # e.g. "converged to Fmax < 1000"
    reached_emtol: bool | None = None    # False when stopped on machine precision / max steps
    steps: int | None = None
    potential: float | None = None       # kJ/mol
    fmax: float | None = None            # kJ/mol/nm
    fmax_atom: int | None = None         # 1-based atom number as printed by GROMACS
    fnorm: float | None = None


@dataclass
class LogInfo:
    path: str
    gromacs_version: str | None = None
    executable: str | None = None
    working_dir: str | None = None
    commands: list[str] = field(default_factory=list)       # one per mdrun session
    started: list[str] = field(default_factory=list)        # ISO timestamps (cluster clock)
    finished: list[str] = field(default_factory=list)
    n_sessions: int = 0
    mdp: dict[str, str] = field(default_factory=dict)       # first Input Parameters block
    mdp_last: dict[str, str] = field(default_factory=dict)  # last block (e.g. after -extend)
    grpopts: dict[str, list[str]] = field(default_factory=dict)
    last_statistics_steps: int | None = None
    performance_ns_day: list[float] = field(default_factory=list)
    wall_s: list[float] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)
    minimization: MinimizationResult | None = None
    trimmed: bool = False                                    # fixture / excerpt marker seen

    # -- convenience ------------------------------------------------------ #
    @property
    def is_minimization(self) -> bool:
        return self.mdp.get("integrator", "") in ("steep", "cg", "l-bfgs")

    def value(self, key: str, default=None):
        """First-block MDP value as float if numeric, else string."""
        v = self.mdp.get(key)
        if v is None:
            return default
        try:
            return float(v)
        except ValueError:
            return v

    @property
    def ref_t(self) -> list[float]:
        return [float(x) for x in self.grpopts.get("ref-t", [])]

    @property
    def simulated_ps(self) -> float | None:
        dt, n = self.value("dt"), self.last_statistics_steps
        if isinstance(dt, float) and n:
            return (n - 1) * dt
        return None

    def to_dict(self) -> dict:
        return asdict(self)


def _parse_param_block(lines: list[str], start: int) -> tuple[dict, dict, int]:
    """Parse an ``Input Parameters:`` dump beginning at ``start``; returns (mdp, grpopts, end)."""
    mdp: dict[str, str] = {}
    grp: dict[str, list[str]] = {}
    in_grp = False
    i = start + 1
    while i < len(lines):
        line = lines[i]
        if line.startswith("grpopts:"):
            in_grp = True
        elif in_grp:
            if not line.strip():
                return mdp, grp, i
            m = re.match(r"\s+([A-Za-z-]+):?\s+(.+?)\s*$", line)
            if m and m.group(1) in ("nrdf", "ref-t", "tau-t", "annealing", "annealing-npoints"):
                grp[m.group(1)] = m.group(2).split()
        else:
            m = re.match(r"\s+([A-Za-z0-9_-]+)\s+=\s+(.+?)\s*$", line)
            if m and m.group(1) not in mdp:
                mdp[m.group(1)] = m.group(2)
            if "ref-p (3x3)" in line or "compressibility (3x3)" in line:
                key = "ref-p_xx" if "ref-p" in line else "compressibility_xx"
                if i + 1 < len(lines):
                    mm = re.search(r"\{\s*([-\d.e+]+)", lines[i + 1])
                    if mm and key not in mdp:
                        mdp[key] = mm.group(1)
        i += 1
    return mdp, grp, i


def _iso(stamp: str) -> str | None:
    try:
        return datetime.strptime(" ".join(stamp.split()), _TIMESTAMP).isoformat()
    except ValueError:
        return None


def parse_log(path: str | Path) -> LogInfo:
    """Parse a GROMACS ``mdrun`` log (full or trimmed)."""
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    info = LogInfo(path=str(path))
    info.trimmed = "lines trimmed by" in text

    m = re.search(r"GROMACS version:\s+(\S+)", text)
    info.gromacs_version = m.group(1) if m else None
    m = re.search(r"Executable:\s+(\S+)", text)
    info.executable = m.group(1) if m else None
    m = re.search(r"Working dir:\s+(\S+)", text)
    info.working_dir = m.group(1) if m else None
    info.commands = [c.strip() for c in re.findall(r"Command line:\n\s+(.+)", text)]
    info.started = [s for s in (_iso(x) for x in re.findall(r"Started mdrun on rank 0 (.+)", text)) if s]
    info.finished = [s for s in (_iso(x) for x in re.findall(r"Finished mdrun on rank 0 (.+)", text)) if s]
    info.n_sessions = len(re.findall(r"Started mdrun on rank 0", text)) or len(info.commands)

    blocks = [i for i, line in enumerate(lines) if line.startswith("Input Parameters:")]
    if blocks:
        info.mdp, info.grpopts, _ = _parse_param_block(lines, blocks[0])
        info.mdp_last, _, _ = _parse_param_block(lines, blocks[-1])

    steps = re.findall(r"Statistics over (\d+) steps", text)
    info.last_statistics_steps = int(steps[-1]) if steps else None
    info.performance_ns_day = [float(x) for x in re.findall(r"Performance:\s+([\d.]+)", text)]
    info.wall_s = [float(x) for x in re.findall(r"Time:\s+[\d.]+\s+([\d.]+)\s+[\d.]+", text)]
    info.counts = {
        "lincs_warnings": len(re.findall(r"LINCS WARNING", text)),
        "wrote_pdb": len(re.findall(r"Wrote pdb", text)),
        "constraint_errors": len(re.findall(r"onstraint error", text)),
        "warnings": len(re.findall(r"(?m)^WARNING", text)),
        "fatal_errors": len(re.findall(r"Fatal error", text)),
    }

    if info.is_minimization or re.search(_MINIMIZERS + r" (converged|did not converge)", text) \
            or "Energy minimization" in text:
        em = MinimizationResult()
        mm = re.search(_MINIMIZERS + r" (converged to .+?) in (\d+) steps", text)
        if mm:
            em.algorithm, em.outcome, em.steps = mm.group(1), mm.group(2), int(mm.group(3))
        elif "reached the maximum number of steps" in text:
            em.outcome = "reached the maximum number of steps"
        em.reached_emtol = bool(em.outcome) and "machine precision" not in em.outcome \
            and "did not reach the requested Fmax" not in text \
            and "maximum number of steps" not in (em.outcome or "")
        for key, pat in (("potential", r"Potential Energy\s+=\s+(\S+)"),
                         ("fmax", r"Maximum force\s+=\s+(\S+)"),
                         ("fnorm", r"Norm of force\s+=\s+(\S+)")):
            found = re.findall(pat, text)
            if found:
                setattr(em, key, float(found[-1]))
        found = re.findall(r"Maximum force\s+=\s+\S+ on atom (\d+)", text)
        em.fmax_atom = int(found[-1]) if found else None
        info.minimization = em
    return info


# --------------------------------------------------------------------------- #
# mdp
# --------------------------------------------------------------------------- #
def parse_mdp(path: str | Path) -> dict[str, str]:
    """Key/value pairs of an ``.mdp`` (comments stripped, ``_`` normalised to ``-``)."""
    out: dict[str, str] = {}
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if "=" in line:
            k, v = (s.strip() for s in line.split("=", 1))
            out[k.replace("_", "-").lower()] = v
    return out


# --------------------------------------------------------------------------- #
# xtc
# --------------------------------------------------------------------------- #
@dataclass
class XtcInfo:
    path: str
    readable: bool
    natoms: int | None = None
    n_frames: int = 0
    t_first_ps: float | None = None
    t_last_ps: float | None = None
    dt_ps: float | None = None            # modal spacing
    uniform_spacing: bool | None = None
    step_first: int | None = None
    step_last: int | None = None
    truncated: bool = False               # last frame incomplete
    error: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


_XTC_MAGIC = 1995


def scan_xtc(path: str | Path, max_frames: int | None = None) -> XtcInfo:
    """
    Walk the frame headers of an XTC file (read-only; no offset cache written).

    Frame layout (XDR, big-endian): magic, natoms, step, time, box[9], natoms,
    then for natoms > 9: precision, minint[3], maxint[3], smallidx, nbytes and
    nbytes of compressed data padded to 4 bytes.
    """
    path = Path(path)
    info = XtcInfo(path=str(path), readable=False)
    try:
        size = path.stat().st_size
        with open(path, "rb") as fh:
            times, steps, pos = [], [], 0
            while pos < size:
                head = fh.read(92)
                if len(head) < 56:
                    info.truncated = len(head) > 0
                    break
                magic, natoms, step = struct.unpack(">iii", head[:12])
                (time,) = struct.unpack(">f", head[12:16])
                if magic != _XTC_MAGIC:
                    info.error = f"bad magic {magic} at byte {pos}"
                    break
                if natoms <= 9:
                    nxt = pos + 56 + natoms * 12
                else:
                    if len(head) < 92:
                        info.truncated = True
                        break
                    (nbytes,) = struct.unpack(">i", head[88:92])
                    nxt = pos + 92 + ((nbytes + 3) // 4) * 4
                if nxt > size:
                    info.truncated = True
                    break
                info.natoms = natoms
                times.append(float(time))
                steps.append(int(step))
                pos = nxt
                fh.seek(pos)
                if max_frames and len(times) >= max_frames:
                    break
        info.n_frames = len(times)
        if times:
            info.readable = True
            info.t_first_ps, info.t_last_ps = times[0], times[-1]
            info.step_first, info.step_last = steps[0], steps[-1]
            if len(times) > 1:
                diffs = [round(b - a, 3) for a, b in zip(times, times[1:])]
                info.dt_ps = max(set(diffs), key=diffs.count)
                info.uniform_spacing = all(abs(d - info.dt_ps) < 1e-3 for d in diffs)
    except OSError as exc:
        info.error = str(exc)
    return info


# --------------------------------------------------------------------------- #
# tpr
# --------------------------------------------------------------------------- #
@dataclass
class TprHeader:
    path: str
    readable: bool
    version: str | None = None
    natoms: int | None = None
    n_tcoupl_groups: int | None = None
    has_velocities: bool | None = None
    has_box: bool | None = None
    error: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def tpr_header(path: str | Path) -> TprHeader:
    """Header of a ``.tpr`` (reads only the first few KB)."""
    path = Path(path)
    out = TprHeader(path=str(path), readable=False)
    try:
        from MDAnalysis.topology.tpr import utils as tpr_utils
        with open(path, "rb") as fh:
            data = tpr_utils.TPXUnpacker(fh.read(8192))
        th = tpr_utils.read_tpxheader(data)
        out.readable = True
        out.version = th.ver_str.decode(errors="replace").replace("VERSION", "").strip() \
            if isinstance(th.ver_str, bytes) else str(th.ver_str)
        out.natoms, out.n_tcoupl_groups = int(th.natoms), int(th.ngtc)
        out.has_velocities, out.has_box = bool(th.bV), bool(th.bBox)
    except Exception as exc:  # unreadable / not a TPR / fake test file
        out.error = f"{type(exc).__name__}: {exc}"
    return out


# --------------------------------------------------------------------------- #
# GROMACS binary (native or WSL)
# --------------------------------------------------------------------------- #
@dataclass
class Gmx:
    kind: str            # "native" | "wsl"
    command: list[str]   # prefix, e.g. ["gmx"] or ["wsl.exe", "-d", "Ubuntu", "--"]
    version: str | None = None


def windows_to_wsl(path: str | Path) -> str:
    """``E:\\a\\b`` -> ``/mnt/e/a/b`` (paths already POSIX are returned unchanged)."""
    s = str(path)
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", s)
    if not m:
        return s.replace("\\", "/")
    return f"/mnt/{m.group(1).lower()}/{m.group(2).replace(chr(92), '/')}"


def find_gmx(prefer_wsl_distro: str | None = None) -> Gmx | None:
    """Locate GROMACS: ``gmx``/``gmx_mpi`` on PATH, else (Windows) inside WSL."""
    for name in ("gmx", "gmx_mpi"):
        exe = shutil.which(name)
        if exe:
            return Gmx("native", [exe])
    if os.name == "nt" and shutil.which("wsl.exe"):
        distro = ["-d", prefer_wsl_distro] if prefer_wsl_distro else []
        try:
            r = subprocess.run(["wsl.exe", *distro, "--", "bash", "-lc", "command -v gmx"],
                               capture_output=True, text=True, timeout=60)
            if r.returncode == 0 and r.stdout.strip():
                return Gmx("wsl", ["wsl.exe", *distro, "--"])
        except (OSError, subprocess.SubprocessError):
            return None
    return None


def _shell_quote(a: str) -> str:
    return "'" + a.replace("'", "'\"'\"'") + "'"


def run_gmx_pipeline(gmx: Gmx, args: list[str], pipe: str, timeout: int = 3600
                     ) -> subprocess.CompletedProcess:
    """``gmx <args> 2>/dev/null | <pipe>`` in bash (native or WSL) -- filters huge dumps at the source."""
    if gmx.kind == "wsl":
        wargs = [windows_to_wsl(a) if re.match(r"^[A-Za-z]:[\\/]", a) else a for a in args]
        exe = "gmx"
    else:
        wargs, exe = args, gmx.command[0]
    script = f"#!/bin/bash\n{_shell_quote(exe)} {' '.join(_shell_quote(a) for a in wargs)} " \
             f"2>/dev/null | {pipe}\n"
    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False, newline="\n") as fh:
        fh.write(script)
        spath = fh.name
    try:
        if gmx.kind == "wsl":
            cmd = gmx.command + ["bash", windows_to_wsl(spath)]
        else:
            bash = shutil.which("bash")
            if bash is None:
                raise RuntimeError("bash is required to filter gmx output")
            cmd = [bash, spath]
        return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    finally:
        try:
            os.unlink(spath)
        except OSError:
            pass


@dataclass
class PositionRestraints:
    tpr: str
    available: bool
    n_restrained: int = 0
    by_force_constant: dict[str, int] = field(default_factory=dict)   # "400" -> atoms
    by_molecule_block: list[int] = field(default_factory=list)       # restrained atoms per block
    error: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def tpr_position_restraints(gmx: Gmx | None, tpr: str | Path) -> PositionRestraints:
    """
    Position restraints stored in a run input (``gmx dump``): restrained atoms per force
    constant (x component, kJ mol⁻¹ nm⁻²) and per molecule block. ``gmx dump`` prints one
    ``Position Rest.`` interaction list per molecule type; its ``nr:`` counts integers
    (2 per restraint: type + atom).
    """
    out = PositionRestraints(tpr=str(tpr), available=False)
    if gmx is None:
        out.error = "GROMACS not available"
        return out
    pattern = r"functype\[[0-9]+\]=POSRES|Position Rest\.:|\(POSRES\)"
    try:
        r = run_gmx_pipeline(gmx, ["dump", "-s", str(tpr)],
                             f"grep -E {_shell_quote(pattern)}", timeout=3600)
    except (OSError, subprocess.SubprocessError, RuntimeError) as exc:
        out.error = str(exc)
        return out
    fc: dict[int, str] = {}
    blocks: list[int] = []
    counts: dict[str, int] = {}
    for line in r.stdout.splitlines():
        m = re.search(r"functype\[(\d+)\]=POSRES.*?fcA=\(\s*([-\d.eE+]+)", line)
        if m:
            fc[int(m.group(1))] = f"{float(m.group(2)):g}"
            continue
        if "Position Rest.:" in line:
            blocks.append(0)
            continue
        m = re.search(r"type=(\d+) \(POSRES\)", line)
        if m and blocks:
            blocks[-1] += 1
            key = fc.get(int(m.group(1)), "?")
            counts[key] = counts.get(key, 0) + 1
    out.available = r.returncode == 0 or bool(blocks)
    out.by_molecule_block = [b for b in blocks if b] if any(blocks) else blocks[:0]
    out.by_force_constant = dict(sorted(counts.items(), key=lambda kv: -float(kv[0])
                                        if kv[0] != "?" else 0))
    out.n_restrained = sum(counts.values())
    if not out.available:
        out.error = (r.stderr or "gmx dump failed")[-300:]
    return out


def run_gmx(gmx: Gmx, args: list[str], stdin: str | None = None,
            timeout: int = 3600) -> subprocess.CompletedProcess:
    """
    Run ``gmx <args>``. For WSL, the command is written to a script file (inline
    ``$variables`` and quoting do not survive ``wsl.exe`` reliably) and Windows
    paths in ``args`` are translated to ``/mnt/<drive>/...``.
    """
    if gmx.kind == "native":
        return subprocess.run(gmx.command + args, input=stdin, capture_output=True,
                              text=True, timeout=timeout)
    wargs = [windows_to_wsl(a) if re.match(r"^[A-Za-z]:[\\/]", a) else a for a in args]
    quoted = " ".join("'" + a.replace("'", "'\"'\"'") + "'" for a in wargs)
    script = f"#!/bin/bash\nset -o pipefail\ngmx {quoted}\n"
    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False, newline="\n") as fh:
        fh.write(script)
        spath = fh.name
    try:
        return subprocess.run(gmx.command + ["bash", windows_to_wsl(spath)], input=stdin,
                              capture_output=True, text=True, timeout=timeout)
    finally:
        try:
            os.unlink(spath)
        except OSError:
            pass
