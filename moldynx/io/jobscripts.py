"""
Chain of custody from job scripts: which file fed which GROMACS step.

Reads shell / SLURM scripts, expands simple ``VAR=value`` assignments and
``${VAR}`` references, and extracts every ``gmx grompp``, ``gmx mdrun`` and
``gmx convert-tpr`` call with its file arguments. Reading only -- nothing runs.
"""

from __future__ import annotations

import re
import shlex
from dataclasses import asdict, dataclass, field
from pathlib import Path

_ASSIGN = re.compile(r'^\s*(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)=["\']?([^"\'\s#]*)["\']?\s*(#.*)?$')
_VAR = re.compile(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?")
_TOOL = re.compile(r"\bgmx(?:_mpi|_d)?\s+(grompp|mdrun|convert-tpr)\b(.*)$")

# the flags that matter for provenance
FLAGS = {
    "grompp": {"-f": "mdp", "-c": "coordinates", "-r": "restraint_reference", "-t": "checkpoint",
               "-p": "topology", "-n": "index", "-o": "output_tpr", "-maxwarn": "maxwarn"},
    "mdrun": {"-s": "tpr", "-deffnm": "deffnm", "-cpi": "checkpoint_in", "-append": "append",
              "-nb": "nb", "-pme": "pme", "-ntmpi": "ntmpi", "-ntomp": "ntomp", "-npme": "npme"},
    "convert-tpr": {"-s": "tpr_in", "-o": "tpr_out", "-extend": "extend_ps", "-until": "until_ps",
                    "-nsteps": "nsteps"},
}
_BOOL_FLAGS = {"-append", "-v"}


@dataclass
class GmxCall:
    script: str
    line_no: int
    tool: str
    args: dict[str, str] = field(default_factory=dict)
    raw: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def parse_job_script(path: str | Path) -> list[GmxCall]:
    path = Path(path)
    env: dict[str, str] = {}
    calls: list[GmxCall] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return calls
    for no, line in enumerate(lines, 1):
        stripped = line.split("#", 1)[0] if not line.lstrip().startswith("#!") else ""
        m = _ASSIGN.match(line)
        if m and "gmx" not in line:
            env[m.group(1)] = _VAR.sub(lambda mm: env.get(mm.group(1), mm.group(0)), m.group(2))
            continue
        expanded = _VAR.sub(lambda mm: env.get(mm.group(1), mm.group(0)), stripped)
        t = _TOOL.search(expanded)
        if not t:
            continue
        tool, rest = t.group(1), t.group(2)
        try:
            toks = shlex.split(rest, posix=True)
        except ValueError:
            toks = rest.split()
        wanted = FLAGS[tool]
        args: dict[str, str] = {}
        i = 0
        while i < len(toks):
            tok = toks[i]
            if tok in ("|", ";", "&&", "||", ">", "2>", "<"):
                break
            if tok in wanted:
                if tok in _BOOL_FLAGS:
                    args[wanted[tok]] = "yes"
                elif i + 1 < len(toks) and not toks[i + 1].startswith("-"):
                    args[wanted[tok]] = toks[i + 1]
                    i += 1
            i += 1
        calls.append(GmxCall(script=path.name, line_no=no, tool=tool, args=args,
                             raw=" ".join(expanded.split())))
    return calls


def chain_of_custody(scripts: list[str | Path]) -> list[dict]:
    """All grompp / mdrun / convert-tpr calls found in the scripts, in order."""
    out = []
    for s in scripts:
        out.extend(c.to_dict() for c in parse_job_script(s))
    return out
