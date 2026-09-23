"""
Recursive discovery of GROMACS (and related) simulation files.

The user only provides a directory; this module finds trajectories, topologies,
structures, energy files, index and topology-parameter files, resolves the best
candidate for each semantic role, and (when several trajectories exist) supports
an interactive choice.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# extension -> semantic category
_TRAJ_EXT = {".xtc", ".trr", ".dcd", ".nc", ".netcdf", ".h5", ".xyz"}
_TOPOLOGY_EXT = {".tpr", ".psf", ".prmtop", ".parm7", ".top"}   # carry bonds/masses
_STRUCTURE_EXT = {".gro", ".pdb", ".g96"}
_ENERGY_EXT = {".edr"}
_INDEX_EXT = {".ndx"}
_GMXTOP_EXT = {".top"}
_LOG_EXT = {".log"}


@dataclass
class FileSet:
    """Resolved semantic file roles for a simulation."""
    directory: Path
    trajectory: Path | None = None
    topology: Path | None = None            # best for MDAnalysis (bonds/masses)
    structure: Path | None = None           # .gro/.pdb single frame
    energy: Path | None = None              # .edr
    index: Path | None = None               # .ndx
    gmx_top: Path | None = None             # .top (#includes)
    log: Path | None = None                 # md .log
    all_trajectories: list[Path] = field(default_factory=list)
    all_files: dict[str, list[Path]] = field(default_factory=dict)

    def available_keys(self) -> set[str]:
        """Semantic keys that are present (used to gate analyses)."""
        keys = set()
        for role in ("trajectory", "topology", "structure", "energy", "index", "gmx_top"):
            if getattr(self, role) is not None:
                keys.add(role)
        # 'topology' also satisfies analyses that only need coordinates.
        if self.topology or self.structure:
            keys.add("topology")
        return keys

    def to_dict(self) -> dict:
        return {role: (str(getattr(self, role)) if getattr(self, role) else None)
                for role in ("directory", "trajectory", "topology", "structure",
                             "energy", "index", "gmx_top", "log")}


def _score_topology(p: Path) -> int:
    """Prefer .tpr (full parameters) > .psf/.prmtop > .pdb/.gro."""
    ext = p.suffix.lower()
    return {".tpr": 5, ".psf": 4, ".prmtop": 4, ".parm7": 4,
            ".pdb": 2, ".gro": 2}.get(ext, 1)


def _largest(paths: list[Path]) -> Path | None:
    return max(paths, key=lambda p: p.stat().st_size) if paths else None


def discover_files(directory: str | Path) -> FileSet:
    """Recursively scan ``directory`` and resolve the best file for each role."""
    directory = Path(directory)
    if not directory.exists():
        raise FileNotFoundError(f"Input directory does not exist: {directory}")

    by_ext: dict[str, list[Path]] = {}
    for p in sorted(directory.rglob("*")):
        if p.is_file():
            by_ext.setdefault(p.suffix.lower(), []).append(p)

    traj = sorted([p for ext in _TRAJ_EXT for p in by_ext.get(ext, [])])
    topo_candidates = [p for ext in _TOPOLOGY_EXT | _STRUCTURE_EXT
                       for p in by_ext.get(ext, [])]
    # best topology: highest score, then largest (bigger .tpr = production run)
    topology = None
    if topo_candidates:
        topology = sorted(topo_candidates,
                          key=lambda p: (_score_topology(p), p.stat().st_size),
                          reverse=True)[0]
    structure = _largest([p for ext in _STRUCTURE_EXT for p in by_ext.get(ext, [])])

    fs = FileSet(
        directory=directory,
        trajectory=_largest(traj),        # default: the largest (usually production)
        topology=topology,
        structure=structure,
        energy=_largest([p for ext in _ENERGY_EXT for p in by_ext.get(ext, [])]),
        index=_largest([p for ext in _INDEX_EXT for p in by_ext.get(ext, [])]),
        gmx_top=_largest([p for ext in _GMXTOP_EXT for p in by_ext.get(ext, [])]),
        log=_largest([p for ext in _LOG_EXT for p in by_ext.get(ext, [])]),
        all_trajectories=traj,
        all_files={k: v for k, v in by_ext.items()},
    )
    return fs


def pick_trajectory(fs: FileSet, interactive: bool = False,
                    ask=None) -> Path | None:
    """
    Choose the trajectory to analyse.

    If several trajectories are found and ``interactive`` is set, ``ask`` (a
    callable returning the selected index) is used; otherwise the largest is
    kept and a note is recorded by the caller.
    """
    if len(fs.all_trajectories) <= 1:
        return fs.trajectory
    if interactive and ask is not None:
        options = [f"{p.name}  ({p.stat().st_size/1e6:.1f} MB)"
                   for p in fs.all_trajectories]
        idx = ask("Multiple trajectories found — which one to analyse?", options)
        return fs.all_trajectories[idx]
    return fs.trajectory  # default: largest
