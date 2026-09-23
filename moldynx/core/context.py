"""
AnalysisContext -- the runtime object every analysis receives.

Provides universe access (full and a cached solute-only trajectory for speed),
frame slicing (start/end/stride), the detected system's selection strings,
per-analysis parameters, RNG seeding, and output-path helpers.

Performance note
----------------
``core_universe()`` streams the full trajectory **once** and writes a
PBC-corrected, solute-only trajectory (protein/nucleic/ligand/cofactor, no
water/ions) that all structural analyses then reuse.  For explicit-solvent
systems this is typically 1-2 orders of magnitude smaller than the full
trajectory, which is what keeps the toolkit tractable on 100 GB+ inputs.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis")
import MDAnalysis as mda  # noqa: E402
from MDAnalysis.transformations import unwrap  # noqa: E402
from tqdm import tqdm  # noqa: E402

from moldynx.core.config import RunConfig
from moldynx.core.system import SystemInfo, ComponentType


class AnalysisContext:
    def __init__(self, config: RunConfig, system: SystemInfo, fileset,
                 provenance=None):
        self.config = config
        self.system = system
        self.fileset = fileset
        self.provenance = provenance
        self.rng = np.random.default_rng(config.seed)
        self._full: mda.Universe | None = None
        self._core: mda.Universe | None = None
        config.ensure_dirs()

    # -- universes -------------------------------------------------------- #
    @property
    def topology(self) -> str:
        return str(self.config.topology or self.fileset.topology or self.fileset.structure)

    @property
    def trajectory(self) -> str:
        return str(self.config.trajectory or self.fileset.trajectory)

    def full_universe(self) -> mda.Universe:
        """The full solvated system (loaded lazily, streamed frame-by-frame)."""
        if self._full is None:
            self._full = mda.Universe(self.topology, self.trajectory)
        return self._full

    def core_selection(self) -> str:
        """Selection string for the 'solute' of interest (no water/ions)."""
        parts = []
        for ct in (ComponentType.PROTEIN, ComponentType.DNA, ComponentType.RNA,
                   ComponentType.LIGAND, ComponentType.COFACTOR):
            if self.system.has(ct):
                parts.append(f"({self.system.components[ct.value].selection})")
        if not parts:  # e.g. pure membrane
            return "not (resname SOL WAT HOH TIP3 SPC NA CL K POT CLA)"
        return " or ".join(parts)

    def core_universe(self) -> mda.Universe:
        """
        A cached, PBC-corrected, solute-only trajectory universe.

        Built once and stored under ``data/core.{pdb,xtc}``; reused on subsequent
        analyses and subsequent runs.
        """
        if self._core is not None:
            return self._core
        core_pdb = self.config.data_dir / "core.pdb"
        core_xtc = self.config.data_dir / "core.xtc"
        if not (core_pdb.exists() and core_xtc.exists()):
            self._extract_core(core_pdb, core_xtc)
        self._core = mda.Universe(str(core_pdb), str(core_xtc))
        return self._core

    def _extract_core(self, core_pdb: Path, core_xtc: Path) -> None:
        u = self.full_universe()
        core = u.select_atoms(self.core_selection())
        if core.n_atoms == 0:
            raise RuntimeError(f"Core selection matched no atoms: {self.core_selection()!r}")
        # make the solute whole across PBC (needs bonds; .tpr provides them)
        try:
            u.trajectory.add_transformations(unwrap(core))
        except Exception:
            pass  # some topologies lack bonds; proceed unwrapped
        sl = self.frame_slice()
        u.trajectory[sl.start or 0]
        core.write(str(core_pdb))
        n = len(range(*sl.indices(len(u.trajectory))))
        with mda.Writer(str(core_xtc), core.n_atoms) as W:
            for _ in tqdm(u.trajectory[sl], total=n,
                          desc="[core-extract] solute", unit="frame"):
                W.write(core)

    # -- frames / time ---------------------------------------------------- #
    def frame_slice(self) -> slice:
        return slice(self.config.start, self.config.end, self.config.stride)

    def times_ns(self, universe: mda.Universe) -> np.ndarray:
        return np.array([ts.time for ts in universe.trajectory]) / 1000.0

    def iter_frames(self, universe: mda.Universe, desc: str = "frames",
                    stride: int = 1):
        n = len(universe.trajectory[::stride])
        return tqdm(universe.trajectory[::stride], total=n, desc=desc, unit="frame")

    # -- selections / params --------------------------------------------- #
    def selection(self, key: str, default: str = "protein") -> str:
        return self.system.selections.get(key, default)

    def params_for(self, name: str) -> dict:
        return self.config.params_for(name)

    # -- output paths ----------------------------------------------------- #
    def csv_path(self, name: str) -> Path: return self.config.results_dir / name
    def fig_path(self, name: str) -> Path: return self.config.figures_dir / name
    def table_path(self, name: str) -> Path: return self.config.tables_dir / name

    def write_csv(self, df: pd.DataFrame, name: str, index: bool = False) -> Path:
        p = self.csv_path(name)
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=index)
        return p
