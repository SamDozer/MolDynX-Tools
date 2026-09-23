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
        A cached, PBC-treated, solute-only trajectory universe.

        Stored under ``data/core.{pdb,xtc}`` with ``data/core_meta.json``. The cache is
        reused only if its key -- input fingerprints, PBC mode, frame slice, solute
        selection -- matches the current run; otherwise it is rebuilt.
        """
        if self._core is not None:
            return self._core
        core_pdb = self.config.data_dir / "core.pdb"
        core_xtc = self.config.data_dir / "core.xtc"
        meta = self._read_core_meta()
        if not (core_pdb.exists() and core_xtc.exists() and meta
                and meta.get("cache_key") == self._core_cache_key()):
            self._extract_core(core_pdb, core_xtc)
        self._core = mda.Universe(str(core_pdb), str(core_xtc))
        return self._core

    # -- cache bookkeeping ------------------------------------------------ #
    _EXTRACT_VERSION = 2   # bump when the extraction algorithm changes

    def _core_cache_key(self) -> dict:
        from moldynx.core.provenance import file_fingerprint
        sl = self.frame_slice()
        return {"version": self._EXTRACT_VERSION,
                "topology": file_fingerprint(self.topology),
                "trajectory": file_fingerprint(self.trajectory),
                "selection": self.core_selection(), "pbc": self.config.pbc,
                "slice": [sl.start, sl.stop, sl.step]}

    def _read_core_meta(self) -> dict | None:
        p = self.config.data_dir / "core_meta.json"
        if not p.exists():
            return None
        try:
            import json
            return json.loads(p.read_text(encoding="utf-8"))
        except ValueError:
            return None

    @property
    def core_meta(self) -> dict:
        """Chain identity, PBC treatment and cache key of the solute trajectory."""
        self.core_universe()
        return self._read_core_meta() or {}

    def chain_groups(self, universe: mda.Universe | None = None) -> list[tuple[dict, "mda.AtomGroup"]]:
        """
        Protein chains of the solute trajectory as ``(record, AtomGroup)`` pairs,
        addressed by atom index ranges (robust to PDB segid truncation).
        """
        u = universe or self.core_universe()
        out = []
        for rec in self.core_meta.get("chains", []):
            ag = u.atoms[rec["core_start"]:rec["core_stop"]]
            if ag.n_atoms:
                out.append((rec, ag))
        return out

    # -- extraction ------------------------------------------------------- #
    def _extract_core(self, core_pdb: Path, core_xtc: Path) -> None:
        import json
        from moldynx.core.pbc import PBCProcessor

        u = self.full_universe()
        core = u.select_atoms(self.core_selection())
        if core.n_atoms == 0:
            raise RuntimeError(f"Core selection matched no atoms: {self.core_selection()!r}")
        core_ix = core.indices

        # units: protein chains (from the full topology) + other bonded molecules
        chains = [c for c in getattr(self.system, "chains", [])
                  if np.isin(np.arange(c.atom_start, c.atom_stop), core_ix).all()]
        units, labels, in_chain = [], [], np.zeros(core.n_atoms, bool)
        for c in chains:
            units.append(u.atoms[c.atom_start:c.atom_stop])
            labels.append(c.segid or f"chain{c.index}")
            lo = int(np.searchsorted(core_ix, c.atom_start))
            in_chain[lo:lo + (c.atom_stop - c.atom_start)] = True
        rest = core[~in_chain]
        if rest.n_atoms:
            try:
                frags = [f.intersection(rest) for f in rest.fragments]
            except Exception:  # no bonds
                frags = [s.atoms.intersection(rest) for s in rest.segments]
            for f in sorted((f for f in frags if f.n_atoms), key=lambda a: int(a.indices[0])):
                units.append(f)
                labels.append(f"{f.residues[0].resname}{f.residues[0].resid}")
        if not units:
            units, labels = [core], ["solute"]
        n_prot = len(chains)
        pairs = [(i, j) for i in range(min(n_prot, 6)) for j in range(i + 1, min(n_prot, 6))]
        pairs += [(0, k) for k in range(n_prot, min(len(units), n_prot + 6))] if n_prot else []

        proc = PBCProcessor(core, units, labels, mode=self.config.pbc, pairs=pairs)
        if not proc.has_bonds and proc.mode != "none":
            print("[core-extract] WARNING: the topology has no bonds -- molecules cannot be "
                  "made whole; this is recorded in data/core_meta.json and results/pbc_summary.json")
        sl = self.frame_slice()
        n = len(range(*sl.indices(len(u.trajectory))))
        wrote_pdb = False
        with mda.Writer(str(core_xtc), core.n_atoms) as W:
            for ts in tqdm(u.trajectory[sl], total=n,
                           desc=f"[core-extract] solute ({proc.mode})", unit="frame"):
                proc.process(ts)
                if not wrote_pdb:
                    core.write(str(core_pdb))
                    wrote_pdb = True
                W.write(core)

        summary = proc.summary()
        self.config.results_dir.mkdir(parents=True, exist_ok=True)
        proc.per_frame().to_csv(self.config.results_dir / "pbc_per_frame.csv", index=False)
        (self.config.results_dir / "pbc_summary.json").write_text(
            json.dumps(summary, indent=2, default=float), encoding="utf-8")
        chain_meta = []
        for c in chains:
            lo = int(np.searchsorted(core_ix, c.atom_start))
            d = c.to_dict()
            d.update(core_start=lo, core_stop=lo + (c.atom_stop - c.atom_start))
            chain_meta.append(d)
        meta = {"cache_key": self._core_cache_key(), "n_atoms": int(core.n_atoms),
                "n_frames": n, "selection": self.core_selection(),
                "pbc_mode": summary["mode"], "made_whole": summary["made_whole"],
                "whole_box_translations_undone": {u_["label"]: u_["whole_box_translations_undone"]
                                                  for u_ in summary["units"]},
                "first_frame_clustered": summary["first_frame_clustered"],
                "checks": summary["checks"], "chains": chain_meta, "units": labels}
        (self.config.data_dir / "core_meta.json").write_text(
            json.dumps(meta, indent=2, default=str), encoding="utf-8")
        if self.provenance is not None:
            self.provenance.data["pbc"] = {k: meta[k] for k in
                                           ("pbc_mode", "made_whole", "checks",
                                            "whole_box_translations_undone")}

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
