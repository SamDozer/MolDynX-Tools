"""
Periodic-boundary handling for the solute trajectory: **diagnose, treat, prove** -- in one pass.

For every frame, before anything is changed, the raw coordinates are measured:
is each molecule split across the box boundary? would its centre of mass need a
whole-box translation to stay continuous with the previous frame? how close are
the partners, measured the PBC-aware way? Then the chosen treatment is applied,
and the result is checked against the raw frame: every atom may only have moved
by an integer combination of box vectors, every bond must be short again, and the
PBC-aware inter-molecular distances must be preserved.

Treatments (``mode``):

* ``none``   -- write raw coordinates (diagnosis and proof still recorded).
* ``whole``  -- make every molecule whole (``unwrap`` over bonded fragments).
* ``nojump`` -- whole, then translate each molecule by whole box vectors so its
  centre of mass follows a continuous path (≙ ``gmx trjconv -pbc nojump``); on the
  first frame each molecule is placed in the image nearest to the largest one
  (≙ ``-pbc cluster``) so a complex starts together.
* ``auto``   -- ``nojump`` for solutes of several molecules, ``whole`` otherwise.

Nothing is repaired silently: if the topology has no bonds the solute cannot be
made whole, and this is recorded (``made_whole = False``), not swallowed.
"""

from __future__ import annotations

import numpy as np

try:  # MDAnalysis >= 2.0
    from MDAnalysis.lib.distances import capped_distance
except ImportError:  # pragma: no cover
    capped_distance = None

SPLIT_TOL_A = 0.1          # raw extent larger than whole extent by more than this -> split
BOND_MAX_A = 2.5           # a whole molecule has no bond longer than this (Å)
DIST_TOL_A = 0.02          # PBC-aware vs processed minimum distance must agree within this
BOX_TRANSLATION_TOL_A = 0.01   # at or below XTC precision (0.001 nm)
CONTACT_A = 4.5            # heavy-atom contact cutoff used for the "still in contact" record
MIN_DIST_CUTOFFS_A = (8.0, 25.0, 60.0)   # widen the search only when nothing is found


def _box_matrix(ts) -> np.ndarray | None:
    m = getattr(ts, "triclinic_dimensions", None)
    if m is None or not np.all(np.isfinite(m)) or np.linalg.det(m) <= 0:
        return None
    return np.asarray(m, dtype=float)


def _extent(pos: np.ndarray) -> float:
    return float((pos.max(0) - pos.min(0)).max()) if len(pos) else 0.0


def _min_dist(a: np.ndarray, b: np.ndarray, box=None) -> float:
    if capped_distance is None or not len(a) or not len(b):
        return float("nan")
    for cutoff in MIN_DIST_CUTOFFS_A:
        _pairs, d = capped_distance(a, b, max_cutoff=cutoff, box=box, return_distances=True)
        if len(d):
            return float(d.min())
    return float("inf")


class PBCProcessor:
    """
    Applies the treatment frame by frame and records diagnosis + proof.

    Parameters
    ----------
    core : the solute AtomGroup that is written out.
    units : list of AtomGroups (subsets of ``core``) treated as rigid-body units for
        continuity -- normally one per molecule (bonded fragment) or chain.
    labels : one label per unit (e.g. chain segids).
    pairs : (i, j) unit pairs whose minimum heavy-atom distance is tracked.
    """

    def __init__(self, core, units, labels, mode: str = "auto", pairs=None):
        self.core = core
        self.units = units
        self.labels = labels
        self.mode = mode if mode != "auto" else ("nojump" if len(units) > 1 else "whole")
        self.pairs = pairs if pairs is not None else (
            [(0, 1)] if len(units) >= 2 else [])
        self.has_bonds = self._has_bonds()
        self.made_whole = self.mode in ("whole", "nojump") and self.has_bonds
        # map each unit to positions inside `core` (core order is what is written)
        pos_in_core = {int(ix): k for k, ix in enumerate(core.indices)}
        self._unit_idx = [np.array([pos_in_core[int(i)] for i in u.indices]) for u in units]
        self._heavy = [np.array([pos_in_core[int(i)] for i in u.select_atoms("not name H*").indices])
                       for u in units]
        self._bonds = self._core_bonds(pos_in_core)
        self._prev_com: list[np.ndarray | None] = [None] * len(units)
        self.rows: list[dict] = []
        # jump *events*: frames where the whole-box translation needed to keep a unit
        # continuous changes; once a unit has crossed, the same translation is applied on
        # every later frame, counted separately as frames_translated
        self.jumps = np.zeros(len(units), dtype=int)
        self.frames_translated = np.zeros(len(units), dtype=int)
        self._applied = [np.zeros(3, dtype=int) for _ in units]
        self.first_frame_cluster = False

    # ------------------------------------------------------------------ #
    def _has_bonds(self) -> bool:
        try:
            return len(self.core.bonds) > 0
        except Exception:  # topology without bond information
            return False

    def _core_bonds(self, pos_in_core) -> np.ndarray:
        if not self.has_bonds:
            return np.empty((0, 2), dtype=int)
        b = self.core.bonds.indices
        keep = [(pos_in_core[int(i)], pos_in_core[int(j)]) for i, j in b
                if int(i) in pos_in_core and int(j) in pos_in_core]
        return np.array(keep, dtype=int) if keep else np.empty((0, 2), dtype=int)

    @staticmethod
    def _com(pos, masses):
        return (pos * masses[:, None]).sum(0) / masses.sum()

    # ------------------------------------------------------------------ #
    def process(self, ts) -> np.ndarray:
        """Treat the current frame of ``core`` in place; return the processed positions."""
        box = _box_matrix(ts)
        dims = ts.dimensions
        raw = self.core.positions.copy()
        masses = self.core.masses
        row: dict = {"frame": int(ts.frame), "time_ps": float(ts.time),
                     "box_a_nm": float(dims[0]) / 10 if dims is not None else np.nan}

        # ---- diagnosis on the raw frame ---------------------------------- #
        raw_ext = [_extent(raw[ix]) for ix in self._unit_idx]
        raw_min = [_min_dist(raw[self._heavy[i]], raw[self._heavy[j]], dims)
                   for i, j in self.pairs]

        # ---- treatment ------------------------------------------------------ #
        if self.made_whole:
            self.core.unwrap(compound="fragments", reference=None, inplace=True)
        pos = self.core.positions
        whole_ext = [_extent(pos[ix]) for ix in self._unit_idx]
        if self.mode == "nojump" and box is not None:
            inv = np.linalg.inv(box)
            coms = [self._com(pos[ix], masses[ix]) for ix in self._unit_idx]
            if self._prev_com[0] is None:            # first frame: cluster around the largest
                big = int(np.argmax([len(ix) for ix in self._unit_idx]))
                for k, ix in enumerate(self._unit_idx):
                    if k == big:
                        continue
                    n = np.round((coms[k] - coms[big]) @ inv)
                    if np.any(n != 0):
                        pos[ix] -= n @ box
                        coms[k] = coms[k] - n @ box
                        self.first_frame_cluster = True
            else:
                for k, ix in enumerate(self._unit_idx):
                    n = np.round((coms[k] - self._prev_com[k]) @ inv).astype(int)
                    if np.any(n != 0):
                        pos[ix] -= n @ box
                        coms[k] = coms[k] - n @ box
                        self.frames_translated[k] += 1
                    if np.any(n != self._applied[k]):
                        self.jumps[k] += 1
                    self._applied[k] = n
            self._prev_com = coms
            self.core.positions = pos

        # ---- proof: processed vs raw ---------------------------------------- #
        pos = self.core.positions
        disp = pos - raw
        if box is not None:
            frac = disp @ np.linalg.inv(box)
            dev = np.abs(disp - np.round(frac) @ box).max() if len(disp) else 0.0
            row["max_dev_from_box_translation_A"] = float(dev)
            for k, ix in enumerate(self._unit_idx):
                n_vec = len({tuple(v) for v in np.round(frac[ix]).astype(int)})
                row[f"translations_{k}"] = n_vec
        for k, ix in enumerate(self._unit_idx):
            row[f"split_raw_{k}"] = bool(raw_ext[k] > whole_ext[k] + SPLIT_TOL_A)
            row[f"extent_whole_nm_{k}"] = whole_ext[k] / 10
        if len(self._bonds):
            bl = np.linalg.norm(pos[self._bonds[:, 0]] - pos[self._bonds[:, 1]], axis=1)
            row["max_bond_A"] = float(bl.max())
        for (i, j), rm in zip(self.pairs, raw_min):
            pm = _min_dist(pos[self._heavy[i]], pos[self._heavy[j]], None)
            row[f"min_dist_raw_pbc_nm_{i}_{j}"] = rm / 10
            row[f"min_dist_processed_nm_{i}_{j}"] = pm / 10
        self.rows.append(row)
        return pos

    # ------------------------------------------------------------------ #
    def summary(self) -> dict:
        import pandas as pd
        df = pd.DataFrame(self.rows)
        n = len(df)
        units = []
        for k, lab in enumerate(self.labels):
            split = df.get(f"split_raw_{k}", pd.Series(dtype=bool)).astype(bool)
            multi = df.get(f"translations_{k}", pd.Series(1, index=df.index)) > 1
            # making a split molecule whole moves only its wrapped atoms, so it shows two
            # translation vectors in exactly the frames where it was split in the raw data
            consistent = (multi == split) if self.made_whole else (~multi | split)
            units.append({
                "label": lab, "n_atoms": int(len(self._unit_idx[k])),
                "frames_split_in_raw": int(split.sum()),
                "first_split_frames": [int(f) for f in df.frame[split].head(10)],
                "frames_with_multiple_translations": int(multi.sum()),
                "multiple_translations_only_where_split": bool(consistent.all()) if n else True,
                "whole_box_translations_undone": int(self.jumps[k]),
                "frames_translated": int(self.frames_translated[k]),
                "extent_whole_nm_mean": float(df.get(f"extent_whole_nm_{k}",
                                                     pd.Series([np.nan])).mean()),
                "extent_whole_nm_max": float(df.get(f"extent_whole_nm_{k}",
                                                    pd.Series([np.nan])).max()),
            })
        pairs = []
        for i, j in self.pairs:
            raw = df.get(f"min_dist_raw_pbc_nm_{i}_{j}")
            proc = df.get(f"min_dist_processed_nm_{i}_{j}")
            if raw is None:
                continue
            diff = (proc - raw).abs()
            pairs.append({
                "pair": [self.labels[i], self.labels[j]],
                "min_dist_raw_pbc_nm": [float(raw.min()), float(raw.mean()), float(raw.max())],
                "frames_without_heavy_atom_contact": int((raw > CONTACT_A / 10).sum()),
                "max_abs_diff_processed_vs_raw_nm": float(diff.max()),
                "frames_separated_in_processed": int((diff > DIST_TOL_A / 10).sum()),
            })
        dev = df.get("max_dev_from_box_translation_A")
        bond = df.get("max_bond_A")
        checks = {
            "only_whole_box_translations": bool(dev is None
                                                or dev.max() <= BOX_TRANSLATION_TOL_A) if n else True,
            "molecules_whole": bool(bond is None or bond.max() < BOND_MAX_A) if n else True,
            "multiple_translations_only_in_split_frames": all(
                u["multiple_translations_only_where_split"] for u in units),
            "interchain_distance_preserved": all(p["frames_separated_in_processed"] == 0
                                                 for p in pairs),
        }
        if not self.made_whole:
            checks["molecules_whole"] = False if self.mode != "none" else checks["molecules_whole"]
        box = df.get("box_a_nm")
        return {
            "mode": self.mode, "made_whole": self.made_whole, "has_bonds": self.has_bonds,
            "first_frame_clustered": self.first_frame_cluster, "n_frames": n,
            "box_a_nm": [float(box.iloc[0]), float(box.iloc[-1]), float(box.min())]
            if box is not None and n else None,
            "half_box_min_nm": float(box.min()) / 2 if box is not None and n else None,
            "max_dev_from_box_translation_A": float(dev.max()) if dev is not None and n else None,
            "max_bond_A": float(bond.max()) if bond is not None and n else None,
            "units": units, "pairs": pairs, "checks": checks,
            "tolerances": {"split_A": SPLIT_TOL_A, "bond_A": BOND_MAX_A,
                           "distance_A": DIST_TOL_A, "box_translation_A": BOX_TRANSLATION_TOL_A},
        }

    def per_frame(self):
        import pandas as pd
        return pd.DataFrame(self.rows)
