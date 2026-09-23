"""
Shrake–Rupley SASA that is correct on multi-frame trajectories.

``mdtraj.shrake_rupley`` (verified with mdtraj 1.11.1) can return wrong values
for some frames of a *multi-frame* call: on a real 100 ns trajectory every other
frame was off by ~0.5 nm², and on a rigid-body test the per-partner area jumped
by 54 nm² -- enough to make buried surface area negative. Single-frame calls were
correct in every test, so this wrapper computes frame by frame.

Cost: single-frame calls run on one core (~4 s per frame for ~13 k atoms);
callers stride long trajectories. A process-pool version is on the roadmap.
"""

from __future__ import annotations

import numpy as np


def shrake_rupley(traj, mode: str = "atom", **kwargs) -> np.ndarray:
    """Drop-in replacement for ``mdtraj.shrake_rupley`` (nm² per atom/residue per frame)."""
    import mdtraj as md
    return np.concatenate([md.shrake_rupley(traj[i], mode=mode, **kwargs)
                           for i in range(traj.n_frames)])
