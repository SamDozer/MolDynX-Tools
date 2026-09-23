"""Multi-frame SASA must equal single-frame SASA (mdtraj 1.11.1 per-thread state bug)."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import FIX


@pytest.fixture
def moving_traj():
    import mdtraj as md
    t = md.load(str(FIX / "structures" / "A8HNE1_ZmBiP2_CA_only.pdb"))
    xyz = np.repeat(t.xyz, 8, axis=0)
    for k in range(8):                           # partner B moves as a rigid body
        xyz[k, 187:, 0] += 0.05 * k * (k + 1) / 2
    return md.Trajectory(xyz, t.topology)


def test_wrapper_matches_single_frame_calls(moving_traj):
    import mdtraj as md
    from moldynx.core.surface import shrake_rupley
    ref = np.array([md.shrake_rupley(moving_traj[k], mode="atom")[0]
                    for k in range(moving_traj.n_frames)])
    got = shrake_rupley(moving_traj, mode="atom")
    np.testing.assert_allclose(got, ref, atol=1e-5)


def test_rigid_partner_area_is_constant_and_buried_area_positive(moving_traj):
    from moldynx.core.surface import shrake_rupley
    b = shrake_rupley(moving_traj.atom_slice(np.arange(187, 850)), mode="atom").sum(1)
    assert np.ptp(b) < 1e-3                       # a rigid body keeps its own area
    a = shrake_rupley(moving_traj.atom_slice(np.arange(187)), mode="atom").sum(1)
    ab = shrake_rupley(moving_traj, mode="atom").sum(1)
    bsa = a + b - ab
    assert (bsa > -1e-3).all() and np.all(np.diff(bsa) <= 1e-3)   # shrinks as B leaves
