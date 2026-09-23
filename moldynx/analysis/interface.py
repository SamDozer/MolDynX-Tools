"""
Interface analysis for complexes (protein–protein / protein–nucleic).

Streaming one frame at a time over heavy atoms, with PBC-aware distances:

* minimum inter-partner heavy-atom distance (is the complex in contact? --
  the centre-of-mass distance is *not* a proximity measure for extended chains);
* residue–residue contacts (heavy atoms < ``contact_cutoff``) per frame and the
  full contact-frequency map; persistent contacts (occupancy ≥ ``persistent``);
* interface residues (heavy atoms < ``interface_cutoff``) per frame and the
  per-residue interface occupancy; core interface residues (occupancy ≥
  ``core_occupancy``) -- reported in full, never truncated;
* residues that come within ``ever_cutoff`` of the partner in *any* frame -- the
  residue set a binding-energy decomposition must cover (late binders included);
* buried surface area (Shrake–Rupley, mdtraj) and interface Cα RMSD;
* the trend of contacts and buried area (slope per 10 ns, p) -- an interface that
  is still growing is reported as such.

Partners: protein vs nucleic acid when present; otherwise the two largest protein
chains, addressed through the chain identity persisted at extraction.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from moldynx import plotting
from moldynx import statistics as st
from moldynx.core.base import BaseAnalysis
from moldynx.core.surface import shrake_rupley
from moldynx.core.system import COMPLEX_SYSTEMS
from moldynx.plotting import PALETTE, SEQ_CMAP

HEAVY = "not name H*"


def _trend(t_ns: np.ndarray, x: np.ndarray) -> dict:
    from scipy import stats
    ok = np.isfinite(x)
    if ok.sum() < 5 or np.ptp(t_ns[ok]) == 0:
        return {"slope_per_10ns": None, "p": None}
    r = stats.linregress(t_ns[ok], x[ok])
    return {"slope_per_10ns": float(r.slope * 10), "p": float(r.pvalue)}


class InterfaceAnalysis(BaseAnalysis):
    name = "interface"
    label = "Interface (contacts, residues, buried area, iRMSD)"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = COMPLEX_SYSTEMS
    order = 40
    default_params = {"contact_cutoff": 4.5, "interface_cutoff": 5.0, "ever_cutoff": 6.0,
                      "core_occupancy": 0.5, "persistent": 0.5, "bsa_stride": 10,
                      "window": 20}
    outputs = ["results/interface_timeseries.csv", "results/interface_residues.csv",
               "results/interface_contact_frequency.csv",
               "results/interface_persistent_contacts.csv",
               "figures/interface.png", "figures/interface_residues.png",
               "figures/interface_contact_map.png"]

    # ------------------------------------------------------------------ #
    def _partners(self, ctx, u):
        """The two interface partners as (label, AtomGroup) pairs, in topology order."""
        if ctx.system.flags.get("has_nucleic"):
            return ("protein", u.select_atoms("protein")), ("nucleic", u.select_atoms("nucleic"))
        chains = ctx.chain_groups(u)
        if len(chains) < 2:
            return None, None
        two = sorted(sorted(chains, key=lambda rc: rc[1].n_atoms, reverse=True)[:2],
                     key=lambda rc: int(rc[1].indices[0]))
        return tuple((rec.get("segid") or f"chain {rec.get('index')}", ag) for rec, ag in two)

    def run(self, ctx) -> dict:
        p = self.params(ctx)
        plotting.set_style()
        u = ctx.core_universe()
        pa, pb = self._partners(ctx, u)
        if pa is None or pa[1].n_atoms == 0 or pb[1].n_atoms == 0:
            return {"status": "skipped",
                    "reason": "could not resolve two interface partners "
                              f"({len(ctx.chain_groups(u))} protein chain(s) recorded)"}
        (nameA, A), (nameB, B) = pa, pb
        hA, hB = A.select_atoms(HEAVY), B.select_atoms(HEAVY)
        resA, resB = A.residues, B.residues
        locA = np.searchsorted(resA.resindices, hA.resindices)
        locB = np.searchsorted(resB.resindices, hB.resindices)
        nA, nB = resA.n_residues, resB.n_residues
        c_cut, i_cut, e_cut = p["contact_cutoff"], p["interface_cutoff"], p["ever_cutoff"]
        search = max(c_cut, i_cut, e_cut, 8.0)

        n = len(u.trajectory)
        times = np.empty(n)
        min_dist = np.full(n, np.nan)
        n_contacts = np.zeros(n, int)
        n_iface = np.zeros(n, int)
        occA, occB = np.zeros(nA), np.zeros(nB)
        everA, everB = np.zeros(nA, bool), np.zeros(nB, bool)
        freq = np.zeros((nA, nB))
        for i, ts in enumerate(ctx.iter_frames(u, desc="[interface]")):
            times[i] = ts.time / 1000.0
            pairs, d = capped_distance(hA.positions, hB.positions, max_cutoff=search,
                                       box=ts.dimensions, return_distances=True)
            if not len(d):
                continue
            min_dist[i] = d.min() / 10.0
            ia, ib = locA[pairs[:, 0]], locB[pairs[:, 1]]
            cm = d < c_cut
            if cm.any():
                rp = np.unique(np.stack([ia[cm], ib[cm]], axis=1), axis=0)
                freq[rp[:, 0], rp[:, 1]] += 1
                n_contacts[i] = len(rp)
            im = d < i_cut
            if im.any():
                ua, ub = np.unique(ia[im]), np.unique(ib[im])
                occA[ua] += 1
                occB[ub] += 1
                n_iface[i] = len(ua) + len(ub)
            em = d < e_cut
            everA[np.unique(ia[em])] = True
            everB[np.unique(ib[em])] = True
        freq /= n
        occA /= n
        occB /= n

        bsa = self._buried_area(ctx, A, B, p["bsa_stride"], n)
        irmsd = self._interface_rmsd(u, A, B, occA, occB, p["core_occupancy"])

        # ---- tables ------------------------------------------------------- #
        ts_df = pd.DataFrame({"time_ns": times, "min_interface_dist_nm": min_dist,
                              "n_inter_contacts": n_contacts, "n_interface_residues": n_iface,
                              "buried_area_nm2": bsa, "interface_ca_rmsd_nm": irmsd,
                              "contacts_movavg": st.moving_average(n_contacts.astype(float),
                                                                   p["window"])})
        ctx.write_csv(ts_df, "interface_timeseries.csv")

        def residue_rows(name, res, occ, ever):
            first = int(res.resids[0])
            return [{"partner": name, "resid": int(r.resid), "partner_residue": int(r.resid) - first + 1,
                     "resname": r.resname, "interface_occupancy": float(o),
                     f"within_{e_cut:g}A_ever": bool(e)}
                    for r, o, e in zip(res, occ, ever)]
        res_df = pd.DataFrame(residue_rows(nameA, resA, occA, everA) +
                              residue_rows(nameB, resB, occB, everB))
        ctx.write_csv(res_df, "interface_residues.csv")
        ctx.write_csv(pd.DataFrame(freq, index=resA.resids, columns=resB.resids),
                      "interface_contact_frequency.csv", index=True)
        fi, fj = np.where(freq > 0)
        contacts = pd.DataFrame({
            "a_resid": resA.resids[fi], "a_resname": resA.resnames[fi],
            "b_resid": resB.resids[fj], "b_resname": resB.resnames[fj],
            "occupancy": freq[fi, fj]}).sort_values("occupancy", ascending=False)
        ctx.write_csv(contacts, "interface_contacts_all.csv")
        persistent = contacts[contacts.occupancy >= p["persistent"]]
        ctx.write_csv(persistent, "interface_persistent_contacts.csv")

        self._figures(ctx, p, nameA, nameB, times, min_dist, n_contacts, bsa, irmsd,
                      resA, resB, occA, occB, freq)

        coreA = [int(r) for r in resA.resids[occA >= p["core_occupancy"]]]
        coreB = [int(r) for r in resB.resids[occB >= p["core_occupancy"]]]
        return {
            "partners": [nameA, nameB],
            "partner_residue_ranges": [[int(resA.resids[0]), int(resA.resids[-1])],
                                       [int(resB.resids[0]), int(resB.resids[-1])]],
            "cutoffs_A": {"contact": c_cut, "interface": i_cut, "ever": e_cut},
            "n_interface_residues": [int((occA > 0).sum()), int((occB > 0).sum())],
            "core_interface_residues": [coreA, coreB],          # full lists, never truncated
            "n_core_interface_residues": [len(coreA), len(coreB)],
            f"residues_within_{e_cut:g}A_ever": [[int(r) for r in resA.resids[everA]],
                                                 [int(r) for r in resB.resids[everB]]],
            "min_interface_dist_nm": st.describe(min_dist[np.isfinite(min_dist)],
                                                 "min_interface_dist_nm"),
            "frames_without_contact": int((~np.isfinite(min_dist)).sum() +
                                          (min_dist[np.isfinite(min_dist)] >= c_cut / 10).sum()),
            "n_inter_contacts": st.describe(n_contacts.astype(float), "n_inter_contacts"),
            "n_persistent_contacts": int(len(persistent)),
            "buried_area_nm2": st.describe(bsa[np.isfinite(bsa)], "buried_area_nm2")
            if np.isfinite(bsa).any() else None,
            "trend": {"contacts": _trend(times, n_contacts.astype(float)),
                      "buried_area": _trend(times, bsa),
                      "interface_residues": _trend(times, n_iface.astype(float))},
            "figure": "interface",
        }

    # ------------------------------------------------------------------ #
    @staticmethod
    def _buried_area(ctx, A, B, stride: int, n: int) -> np.ndarray:
        """Buried surface area (nm²) = SASA(A) + SASA(B) − SASA(A∪B), heavy atoms + H (mdtraj)."""
        out = np.full(n, np.nan)
        try:
            import mdtraj as md
            traj = md.load(str(ctx.config.data_dir / "core.xtc"),
                           top=str(ctx.config.data_dir / "core.pdb"), stride=stride)
        except Exception:
            return out
        ia, ib = A.indices, B.indices
        both = np.concatenate([ia, ib])
        sa = shrake_rupley(traj.atom_slice(ia), mode="atom").sum(1)
        sb = shrake_rupley(traj.atom_slice(ib), mode="atom").sum(1)
        sab = shrake_rupley(traj.atom_slice(both), mode="atom").sum(1)
        out[::stride] = (sa + sb - sab)[: len(out[::stride])]
        return out

    @staticmethod
    def _interface_rmsd(u, A, B, occA, occB, core: float) -> np.ndarray:
        """Cα RMSD of the core interface (both partners), superposed on itself, vs frame 0."""
        from MDAnalysis.analysis import rms
        n = len(u.trajectory)
        sel_res = list(A.residues[occA >= core]) + list(B.residues[occB >= core])
        if len(sel_res) < 3:
            return np.full(n, np.nan)
        ca = sum((r.atoms.select_atoms("name CA") for r in sel_res[1:]),
                 sel_res[0].atoms.select_atoms("name CA"))
        if ca.n_atoms < 3:
            return np.full(n, np.nan)
        r = rms.RMSD(ca, ca, ref_frame=0).run()
        return r.results.rmsd[:, 2] / 10.0

    def _figures(self, ctx, p, nameA, nameB, times, min_dist, n_contacts, bsa, irmsd,
                 resA, resB, occA, occB, freq) -> None:
        plt = plotting.style.plt
        rows = [("Min. heavy-atom\ndistance (nm)", min_dist, PALETTE["primary"]),
                ("Residue–residue\ncontacts", n_contacts, PALETTE["secondary"])]
        if np.isfinite(bsa).any():
            rows.append(("Buried area (nm²)", bsa, PALETTE["green"]))
        if np.isfinite(irmsd).any():
            rows.append(("Interface Cα\nRMSD (nm)", irmsd, PALETTE["purple"]))
        fig, axes = plt.subplots(len(rows), 1, figsize=(7.6, 2.2 * len(rows) + 0.6), sharex=True)
        for ax, (lab, y, col) in zip(np.atleast_1d(axes), rows):
            ok = np.isfinite(y)
            ax.plot(times[ok], y[ok], color=col, lw=0.9, alpha=0.5)
            ma = st.moving_average(np.where(ok, y, np.nan).astype(float), p["window"])
            ax.plot(times, ma, color=col, lw=2.0)
            ax.set_ylabel(lab, fontsize=10)
        np.atleast_1d(axes)[0].set_title(f"Interface: {nameA} – {nameB}")
        np.atleast_1d(axes)[-1].set_xlabel("Time (ns)")
        plotting.save_figure(fig, ctx.fig_path("interface"), dpi=ctx.config.dpi)

        fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.4, 6.0))
        a1.bar(resA.resids, occA * 100, color=PALETTE["primary"], width=1.0)
        a1.set_ylabel("Interface\noccupancy (%)")
        a1.set_title(f"{nameA} interface residues")
        a2.bar(resB.resids, occB * 100, color=PALETTE["green"], width=1.0)
        a2.set_ylabel("Interface\noccupancy (%)")
        a2.set_xlabel("Residue number")
        a2.set_title(f"{nameB} interface residues")
        plotting.save_figure(fig, ctx.fig_path("interface_residues"), dpi=ctx.config.dpi)

        fig, ax = plotting.new_axes(figsize=(7.2, 5.6))
        im = ax.imshow(freq, cmap=SEQ_CMAP, origin="lower", aspect="auto", vmin=0, vmax=1,
                       extent=[resB.resids[0], resB.resids[-1], resA.resids[0], resA.resids[-1]])
        ax.set_xlabel(f"{nameB} residue")
        ax.set_ylabel(f"{nameA} residue")
        ax.set_title("Inter-partner contact frequency")
        fig.colorbar(im, ax=ax, shrink=0.85, label="Occupancy (fraction of frames)")
        plotting.save_figure(fig, ctx.fig_path("interface_contact_map"), dpi=ctx.config.dpi)
