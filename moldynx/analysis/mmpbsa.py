"""
MM-GBSA / MM-PBSA binding free energy with gmx_MMPBSA
(https://github.com/Valdes-Tresanca-MS/gmx_MMPBSA; see :mod:`moldynx.binding`).

* If gmx_MMPBSA outputs exist in ``<output>/binding_energy/{gb,pb}`` they are
  analysed: autocorrelation-corrected means per window, drift, GB vs PB, hotspots in
  biological numbering, closure of the per-residue decomposition, entropy validity.
* Otherwise a complete, ready-to-run package is prepared (protein-only system,
  derived ionic strength, explicit decomposition residues, probe-gated run script)
  and the module reports ``prepared`` with the command -- never invented numbers.
  ``moldynx binding-energy --execute`` runs it (Linux/WSL with gmx_MMPBSA).

Results are approximate end-point estimates, not experimental affinities.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from moldynx import plotting
from moldynx.binding import CITATIONS, GMX_MMPBSA_URL
from moldynx.core import annotations as ann
from moldynx.core.base import BaseAnalysis
from moldynx.core.system import COMPLEX_SYSTEMS
from moldynx.plotting import PALETTE


def _windows(ctx, t_end: float) -> dict:
    spec = (ctx.config.binding_energy or {}).get("primary_window")
    win = {f"0-{t_end:g} ns": (0.0, t_end)}
    if isinstance(spec, (list, tuple)) and len(spec) == 2:
        win[f"{spec[0]:g}-{spec[1]:g} ns"] = (float(spec[0]), float(spec[1]))
    elif isinstance(spec, str) and re.match(r"^\s*[\d.]+\s*-\s*[\d.]+", spec):
        a, b = (float(x) for x in re.findall(r"[\d.]+", spec)[:2])
        win[f"{a:g}-{b:g} ns"] = (a, b)
    else:
        win[f"{0.8 * t_end:g}-{t_end:g} ns"] = (0.8 * t_end, t_end)
    return win


class MMPBSA(BaseAnalysis):
    name = "mmpbsa"
    label = "MM-GBSA / MM-PBSA binding free energy (gmx_MMPBSA)"
    category = "interactions"
    required_files = {"trajectory", "topology"}
    supported_systems = COMPLEX_SYSTEMS
    order = 170                  # after interface (residue set) and analysis_window
    default_params = {"target_spacing_ns": 1.0, "ever_cutoff_column": "within_6A_ever"}
    outputs = ["binding_energy/", "results/binding_energy_summary.json",
               "figures/binding_energy_timeseries.png", "figures/binding_energy_hotspots.png"]

    def run(self, ctx) -> dict:
        be_dir = ctx.config.output_dir / "binding_energy"
        if (be_dir / "gb" / "FINAL_RESULTS_MMGBSA.csv").exists():
            return self.analyse_outputs(ctx, be_dir)
        return self.prepare(ctx, be_dir)

    # ------------------------------------------------------------------ #
    def _partners(self, ctx):
        chains = sorted(ctx.core_meta.get("chains", []), key=lambda c: c["core_stop"] - c["core_start"],
                        reverse=True)[:2]
        if len(chains) < 2:
            return None
        roles = {c.get("segid"): c.get("role") for c in (ctx.config.annotations or {}).get("chains", [])}
        rec = next((c for c in chains if roles.get(c["segid"]) == "receptor"), chains[0])
        lig = next(c for c in chains if c is not rec)
        order = sorted(chains, key=lambda c: c["core_start"])
        return rec, lig, {c["segid"]: "ABCDEFGH"[k] for k, c in enumerate(order)}

    def prepare(self, ctx, be_dir: Path) -> dict:
        from moldynx.binding import prepare as prep
        from moldynx.io import gromacs
        fs = ctx.fileset
        missing = [k for k in ("gmx_top", "topology") if getattr(fs, k, None) is None]
        if missing or Path(fs.topology).suffix.lower() != ".tpr":
            return {"status": "skipped",
                    "reason": f"binding energy needs the run input (.tpr) and topol.top + toppar/ "
                              f"(missing: {', '.join(missing) or 'a .tpr run input'})"}
        partners = self._partners(ctx)
        if partners is None:
            return {"status": "skipped", "reason": "fewer than two protein chains"}
        rec, lig, letters = partners
        u = ctx.core_universe()
        n = len(u.trajectory)
        dt = float(u.trajectory[1].time - u.trajectory[0].time) / 1000.0 if n > 1 else 1.0
        interval = max(1, int(round(self.params(ctx)["target_spacing_ns"] / dt))) if dt else 1
        u.trajectory[0]
        volume = abs(float(np.linalg.det(u.trajectory.ts.triclinic_dimensions))) / 1000.0
        ionic = prep.derive_ionic_strength(getattr(ctx.system, "ion_counts", {}) or {}, volume)
        if not ionic["ions"]:
            ionic = {**ionic, "ionic_strength_M": 0.0,
                     "note": "no ions found in the run input: ionic strength 0 (set "
                             "binding_energy.ionic_strength to override)"}
        override = (ctx.config.binding_energy or {}).get("ionic_strength")
        if isinstance(override, (int, float)):
            ionic = {**ionic, "ionic_strength_M": float(override),
                     "note": f"user-specified ({override} M); derived value kept for reference"}
        temp = 300.0
        if fs.log is not None:
            ref = gromacs.parse_log(fs.log).ref_t
            temp = ref[0] if ref else temp
        local = {"R": [], "L": []}
        res_csv = ctx.csv_path("interface_residues.csv")
        if res_csv.exists():
            df = pd.read_csv(res_csv)
            col = next((c for c in df.columns if c.startswith("within_") and c.endswith("_ever")), None)
            if col:
                for side, c in (("R", rec), ("L", lig)):
                    sel = df[(df.partner == c["segid"]) & df[col].astype(bool)]
                    local[side] = sel.partner_residue.astype(int).tolist()
        solute_xtc = ctx.config.data_dir / "core.xtc"
        if u.atoms.n_atoms != sum(c["core_stop"] - c["core_start"] for c in (rec, lig)):
            return {"status": "skipped",
                    "reason": "the solute trajectory holds more than the two partners; "
                              "binding energy for mixed solutes is not automated yet"}
        top_text = Path(fs.gmx_top).read_text(encoding="utf-8", errors="replace")
        g = lambda c: {"name": c.get("segid"), "letter": letters[c["segid"]],  # noqa: E731
                       "atom_start1": c["core_start"] + 1, "atom_stop1": c["core_stop"]}
        info = prep.write_package(
            be_dir, tpr_src=Path(fs.topology), top_text=top_text,
            toppar_src=Path(fs.toppar) if fs.toppar else None, solute_xtc=solute_xtc,
            receptor=g(rec), ligand=g(lig), temperature=temp, ionic=ionic,
            frames={"start": 1, "end": n, "interval": interval, "frame_dt_ns": dt},
            print_res_local=local, meta={"citations": CITATIONS})
        gmx = gromacs.find_gmx()
        tpr_ok, why = (prep.make_complex_tpr(gmx, Path(fs.topology), be_dir) if gmx
                       else (False, "GROMACS not found: run make_ndx/convert-tpr as described"))
        return {"status": "prepared", "directory": str(be_dir), "complex_tpr": tpr_ok,
                "complex_tpr_note": why, "receptor": rec["segid"], "ligand": lig["segid"],
                "ionic_strength_M": ionic["ionic_strength_M"], "temperature_K": temp,
                "frames": info["frames"], "print_res": info["print_res"],
                "how_to_run": f"on Linux/WSL with gmx_MMPBSA ({GMX_MMPBSA_URL}): "
                              f"cd {be_dir} && bash run_mmpbsa.sh  (or: moldynx binding-energy "
                              f"--run {ctx.config.output_dir} --execute)"}

    # ------------------------------------------------------------------ #
    def analyse_outputs(self, ctx, be_dir: Path) -> dict:
        from moldynx.binding.analyse import analyse
        meta_p = be_dir / "prepare_meta.json"
        meta = json.loads(meta_p.read_text(encoding="utf-8")) if meta_p.exists() else {}
        dt = meta.get("frames", {}).get("frame_dt_ns") or (ctx.config.binding_energy or {}).get(
            "frame_dt_ns", 0.1)
        chains = ann.chain_annotations(ctx.core_meta.get("chains", []), ctx.config.annotations or {})
        rec = chains.get(meta.get("receptor", {}).get("name"))
        lig = chains.get(meta.get("ligand", {}).get("name"))
        mapper = lambda c: (lambda k: (c.display, c.bio(c.resid_first + k - 1))) if c else None  # noqa: E731
        gb = be_dir / "gb"
        pb_csv = be_dir / "pb" / "FINAL_RESULTS_MMPBSA.csv"
        dec = gb / "FINAL_DECOMP_MMGBSA.csv"
        from moldynx.binding.analyse import read_delta
        t_end = float(read_delta(gb / "FINAL_RESULTS_MMGBSA.csv", dt).time_ns.max())
        dat = gb / "FINAL_RESULTS_MMGBSA.dat"
        res = analyse(gb / "FINAL_RESULTS_MMGBSA.csv", pb_csv if pb_csv.exists() else None,
                      frame_dt_ns=dt, windows=_windows(ctx, t_end),
                      decomp_csv=dec if dec.exists() else None,
                      gb_dat=dat if dat.exists() else None,
                      receptor=mapper(rec), ligand=mapper(lig))
        res["components"].to_csv(ctx.csv_path("binding_energy_components.csv"), index=False)
        frames = res["gb"][["frame", "time_ns", "TOTAL"]].rename(columns={"TOTAL": "dG_GB"})
        if res["pb"] is not None:
            frames["dG_PB"] = res["pb"].TOTAL.values
        frames.to_csv(ctx.csv_path("binding_energy_per_frame.csv"), index=False)
        for w, df in res.get("decomposition", {}).items():
            df.to_csv(ctx.csv_path(f"binding_energy_residues_{w.replace(' ', '')}.csv"), index=False)
        self._figures(ctx, res)
        summary = {k: res[k] for k in ("checks", "windows", "headline", "trend_per_ns", "gb_vs_pb")}
        summary.update({"entropy": res.get("entropy"), "closure": res.get("closure"),
                        "residues_in_decomposition": res.get("residues_in_decomposition"),
                        "hotspots": {w: df.head(20)[["partner", "residue", "resname", "TOTAL"]]
                                     .to_dict("records")
                                     for w, df in res.get("decomposition", {}).items()},
                        "method": {"tool": GMX_MMPBSA_URL, "citations": CITATIONS,
                                   "prepared": meta},
                        "caveat": "approximate end-point estimates, not experimental affinities"})
        ctx.csv_path("binding_energy_summary.json").write_text(
            json.dumps(summary, indent=2, default=float), encoding="utf-8")
        summary["status"] = "analysed"
        summary["figure"] = "binding_energy_timeseries"
        return summary

    def _figures(self, ctx, res) -> None:
        plotting.set_style()
        plt = plotting.style.plt
        gb, pb = res["gb"], res["pb"]
        fig, ax = plotting.new_axes(figsize=(8, 4.4))
        ax.plot(gb.time_ns, gb.TOTAL, color=PALETTE["primary"], lw=1.2, label="MM-GBSA")
        if pb is not None:
            ax.plot(pb.time_ns, pb.TOTAL, color=PALETTE["secondary"], lw=1.2, label="MM-PBSA")
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("ΔG_bind (kcal/mol)")
        ax.set_title("Binding free energy over time (end-point estimate)")
        ax.legend()
        plotting.save_figure(fig, ctx.fig_path("binding_energy_timeseries"), dpi=ctx.config.dpi)
        dec = res.get("decomposition", {})
        if dec:
            w = list(dec)[-1]
            top = dec[w].head(20).iloc[::-1]
            fig, ax = plt.subplots(figsize=(7, 6))
            lab = [f"{p or s} {r}{n}" for p, s, r, n in zip(top.partner, top.side, top.resname,
                                                              top.residue)]
            ax.barh(lab, top.TOTAL, color=PALETTE["primary"])
            ax.set_xlabel("Per-residue ΔG (kcal/mol)")
            ax.set_title(f"Top 20 residues, {w}")
            plotting.save_figure(fig, ctx.fig_path("binding_energy_hotspots"), dpi=ctx.config.dpi)
