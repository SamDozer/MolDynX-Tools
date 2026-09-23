"""
Preparation and equilibration audit: minimisation, NVT, NPT -- and how they fed production.

From the stage logs, energy files, run inputs and job scripts found at intake:

* stage parameters (integrator, dt, steps, thermostat/barostat and targets, restraints
  reference scaling, continuation), sessions, timings, warnings;
* the minimisation outcome -- whether the requested force tolerance was reached --
  and where the largest residual force sits (chain / residue / atom);
* crash-dump accounting (``stepN{b,c}.pdb`` vs ``Wrote pdb`` events);
* energy-file statistics per stage: whole-stage and tail means, residual drift,
  time for temperature and density to settle;
* position restraints stored in each run input (and none in production);
* protonation states of titratable residues in the simulated system;
* the chain of custody (grompp -c/-r/-t, convert-tpr) and the stage timeline.

Missing inputs are reported, never invented.
"""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from moldynx import plotting
from moldynx.core.base import BaseAnalysis
from moldynx.io import gromacs
from moldynx.io.jobscripts import chain_of_custody
from moldynx.plotting import PALETTE

STAGES = ("em", "nvt", "npt")
LABEL = {"em": "Energy minimisation", "nvt": "NVT equilibration", "npt": "NPT equilibration",
         "production": "Production"}
TERMS = ("Temperature", "Pressure", "Density", "Volume", "Box-X", "Potential", "Total Energy")

# residue/atom signatures of non-default protonation states (CHARMM and AMBER names)
PROTONATION = [
    ("GLU", "HE2", "protonated glutamate"), ("GLUP", None, "protonated glutamate"),
    ("GLH", None, "protonated glutamate"), ("ASP", "HD2", "protonated aspartate"),
    ("ASPP", None, "protonated aspartate"), ("ASH", None, "protonated aspartate"),
    ("HSP", None, "doubly protonated histidine"), ("HIP", None, "doubly protonated histidine"),
    ("HSE", None, "histidine, Nε-H tautomer"), ("HIE", None, "histidine, Nε-H tautomer"),
    ("LSN", None, "neutral lysine"), ("LYN", None, "neutral lysine"),
    ("CYM", None, "deprotonated cysteine"),
]


# --------------------------------------------------------------------------- #
# statistics helpers
# --------------------------------------------------------------------------- #
def _drift(t: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    """Linear slope per ns over the second half, and its p-value."""
    from scipy import stats
    h = len(x) // 2
    if len(x) - h < 5 or np.ptp(t[h:]) == 0:
        return float("nan"), float("nan")
    r = stats.linregress(t[h:], x[h:])
    return float(r.slope * 1000.0), float(r.pvalue)


def energy_statistics(df: pd.DataFrame, tail_ps: float = 500.0) -> dict:
    """Whole-stage and tail statistics of the standard terms of one stage."""
    t = df["Time"].to_numpy(float)
    dur = float(t[-1] - t[0]) if len(t) > 1 else 0.0
    tail = min(tail_ps, 0.5 * dur) if dur > 0 else 0.0
    tail_mask = t >= t[-1] - tail
    out: dict = {"n_records": int(len(df)), "t_first_ps": float(t[0]), "t_last_ps": float(t[-1]),
                 "tail_ps": tail, "terms": {}}
    for term in TERMS:
        if term not in df:
            continue
        x = df[term].to_numpy(float)
        slope, p = _drift(t, x)
        out["terms"][term] = {
            "mean": float(x.mean()), "sd": float(x.std(ddof=1)) if len(x) > 1 else 0.0,
            "first": float(x[0]), "last": float(x[-1]),
            "tail_mean": float(x[tail_mask].mean()),
            "tail_sd": float(x[tail_mask].std(ddof=1)) if tail_mask.sum() > 1 else 0.0,
            "drift_per_ns": slope, "drift_p": p,
        }
    return out


def settle_time(t: np.ndarray, x: np.ndarray, target: float, tol: float) -> float | None:
    """First time after which ``x`` stays within ``target ± tol`` until the end."""
    inside = np.abs(x - target) <= tol
    if not inside[-1]:
        return None
    bad = np.where(~inside)[0]
    return float(t[0] if not len(bad) else t[min(bad[-1] + 1, len(t) - 1)])


def fraction_time(t: np.ndarray, x: np.ndarray, frac: float = 0.999) -> float | None:
    """First time ``x`` comes within ``1 - frac`` of its final value (mean of the last 10 %)."""
    if len(x) < 10:
        return None
    final = x[-max(5, len(x) // 10):].mean()
    ok = np.abs(x - final) <= abs(final) * (1 - frac)
    idx = np.argmax(ok) if ok.any() else None
    return float(t[idx]) if idx is not None else None


class EquilibrationAudit(BaseAnalysis):
    name = "equilibration"
    label = "Preparation and equilibration audit"
    category = "preparation"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    order = 2
    default_params = {"tail_ps": 500.0, "nvt_settle_tol_K": 2.0}
    outputs = ["results/equilibration_summary.json", "results/equilibration_stats.csv",
               "figures/equilibration_overview.png"]

    # ------------------------------------------------------------------ #
    def run(self, ctx) -> dict:
        p = self.params(ctx)
        fs = ctx.fileset
        present = {s: {k: [str(x) for x in v] for k, v in fs.stages.get(s, {}).items()}
                   for s in STAGES}
        if not any(fs.stage_file(s, k) for s in STAGES for k in ("log", "edr", "tpr")):
            return {"status": "skipped",
                    "reason": "no minimisation / NVT / NPT files were found with the production run"}
        summary: dict = {"stages": {}, "missing": [], "files": present}
        logs: dict[str, gromacs.LogInfo] = {}
        for stage in STAGES + ("production",):
            log = fs.log if stage == "production" else fs.stage_file(stage, "log")
            if log is None:
                if stage != "production":
                    summary["missing"].append(f"{LABEL[stage]}: no .log")
                continue
            info = gromacs.parse_log(log)
            logs[stage] = info
            summary["stages"][stage] = self._stage_record(info)
        for stage in STAGES:
            if fs.stage_file(stage, "edr") is None:
                summary["missing"].append(f"{LABEL[stage]}: no .edr")

        # -- energies ----------------------------------------------------- #
        series = self._energies(ctx, p, summary)
        # -- minimisation: Fmax location + crash dumps ---------------------- #
        if "em" in logs and logs["em"].minimization:
            summary["stages"]["em"]["fmax_location"] = self._locate_atom(
                ctx, logs["em"].minimization.fmax_atom)
        summary["crash_dumps"] = self._crash_dumps(fs, logs.get("em"))
        # -- restraints ---------------------------------------------------- #
        summary["position_restraints"] = self._restraints(ctx, fs)
        # -- protonation, custody, timeline -------------------------------- #
        summary["protonation"] = self._protonation(ctx)
        summary["chain_of_custody"] = chain_of_custody(fs.job_scripts)
        summary["timeline"] = self._timeline(logs)
        summary["mdp_files"] = {k: str(v) for k, v in fs.mdp.items()}
        # thermostat targets of the dynamics stages (minimisation has no temperature)
        summary["temperatures_K"] = {s: logs[s].ref_t[:1] for s in logs
                                     if logs[s].ref_t and not logs[s].is_minimization
                                     and any(t > 0 for t in logs[s].ref_t)}
        summary["verdicts"] = self._verdicts(summary)

        ctx.csv_path("equilibration_summary.json").write_text(
            json.dumps(summary, indent=2, default=str), encoding="utf-8")
        fig = self._figure(ctx, series, summary)
        brief = {
            "stages": {s: {k: v for k, v in r.items()
                           if k in ("integrator", "dt_ps", "nsteps", "simulated_ps", "tcoupl",
                                    "ref_t_K", "pcoupl", "minimization")}
                       for s, r in summary["stages"].items()},
            "missing": summary["missing"], "verdicts": summary["verdicts"],
            "restrained_atoms": {s: r.get("n_restrained") for s, r in
                                 summary["position_restraints"].items()},
            "protonation": summary["protonation"].get("non_default", []),
        }
        if fig:
            brief["figure"] = fig
        return brief

    # ------------------------------------------------------------------ #
    @staticmethod
    def _stage_record(info: gromacs.LogInfo) -> dict:
        m = info.mdp
        rec = {"log": info.path, "gromacs": info.gromacs_version, "commands": info.commands,
               "sessions": info.n_sessions, "started": info.started[:1],
               "finished": info.finished[-1:], "ns_per_day": info.performance_ns_day[-1:],
               "wall_h": round(sum(info.wall_s) / 3600, 3) if info.wall_s else None,
               "integrator": m.get("integrator"), "dt_ps": m.get("dt"), "nsteps": m.get("nsteps"),
               "simulated_ps": info.simulated_ps, "tcoupl": m.get("tcoupl"),
               "ref_t_K": info.ref_t, "tau_t_ps": info.grpopts.get("tau-t"),
               "nrdf": info.grpopts.get("nrdf"), "pcoupl": m.get("pcoupl"),
               "pcoupltype": m.get("pcoupltype"), "tau_p_ps": m.get("tau-p"),
               "ref_p_bar": m.get("ref-p_xx"), "compressibility": m.get("compressibility_xx"),
               "refcoord_scaling": m.get("refcoord-scaling"), "continuation": m.get("continuation"),
               "emtol": m.get("emtol"), "coulombtype": m.get("coulombtype"),
               "rcoulomb": m.get("rcoulomb"), "rvdw": m.get("rvdw"),
               "vdw_modifier": m.get("vdw-modifier"), "fourier_nx": m.get("fourier-nx"),
               "constraint_algorithm": m.get("constraint-algorithm"),
               "nstlist": m.get("nstlist"), "counts": info.counts}
        if info.minimization:
            rec["minimization"] = info.minimization.__dict__
        return rec

    def _energies(self, ctx, p, summary) -> dict:
        try:
            import panedr
        except ImportError:
            summary["missing"].append("panedr not installed — energy-file statistics skipped "
                                      "(pip install panedr)")
            return {}
        series = {}
        for stage in STAGES:
            edr = ctx.fileset.stage_file(stage, "edr")
            if edr is None:
                continue
            try:
                df = panedr.edr_to_df(str(edr))
            except Exception as exc:  # corrupt / unreadable energy file
                summary["missing"].append(f"{LABEL[stage]}: energy file unreadable ({exc})")
                continue
            keep = ["Time"] + [c for c in TERMS if c in df]
            series[stage] = df[keep]
            df[keep].to_csv(ctx.csv_path(f"equilibration_timeseries_{stage}.csv"), index=False)
            st = energy_statistics(df, p["tail_ps"])
            rec = summary["stages"].setdefault(stage, {})
            rec["energy"] = st
            t = df["Time"].to_numpy(float)
            ref = rec.get("ref_t_K") or []
            if stage == "nvt" and "Temperature" in df and ref:
                rec["temperature_settled_ps"] = settle_time(
                    t, df["Temperature"].to_numpy(float), ref[0], p["nvt_settle_tol_K"])
                if rec["temperature_settled_ps"] is not None:
                    m = t >= rec["temperature_settled_ps"]
                    x = df["Temperature"].to_numpy(float)[m]
                    rec["temperature_after_settling"] = {"mean": float(x.mean()),
                                                         "sd": float(x.std(ddof=1))}
            if stage == "npt" and "Density" in df:
                rec["density_reached_99_9pct_ps"] = fraction_time(t, df["Density"].to_numpy(float))
            if "Box-X" in df:
                rec["box_x_nm"] = [float(df["Box-X"].iloc[0]), float(df["Box-X"].iloc[-1])]
            if "Volume" in df:
                v0, v1 = float(df["Volume"].iloc[0]), float(df["Volume"].iloc[-1])
                rec["volume_nm3"] = [v0, v1, (v1 / v0 - 1) * 100 if v0 else None]
        rows = []
        for stage, rec in summary["stages"].items():
            for term, s in rec.get("energy", {}).get("terms", {}).items():
                rows.append({"stage": stage, "term": term, **s})
        if rows:
            pd.DataFrame(rows).to_csv(ctx.csv_path("equilibration_stats.csv"), index=False)
        return series

    @staticmethod
    def _locate_atom(ctx, atom_number: int | None) -> dict | None:
        """GROMACS 1-based atom number -> chain / residue / atom (from the run input)."""
        if not atom_number:
            return None
        try:
            u = ctx.full_universe()
            a = u.atoms[atom_number - 1]
        except Exception:
            return None
        loc = {"atom_number": atom_number, "atom": a.name, "resname": a.resname,
               "resid": int(a.resid), "segid": a.segid}
        for c in getattr(ctx.system, "chains", []):
            if c.atom_start <= atom_number - 1 < c.atom_stop:
                loc["chain"] = c.segid
                loc["chain_residue"] = int(a.resid) - c.resid_first + 1
        return loc

    @staticmethod
    def _crash_dumps(fs, em_log) -> dict:
        steps = sorted({int(m.group(1)) for p in fs.crash_dumps
                        if (m := re.match(r"step-?(\d+)[bc]\.pdb$", p.name, re.I))})
        pairs = len(steps)
        events = em_log.counts.get("wrote_pdb", 0) if em_log else None
        return {"files": len(fs.crash_dumps), "steps": steps, "pairs": pairs,
                "wrote_pdb_events_in_em_log": events,
                "all_accounted_for": (events is not None and events == pairs) if pairs else True}

    @staticmethod
    def _restraints(ctx, fs) -> dict:
        gmx = gromacs.find_gmx()
        out = {}
        for stage in STAGES:
            tpr = fs.stage_file(stage, "tpr")
            if tpr is not None:
                out[stage] = gromacs.tpr_position_restraints(gmx, tpr).to_dict()
        if fs.topology is not None and Path(fs.topology).suffix.lower() == ".tpr":
            out["production"] = gromacs.tpr_position_restraints(gmx, fs.topology).to_dict()
        return out

    @staticmethod
    def _protonation(ctx) -> dict:
        try:
            u = ctx.full_universe()
            prot = u.select_atoms("protein")
        except Exception as exc:
            return {"error": str(exc)}
        found = []
        chains = getattr(ctx.system, "chains", [])
        for res in prot.residues:
            rn = res.resname.upper()
            names = set(res.atoms.names)
            for resname, atom, what in PROTONATION:
                if rn == resname and (atom is None or atom in names):
                    rec = {"resname": rn, "resid": int(res.resid), "state": what}
                    for c in chains:
                        if c.atom_start <= int(res.atoms.indices[0]) < c.atom_stop:
                            rec["chain"] = c.segid
                            rec["chain_residue"] = int(res.resid) - c.resid_first + 1
                    found.append(rec)
                    break
        his = {}
        for res in prot.residues:
            if res.resname.upper() in ("HSD", "HSE", "HSP", "HID", "HIE", "HIP", "HIS"):
                his[res.resname.upper()] = his.get(res.resname.upper(), 0) + 1
        return {"non_default": [r for r in found if "tautomer" not in r["state"]],
                "histidine_tautomers": his,
                "note": "HSD/HID (Nδ-H) is the CHARMM/AMBER default neutral histidine"}

    @staticmethod
    def _timeline(logs: dict) -> list[dict]:
        rows = []
        for stage in STAGES + ("production",):
            info = logs.get(stage)
            if info is None:
                continue
            start = info.started[0] if info.started else None
            end = info.finished[-1] if info.finished else None
            rows.append({"stage": stage, "started": start, "finished": end})
        for prev, cur in zip(rows, rows[1:]):
            if prev["finished"] and cur["started"]:
                gap = (datetime.fromisoformat(cur["started"]) -
                       datetime.fromisoformat(prev["finished"])).total_seconds()
                cur["gap_after_previous_s"] = gap
        return rows

    @staticmethod
    def _verdicts(summary: dict) -> list[str]:
        """Plain statements of what the evidence shows (documents word them further)."""
        v = []
        em = summary["stages"].get("em", {}).get("minimization")
        if em:
            if em.get("reached_emtol") is False:
                v.append(f"minimisation stopped ({em.get('outcome')}) without reaching the "
                         f"requested force tolerance; Fmax = {em.get('fmax'):.4g} kJ/mol/nm")
            elif em.get("reached_emtol"):
                v.append(f"minimisation {em.get('outcome')} in {em.get('steps')} steps")
        warn = sum(r.get("counts", {}).get("lincs_warnings", 0) +
                   r.get("counts", {}).get("fatal_errors", 0)
                   for s, r in summary["stages"].items() if s != "em")
        v.append("no LINCS warnings or fatal errors in the MD stages" if warn == 0 else
                 f"{warn} LINCS warnings / fatal errors in the MD stages")
        npt = summary["stages"].get("npt", {}).get("energy", {}).get("terms", {})
        for term, unit in (("Density", "kg m⁻³"), ("Volume", "nm³"), ("Potential", "kJ mol⁻¹")):
            s = npt.get(term)
            if s and s["drift_p"] == s["drift_p"] and s["drift_p"] < 0.01:
                v.append(f"NPT {term.lower()} still drifting at the end of the stage "
                         f"({s['drift_per_ns']:+.3g} {unit} per ns, p = {s['drift_p']:.1g})")
        temps = summary.get("temperatures_K", {})
        distinct = {tuple(x) for x in temps.values() if x}
        if len(distinct) > 1:
            v.append("reference temperature changes between stages: "
                     + ", ".join(f"{LABEL[s]} {x[0]:g} K" for s, x in temps.items() if x))
        pr = summary.get("position_restraints", {})
        if pr.get("production", {}).get("available"):
            n = pr["production"]["n_restrained"]
            v.append("production run input has no position restraints" if n == 0 else
                     f"production run input still contains {n} position restraints")
        cd = summary.get("crash_dumps", {})
        if cd.get("files"):
            v.append(f"{cd['pairs']} minimisation crash-dump pairs; "
                     + ("all matched to 'Wrote pdb' events in the minimisation log"
                        if cd.get("all_accounted_for") else
                        f"{cd.get('wrote_pdb_events_in_em_log')} 'Wrote pdb' events in the "
                        "surviving minimisation log — the rest cannot be matched to a log"))
        if not summary.get("mdp_files"):
            v.append("stage .mdp files are not present: MDP-only settings (gen-vel, define) "
                     "are not recorded")
        return v

    # ------------------------------------------------------------------ #
    def _figure(self, ctx, series: dict, summary: dict) -> str | None:
        if not series:
            return None
        plotting.set_style()
        plt = plotting.style.plt
        panels = []
        if "em" in series and "Potential" in series["em"]:
            panels.append(("em", "Potential", "Potential energy (10⁶ kJ mol⁻¹)", 1e-6, "step"))
        if "nvt" in series and "Temperature" in series["nvt"]:
            panels.append(("nvt", "Temperature", "Temperature (K)", 1, "ps"))
        for term, lab in (("Temperature", "Temperature (K)"), ("Pressure", "Pressure (bar)"),
                          ("Density", "Density (kg m⁻³)"), ("Volume", "Volume (nm³)")):
            if "npt" in series and term in series["npt"]:
                panels.append(("npt", term, lab, 1, "ns"))
        if not panels:
            return None
        ncol = 3 if len(panels) > 4 else 2
        nrow = int(np.ceil(len(panels) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 3.3 * nrow), squeeze=False)
        for ax in axes.flat[len(panels):]:
            ax.axis("off")
        for ax, (stage, term, lab, scale, tunit) in zip(axes.flat, panels):
            df = series[stage]
            t = df["Time"].to_numpy(float)
            x = df[term].to_numpy(float) * scale
            tt = t / 1000.0 if tunit == "ns" else t
            ax.plot(tt, x, lw=0.8, color=PALETTE["primary"])
            ref = summary["stages"].get(stage, {}).get("ref_t_K") or []
            if term == "Temperature" and ref:
                ax.axhline(ref[0], ls="--", lw=1, color=PALETTE["accent"])
            if term == "Pressure":
                run = pd.Series(x).rolling(50, center=True, min_periods=10).mean()
                ax.plot(tt, run, lw=1.4, color=PALETTE["secondary"])
                lo, hi = np.nanpercentile(x, [1, 99])
                ax.set_ylim(lo - 0.1 * abs(hi - lo), hi + 0.1 * abs(hi - lo))
            ax.set_title(f"{LABEL[stage]} — {term.lower()}", fontsize=10)
            ax.set_xlabel({"step": "Step", "ps": "Time (ps)", "ns": "Time (ns)"}[tunit],
                          fontsize=9)
            ax.set_ylabel(lab, fontsize=9)
            ax.tick_params(labelsize=8)
        plotting.save_figure(fig, ctx.fig_path("equilibration_overview"), dpi=ctx.config.dpi)
        return "equilibration_overview"
