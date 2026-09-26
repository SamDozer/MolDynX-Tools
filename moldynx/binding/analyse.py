"""
Analysis of gmx_MMPBSA outputs (``FINAL_RESULTS_*.csv``, ``FINAL_DECOMP_MMGBSA.csv``, ``*.dat``).

All means carry autocorrelation-corrected standard errors (statistical inefficiency,
N_eff = N/g). Values are approximate end-point estimates, **not** experimental
affinities; entropy terms are reported only if σ(interaction energy) passes the
validity limits gmx_MMPBSA itself warns about (IE < 3.6, C2 < 6.0 kcal/mol).
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np
import pandas as pd

from moldynx.statistics import describe_correlated, drift

BONDED = ["BOND", "ANGLE", "DIHED", "UB", "IMP", "CMAP", "1-4 VDW", "1-4 EEL"]
COMPONENTS = [("ΔE_vdW", "VDWAALS", "VDWAALS"), ("ΔE_elec", "EEL", "EEL"),
              ("ΔG_polar", "EGB", "EPB"), ("ΔG_nonpolar", "ESURF", "ENPOLAR"),
              ("ΔG_gas", "GGAS", "GGAS"), ("ΔG_solv", "GSOLV", "GSOLV"),
              ("ΔG_bind", "TOTAL", "TOTAL")]
ENTROPY_LIMITS = {"IE": 3.6, "C2": 6.0}          # kcal/mol, as gmx_MMPBSA warns


def read_delta(csv: str | Path, frame_dt_ns: float) -> pd.DataFrame:
    """The ``Delta Energy Terms`` block, with time from the original frame numbers."""
    lines = Path(csv).read_text(encoding="utf-8").splitlines()
    i = next(k for k, line in enumerate(lines) if line.strip().startswith("Delta Energy Terms"))
    rows = []
    for line in lines[i + 1:]:
        if line.startswith("Frame") or (line[:1].isdigit()):
            rows.append(line)
        elif rows:
            break
    df = pd.read_csv(io.StringIO("\n".join(rows))).rename(columns={"Frame #": "frame"})
    df["time_ns"] = (df.frame - 1) * frame_dt_ns
    return df


def read_decomp(csv: str | Path, frame_dt_ns: float, section: str = "DELTAS",
                kind: str = "Total Decomposition Contribution") -> pd.DataFrame:
    """Per-residue decomposition rows (``L::THR:1`` / ``R::GLU:205`` labels split out)."""
    lines = Path(csv).read_text(encoding="utf-8").splitlines()
    i = next(k for k, line in enumerate(lines) if line.strip().startswith(section))
    j = next(k for k in range(i, len(lines)) if lines[k].strip().startswith(kind))
    rows = [lines[j + 1]]
    for line in lines[j + 2:]:
        if not line or not line[:1].isdigit():
            break
        rows.append(line)
    df = pd.read_csv(io.StringIO("\n".join(rows))).rename(columns={"Frame #": "frame",
                                                                    "Residue": "label"})
    m = df.label.str.extract(r"^([LR])::([A-Z0-9]+):(-?\d+)$")
    df["side"], df["resname"], df["resnum"] = m[0], m[1], m[2].astype(int)
    df["time_ns"] = (df.frame - 1) * frame_dt_ns
    return df


def read_entropy(dat: str | Path) -> dict:
    """σ(interaction energy) and −TΔS for IE and C2 from a ``FINAL_RESULTS_*.dat``."""
    text = Path(dat).read_text(encoding="utf-8", errors="replace")
    out = {}
    for kind in ("IE", "C2"):
        m = re.search(rf"^\S+\s+{kind}\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", text, re.M)
        if m:
            sigma, value = float(m.group(1)), float(m.group(2))
            out[kind] = {"sigma_int_kcal": sigma, "minus_TdS_kcal": value,
                         "sem_kcal": float(m.group(4)), "limit_kcal": ENTROPY_LIMITS[kind],
                         "valid": sigma < ENTROPY_LIMITS[kind]}
    return out


def _window_mask(t: np.ndarray, window) -> np.ndarray:
    lo, hi = window
    return (t >= lo - 1e-9) & (t <= hi + 1e-9)


def analyse(gb_csv, pb_csv=None, frame_dt_ns: float = 0.1, windows: dict | None = None,
            decomp_csv=None, gb_dat=None, receptor=None, ligand=None) -> dict:
    """
    Parameters
    ----------
    windows : {"0-100 ns": (0, 100), ...}; defaults to the full run and its final 20 %.
    receptor / ligand : callables mapping a gmx_MMPBSA residue number (numbered from 1
        within each partner) to a (display, biological number) tuple; identity if None.
    """
    gb = read_delta(gb_csv, frame_dt_ns)
    pb = read_delta(pb_csv, frame_dt_ns) if pb_csv else None
    t = gb.time_ns.to_numpy()
    if windows is None:
        end = float(t[-1])
        windows = {f"{t[0]:g}-{end:g} ns": (float(t[0]), end),
                   f"{0.8 * end:g}-{end:g} ns": (0.8 * end, end)}
    checks = {"frames": int(len(gb)),
              "max_abs_bonded_delta": float(max(np.abs(gb[c]).max() for c in BONDED if c in gb)),
              "same_frames_gb_pb": bool(pb is None or (gb.frame.values == pb.frame.values).all())}
    if pb is not None:
        checks["max_abs_ggas_gb_minus_pb"] = float(np.abs(gb.GGAS - pb.GGAS).max())
    checks["single_trajectory_consistent"] = checks["max_abs_bonded_delta"] < 1e-6

    rows = []
    for method, df in (("GB", gb), ("PB", pb)):
        if df is None:
            continue
        for wname, w in windows.items():
            mask = _window_mask(df.time_ns.to_numpy(), w)
            for label, cg, cp in COMPONENTS:
                col = cg if method == "GB" else cp
                if col not in df:
                    continue
                d = describe_correlated(df[col].to_numpy()[mask])
                rows.append({"method": method, "window": wname, "term": label, **d})
    comp = pd.DataFrame(rows)
    headline = comp[comp.term == "ΔG_bind"][["method", "window", "mean", "sem", "sd", "n",
                                              "n_eff", "stat_ineff"]].reset_index(drop=True)
    trend = {m: drift(df.time_ns.to_numpy(), df.TOTAL.to_numpy())
             for m, df in (("GB", gb), ("PB", pb)) if df is not None}
    gbpb = None
    if pb is not None:
        r = float(np.corrcoef(gb.TOTAL, pb.TOTAL)[0, 1])
        gbpb = {"pearson_r": r, "mean_offset_pb_minus_gb": float((pb.TOTAL - gb.TOTAL).mean())}

    out = {"checks": checks, "windows": {k: list(v) for k, v in windows.items()},
           "headline": headline.to_dict("records"), "components": comp,
           "trend_per_ns": {m: {"slope": d["slope"], "p": d["p_slope"], "half_p": d["half_p"]}
                            for m, d in trend.items()},
           "gb_vs_pb": gbpb, "gb": gb, "pb": pb}
    if gb_dat:
        out["entropy"] = read_entropy(gb_dat)
    if decomp_csv:
        out.update(_decomposition(decomp_csv, frame_dt_ns, windows, gb, receptor, ligand))
    return out


def _decomposition(decomp_csv, frame_dt_ns, windows, gb, receptor, ligand) -> dict:
    dec = read_decomp(decomp_csv, frame_dt_ns)
    ident = lambda n: (None, n)  # noqa: E731
    per_window, closure = {}, {}
    for wname, w in windows.items():
        sub = dec[_window_mask(dec.time_ns.to_numpy(), w)]
        res = sub.groupby(["side", "resname", "resnum"]).TOTAL.mean().reset_index()
        mapped = [(receptor or ident)(n) if s == "R" else (ligand or ident)(n)
                  for s, n in zip(res.side, res.resnum)]
        res["partner"] = [m[0] for m in mapped]
        res["residue"] = [m[1] for m in mapped]
        res = res.sort_values("TOTAL").reset_index(drop=True)
        per_window[wname] = res
        total = gb.TOTAL[_window_mask(gb.time_ns.to_numpy(), w)].mean()
        s = float(res.TOTAL.sum())
        closure[wname] = {"sum_residues": s, "total": float(total),
                          "unattributed": float(total - s),
                          "closes": bool(abs(total - s) <= 1.0)}
    return {"decomposition": per_window, "closure": closure,
            "residues_in_decomposition": int(dec.drop_duplicates("label").shape[0])}
