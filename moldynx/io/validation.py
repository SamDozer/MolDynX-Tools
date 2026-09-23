"""
Validation of a discovered :class:`FileSet` and the **capability matrix**.

Produces clear, actionable messages: what is missing, why it is needed, what it
disables, and how to obtain it -- so a researcher is never left guessing, and no
analysis silently runs on the wrong or missing inputs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from moldynx.io.discovery import FileSet

# role -> (why it is needed, how to obtain it)
_GUIDANCE = {
    "trajectory": ("the time-series coordinates every dynamic analysis needs",
                   "produce one with `gmx mdrun` (writes .xtc/.trr) or `gmx trjconv`"),
    "topology": ("atom names, bonds, masses and charges used for selections/analyses",
                 "use the run input `.tpr` (best) or a `.gro`/`.pdb` of the system"),
    "energy": ("thermodynamic observables (T, P, density, energies)",
               "the `.edr` written by `gmx mdrun`; required only for energy analysis"),
    "index": ("custom atom groups",
              "create with `gmx make_ndx`; optional -- MolDynX builds its own groups"),
    "gmx_top": ("the processed topology with #include's",
                "the `topol.top` used by grompp (with its `toppar/`); required for MM-GBSA/PBSA"),
}

# capability -> (required keys (all), description, what happens without it)
CAPABILITIES: list[tuple[str, tuple[str, ...], str, str]] = [
    ("core analyses", ("trajectory", "topology"),
     "every structural / dynamic analysis", "nothing can run"),
    ("energies", ("energy",),
     "thermodynamic stability observables", "energy analysis is skipped"),
    ("methods & completion proof", ("production_log",),
     "Methods table, mdrun sessions, extension audit, proof the run finished",
     "Methods come from the run input only; completion is reported as not proven"),
    ("binding energy (MM-GBSA/PBSA)", ("gmx_top", "toppar"),
     "a protein-only complex system for gmx_MMPBSA",
     "skipped; a ready-to-run package is still written"),
    ("minimisation audit", ("em",), "outcome, Fmax, crash-dump accounting",
     "the equilibration report lists it as missing"),
    ("NVT audit", ("nvt",), "temperature equilibration", "listed as missing"),
    ("NPT audit", ("npt",), "pressure/density equilibration", "listed as missing"),
    ("restraint & stage-parameter audit", ("stage_tprs",),
     "position restraints and exact stage parameters from the run inputs",
     "reported as not recoverable"),
    ("MDP-only settings", ("mdp",), "gen-vel, gen-temp, define",
     "reported as not recorded (gen-vel inferred only from evidence)"),
    ("chain of custody", ("job_scripts",), "grompp -c/-r/-t links between stages",
     "inferred from timestamps and continuation flags, labelled inferred"),
]


@dataclass
class ValidationResult:
    ok: bool
    errors: list[str]
    warnings: list[str]
    capabilities: list[dict] = field(default_factory=list)

    def report(self) -> str:
        lines = [f"  [ERROR] {e}" for e in self.errors]
        lines += [f"  [warn ] {w}" for w in self.warnings]
        return "\n".join(lines) if lines else "  all required inputs present"

    def capability_table(self) -> str:
        rows = [f"  {'capability':36s} {'status':11s} note"]
        for c in self.capabilities:
            rows.append(f"  {c['capability']:36s} {c['status']:11s} "
                        f"{c['enables'] if c['status'] == 'available' else c['without']}")
        return "\n".join(rows)


def capability_matrix(fs: FileSet) -> list[dict]:
    keys = fs.available_keys()
    out = []
    for name, req, enables, without in CAPABILITIES:
        missing = [k for k in req if k not in keys]
        out.append({"capability": name, "status": "available" if not missing else "unavailable",
                    "requires": list(req), "missing": missing, "enables": enables,
                    "without": without})
    return out


def validate_fileset(fs: FileSet, allow_ambiguous: bool = False) -> ValidationResult:
    """Check the minimum inputs and the consistency of the chosen production run."""
    errors, warnings = [], []

    def explain(role: str) -> str:
        why, how = _GUIDANCE.get(role, ("", ""))
        return f"missing '{role}' — needed for {why}. To obtain it: {how}."

    if fs.trajectory is None:
        errors.append(explain("trajectory"))
    if fs.topology is None and fs.structure is None:
        errors.append(explain("topology"))

    ev = fs.evidence
    t = (ev.get("trajectory_candidates") or {}).get(str(fs.trajectory)) if fs.trajectory else None
    h = (ev.get("tpr_candidates") or {}).get(str(fs.topology)) if fs.topology else None
    if t and h and t.get("readable") and h.get("readable") and t["natoms"] != h["natoms"]:
        errors.append(f"atom count mismatch: {fs.trajectory.name} has {t['natoms']} atoms but "
                      f"{fs.topology.name} has {h['natoms']} — they are not from the same run")
    for a in fs.ambiguities:
        (warnings if allow_ambiguous else errors).append(
            f"ambiguous input: {a} — choose explicitly (--traj/--top, or --interactive)")

    if fs.energy is None:
        warnings.append("no production '.edr' found — energy analysis will be skipped.")
    if fs.gmx_top is None:
        warnings.append("no 'topol.top' found — MM-GBSA/PBSA cannot be prepared.")
    if not ev.get("verified", True) and fs.trajectory is not None:
        warnings.append("production run chosen by stage/name only (trajectory headers could not "
                        "be read) — the choice is unverified")
    warnings.extend(fs.warnings)
    missing_refs = ev.get("job_script_inputs_not_found")
    if missing_refs:
        warnings.append("inputs referenced by job scripts are not in this folder: "
                        + ", ".join(sorted(missing_refs)))

    return ValidationResult(ok=not errors, errors=errors, warnings=warnings,
                            capabilities=capability_matrix(fs))
