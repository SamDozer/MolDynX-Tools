"""
Validation of a discovered :class:`FileSet`.

Produces clear, actionable messages: what is missing, why it is needed, and how
to obtain it -- so a researcher is never left guessing.
"""

from __future__ import annotations

from dataclasses import dataclass

from moldynx.io.discovery import FileSet

# role -> (why it is needed, how to obtain it)
_GUIDANCE = {
    "trajectory": ("the time-series coordinates every dynamic analysis needs",
                   "produce one with `gmx mdrun` (writes .xtc/.trr) or `gmx trjconv`"),
    "topology": ("atom names, bonds, masses and charges used for selections/analyses",
                 "use the run input `.tpr` (best) or a `.gro`/`.pdb` of the system"),
    "energy": ("thermodynamic observables (T, P, density, energies)",
               "the `.edr` written by `gmx mdrun`; required only for energy analysis"),
    "index": ("custom atom groups (e.g. for MM/PBSA group definitions)",
              "create with `gmx make_ndx`; optional for most analyses"),
    "gmx_top": ("the processed topology with #include's",
                "the `topol.top` used for the run; required only for MM/PBSA"),
}


@dataclass
class ValidationResult:
    ok: bool
    errors: list[str]
    warnings: list[str]

    def report(self) -> str:
        lines = []
        for e in self.errors:
            lines.append(f"  [ERROR] {e}")
        for w in self.warnings:
            lines.append(f"  [warn ] {w}")
        return "\n".join(lines) if lines else "  all required inputs present"


def validate_fileset(fs: FileSet) -> ValidationResult:
    """Check that the minimum inputs (trajectory + topology) are present."""
    errors, warnings = [], []

    def explain(role: str) -> str:
        why, how = _GUIDANCE.get(role, ("", ""))
        return f"missing '{role}' — needed for {why}. To obtain it: {how}."

    if fs.trajectory is None:
        errors.append(explain("trajectory"))
    if fs.topology is None and fs.structure is None:
        errors.append(explain("topology"))

    # Non-fatal: features that will simply be skipped.
    if fs.energy is None:
        warnings.append("no '.edr' found — energy/stability analysis will be skipped.")
    if fs.gmx_top is None:
        warnings.append("no '.top' found — MM/PBSA workflow generation will be limited.")
    if len(fs.all_trajectories) > 1:
        warnings.append(f"{len(fs.all_trajectories)} trajectories found — "
                        f"defaulting to the largest ({fs.trajectory.name}); "
                        f"use --interactive to choose.")

    return ValidationResult(ok=not errors, errors=errors, warnings=warnings)
