"""
MolDynX Tools (``moldynx``; formerly ``mdforge``) -- a reusable, extensible
framework for analysing GROMACS molecular dynamics simulations of arbitrary
biomolecular systems.

Point it at a simulation directory and it will discover the files, detect the
system composition (protein / ligand / nucleic acid / membrane / ions / ...),
select the appropriate analyses, run them with streaming-friendly performance,
and produce publication-quality figures plus a fully reproducible report.

Public API
----------
    from moldynx import analyze, detect_system, __version__
"""

from __future__ import annotations

__version__ = "0.3.0"
__display_name__ = "MolDynX Tools"
__author__ = "Hossam Mahmoud"

# Re-export the most commonly used entry points (kept import-light).
from moldynx.core.system import detect_system, SystemInfo, SystemType, ComponentType  # noqa: E402
from moldynx.core.registry import registry  # noqa: E402

__all__ = [
    "__version__",
    "detect_system",
    "SystemInfo",
    "SystemType",
    "ComponentType",
    "registry",
    "analyze",
]


def analyze(*args, **kwargs):
    """Convenience wrapper around :func:`moldynx.core.pipeline.run_pipeline`."""
    from moldynx.core.pipeline import run_pipeline
    return run_pipeline(*args, **kwargs)
