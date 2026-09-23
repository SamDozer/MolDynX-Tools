"""
BaseAnalysis -- the contract every analysis (built-in or plugin) implements.

Each analysis declares, as class attributes:

* ``name``              -- unique short key (used in the CLI and config).
* ``label``             -- human-readable title.
* ``category``          -- grouping ("structure", "dynamics", "interactions"...).
* ``required_files``    -- semantic input keys it needs ({"trajectory","topology",
                           "energy",...}); the pipeline verifies availability.
* ``supported_systems`` -- set of SystemType it applies to, or {"*"} for all.
* ``outputs``           -- declared output basenames (for documentation/manifest).
* ``default_params``    -- parameters with defaults (overridable via config/CLI).

Subclasses are **auto-registered** on definition, so adding an analysis (or a
plugin) is simply a matter of subclassing ``BaseAnalysis`` and being imported.
"""

from __future__ import annotations

import abc
from typing import TYPE_CHECKING

from moldynx.core.system import SystemType

if TYPE_CHECKING:  # avoid import cycle at runtime
    from moldynx.core.context import AnalysisContext
    from moldynx.core.system import SystemInfo


class BaseAnalysis(abc.ABC):
    # -- declaration (override in subclasses) ----------------------------- #
    name: str = ""
    label: str = ""
    category: str = "general"
    required_files: set[str] = {"trajectory", "topology"}
    supported_systems: set = {"*"}
    outputs: list[str] = []
    default_params: dict = {}
    # Execution order (lower runs first). Analyses that read other analyses'
    # CSVs (convergence, statistics) use a higher value so they run last.
    order: int = 100

    # -- auto-registration ------------------------------------------------ #
    def __init_subclass__(cls, register: bool = True, **kwargs):
        super().__init_subclass__(**kwargs)
        # Only register concrete, named analyses (skip intermediate bases).
        if register and getattr(cls, "name", ""):
            from moldynx.core.registry import registry
            registry.register(cls)

    # -- applicability ---------------------------------------------------- #
    @classmethod
    def is_applicable(cls, system: "SystemInfo") -> bool:
        """Whether this analysis applies to a detected system."""
        if "*" in cls.supported_systems:
            return True
        return system.system_type in cls.supported_systems

    @classmethod
    def missing_files(cls, available: set[str]) -> set[str]:
        """Required semantic file keys that are not available."""
        return set(cls.required_files) - set(available)

    # -- execution -------------------------------------------------------- #
    @abc.abstractmethod
    def run(self, ctx: "AnalysisContext") -> dict:
        """
        Execute the analysis and return a JSON-serialisable summary dict.

        Implementations should read inputs via ``ctx`` (universe, selections,
        params), write CSV/figure outputs under ``ctx.results_dir`` /
        ``ctx.figures_dir``, and return a compact summary used by the report.
        """
        raise NotImplementedError

    # -- helpers ---------------------------------------------------------- #
    def params(self, ctx: "AnalysisContext") -> dict:
        """Merge default_params with any user overrides for this analysis."""
        merged = dict(self.default_params)
        merged.update(ctx.params_for(self.name))
        return merged

    def describe(self) -> dict:
        return {
            "name": self.name,
            "label": self.label,
            "category": self.category,
            "required_files": sorted(self.required_files),
            "supported_systems": ["*"] if "*" in self.supported_systems
            else sorted(s.value if isinstance(s, SystemType) else s
                        for s in self.supported_systems),
            "outputs": list(self.outputs),
            "default_params": dict(self.default_params),
        }
