"""
Analysis registry + plugin discovery.

The global ``registry`` collects every :class:`BaseAnalysis` subclass (built-in
or plugin).  It supports:

* automatic selection of analyses for a detected system,
* discovery of plugins from (a) the built-in ``moldynx.analysis.plugins`` package,
  (b) any user directory passed via ``--plugin-dir`` / config, and
  (c) installed packages advertising a ``moldynx.plugins`` entry-point.

Adding an analysis requires no change to the core: define a ``BaseAnalysis``
subclass and make sure its module is imported.
"""

from __future__ import annotations

import importlib
import importlib.util
import pkgutil
import sys
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from moldynx.core.base import BaseAnalysis
    from moldynx.core.system import SystemInfo


class AnalysisRegistry:
    def __init__(self) -> None:
        self._analyses: dict[str, type] = {}

    # -- registration ----------------------------------------------------- #
    def register(self, cls: type) -> type:
        key = getattr(cls, "name", "")
        if not key:
            raise ValueError(f"Analysis {cls!r} has no 'name'.")
        self._analyses[key] = cls
        return cls

    def get(self, name: str) -> type:
        if name not in self._analyses:
            raise KeyError(f"Unknown analysis '{name}'. "
                           f"Available: {', '.join(sorted(self._analyses))}")
        return self._analyses[name]

    def all(self) -> dict[str, type]:
        return dict(self._analyses)

    def names(self) -> list[str]:
        return sorted(self._analyses)

    # -- selection -------------------------------------------------------- #
    def select(self, system: "SystemInfo", available_files: set[str],
               requested: list[str] | None = None,
               run_all: bool = False) -> tuple[list[type], list[tuple[str, str]]]:
        """
        Choose analyses for a system.

        Returns ``(selected, skipped)`` where ``skipped`` is a list of
        ``(name, reason)`` explaining why an analysis was not selected -- this
        powers the ``--plan`` dry-run.
        """
        self.ensure_builtins_loaded()
        selected, skipped = [], []
        for name in (requested or sorted(self._analyses)):
            cls = self._analyses.get(name)
            if cls is None:
                skipped.append((name, "not registered"))
                continue
            if not run_all and requested is None and not cls.is_applicable(system):
                skipped.append((name, f"not applicable to {system.system_type.value}"))
                continue
            missing = cls.missing_files(available_files)
            if missing:
                skipped.append((name, f"missing required files: {', '.join(sorted(missing))}"))
                continue
            selected.append(cls)
        # Run in declared order (convergence/statistics last), then alphabetically.
        selected.sort(key=lambda c: (getattr(c, "order", 100), c.name))
        return selected, skipped

    # -- plugin / builtin loading ---------------------------------------- #
    _builtins_loaded = False

    def ensure_builtins_loaded(self) -> None:
        if self._builtins_loaded:
            return
        import moldynx.analysis as analysis_pkg
        for mod in pkgutil.iter_modules(analysis_pkg.__path__):
            if mod.name.startswith("_") or mod.name == "plugins":
                continue
            importlib.import_module(f"moldynx.analysis.{mod.name}")
        self.discover_plugins()  # built-in plugins package
        self._builtins_loaded = True

    def discover_plugins(self, plugin_dirs: list[str | Path] | None = None) -> list[str]:
        """Import plugin modules from the plugins package, dirs, and entry-points."""
        loaded: list[str] = []
        # (a) built-in plugins package
        try:
            import moldynx.analysis.plugins as plug_pkg
            for mod in pkgutil.iter_modules(plug_pkg.__path__):
                if not mod.name.startswith("_"):
                    importlib.import_module(f"moldynx.analysis.plugins.{mod.name}")
                    loaded.append(mod.name)
        except Exception:
            pass
        # (b) user directories
        for d in (plugin_dirs or []):
            d = Path(d)
            if not d.is_dir():
                continue
            for pyfile in sorted(d.glob("*.py")):
                if pyfile.name.startswith("_"):
                    continue
                spec = importlib.util.spec_from_file_location(
                    f"moldynx_plugin_{pyfile.stem}", pyfile)
                if spec and spec.loader:
                    module = importlib.util.module_from_spec(spec)
                    sys.modules[spec.name] = module
                    spec.loader.exec_module(module)
                    loaded.append(pyfile.stem)
        # (c) installed entry-points; "mdforge.plugins" is the pre-rename group name
        try:
            from importlib.metadata import entry_points
            eps = entry_points()
            for group_name in ("moldynx.plugins", "mdforge.plugins"):
                group = eps.select(group=group_name) if hasattr(eps, "select") \
                    else eps.get(group_name, [])
                for ep in group:
                    if group_name == "mdforge.plugins":
                        warnings.warn(
                            f"plugin {ep.name!r} uses the deprecated 'mdforge.plugins' "
                            f"entry-point group; advertise it under 'moldynx.plugins'",
                            DeprecationWarning, stacklevel=2)
                    ep.load()
                    loaded.append(ep.name)
        except Exception:
            pass
        return loaded


# module-level singleton
registry = AnalysisRegistry()
