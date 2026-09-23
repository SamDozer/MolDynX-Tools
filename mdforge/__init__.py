"""
Deprecated import name for MolDynX Tools.

``mdforge`` was renamed to ``moldynx`` in 0.3.0. This shim keeps existing code and
third-party plugins working for one minor version: ``import mdforge`` and any
``mdforge.<submodule>`` import resolve to the same module objects as
``moldynx.<submodule>`` (so analyses register once, in one registry), and a
``DeprecationWarning`` is emitted. It will be removed in 0.4.0.
"""

from __future__ import annotations

import importlib
import importlib.abc
import importlib.util
import sys
import warnings

warnings.warn(
    "the 'mdforge' package has been renamed to 'moldynx' (MolDynX Tools); "
    "'import mdforge' will stop working in 0.4.0",
    DeprecationWarning,
    stacklevel=2,
)

import moldynx as _moldynx  # noqa: E402
from moldynx import *  # noqa: E402,F401,F403
from moldynx import __version__  # noqa: E402,F401


class _AliasFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Resolve ``mdforge.x.y`` to the already-importable ``moldynx.x.y``."""

    prefix = "mdforge."

    def find_spec(self, fullname, path=None, target=None):
        if fullname.startswith(self.prefix):
            real = "moldynx." + fullname[len(self.prefix):]
            if importlib.util.find_spec(real) is None:
                return None
            return importlib.util.spec_from_loader(fullname, self, is_package=True)
        return None

    def create_module(self, spec):
        return importlib.import_module("moldynx." + spec.name[len(self.prefix):])

    def exec_module(self, module):  # the real module is already executed
        pass


if not any(isinstance(f, _AliasFinder) for f in sys.meta_path):
    sys.meta_path.insert(0, _AliasFinder())

__path__ = []  # make this a package so submodule imports reach the finder
