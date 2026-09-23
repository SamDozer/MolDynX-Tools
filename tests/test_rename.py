"""The mdforge -> moldynx (MolDynX Tools) rename and its one-version compatibility shim."""

from __future__ import annotations

import importlib
import sys
import warnings


def test_new_package_identity():
    import moldynx

    assert moldynx.__display_name__ == "MolDynX Tools"
    assert moldynx.__version__.startswith("0.3")


def test_legacy_import_warns_and_aliases():
    for name in [m for m in sys.modules if m == "mdforge" or m.startswith("mdforge.")]:
        del sys.modules[name]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mdforge = importlib.import_module("mdforge")
    assert any(issubclass(w.category, DeprecationWarning) and "moldynx" in str(w.message)
               for w in caught)

    import moldynx
    assert mdforge.__version__ == moldynx.__version__


def test_legacy_submodules_are_the_same_objects():
    importlib.import_module("mdforge")
    from mdforge.core.base import BaseAnalysis as Old
    from moldynx.core.base import BaseAnalysis as New
    assert Old is New

    from mdforge.core.registry import registry as old_registry
    from moldynx.core.registry import registry as new_registry
    assert old_registry is new_registry        # one registry, analyses register once
