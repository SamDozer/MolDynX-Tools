"""Unit tests for the MolDynX Tools core (no large trajectory data required)."""

from __future__ import annotations

import numpy as np

from moldynx.core.system import (
    _classify_resname, _build_flags, _classify_system,
    ComponentType, SystemType, SystemInfo, Component,
)
from moldynx.core.registry import registry
from moldynx import statistics as st


# --------------------------------------------------------------------------- #
# System detection logic
# --------------------------------------------------------------------------- #
def test_classify_resname():
    assert _classify_resname("ALA") is ComponentType.PROTEIN
    assert _classify_resname("HSD") is ComponentType.PROTEIN     # CHARMM variant
    assert _classify_resname("DA") is ComponentType.DNA
    assert _classify_resname("U") is ComponentType.RNA
    assert _classify_resname("TIP3") is ComponentType.WATER
    assert _classify_resname("POT") is ComponentType.ION
    assert _classify_resname("POPC") is ComponentType.LIPID
    assert _classify_resname("HEM") is ComponentType.COFACTOR
    assert _classify_resname("LIG") is ComponentType.OTHER       # -> ligand


def _make_info(**components) -> SystemInfo:
    comps = {}
    for ctype, n in components.items():
        ct = ComponentType(ctype)
        comps[ct.value] = Component(ct, [ctype.upper()], n, n, "sel")
    info = SystemInfo(n_atoms=1000, system_type=SystemType.UNKNOWN, components=comps)
    return info


def test_system_type_classification():
    def classify(chains=1, **comp):
        info = _make_info(**comp)
        info.n_protein_chains = chains
        info.protein_chain_lengths = [187] * chains
        info.flags = _build_flags(info)
        return _classify_system(info, peptide_cutoff=30)

    assert classify(chains=1, protein=187) is SystemType.PROTEIN_ONLY
    assert classify(chains=2, protein=374) is SystemType.PROTEIN_PROTEIN
    assert classify(chains=1, protein=187, ligand=1) is SystemType.PROTEIN_LIGAND
    assert classify(chains=1, protein=187, dna=20) is SystemType.PROTEIN_DNA
    assert classify(chains=1, protein=187, rna=20) is SystemType.PROTEIN_RNA
    assert classify(chains=1, protein=187, lipid=200) is SystemType.PROTEIN_MEMBRANE


# --------------------------------------------------------------------------- #
# Registry / plugin discovery
# --------------------------------------------------------------------------- #
def test_registry_loads_builtins_and_plugin():
    registry.ensure_builtins_loaded()
    names = registry.names()
    for expected in ("rmsd", "rmsf", "rog", "sasa"):
        assert expected in names
    # the bundled example plugin must be auto-discovered
    assert "end_to_end" in names


def test_registry_selection_respects_system():
    registry.ensure_builtins_loaded()
    info = _make_info(protein=187, water=1000)
    info.n_protein_chains = 1
    info.flags = _build_flags(info)
    info.system_type = SystemType.PROTEIN_ONLY
    selected, _ = registry.select(info, available_files={"trajectory", "topology"})
    assert any(c.name == "rmsd" for c in selected)


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
def test_describe_and_ci():
    x = np.arange(100.0)
    d = st.describe(x, "x")
    assert d["n"] == 100
    assert abs(d["mean"] - 49.5) < 1e-9
    assert d["ci95_low"] < d["mean"] < d["ci95_high"]


def test_moving_average_and_plateau():
    x = np.concatenate([np.linspace(0, 1, 50), np.ones(50)])
    ma = st.moving_average(x, 10)
    assert len(ma) == len(x)
    conv = st.plateau_detection(np.arange(len(x)), x)
    assert conv["converged"] is True     # flat tail


def test_bootstrap_ci():
    rng = np.random.default_rng(0)
    x = rng.normal(5.0, 1.0, 500)
    b = st.bootstrap_ci(x, seed=0)
    assert b["ci_low"] < b["estimate"] < b["ci_high"]
    assert abs(b["estimate"] - 5.0) < 0.3
