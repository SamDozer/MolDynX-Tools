"""
Automatic biomolecular system detection.

Given a topology (and optionally coordinates), classify every residue into a
component type (protein, DNA, RNA, lipid, water, ion, ligand, cofactor, ...),
count chains, and infer the overall :class:`SystemType`.  The result carries
ready-to-use MDAnalysis selection strings so downstream analyses never have to
re-derive selections.

The detector is deliberately conservative and data-driven: anything that is not
a recognised polymer/solvent/ion is treated as a (potentially interesting)
"ligand/other" component and surfaced to the user rather than silently ignored.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis")
import MDAnalysis as mda  # noqa: E402


# --------------------------------------------------------------------------- #
# Residue-name dictionaries (extend freely; matching is case-insensitive)
# --------------------------------------------------------------------------- #
_AMINO = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # protonation / naming variants (CHARMM / AMBER)
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "CYX", "CYM", "LYN", "ASH",
    "GLH", "MSE", "SEC", "PYL", "ACE", "NME", "NMA", "HISE", "HISD",
}
_DNA = {"DA", "DT", "DG", "DC", "DA5", "DA3", "DT5", "DT3", "DG5", "DG3",
        "DC5", "DC3", "DI", "THY", "ADE", "GUA", "CYT"}
_RNA = {"A", "U", "G", "C", "RA", "RU", "RG", "RC", "RA5", "RU5", "RG5", "RC5",
        "RA3", "RU3", "RG3", "RC3", "URA"}
_WATER = {"SOL", "WAT", "HOH", "TIP3", "TIP4", "TIP5", "TIP3P", "SPC", "SPCE",
          "T3P", "T4P", "OPC", "TP3", "TP4"}
_IONS = {"NA", "CL", "K", "POT", "CLA", "SOD", "MG", "CA", "ZN", "FE", "MN",
         "CU", "CO", "NI", "LI", "RB", "CS", "BR", "IOD", "I", "CAL", "MG2",
         "ZN2", "NA+", "CL-", "K+", "CA2", "FE2", "FE3", "F", "CD"}
_LIPIDS = {"POPC", "POPE", "POPS", "POPG", "POPI", "POPA", "DOPC", "DOPE",
           "DPPC", "DMPC", "DSPC", "DLPC", "DPPE", "DPPG", "CHL1", "CHOL",
           "SAPI", "PSM", "DPSM", "LPPC", "PLPC", "PLPE", "OANL", "SDPC",
           "POP2", "PIP2", "CER", "DAG", "TAG"}
# A curated subset of biologically important non-polymer molecules that are
# usually *cofactors* rather than drug-like ligands.
_COFACTORS = {"HEM", "HEME", "HEC", "HEB", "NAD", "NAI", "NAP", "NDP", "NAJ",
              "FAD", "FMN", "FDA", "ATP", "ADP", "AMP", "GTP", "GDP", "GNP",
              "ANP", "SAM", "SAH", "PLP", "TPP", "COA", "COO", "COENZYME",
              "BTN", "MTE", "MGD", "F43", "B12", "CLA", "BCL", "PQQ"}
# Common crystallographic additives / buffer that are not the ligand of interest
_BUFFER_LIKE = {"GOL", "EDO", "PEG", "PG4", "SO4", "PO4", "ACT", "DMS", "MPD",
                "FMT", "TRS", "BME", "IPA", "MES", "EPE", "CIT", "TLA"}


class ComponentType(str, Enum):
    PROTEIN = "protein"
    DNA = "dna"
    RNA = "rna"
    LIPID = "lipid"
    WATER = "water"
    ION = "ion"
    COFACTOR = "cofactor"
    LIGAND = "ligand"
    OTHER = "other"


class SystemType(str, Enum):
    PROTEIN_ONLY = "protein_only"
    PROTEIN_PROTEIN = "protein_protein"
    PROTEIN_PEPTIDE = "protein_peptide"
    PROTEIN_LIGAND = "protein_ligand"
    PROTEIN_DNA = "protein_dna"
    PROTEIN_RNA = "protein_rna"
    PROTEIN_MEMBRANE = "protein_membrane"
    NUCLEIC_ONLY = "nucleic_only"
    MEMBRANE_ONLY = "membrane_only"
    MIXED = "mixed"
    UNKNOWN = "unknown"


# Convenience set: every system type that contains a protein (analyses needing a
# protein declare ``supported_systems = PROTEIN_SYSTEMS``).
PROTEIN_SYSTEMS = {
    SystemType.PROTEIN_ONLY, SystemType.PROTEIN_PROTEIN, SystemType.PROTEIN_PEPTIDE,
    SystemType.PROTEIN_LIGAND, SystemType.PROTEIN_DNA, SystemType.PROTEIN_RNA,
    SystemType.PROTEIN_MEMBRANE, SystemType.MIXED,
}
# Systems with >= 2 macromolecular partners (interface analyses apply).
COMPLEX_SYSTEMS = {
    SystemType.PROTEIN_PROTEIN, SystemType.PROTEIN_PEPTIDE, SystemType.PROTEIN_DNA,
    SystemType.PROTEIN_RNA,
}


def _classify_resname(resname: str) -> ComponentType:
    r = resname.strip().upper()
    if r in _AMINO:
        return ComponentType.PROTEIN
    if r in _DNA:
        return ComponentType.DNA
    if r in _RNA:
        return ComponentType.RNA
    if r in _WATER:
        return ComponentType.WATER
    if r in _IONS:
        return ComponentType.ION
    if r in _LIPIDS:
        return ComponentType.LIPID
    if r in _COFACTORS:
        return ComponentType.COFACTOR
    return ComponentType.OTHER  # candidate ligand / buffer — refined below


@dataclass
class Component:
    """A distinct molecular component and how to select it."""
    ctype: ComponentType
    resnames: list[str]
    n_residues: int
    n_molecules: int
    selection: str


@dataclass
class SystemInfo:
    """Structured description of a detected simulation system."""
    n_atoms: int
    system_type: SystemType
    components: dict[str, Component] = field(default_factory=dict)
    n_protein_chains: int = 0
    protein_chain_lengths: list[int] = field(default_factory=list)
    ligand_resnames: list[str] = field(default_factory=list)
    selections: dict[str, str] = field(default_factory=dict)
    flags: dict[str, bool] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    # -- convenience ---------------------------------------------------- #
    def has(self, ctype: ComponentType | str) -> bool:
        key = ctype.value if isinstance(ctype, ComponentType) else ctype
        return key in self.components and self.components[key].n_residues > 0

    def to_dict(self) -> dict:
        return {
            "n_atoms": self.n_atoms,
            "system_type": self.system_type.value,
            "n_protein_chains": self.n_protein_chains,
            "protein_chain_lengths": self.protein_chain_lengths,
            "ligand_resnames": self.ligand_resnames,
            "components": {
                k: {"type": c.ctype.value, "resnames": c.resnames,
                    "n_residues": c.n_residues, "n_molecules": c.n_molecules}
                for k, c in self.components.items()
            },
            "selections": self.selections,
            "flags": self.flags,
            "notes": self.notes,
        }

    def summary(self) -> str:
        lines = [f"System type : {self.system_type.value}",
                 f"Atoms       : {self.n_atoms:,}",
                 f"Protein     : {self.n_protein_chains} chain(s), "
                 f"lengths {self.protein_chain_lengths}"]
        for key, c in self.components.items():
            if c.n_residues:
                lines.append(f"  {c.ctype.value:9s}: {c.n_residues} residues "
                             f"({', '.join(c.resnames[:6])}"
                             f"{' …' if len(c.resnames) > 6 else ''})")
        return "\n".join(lines)


def _count_protein_chains(u: mda.Universe) -> list[int]:
    """Return the residue count of each protein chain (by segment/fragment)."""
    protein = u.select_atoms("protein")
    if protein.n_atoms == 0:
        # fall back to our amino-acid dictionary if MDAnalysis 'protein' misses
        sel = " or ".join(f"resname {r}" for r in sorted(_AMINO))
        protein = u.select_atoms(sel)
    if protein.n_atoms == 0:
        return []
    lengths = []
    try:
        for seg in protein.segments:
            n = seg.atoms.select_atoms("name CA").n_residues
            if n:
                lengths.append(int(n))
    except Exception:
        pass
    if not lengths:
        # single chain fallback
        lengths = [int(protein.select_atoms("name CA").n_residues)]
    return [n for n in lengths if n > 0]


def detect_system(topology: str | Path, coordinates: str | Path | None = None,
                  peptide_cutoff: int = 30) -> SystemInfo:
    """
    Detect the composition and type of an MD system from its topology.

    Parameters
    ----------
    topology : path to a .tpr/.gro/.pdb/.psf that carries residue names.
    coordinates : optional coordinate file (only needed if the topology lacks it).
    peptide_cutoff : a protein chain shorter than this many residues attached to a
        larger chain marks a protein-peptide system.

    Returns
    -------
    SystemInfo
    """
    topology = str(topology)
    u = mda.Universe(topology, str(coordinates)) if coordinates else mda.Universe(topology)

    # -- classify residues ------------------------------------------------ #
    buckets: dict[ComponentType, dict] = {}
    for res in u.residues:
        ct = _classify_resname(res.resname)
        b = buckets.setdefault(ct, {"resnames": set(), "n_res": 0})
        b["resnames"].add(res.resname.strip().upper())
        b["n_res"] += 1

    # Refine OTHER into ligand vs buffer-like; keep both under LIGAND but note buffers.
    components: dict[str, Component] = {}
    ligand_resnames: list[str] = []
    notes: list[str] = []
    for ct, b in buckets.items():
        resnames = sorted(b["resnames"])
        real_ct = ct
        if ct == ComponentType.OTHER:
            real_ct = ComponentType.LIGAND
            ligs = [r for r in resnames if r not in _BUFFER_LIKE]
            buf = [r for r in resnames if r in _BUFFER_LIKE]
            ligand_resnames.extend(ligs)
            if buf:
                notes.append(f"Ignoring buffer/crystallisation additives: {', '.join(buf)}")
        sel = _selection_for(real_ct, resnames)
        # merge into existing (e.g. OTHER folded into LIGAND)
        if real_ct.value in components:
            existing = components[real_ct.value]
            existing.resnames = sorted(set(existing.resnames) | set(resnames))
            existing.n_residues += b["n_res"]
        else:
            components[real_ct.value] = Component(
                ctype=real_ct, resnames=resnames, n_residues=b["n_res"],
                n_molecules=b["n_res"], selection=sel)

    chain_lengths = _count_protein_chains(u)
    info = SystemInfo(
        n_atoms=int(u.atoms.n_atoms),
        system_type=SystemType.UNKNOWN,
        components=components,
        n_protein_chains=len(chain_lengths),
        protein_chain_lengths=sorted(chain_lengths, reverse=True),
        ligand_resnames=sorted(set(ligand_resnames)),
        notes=notes,
    )
    info.selections = _build_selections(info)
    info.flags = _build_flags(info)
    info.system_type = _classify_system(info, peptide_cutoff)
    return info


def _selection_for(ct: ComponentType, resnames: list[str]) -> str:
    if ct == ComponentType.PROTEIN:
        return "protein"
    if ct in (ComponentType.DNA, ComponentType.RNA):
        return "nucleic"
    joined = " ".join(resnames)
    return f"resname {joined}" if resnames else "not all"


def _build_selections(info: SystemInfo) -> dict[str, str]:
    sel = {}
    if info.has(ComponentType.PROTEIN):
        sel["protein"] = "protein"
        sel["backbone"] = "protein and backbone"
        sel["ca"] = "protein and name CA"
    for ct in (ComponentType.DNA, ComponentType.RNA):
        if info.has(ct):
            sel["nucleic"] = "nucleic"
    for ct in (ComponentType.LIGAND, ComponentType.LIPID, ComponentType.ION,
               ComponentType.WATER, ComponentType.COFACTOR):
        if info.has(ct):
            sel[ct.value] = info.components[ct.value].selection
    return sel


def _build_flags(info: SystemInfo) -> dict[str, bool]:
    return {
        "has_protein": info.has(ComponentType.PROTEIN),
        "has_dna": info.has(ComponentType.DNA),
        "has_rna": info.has(ComponentType.RNA),
        "has_nucleic": info.has(ComponentType.DNA) or info.has(ComponentType.RNA),
        "has_lipid": info.has(ComponentType.LIPID),
        "has_ligand": info.has(ComponentType.LIGAND),
        "has_cofactor": info.has(ComponentType.COFACTOR),
        "has_ions": info.has(ComponentType.ION),
        "has_water": info.has(ComponentType.WATER),
        "multi_chain": info.n_protein_chains >= 2,
    }


def _classify_system(info: SystemInfo, peptide_cutoff: int) -> SystemType:
    f = info.flags
    if not f["has_protein"]:
        if f["has_nucleic"]:
            return SystemType.NUCLEIC_ONLY
        if f["has_lipid"]:
            return SystemType.MEMBRANE_ONLY
        return SystemType.UNKNOWN
    # protein present
    if f["has_lipid"]:
        return SystemType.PROTEIN_MEMBRANE
    if f["has_dna"]:
        return SystemType.PROTEIN_DNA
    if f["has_rna"]:
        return SystemType.PROTEIN_RNA
    if f["has_ligand"]:
        return SystemType.PROTEIN_LIGAND
    if info.n_protein_chains >= 2:
        lengths = sorted(info.protein_chain_lengths, reverse=True)
        if len(lengths) >= 2 and lengths[-1] < peptide_cutoff <= lengths[0]:
            return SystemType.PROTEIN_PEPTIDE
        return SystemType.PROTEIN_PROTEIN
    return SystemType.PROTEIN_ONLY
