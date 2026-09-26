"""
Regions of the system from user annotations: domains (transferred by homology),
motifs, docking-site residues -- and how much of each region is at the interface.

Runs only when the configuration has an ``annotations:`` block; otherwise it says so.
All numbers are reported in biological numbering with the trajectory offset stated.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from moldynx.core import annotations as ann
from moldynx.core.base import BaseAnalysis


def local_to_bio(chain: ann.ChainAnnotation, position: int) -> int:
    """1-based position in the chain's sequence -> biological residue number."""
    return chain.bio(chain.resid_first + int(position) - 1)


class Annotation(BaseAnalysis):
    name = "annotation"
    label = "Domains, motifs and docking site (user annotations)"
    category = "annotation"
    required_files = {"trajectory", "topology"}
    supported_systems = {"*"}
    order = 150          # after interface (40) and rmsf
    outputs = ["results/annotation.json", "results/region_summary.csv"]

    def run(self, ctx) -> dict:
        a = ctx.config.annotations or {}
        chains = ctx.core_meta.get("chains", [])
        ca = ann.chain_annotations(chains, a)
        doc: dict = {"chains": {s: {"display": c.display, "role": c.role,
                                    "trajectory_resids": [c.resid_first, c.resid_last],
                                    "biological_range": list(c.bio_range),
                                    "numbering_offset": c.offset} for s, c in ca.items()},
                     "domains": {}, "motifs": [], "docking_site": {},
                     "unresolved_metadata": a.get("unresolved_metadata", {})}
        regions = []                                        # (segid, kind, name, bio_start, bio_end)
        for seg, spec in (a.get("domains") or {}).items():
            if seg not in ca:
                doc["domains"][seg] = {"error": f"chain {seg} not in this system"}
                continue
            ref = spec.get("reference_sequence") or (
                ann.read_fasta_sequence(Path(spec["reference_fasta"]))
                if spec.get("reference_fasta") else "")
            if not ref:
                doc["domains"][seg] = {"error": "no reference_sequence / reference_fasta"}
                continue
            try:
                mapped = ann.transfer_regions(ca[seg].sequence, ref, spec.get("regions", {}))
            except ImportError:
                doc["domains"][seg] = {"error": "biopython is required for domain transfer"}
                continue
            for m in mapped:                              # sequence position -> biological
                for k in ("start", "end"):
                    if m[k]:
                        m[k] = local_to_bio(ca[seg], m[k])
            doc["domains"][seg] = {"reference": spec.get("reference_name", "reference"),
                                   "regions": mapped}
            regions += [(seg, "domain", m["region"], m["start"], m["end"]) for m in mapped
                        if m["start"] and m["end"]]
        for m in a.get("motifs") or []:
            segs = [m["chain"]] if m.get("chain") else list(ca)
            for seg in segs:
                if seg not in ca:
                    continue
                try:
                    loc = ann.locate_motif(ca[seg].sequence, m["sequence"])
                except ImportError:
                    loc = {"match": "exact-only search (biopython missing)"}
                for k in ("start", "end"):
                    if loc.get(k):
                        loc[k] = local_to_bio(ca[seg], loc[k])
                doc["motifs"].append({"name": m["name"], "chain": seg, **loc})
                if loc.get("start"):
                    regions.append((seg, "motif", m["name"], loc["start"], loc["end"]))
        for seg, residues in (a.get("docking_site") or {}).items():
            doc["docking_site"][seg] = sorted(int(r) for r in residues)

        # interface involvement per region (from the interface analysis, if it ran)
        rows = []
        res_csv = ctx.csv_path("interface_residues.csv")
        iface = pd.read_csv(res_csv) if res_csv.exists() else None
        for seg, kind, name, b0, b1 in regions:
            row = {"chain": seg, "display": ca[seg].display, "kind": kind, "region": name,
                   "bio_start": b0, "bio_end": b1, "n_residues": b1 - b0 + 1}
            if iface is not None:
                sub = iface[(iface.partner == seg) &
                            iface.resid.between(ca[seg].md(b0), ca[seg].md(b1))]
                row.update({"interface_residues": int((sub.interface_occupancy > 0).sum()),
                            "core_residues": int((sub.interface_occupancy >= 0.5).sum()),
                            "occupancy_sum": float(sub.interface_occupancy.sum())})
            rows.append(row)
        for seg, residues in doc["docking_site"].items():
            if iface is not None and seg in ca:
                md_res = [ca[seg].md(r) for r in residues]
                sub = iface[(iface.partner == seg) & iface.resid.isin(md_res)]
                rows.append({"chain": seg, "display": ca[seg].display, "kind": "docking site",
                             "region": "docking site", "n_residues": len(residues),
                             "interface_residues": int((sub.interface_occupancy > 0).sum()),
                             "core_residues": int((sub.interface_occupancy >= 0.5).sum()),
                             "occupancy_sum": float(sub.interface_occupancy.sum())})
        pd.DataFrame(rows).to_csv(ctx.csv_path("region_summary.csv"), index=False)
        ctx.csv_path("annotation.json").write_text(json.dumps(doc, indent=2, default=str),
                                                    encoding="utf-8")
        return {"annotated": bool(a), "chains": doc["chains"],
                "n_domains": sum(len(d.get("regions", [])) for d in doc["domains"].values()),
                "uncertain_domain_edges": [f"{s}:{m['region']}" for s, d in doc["domains"].items()
                                           for m in d.get("regions", []) if m["uncertain"]],
                "motifs": [{k: m.get(k) for k in ("name", "chain", "match", "start", "end")}
                           for m in doc["motifs"]],
                "regions": rows,
                "note": None if a else "no annotations in the configuration: display names "
                                       "default to segids and no regions are defined"}
