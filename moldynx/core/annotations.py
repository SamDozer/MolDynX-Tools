"""
User-supplied annotations: display names, biological numbering, domains by homology,
motifs, docking-site residues. Nothing here is invented -- every region comes from
the ``annotations:`` block of the run configuration.

Schema (YAML)::

    annotations:
      chains:
        - {segid: seg_0_PROA, display: "α-zein Q946V6", role: ligand}
        - {segid: seg_1_PROB, display: "ZmBiP2", role: receptor, numbering_offset: 213}
      domains:                      # transferred by alignment, never by copying numbers
        seg_1_PROB:
          reference_name: "UniProt P11021 (human BiP)"
          reference_sequence: "MKLSLVAAMLLLLSAARA..."   # or reference_fasta: path
          regions: {NBD: [26, 405], SBDbeta: [418, 507]}   # reference numbering
      motifs:
        - {name: "Motif 1", sequence: "CSQAPIASLLPPYLSPAVSSVC", chain: seg_0_PROA}
      docking_site: {seg_1_PROB: [405, 434, 435, 438]}      # biological numbering
      unresolved_metadata: {force_field_variant: null, salt_concentration_M: null}

Biological numbering = trajectory resid − ``numbering_offset``; the default offset
makes each chain start at 1.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ChainAnnotation:
    segid: str
    display: str
    role: str | None
    resid_first: int
    resid_last: int
    offset: int
    sequence: str

    def bio(self, resid: int) -> int:
        return int(resid) - self.offset

    def md(self, bio: int) -> int:
        return int(bio) + self.offset

    @property
    def bio_range(self) -> tuple[int, int]:
        return self.bio(self.resid_first), self.bio(self.resid_last)


def chain_annotations(chains: list[dict], annotations: dict) -> dict[str, ChainAnnotation]:
    """Merge persisted chain records with the user's chain annotations."""
    user = {c.get("segid"): c for c in (annotations or {}).get("chains", []) if c.get("segid")}
    out = {}
    for rec in chains:
        seg = rec.get("segid") or f"chain{rec.get('index')}"
        u = user.get(seg, {})
        out[seg] = ChainAnnotation(
            segid=seg, display=u.get("display") or seg, role=u.get("role"),
            resid_first=int(rec["resid_first"]), resid_last=int(rec["resid_last"]),
            offset=int(u.get("numbering_offset", int(rec["resid_first"]) - 1)),
            sequence=rec.get("sequence", ""))
    return out


# --------------------------------------------------------------------------- #
# sequence alignment
# --------------------------------------------------------------------------- #
def _aligner(mode: str = "global"):
    from Bio import Align
    from Bio.Align import substitution_matrices
    a = Align.PairwiseAligner()
    a.mode = mode
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -10.0, -0.5
    return a


def residue_map(query: str, reference: str) -> tuple[dict[int, int], float]:
    """1-based reference position -> 1-based query position, and % identity over aligned pairs."""
    aln = _aligner("global").align(query, reference)[0]
    ref_to_q: dict[int, int] = {}
    same = pairs = 0
    for (qs, qe), (rs, re_) in zip(*aln.aligned):
        for k in range(qe - qs):
            ref_to_q[rs + k + 1] = qs + k + 1
            pairs += 1
            same += query[qs + k] == reference[rs + k]
    return ref_to_q, (100.0 * same / pairs if pairs else 0.0)


def transfer_regions(query: str, reference: str, regions: dict[str, list[int]],
                     search: int = 15) -> list[dict]:
    """
    Map region edges from reference to query numbering through a global alignment
    (BLOSUM62, gap −10/−0.5). An edge that falls in an alignment gap is moved to
    the nearest aligned reference position and flagged ``uncertain`` with its shift.
    """
    ref_to_q, ident = residue_map(query, reference)
    out = []
    for name, (start, end) in regions.items():
        edges, notes = [], []
        for pos, direction in ((int(start), 1), (int(end), -1)):
            if pos in ref_to_q:
                edges.append(ref_to_q[pos])
                continue
            found = None
            for d in range(1, search + 1):
                for cand in (pos + direction * d, pos - direction * d):
                    if cand in ref_to_q:
                        found = (cand, d)
                        break
                if found:
                    break
            if found:
                edges.append(ref_to_q[found[0]])
                notes.append(f"reference {pos} is in an alignment gap; used {found[0]} (±{found[1]})")
            else:
                edges.append(None)
                notes.append(f"reference {pos} has no aligned residue within ±{search}")
        out.append({"region": name, "reference": [int(start), int(end)],
                    "start": edges[0], "end": edges[1], "uncertain": bool(notes),
                    "notes": notes, "identity_pct": round(ident, 1)})
    return out


def locate_motif(sequence: str, motif: str) -> dict:
    """Exact match first; else the best local alignment, with partial-match bookkeeping."""
    i = sequence.find(motif)
    if i >= 0:
        return {"start": i + 1, "end": i + len(motif), "match": "exact", "identity_pct": 100.0,
                "covered": len(motif), "length": len(motif)}
    try:
        aln = _aligner("local").align(sequence, motif)[0]
    except Exception:
        return {"match": "none", "length": len(motif)}
    blocks = list(zip(*aln.aligned))
    if not blocks:
        return {"match": "none", "length": len(motif)}
    s0, s1 = blocks[0][0][0], blocks[-1][0][1]
    m0, m1 = blocks[0][1][0], blocks[-1][1][1]
    same = sum(sequence[qs + k] == motif[ms + k]
               for (qs, qe), (ms, _me) in blocks for k in range(qe - qs))
    covered = sum(qe - qs for (qs, qe), _ in blocks)
    return {"start": s0 + 1, "end": s1, "match": "partial" if covered < len(motif) else "similar",
            "identity_pct": round(100.0 * same / max(covered, 1), 1),
            "motif_positions_covered": [m0 + 1, m1], "covered": covered, "length": len(motif)}


def read_fasta_sequence(path: str | Path) -> str:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    return "".join(x.strip() for x in lines if x and not x.startswith(">"))
