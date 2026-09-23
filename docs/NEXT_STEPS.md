# MolDynX Tools 0.3 — remaining work

Work paused after Phase 5 (part 1) on branch `feature/moldynx-dataset-pipeline`. The goal of
0.3: running MolDynX on any GROMACS/CHARMM-GUI simulation folder always yields the same audited,
reproducible dataset. The quality baseline is two hand-made analyses of protein–protein
complexes (α-zein A8HNE1–ZmBiP2, α-zein Q946V6–ZmBiP2). They are **test cases only**: no code
may name a system, chain, residue range or folder. The full specification is in the handoff
package `MolDynX_handoff_2026-09-23/PROMPT.md` (kept outside the repo); the essentials are below.

## Done

| Phase | Result |
|---|---|
| 1 Rename | mdforge → MolDynX Tools (`moldynx`), `mdforge` shim for one version |
| 2 Intake | evidence-based discovery, capability matrix, `moldynx intake`, INTAKE_REPORT |
| 3 PBC | one-pass diagnose → treat → prove (`core/pbc.py`), chain identity, cache key, `pbc_validation` |
| 4 Equilibration | `equilibration` analysis: stage parameters, minimisation outcome + Fmax location, crash dumps, energy statistics, position restraints (`gmx dump`), protonation, chain of custody, timeline, figure |
| 5 (part) | `interface` rewritten (full core lists, residues ever within 6 Å, buried area, iRMSD, trends); `core/surface.py` works around an mdtraj multi-frame SASA bug |

Phases 2–4 reproduce the manual audit of both reference datasets exactly (see CHANGELOG).

## Remaining

### Phase 5 — finish complex analyses
- Run `interface` end to end on a real 100 ns complex; compare with the legacy values
  (A8HNE1: 27.9 nm² buried, mean minimum distance 0.266 nm, 34 persistent contacts).
- Port from the legacy pipeline (`alpha-zein-md-analysis/scripts/`), generalised, as
  `BaseAnalysis` classes: `contact_lifetime`, `dynamic_network`, `water_bridges`, `porcupine`,
  eigenvalue spectrum into `pca`. Do **not** port `report_advanced.py` or legacy `mmgbsa.py`.
- SASA speed: single-frame calls are correct but serial (~4 s/frame for 13 k atoms);
  `sasa` and `interface` stride by 10. Parallelise with a process pool, then restore stride 1.

### Phase 6 — annotation and analysis window
- `annotations:` in the YAML config (template written by `moldynx intake --detect`):
  chain display names and roles, biological numbering offsets, homology reference, motifs file,
  docking-site residues, unresolved metadata (force-field variant, salt, protonation rationale —
  never guessed).
- Domains by homology only (Needleman–Wunsch, BLOSUM62, gap −10/−0.5; CATH edges of a reference
  structure through SIFTS; map author numbering to UniProt by alignment; DSSP cross-check; flag
  uncertain edges). Motifs only from a user-supplied file; report partial matches honestly.
- Window selection: Chodera t₀ per binding-relevant observable, slope + Welch half-vs-half;
  if nothing is stationary, report time-resolved values and ask for the primary window.

### Phase 7 — MM-GBSA / MM-PBSA (rewrite `analysis/mmpbsa.py`)
- prepare → probe → run → analyse. Protein-only `complex.tpr` (convert-tpr on the Protein group),
  `complex.top` with protein molecule types only, `complex.ndx` from chain identity
  (0 = receptor, 1 = ligand); **no Amber `forcefields=` line**.
- Probe gate: 11 frames GB; abort unless every bonded Δ term is exactly 0.
- GB: igb 5, PBRadii 3, idecomp 2, IE + C2 entropy (reported only if σ < 3.6 / 6.0 kcal/mol).
  PB: PBRadii 7, inp 1, fillratio 4. Same frames (1 ns spacing). Ionic strength **derived** from
  ion counts (`SystemInfo.ion_counts`) and box volume, or labelled assumed.
- Decomposition: explicit `print_res` = `interface` residues within 6 Å ever
  (`residues_within_6A_ever`), closure check Σ residues vs ΔG_bind (Q946V6 legacy gap +3.73
  kcal/mol at 80–100 ns should close).
- MPI ranks ≤ physical cores; PB ≈ 3.3 GB per rank. Analysis: autocorrelation-corrected SEM,
  windows, drift, GB vs PB, hotspots, region sums. Test data: `tests/fixtures/gmx_mmpbsa/`.

### Phase 8 — documents
- Data-driven README, AUDIT_PROVENANCE, PBC_VALIDATION, EQUILIBRATION, BINDING_ENERGY,
  analysis report; Observations / Interpretation / Limitations / Conclusions; no hard-coded
  system names or conclusions; never call a system "stable" from an RMSD plateau (current
  `report/generator.py` `_interpretation` does — fix); render Markdown properly to
  self-contained HTML; validator (pair `**` across the whole document).

### Phase 9 — dataset, verify, package
- `moldynx dataset` (numbered layout 0_raw_production … 7_binding_energy with SHA-256 manifests),
  `moldynx verify` (layout, references, manifests, Rg smoke test), `moldynx package`
  (analysis-only zip when raw data > ~1 GB; stub injected into the zip only; extract + verify).

### Phase 10 — wrap-up
- README rewrite (input-file contract, annotations schema, new commands), QUICKSTART,
  CITATION/.zenodo titles (done), ROADMAP ticks, fix README claim about a `legacy/` folder.
- The user renames the GitHub repository and creates the release (Zenodo DOI); not automated.

## Acceptance values (manual audit, for the desktop data)

A8HNE1 — MM-GBSA 0–100 ns −55.2 ± 4.0, 80–100 −69.8 ± 3.7; MM-PBSA −57.2 ± 7.4 / −76.6 ± 4.6
kcal/mol; r(GB,PB) 0.86. Q946V6 — MM-GBSA −33.4 ± 2.1 / −37.8 ± 2.4; MM-PBSA −40.7 ± 2.5 /
−47.1 ± 3.3; r 0.73; entropy σ 86.7 → invalid.
