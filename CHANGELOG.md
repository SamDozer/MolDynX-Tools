# Changelog

## 0.3.0 — unreleased

### Renamed: mdforge → **MolDynX Tools**

- Import package `mdforge` → `moldynx`; CLI `mdforge` → `moldynx`; distribution
  `molecular-dynamics-forge` → `moldynx-tools`.
- Default output folder `mdforge_results` → `moldynx_results`; manifest key
  `mdforge_version` → `moldynx_version`; plugin entry-point group `mdforge.plugins` →
  `moldynx.plugins`.
- **Compatibility for one minor version (removed in 0.4.0):** `import mdforge` and every
  `mdforge.<submodule>` import still work (same module objects, one registry) with a
  `DeprecationWarning`; the `mdforge` console script is kept as an alias; plugins advertised
  under `mdforge.plugins` are still loaded, with a warning.
- The Zenodo concept DOI (10.5281/zenodo.21265946) is unchanged; citation metadata now reads
  "MolDynX Tools (formerly mdforge)".

### Intake: the right files, chosen by evidence

- **Fixed:** discovery picked the *largest* file per role. On real CHARMM-GUI folders that selects
  an equilibration run input as the topology (equilibration `.tpr` files can be larger than the
  production one) and a minimisation crash dump (`stepNc.pdb`) as the structure.
- Files are now classified by simulation **stage** (setup / minimisation / NVT / NPT /
  production); folders of earlier analysis output, GROMACS backups and caches are ignored; the
  production trajectory is the longest complete one whose atom count equals the run input's,
  verified from file headers; ties and mismatches are reported as ambiguities (errors unless
  `--allow-ambiguous`, `--traj/--top` or `--interactive`).
- New read-only evidence readers (`moldynx.io.gromacs`): mdrun logs (sessions, parameters,
  coupling groups, minimisation outcome, warnings), `.mdp`, XTC frame headers (no offset cache is
  written into the simulation folder), `.tpr` headers, and a GROMACS runner that falls back to WSL.
- New `moldynx intake` command and `INTAKE_REPORT.md` / `intake_manifest.json` (also written by
  every `analyze` run): canonical run and why, per-stage summary, temperature changes between
  stages, a **capability matrix** (what each missing file disables), cluster-clock offset, run
  extensions, and inputs referenced by job scripts that are absent.
- `SystemInfo` now carries per-chain records (segid, atom range, residue range, sequence) and ion
  counts.
- The run manifest records the evidence chain and fingerprints every file the run depends on.
- `analyze` output now defaults to `./moldynx_results/<input folder name>` (was
  `./moldynx_results`), so runs of different simulations do not overwrite each other.

### Periodic boundaries: diagnose, treat, prove — in one pass

- **Fixed:** solute extraction applied `unwrap` inside `try/except: pass`, so a topology without
  bonds silently produced a trajectory that was never made whole; the cached solute trajectory was
  reused whenever the files existed, whatever the inputs, frame slice or treatment.
- New `moldynx.core.pbc`: for every frame the raw coordinates are diagnosed first (molecules split
  across the boundary, centre-of-mass continuity, PBC-aware minimum distance between partners),
  then treated (`--pbc auto|none|whole|nojump`; `auto` = nojump for multi-molecule solutes, with
  first-frame clustering), then **proven** against the raw frame: only whole-box translations,
  all bonds short, partner distances preserved, and a split molecule shows two translation vectors
  exactly in the frames where it was split. Nothing is repaired silently.
- New `pbc_validation` analysis (runs first) with `results/pbc_summary.json`,
  `results/pbc_per_frame.csv` and a figure.
- The solute cache (`data/core_meta.json`) is keyed on input fingerprints, PBC mode, frame slice
  and selection.
- **Fixed:** `interface` could not resolve the two partners of CHARMM-GUI complexes: the PDB
  format truncates `seg_0_PROA` / `seg_1_PROB` to `seg_`, so both chains read back from the cached
  structure as one segment. Chain identity is now persisted as atom-index ranges
  (`AnalysisContext.chain_groups()`); `interface` returns an explicit `skipped` reason instead of a
  silent note when partners cannot be resolved.
- Validated on two 100 ns protein–protein trajectories (1.6 M and 2.0 M atoms): split-frame
  counts, zero whole-box translations, minimum inter-chain distances and molecular extents
  reproduce an independent manual audit exactly.

### Preparation and equilibration audit

- New `equilibration` analysis: per-stage parameters from the logs, minimisation outcome (states
  plainly when the force tolerance was *not* reached) and the chain/residue/atom carrying the
  largest residual force, crash-dump accounting, energy-file statistics (tail means, residual
  drift, settling times), position restraints per force constant and molecule type (`gmx dump`,
  native or WSL; production must have none), protonation states, chain of custody from job
  scripts (`moldynx.io.jobscripts`), stage timeline, figure. Missing inputs are reported.
- Reproduces the manual audit of both reference datasets and additionally located the largest
  residual force of one system on a protonated aspartate.

### Interface and surface area

- `interface` rewritten from the validated legacy suite: residue–residue contacts, interface
  occupancy with **full** core lists, persistent contacts, residues ever within 6 Å (the set a
  binding-energy decomposition must cover), buried area, interface Cα RMSD, trends.
- **Fixed:** `mdtraj.shrake_rupley` (1.11.1) returns wrong values for some frames of a
  multi-frame call (±0.5 nm² on a real trajectory; negative buried areas in a rigid-body test).
  `moldynx.core.surface.shrake_rupley` computes frame by frame; `sasa` and `interface` buried
  area now stride by 10 frames by default (parameters `stride` / `bsa_stride`).

### Complex analyses, annotations, stationarity

- `contact_lifetime`, `water_bridges`, `porcupine` (ported from a validated legacy pipeline,
  generalised; porcupine uses an SVD and its own aligned copy of the trajectory).
- `annotation`: biological numbering, domains transferred by global alignment (BLOSUM62,
  gap −10/−0.5) with uncertain edges flagged, motif location with partial matches, docking-site
  involvement — only from the `annotations:` block of the configuration.
- `analysis_window`: Chodera equilibration detection, drift and half-vs-half tests; the
  averaging window is a user decision (`binding_energy.primary_window`) when nothing is
  stationary.
- `moldynx.statistics`: statistical inefficiency, corrected SEM, equilibration detection, drift.

### MM-GBSA / MM-PBSA with gmx_MMPBSA

- **Replaced** the placeholder `mmpbsa` (fixed groups `r 1-100`/`r 101-9999`, an Amber
  force-field line in a CHARMM/GROMACS workflow, the full solvated system, never run; the
  contact-count proxy is dropped — interface occupancy covers it) with a complete workflow built
  on [gmx_MMPBSA](https://github.com/Valdes-Tresanca-MS/gmx_MMPBSA): protein-only system,
  groups from chain identity, derived ionic strength, explicit decomposition residues, probe gate,
  GB + PB, analysis with autocorrelation-corrected errors, drift, GB vs PB, hotspots, closure and
  an entropy validity gate. New `moldynx binding-energy [--execute]`. Reproduces a manual
  analysis of real outputs exactly.

### Documents and datasets

- `PBC_VALIDATION`, `EQUILIBRATION`, `BINDING_ENERGY` documents written from the results files;
  Markdown rendered to self-contained HTML (the old renderer only deleted `**`); report wording
  no longer claims stability from an RMSD plateau.
- `moldynx dataset` / `verify` / `package`: numbered dataset, raw-file manifest,
  `verify_dataset.py` shipped inside, zip extracted and re-verified.
- New dependency: `markdown`; optional `biopython` (`[annotation]`).

See `docs/WHATS_NEW_0.3.md` for the comparison with mdforge 0.2 and `docs/NEXT_STEPS.md` for
what is left.
