# MolDynX Tools 0.3 — what changed compared with mdforge 0.2

mdforge 0.2 was a good general framework: plugin registry, automatic system detection and
analysis selection, a run manifest, figures and a report. MolDynX Tools 0.3 keeps all of that
and adds what was missing to *trust* the results: evidence for every input choice, proof for
every trajectory transformation, an audit of how the system was prepared, binding energies
that are actually computed, and documents in which every number traces back to a file.

The quality bar was two hand-made, fully audited analyses of 100 ns protein–protein complexes.
Where MolDynX now automates a step of those analyses, the tests check it reproduces them.

## 1. Choosing the input files

| mdforge 0.2 | MolDynX Tools 0.3 |
|---|---|
| Picked the **largest** file per role. On a real CHARMM-GUI folder this took an *equilibration* run input as the topology (it was larger than the production one) and a minimisation crash dump (`stepNc.pdb`) as the structure. | Files are classified by simulation stage. The production trajectory is the longest complete one whose atom count equals the run input's, read from file headers. Ties are reported as ambiguities, never resolved silently. |
| Scanned every subfolder, including earlier analysis outputs. | Ignores earlier analysis folders, GROMACS backups, caches and crash dumps, and says so. |
| Presence checks only. | Capability matrix: each missing file and what it disables. |
| — | `moldynx intake` / `INTAKE_REPORT.md`: extended runs (declared vs. executed steps), temperature changes between stages, the cluster-clock offset, inputs that job scripts reference but that are absent. |
| — | Read-only readers for mdrun logs, `.mdp`, XTC frame headers (no cache files written into the data folder) and `.tpr` headers; a GROMACS runner that falls back to WSL on Windows. |

## 2. Periodic boundaries

| mdforge 0.2 | MolDynX Tools 0.3 |
|---|---|
| `unwrap` inside `try/except: pass`: a topology without bonds silently gave a trajectory that was never made whole. | One pass per frame: **diagnose** the raw coordinates, **treat** them (`none`/`whole`/`nojump`, first-frame clustering for complexes), **prove** the result (only whole-box translations, bonds short, partner distances preserved). Failures are recorded, never swallowed. |
| Cached solute trajectory reused whenever the files existed. | Cache keyed on input fingerprints, PBC mode, frame slice and selection. |
| Chain identity lost: PDB truncates CHARMM-GUI segment IDs (`seg_0_PROA`, `seg_1_PROB` → `seg_`), so `interface` could not find two partners. | Chain identity persisted as atom-index ranges; every complex analysis uses it. |

## 3. Preparation and equilibration (new)

Minimisation outcome (states plainly when the force tolerance was *not* reached) and the residue
carrying the largest residual force; NVT/NPT statistics, settling times and residual drift;
position restraints read from each run input (and none in production); protonation states;
crash-dump accounting; chain of custody from job scripts; stage timeline — as
`EQUILIBRATION.md`. On one reference system it located the largest residual force on a
protonated aspartate, which the manual audit had not connected.

## 4. Complex analyses

| mdforge 0.2 | MolDynX Tools 0.3 |
|---|---|
| `interface` "implemented but not validated"; atom-pair contacts including hydrogens. | Rewritten from a validated suite: residue–residue heavy-atom contacts, full (never truncated) core interface lists, persistent contacts, residues within 6 Å at any time, buried area, interface Cα RMSD, trends. |
| — | `contact_lifetime`, `water_bridges`, `porcupine` (SVD, no 3N × 3N matrix). |
| — | `annotation`: biological numbering, domains transferred by sequence alignment (uncertain edges flagged), motif location with partial matches, docking-site involvement — only from user input. |
| — | `analysis_window`: Chodera equilibration detection and drift tests; the averaging window is a user decision when nothing is stationary. |
| SASA via one multi-frame `mdtraj.shrake_rupley` call. | **Bug found:** mdtraj 1.11.1 returns wrong values for some frames of a multi-frame call (±0.5 nm² on a real trajectory; negative buried areas in a rigid-body test). MolDynX computes frame by frame. |

## 5. Binding energy

| mdforge 0.2 | MolDynX Tools 0.3 |
|---|---|
| Wrote a placeholder script: groups `r 1-100` / `r 101-9999`, an Amber `forcefields = "oldff/leaprc.ff99SB, leaprc.gaff"` line in a CHARMM/GROMACS workflow, the full solvated system, never ran. A contact-count "proxy" was shown under an MM/PBSA label. | Full [gmx_MMPBSA](https://github.com/Valdes-Tresanca-MS/gmx_MMPBSA) workflow: protein-only system cut from the run input, groups from chain identity, **ionic strength derived** from the ions in the system, **explicit decomposition residues** (every residue that contacted the partner in any frame — `within 6` misses late binders), a **probe gate** (bonded Δ terms must be zero), GB and PB runs, and an analysis with autocorrelation-corrected errors, drift, GB vs PB, hotspots in biological numbering, decomposition closure and an entropy validity gate. Reproduces a manual analysis of real outputs exactly. gmx_MMPBSA and MMPBSA.py are credited and cited in every output. |

## 6. Documents and datasets

| mdforge 0.2 | MolDynX Tools 0.3 |
|---|---|
| One generic report; the HTML renderer deleted `**` and ignored inline code, links and subscripts. | Report plus `PBC_VALIDATION`, `EQUILIBRATION` and `BINDING_ENERGY` documents, written from the results files, with Observations / Interpretation / Limitations; proper Markdown → self-contained HTML. |
| Interpretation said "**Stability.** Backbone RMSD reached a plateau". | Wording rules: no stability claim from an RMSD plateau; binding energies are end-point estimates, not affinities; entropy only when valid; unknown metadata stated as not recorded. |
| — | `moldynx dataset` / `verify` / `package`: numbered dataset, raw-file manifest, `verify_dataset.py` shipped inside, zip extracted and re-verified. |

## 7. Engineering

- Renamed to **MolDynX Tools** (`moldynx`); `import mdforge`, the `mdforge` command and the
  `mdforge.plugins` entry-point group keep working until 0.4 with a deprecation warning.
- Autocorrelation statistics (`statistical_inefficiency`, corrected SEM, equilibration
  detection, drift) in `moldynx.statistics`.
- Tests: from 7 unit tests to a suite running on real, trimmed data from two production
  simulations (logs, energy file, gmx_MMPBSA outputs, raw-folder listings).
- Default output `./moldynx_results/<input folder>`; MolDynX never writes into the simulation
  folder.
