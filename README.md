# MolDynX Tools — audited, reproducible GROMACS MD analysis

[![CI](https://github.com/SamDozer/molecular-dynamics-forge/actions/workflows/ci.yml/badge.svg)](https://github.com/SamDozer/molecular-dynamics-forge/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21265946.svg)](https://doi.org/10.5281/zenodo.21265946)

**MolDynX Tools** (`moldynx`; formerly *mdforge*) turns a GROMACS simulation folder into a
complete, **audited** analysis: it finds the right files by evidence, proves that trajectory
processing changed nothing it should not, audits how the system was prepared and equilibrated,
analyses the structure, dynamics and interfaces, computes MM-GBSA/MM-PBSA binding energies with
gmx_MMPBSA, and writes documents and a self-verifying dataset in which every number traces back
to a file.

> What is new compared with mdforge 0.2: [docs/WHATS_NEW_0.3.md](docs/WHATS_NEW_0.3.md).

---

## Quick start

```bash
git clone https://github.com/SamDozer/molecular-dynamics-forge
cd molecular-dynamics-forge
python -m pip install -e ".[all]"            # or ".[dev]" for tests

moldynx intake  --input /path/to/sim_dir     # which files form the run, what is missing, why
moldynx analyze --input /path/to/sim_dir     # everything applicable -> ./moldynx_results/<folder>
moldynx binding-energy --input /path/to/sim_dir [--execute]   # MM-GBSA/PBSA (gmx_MMPBSA)
moldynx dataset --run moldynx_results/<folder> --out my_dataset
moldynx package my_dataset                   # zip, extract, verify
```

MolDynX never writes into the simulation folder.

## What a run does

| Stage | What happens | Output |
|---|---|---|
| **Intake** | files classified by simulation stage (setup / minimisation / NVT / NPT / production); the production run chosen by evidence (atom counts, time span, log completion) — never by size or name; earlier analysis folders, backups and crash dumps ignored; extensions, clock offsets and inputs referenced by job scripts but absent are reported | `intake/INTAKE_REPORT.md`, capability matrix |
| **PBC: diagnose → treat → prove** | each raw frame is measured (split molecules, continuity, PBC-aware partner distance), treated (`--pbc auto\|none\|whole\|nojump`), then checked against the raw frame | `PBC_VALIDATION.md`, `pbc_summary.json` |
| **Preparation audit** | minimisation outcome (and where the largest force sits), NVT/NPT statistics and residual drift, position restraints, protonation states, chain of custody, timeline | `EQUILIBRATION.md` |
| **Analyses** | RMSD, RMSF, Rg, SASA, PCA, DCCM, clustering, secondary structure, H-bonds, salt bridges, RIN, ProLIF, … plus for complexes: interface suite, contact lifetimes, water bridges, porcupine | figures, CSVs |
| **Annotations** | display names, biological numbering, domains transferred by alignment, motifs, docking site — only from your config | `annotation.json`, `region_summary.csv` |
| **Stationarity** | Chodera equilibration detection and drift tests; no silent choice of averaging window | `window_selection.json` |
| **Binding energy** | gmx_MMPBSA package (protein-only system, derived ionic strength, explicit decomposition residues, probe gate) → analysis with autocorrelation-corrected errors | `BINDING_ENERGY.md` |
| **Documents & dataset** | report + companion documents (Markdown and self-contained HTML), numbered dataset with `verify_dataset.py` | `report/`, dataset folder, zip |

## The input-file contract

| If this is present… | …this becomes possible | Without it |
|---|---|---|
| production run input (`.tpr`) + trajectory with equal atom counts | every analysis | **stop** |
| production `.edr` | energy analysis | skipped, with reason |
| production `.log` | Methods, sessions, extension audit, completion proof | completion reported as not proven |
| `topol.top` + `toppar/` | MM-GBSA/PBSA | skipped; package still prepared |
| minimisation / NVT / NPT `.log` + `.edr` | preparation audit | listed as missing |
| stage `.tpr` files | restraint audit | "not recoverable" |
| stage `.mdp` files | MDP-only settings (`gen-vel`, `define`) | "not recorded" |
| job scripts | explicit chain of custody | inferred, labelled inferred |

## Configuration

Everything can be driven from one YAML file (`moldynx analyze --config config.yaml`); see
[`examples/alpha_zein_A8HNE1/config.yaml`](examples/alpha_zein_A8HNE1/config.yaml). Facts that
cannot be read from files are supplied — never guessed:

```yaml
annotations:
  chains:
    - {segid: seg_0_PROA, display: "Partner A", role: ligand}
    - {segid: seg_1_PROB, display: "Partner B", role: receptor, numbering_offset: 187}
  domains:
    seg_1_PROB: {reference_name: "UniProt P11021", reference_sequence: "MKLS...",
                 regions: {NBD: [26, 405]}}
  motifs:
    - {name: "Motif 1", sequence: "CSQAPIASLLPPYLSPAVSSVC", chain: seg_0_PROA}
  docking_site: {seg_1_PROB: [405, 434, 435]}
  unresolved_metadata: {force_field_variant: null, salt_concentration_M: null}
binding_energy:
  primary_window: [80, 100]      # ns; omit and MolDynX shows the full run + final 20 %
```

## Binding energy (gmx_MMPBSA)

MM-GBSA/MM-PBSA calculations are performed with
[**gmx_MMPBSA**](https://github.com/Valdes-Tresanca-MS/gmx_MMPBSA) (Valdés-Tresanca et al.,
*J. Chem. Theory Comput.* 2021, 17, 6281) on top of AmberTools MMPBSA.py (Miller et al., 2012);
please cite both. MolDynX prepares the calculation (protein-only system; ionic strength derived
from the ions in the run input; decomposition of every residue that contacted the partner in any
frame), gates it (a short probe run must show zero bonded Δ terms), runs it on Linux/WSL, and
analyses it (windows, drift, GB vs PB, hotspots, decomposition closure, entropy validity).
Results are end-point estimates, not experimental affinities.

## Reproducibility and verification

Every run writes `manifest.json/yaml` (versions, git commit, parameters, input fingerprints and
the evidence for every chosen file). Every dataset ships `5_validation/verify_dataset.py`
(layout, links, raw-file manifest, PBC proof, Rg smoke test):

```bash
moldynx verify my_dataset [--full]
```

## Architecture

```
moldynx/
  core/      system (detection, chain records) · pbc · surface · annotations · context · pipeline
             registry (+plugins) · config (YAML+CLI) · provenance
  io/        discovery (stage classification, evidence) · validation (capability matrix)
             gromacs (log/mdp/xtc/tpr readers, gmx runner with WSL fallback) · jobscripts · intake
  analysis/  32 analyses (incl. the bundled plugin) + plugins/ (auto-discovered)
  binding/   gmx_MMPBSA prepare + analyse
  report/    generator + documents (Markdown → self-contained HTML)
  dataset.py · verify_template.py
mdforge/     deprecated import shim (removed in 0.4)
```

Every analysis subclasses `BaseAnalysis` and declares `required_files`, `supported_systems`,
`outputs` and `default_params`; `moldynx analyze --plan` shows what runs and why.

## Tested

66 tests run in CI on real, trimmed data from two 100 ns protein–protein simulations (GROMACS
logs, an energy file, gmx_MMPBSA outputs, raw-folder listings) — see
[`tests/fixtures/README.md`](tests/fixtures/README.md). Intake, PBC, preparation audit and
MM-GBSA/PBSA analysis reproduce an independent manual analysis of those datasets exactly.

## Container

```bash
docker build -t moldynx .
docker run --rm -v /data/sim:/sim moldynx analyze --input /sim --output /sim/results
```

## Roadmap

See [ROADMAP.md](ROADMAP.md) — next: comparison mode (`moldynx compare`), parallel execution,
a process-pool SASA backend, membrane and multi-engine support.

## Citation

> Mahmoud, H. *MolDynX Tools: a reusable, reproducible analysis framework for GROMACS molecular
> dynamics simulations.* Zenodo. https://doi.org/10.5281/zenodo.21265946

A machine-readable [`CITATION.cff`](CITATION.cff) is included.

## License

MIT — see [LICENSE](LICENSE).
