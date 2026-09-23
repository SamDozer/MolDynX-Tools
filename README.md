# MolDynX Tools — reusable, reproducible, audited GROMACS MD analysis

[![CI](https://github.com/SamDozer/molecular-dynamics-forge/actions/workflows/ci.yml/badge.svg)](https://github.com/SamDozer/molecular-dynamics-forge/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21265946.svg)](https://doi.org/10.5281/zenodo.21265946)

**MolDynX Tools** (`moldynx`; formerly *mdforge*) turns a GROMACS simulation directory into a complete, publication-quality,
fully reproducible analysis — with minimal input. Point it at a folder; it discovers
the files, **detects the system** (protein / ligand / DNA / RNA / membrane / ions /
multi-chain / …), **auto-selects the right analyses**, runs them with streaming-friendly
performance, and produces figures, tables, a provenance manifest, and a report.

> This repository began as a project-specific pipeline for an α-zein simulation
> (now preserved in [`legacy/`](legacy) and [`examples/`](examples)) and was
> refactored into this general toolkit.

---

## Highlights

- **Zero-config detection** — recursively finds `*.xtc/.trr/.tpr/.gro/.edr/.ndx/.top`,
  validates them (clear messages for anything missing), and classifies the system.
- **Automatic module selection** — each analysis declares the system types and files
  it supports; the pipeline runs exactly what applies (`--plan` shows *why*).
- **Extensible via plugins** — drop a `BaseAnalysis` subclass into
  `moldynx/analysis/plugins/` (or a `--plugin-dir`) and it is auto-discovered.
- **Config-driven** — describe a whole run in `config.yaml` and re-run with one command
  (ideal for HPC/batch).
- **Reproducible by construction** — every run writes `manifest.json/yaml` with library
  versions, git commit, input-file fingerprints, seeds, parameters and runtimes.
- **Publication-quality output** — 300-dpi PNG **and** vector PDF, consistent
  Nature-like style, plus Markdown + self-contained HTML (+ optional PDF) reports.
- **Scales** — streams frame-by-frame and caches a solute-only trajectory, so large
  (100 GB+) explicit-solvent runs stay tractable.

## Supported systems

protein-only · protein–protein · protein–peptide · protein–ligand · protein–DNA ·
protein–RNA · protein–membrane · multi-chain · ions · cofactors · mixed biomolecular
systems (detected automatically; override with `--system-type`).

## Install

```bash
git clone https://github.com/SamDozer/molecular-dynamics-forge
cd molecular-dynamics-forge
python -m pip install -e ".[all]"      # or ".[dev]" for tests
```

## Usage

```bash
moldynx detect  --input /path/to/sim_dir            # what's in my system?
moldynx analyze --input /path/to/sim_dir --plan     # what would run, and why?
moldynx analyze --input /path/to/sim_dir -o results # run everything applicable
moldynx analyze --config examples/alpha_zein_A8HNE1/config.yaml   # reproducible
moldynx list-analyses                                # registered analyses (incl. plugins)
```

See [`docs/QUICKSTART.md`](docs/QUICKSTART.md) for all options and the plugin template.

## Architecture

```
moldynx/
  core/       system.py (detection) · base.py (BaseAnalysis) · registry.py (+plugins)
              context.py · config.py (YAML+CLI) · provenance.py · pipeline.py
  io/         discovery.py · validation.py
  statistics/ descriptive · timeseries · correlation · bootstrap
  plotting/   style · figures (PNG+PDF, 300 dpi)
  analysis/   rmsd · rmsf · rog · sasa · … + plugins/ (auto-discovered)
  report/     generator (Markdown + HTML + PDF)
  cli/        main (analyze/detect/list) · interactive
tests/  ·  examples/  ·  docs/  ·  legacy/  ·  Dockerfile  ·  .github/workflows/
```

Every analysis subclasses `BaseAnalysis`, declaring `required_files`,
`supported_systems`, `outputs` and `default_params`; subclasses auto-register, so
**detection → selection → run → report** is entirely data-driven.

## Analyses (24 built-in)

| System scope | Analyses |
|---|---|
| **Any system** | RMSD, radius of gyration, SASA, COM, H-bonds, energies (`.edr`), ProLIF, MM/PBSA workflow, statistics |
| **Protein** | RMSF, structural descriptors (Dmax/κ²/volume), native contacts (Q), contact map, DSSP secondary structure, RIN, PCA + free-energy landscape, DCCM, clustering, salt bridges, convergence (RMSIP/block-avg), end-to-end *(plugin)* |
| **Complex** (protein–protein/–nucleic) | interface (BSA, contacts, iRMSD) |
| **Protein–ligand** | ligand RMSD, ligand contacts, binding pocket |
| **Protein–DNA/RNA** | protein–nucleic contacts, nucleic RMSD |

Each is a drop-in `BaseAnalysis`; the pipeline runs only those applicable to the
detected system (`moldynx list-analyses` shows all; `--plan` shows what runs and why).
The complex/ligand/nucleic modules are implemented and gate correctly but await
validation on a matching test trajectory.

## Reproducibility

```bash
cat results/manifest.json     # versions, git commit, seeds, params, input hashes, runtimes
```

## Container

```bash
docker build -t moldynx .
docker run --rm -v /data/sim:/sim moldynx analyze --input /sim --output /sim/results
```

## Roadmap

See **[ROADMAP.md](ROADMAP.md)** for the plan — the flagship being a
**comparison mode** (`moldynx compare`) that overlays control vs. protein–ligand /
protein–protein systems on shared axes (ΔRMSF maps, common-subspace PCA, ensemble
similarity), plus parallel execution, a functional API, membrane and multi-engine
support — drawing design influence from
[MDAnalysis](https://github.com/MDAnalysis/mdanalysis) and
[mdtraj](https://github.com/mdtraj/mdtraj).

## Citation

If you use MolDynX Tools, please cite it (concept DOI — always resolves to the latest version):

> Mahmoud, H. *MolDynX Tools: a reusable, reproducible analysis framework for GROMACS
> molecular dynamics simulations.* Zenodo. https://doi.org/10.5281/zenodo.21265946

A machine-readable [`CITATION.cff`](CITATION.cff) is included (GitHub shows a
"Cite this repository" button).

## License

MIT — see [LICENSE](LICENSE).
