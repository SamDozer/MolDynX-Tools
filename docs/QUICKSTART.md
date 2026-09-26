# MolDynX Tools — Quick Start

## Install

```bash
git clone https://github.com/SamDozer/MolDynX-Tools
cd MolDynX-Tools
python -m pip install -e ".[all]"     # core + energy + fingerprints + pdf + ui
```
or with conda:
```bash
conda env create -f environment.yml && conda activate moldynx
python -m pip install -e .
```

## 0. Check the inputs

```bash
moldynx intake --input /path/to/simulation_dir [--detect]
```
Writes `INTAKE_REPORT.md`: which files form the production run and the evidence for it, each
simulation stage, what every missing file disables, run extensions, temperature changes between
stages and inputs your job scripts reference but that are absent. Nothing is written into the
simulation folder.

## 1. Detect your system

```bash
moldynx detect --input /path/to/simulation_dir
```
Prints the discovered files, validation status, and the detected composition
(protein chains, ligands, nucleic acids, lipids, ions, water) and system type.

## 2. See what would run (dry-run)

```bash
moldynx analyze --input /path/to/simulation_dir --plan
```
Shows which analyses are selected for the detected system — and *why* each other
one was skipped.

## 3. Run the full analysis

```bash
moldynx analyze --input /path/to/simulation_dir --output results/
```
Discovers files → detects the system → auto-selects analyses → runs them →
writes CSVs, 300-dpi PNG+PDF figures, a `manifest.json/yaml` provenance record,
and a Markdown+HTML report.

### Handy options
```bash
--stride 10                 # analyse every 10th frame (quick look)
--analyses rmsd,rog,sasa    # run a specific subset
--exclude sasa              # skip an analysis
--all                       # run every applicable analysis
--interactive               # ask when the choice is ambiguous
--report md,html,pdf        # report formats
--plugin-dir ./my_plugins   # load extra drop-in analyses
--config config.yaml        # drive everything from a YAML file
```

## 4. Reproducible, config-driven workflow

```bash
moldynx analyze --config examples/alpha_zein_A8HNE1/config.yaml
```
Every run writes a `manifest.json` capturing library versions, git commit, input
fingerprints, seeds, parameters and per-analysis runtimes — enough for another
researcher to reproduce the analysis exactly.

## 5. Binding energy (MM-GBSA / MM-PBSA with gmx_MMPBSA)

```bash
moldynx binding-energy --input /path/to/simulation_dir            # prepare the package
moldynx binding-energy --input /path/to/simulation_dir --execute  # run it (Linux/WSL, hours)
```
Needs [gmx_MMPBSA](https://github.com/Valdes-Tresanca-MS/gmx_MMPBSA) in a conda environment
(default name `gmxMMPBSA`, override with `MOLDYNX_MMPBSA_ENV`), e.g.
`conda create -n gmxMMPBSA -c conda-forge --override-channels python=3.12 "ambertools>=24.8,<27"
"mpi4py>=4.0.1,<5" "numpy<2"` then `pip install gmx_MMPBSA`. When the outputs exist, `analyze`
reads them and writes `BINDING_ENERGY.md`.

## 6. Share a dataset

```bash
moldynx dataset --run moldynx_results/<folder> --out my_dataset [--include-raw]
moldynx package my_dataset            # zip -> extract -> run verify_dataset.py
moldynx verify  my_dataset --full     # re-check at any time
```

## 5. Add your own analysis (plugin)

Drop a file into `moldynx/analysis/plugins/` (or any `--plugin-dir`):

```python
from moldynx.core.base import BaseAnalysis

class MyAnalysis(BaseAnalysis):
    name = "my_analysis"
    label = "My analysis"
    supported_systems = {"*"}
    def run(self, ctx):
        u = ctx.core_universe()
        ...                       # compute
        ctx.write_csv(df, "my_analysis.csv")
        return {"figure": "my_analysis"}
```
It is auto-discovered and appears in `moldynx list-analyses` immediately.
See `moldynx/analysis/plugins/example_end_to_end.py` for a complete template.
