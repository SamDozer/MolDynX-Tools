# MolDynX Tools — Quick Start

## Install

```bash
git clone https://github.com/SamDozer/molecular-dynamics-forge
cd molecular-dynamics-forge
python -m pip install -e ".[all]"     # core + energy + fingerprints + pdf + ui
```
or with conda:
```bash
conda env create -f environment.yml && conda activate zein-md
python -m pip install -e .
```

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
