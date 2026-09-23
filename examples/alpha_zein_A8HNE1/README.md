# Example: α-zein (A8HNE1)

This directory holds the **reference example** that seeded the toolkit: a 100 ns
GROMACS simulation of an AlphaFold 3 α-zein model (187 residues, explicit
solvent, CHARMM36).

- `config.yaml` — a complete, reproducible MolDynX Tools workflow for this system.
- `figures/`, `results/`, `tables/`, `report/` — outputs produced by the original
  project-specific pipeline (now preserved in [`../../legacy/`](../../legacy)),
  kept here as a worked example of what a full analysis looks like.

## Reproduce with MolDynX Tools

```bash
# edit input_dir in config.yaml to point at your copy of the raw files, then:
moldynx analyze --config examples/alpha_zein_A8HNE1/config.yaml
```

The raw simulation files (`step5_production.xtc/.tpr/.edr/...`) are **not**
included (multi-GB); point `input_dir` at your own copy.

## Key finding (reference)

The extended AlphaFold model **compacts** (Rg 4.2 → 3.5 nm) while retaining its
helical content (~55 %) and most native contacts (Q ≈ 0.89) — characteristic of
an elongated, flexible seed-storage protein. See `report/analysis_report.md`.
