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
