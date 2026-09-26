# MolDynX Tools — next steps after 0.3

All ten phases of the 0.3 plan are implemented (see [WHATS_NEW_0.3.md](WHATS_NEW_0.3.md) and
the CHANGELOG). The reference datasets used as the quality bar (two 100 ns protein–protein
complexes) are test cases only; nothing in the code names a system.

## Left for the maintainer (cannot be automated)

1. **Rename the GitHub repository** (Settings → General → Repository name), e.g. to
   `moldynx-tools`; GitHub redirects the old URL. Then replace `molecular-dynamics-forge` with
   the new name in `README.md` (badges, clone commands), `docs/QUICKSTART.md`, `CITATION.cff`
   and `pyproject.toml` (7 URLs).
2. **Release 0.3.0** on GitHub when satisfied; Zenodo mints a new version DOI under the existing
   concept DOI 10.5281/zenodo.21265946.
3. Decide whether the real test fixtures (trimmed logs and gmx_MMPBSA outputs from the reference
   simulations) should stay public or be replaced by synthetic files.

## Not yet validated end to end on real data

- `interface`, `contact_lifetime`, `water_bridges`, `porcupine`, `annotation` and the full
  `analyze` → `dataset` → `package` chain were tested on real geometry and synthetic
  trajectories, not yet on a full 100 ns run. Reference values from the manual analysis of one
  complex: 27.9 nm² buried area, 0.266 nm mean minimum distance, 34 persistent contacts.
- `moldynx binding-energy --execute` was prepared and gated on real data (complex.tpr, topology,
  index, derived ionic strength 0.158 M); a full GB + PB run (≈ 2.5 h) through the new script is
  still to be done. The explicit decomposition residue set should close the +3.7 kcal/mol gap
  seen with `print_res = "within 6"` in one reference system.

## Engineering backlog

- SASA is computed frame by frame (mdtraj 1.11.1 multi-frame bug) and runs on one core;
  `sasa` and the buried area stride by 10. A process-pool backend would allow stride 1.
- Comparison mode (`moldynx compare`, ROADMAP v0.3) should open with a table of preparation
  differences (equilibration temperature, barostat, minimisation outcome, protonation) — the
  equilibration audit already produces every field.
- Domain boundaries from CATH/SIFTS for a reference structure are supplied by the user today;
  fetching and caching them would remove a manual step.
- Remove the `mdforge` compatibility shim in 0.4.
