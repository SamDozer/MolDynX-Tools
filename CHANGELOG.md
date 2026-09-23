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
