# Test fixtures (real data, trimmed — 0.96 MB total)

All files come from the two finished datasets. Logs are trimmed: the head (banner, command line,
the full `Input Parameters:` dump and the `grpopts` block) and the tail (averages, timings,
minimisation result, `Finished mdrun`) are kept verbatim, and every marker line in between
(`Command line`, `Started mdrun`, `Finished mdrun`, `Statistics over`, checkpoint reads, `Wrote pdb`,
`Maximum force`, …) is kept with a line of context. Trimmed stretches are marked with
`[... N lines trimmed by the mdforge fixture builder ...]` — parsers must tolerate that line.

## gromacs_logs/ — expected parser output

| Fixture | Expected |
|---|---|
| `A8HNE1_step4.0_minimization_machine_precision.log` | integrator `steep`, `emtol` 100, `nsteps` 50000; outcome **machine precision, tolerance not reached**, 3056 steps; Epot −2.6774416e+07; Fmax 7.2325469e+03 on atom 6116; ‖F‖ 11.15; 4 `Wrote pdb` events; command `-ntmpi 1 -ntomp 8` |
| `A8HNE1_step4.1_equilibration_NVT.log` | md, dt 0.001, nsteps 125000; V-rescale; ref-t 303.15 ×2; tau-t 1 ×2; nrdf 33024 / 3.11729e+06; pcoupl No; continuation false; 1 session; 5.418 ns/day; started Thu Jun 4 04:11:36 2026, finished 04:44:50 |
| `A8HNE1_step4.2_equilibration_NPT_303K.log` | nsteps 2000000; ref-t **303.15**; pcoupl Berendsen, isotropic, tau-p 1, ref-p 1.0, compressibility 4.5e-05; refcoord-scaling COM; continuation true; 5.397 ns/day; finished Thu Jun 4 13:38:39 2026 |
| `A8HNE1_step5_production_11_sessions.log` | dt 0.002; first-block nsteps 500000 (1-ns TPR, later extended); V-rescale 303.15; C-rescale tau-p 5; **11** `Started mdrun` sessions; last `Statistics over 50000001 steps` (= 100 ns); finished Tue Jun 16 13:40:25 2026 |
| `Q946V6_step4.0_minimization_converged.log` | `emtol` 1000, `nsteps` 5000; outcome **converged** (Fmax < 1000) in 2027 steps; Fmax 7.1497339e+02 on atom 6883; 17 `Wrote pdb` events |
| `Q946V6_step4.2_equilibration_NPT_310K.log` | ref-t **310** ×2 (differs from production 303.15); Berendsen; 3.349 ns/day |
| `Q946V6_step5_production_12_sessions.log` | **12** sessions; last `Statistics over 50000001 steps`; 9.09 ns/day in the final session |

The two NPT logs together are the test for "surface preparation differences between datasets"
(303.15 K vs 310 K).

## edr/

- `A8HNE1_step4.1_equilibration.edr` — a real GROMACS energy file (NVT, 126 records, 0–125 ps).
  With panedr: T[0] = 304.80 K, min T = 211.20 K at 1 ps, mean T over t ≥ 25 ps = 303.187 ± 0.249 K,
  stays within ±2 K of 303.15 from 13 ps on.
- `A8HNE1_timeseries_step4.*.csv` — the series used for the equilibration figure. From the NPT CSV
  the last-500-ps statistics are T 303.145 ± 0.261 K, P 1.65 ± 17.70 bar, density
  1014.53 ± 0.27 kg m⁻³; Box-X 25.500 → 24.989 nm; volume 16,581 → 15,605 nm³.
  (Complete series: 2409 EM records, 126 NVT, 2001 NPT.)

## gmx_mmpbsa/ — α-zein Q946V6–ZmBiP2, 101 frames (0–100 ns at 1 ns)

- `gb/FINAL_RESULTS_MMGBSA.csv`, `pb/FINAL_RESULTS_MMPBSA.csv` — complete files. Sections:
  `Complex / Receptor / Ligand / Delta Energy Terms`, each with a `Frame #,…` header and 101 rows.
  Expected from the Delta section (kcal/mol, mean ± autocorrelation-corrected SEM, frame f ↦
  t = (f − 1) × 0.1 ns): GB 0–100 ns −33.4 ± 2.1, 80–100 ns −37.8 ± 2.4; PB −40.7 ± 2.5 /
  −47.1 ± 3.3; Pearson r(GB, PB) = 0.73; every bonded Δ term (BOND, ANGLE, DIHED, UB, IMP, CMAP,
  1-4 VDW, 1-4 EEL) is exactly 0; GGAS identical in GB and PB.
- `gb/FINAL_DECOMP_MMGBSA_frames1-3.csv` — every section of the decomposition file
  (`Complex / Receptor / Ligand / DELTAS` × `TDC / SDC / BDC`), frames 1–3 only. Residue labels
  look like `L::THR:1` / `R::<RES>:<n>` (R = receptor, renumbered from 1 by gmx_MMPBSA — verify
  against the sequence; L = ligand). Parse with `^([LR])::([A-Z0-9]+):(-?\d+)$`.
- `*/gmx_MMPBSA.log`, `*.dat`, `mmgbsa.in`, `mmpbsa.in` — the exact inputs and logs. The GB
  `.dat` contains `ENTROPY RESULTS (INTERACTION ENTROPY)` and `(C2 ENTROPY)` blocks, each followed by
  gmx_MMPBSA's own warning that σ(Int. Energy) exceeds 3.6 kcal/mol (IE) / 6.0 kcal/mol (C2);
  σ = 86.7 kcal/mol here, so the entropy gate must reject both.

## raw_trees/ — for discovery and canonical-run tests

- `A8HNE1_raw_tree.json` — every file of the A8HNE1 raw folder (464 files: path, bytes, mtime), plus
  an `expected` block: the canonical production TPR/XTC/EDR/log, `topol.top` + `toppar/`, the three
  equilibration stages, what a size-based heuristic would wrongly pick (the NPT TPR as topology, a
  minimisation crash dump as structure), and which files are missing from the archive. Earlier
  analysis outputs under `Production/analysis*` and `Production/mmpbsa_BiP2/` are marked — discovery
  must ignore them. Build a fake tree from this listing (zero-byte or sparse files, sizes via
  metadata) — do not commit multi-GB files.
- `Q946V6_raw_archive_inventory.csv` — all 113 files of the Q946V6 archive with SHA-256, role and
  retained/excluded status: the answer key for classification, duplicate detection (3 byte-identical
  pairs among the excluded `stepN{b,c}.pdb` dumps) and canonical-run identification
  (`step5_production-001.xtc` + `step5_production_ext.tpr`).

## structures/

- `A8HNE1_ZmBiP2_CA_only.pdb` — Cα atoms of the processed two-chain complex (850 residues:
  α-zein A8HNE1 resid 1–187, ZmBiP2 resid 188–850 = ZmBiP2 1–663 + 187). For numbering-offset,
  chain-split and selection tests. Segment IDs are the CHARMM-GUI style `seg_0_PROA`/`seg_1_PROB`
  truncated to 4 characters by the PDB format (`seg_`) — a realistic trap for segid-based
  selections.
