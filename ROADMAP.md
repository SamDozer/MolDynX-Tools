# MolDynX Tools Roadmap

> Turning MolDynX Tools from a solid single-run analyzer into a **comparative, parallel,
> multi-engine** MD analysis platform — inspired by the best of
> [MDAnalysis](https://github.com/MDAnalysis/mdanalysis) and
> [mdtraj](https://github.com/mdtraj/mdtraj), built to give the user *more options*
> and to **portray many aspects of a system visually in one place**.

This document is a living plan: vision → design influences → the flagship
**comparison mode** → a versioned feature roadmap with pseudocode and task lists →
architecture upgrades → testing/docs/release plan.

---

## 0. Guiding principles

1. **Compare, don't just describe.** The scientific payoff of MD is usually *"what
   changed?"* — apo vs. holo, wild-type vs. mutant, monomer vs. complex. Every
   analysis should be overlay-able across systems on one figure.
2. **More options, sane defaults.** Expose reference choice, alignment selection,
   units, PBC treatment, per-analysis parameters — but never require them.
3. **Scale honestly.** Stream and parallelize so "100 GB" is benchmarked, not
   claimed.
4. **Reproducible & citable.** Manifest + method citations for everything.
5. **Don't reinvent kernels.** Lean on MDAnalysis/mdtraj C code; add value in
   orchestration, comparison, detection, reporting.

---

## 1. Design influences — what MDAnalysis & mdtraj do better, and what we adopt

| Capability | Their approach | MolDynX today | Planned adoption |
|---|---|---|---|
| **Per-frame analysis + parallelism** | MDAnalysis `AnalysisBase`: `_prepare/_single_frame/_conclude` + `run(backend="multiprocessing"/"dask", n_workers, n_parts)` with split-apply-combine (`_get_aggregator()` → `ResultsGroup`, `ndarray_vstack`) | `run(ctx)` monolithic, serial | **`FrameAnalysis` base** with the same lifecycle + a parallel executor (§4.1) |
| **Ensemble similarity** | `encore`/`mdaencore`: `hes()`, `ces()`, `dres()`, `*_convergence()` (JS divergence between ensembles) | none | Wrap `mdaencore` for **comparison mode + convergence** (§3, §4.4) |
| **Functional API** | mdtraj `compute_rmsd/dssp/contacts/gyration_tensor/…` — composable, notebook-friendly | class + CLI only | Add `moldynx.compute.*` thin functions (§4.2) |
| **Chunked streaming** | mdtraj `md.iterload(traj, chunk=100, top=…)` | pre-extract a solute trajectory once | Add real **chunked iteration** + keep extraction as an optimization (§4.3) |
| **C-kernel speed** | mdtraj SIMD RMSD (QCP), contacts, SASA, DSSP, H-bonds | mixed (MDAnalysis Python loops in places) | Route hot paths (pairwise-RMSD, contacts, hbonds) through mdtraj kernels (§4.5) |
| **On-the-fly transforms** | MDAnalysis `transformations` (unwrap/center/fit) composed on the trajectory | hard-coded unwrap in core extraction | Expose a **transform pipeline** via config (§4.6) |
| **Citation duty** | MDAnalysis prints per-method citations (`Duecredit`) | none | Collect per-analysis **references** into report + manifest (§4.7) |
| **Rich selections/options** | powerful selection DSL, updating selections, reference choices | fixed `protein`/`name CA` | Per-analysis **selection/reference/alignment/unit** options (§5) |
| **Format breadth** | both read AMBER/NAMD/CHARMM/… by extension | GROMACS-centric discovery | Broaden discovery + detection to any MDAnalysis-readable engine (§v0.4) |

**Principle:** *influence, not copy.* We reimplement the patterns against MolDynX's
own `ctx`/registry design and add the comparison/report/detection layers those
libraries deliberately leave out.

---

## 2. Feature themes (backlog, grouped)

- **Comparison & statistics** — overlays, difference maps, ensemble similarity, per-metric hypothesis tests, replica averaging.
- **Performance** — parallel backends, chunked streaming, content-addressed cache/resume, benchmarks.
- **Coverage** — membrane suite, multi-engine, validated complex modules, advanced dynamics (dPCA, LMI-DCCM, MSM).
- **UX & ecosystem** — functional API, `rich` UI, interactive HTML, docs site, PyPI/conda-forge, `moldynx reproduce`.

---

## 3. FLAGSHIP — Comparison mode (`moldynx compare`)  ⭐ v0.3

**Goal:** given several labeled systems — e.g. `control` (protein alone), `+ligand`,
`+partner` (protein–protein) — analyze the **common entity** (the protein) in each
and **overlay every metric on shared axes**, quantify the differences statistically,
and answer *"what does binding/partnering do to the protein?"*.

### 3.1 Command-line interface

```bash
# explicit labels -> directories
moldynx compare \
    --system control=/runs/apo \
    --system holo=/runs/with_ligand \
    --system dimer=/runs/complex \
    --common "protein and name CA" \      # entity compared across systems
    --reference control \                 # align/compare against this system's frames
    --analyses rmsd,rmsf,rog,sasa,pca,dccm,secondary_structure,hbonds \
    --output compare_results/

# or drive it from YAML for reproducibility / HPC
moldynx compare --config compare.yaml
```

```yaml
# compare.yaml
systems:
  control:  { input: /runs/apo }
  holo:     { input: /runs/with_ligand, ligand: LIG }
  dimer:    { input: /runs/complex, chains: "A" }   # compare chain A only
common_selection: "protein and name CA"
reference: control            # common subspace / alignment reference
stride: 5
analyses: [rmsd, rmsf, rog, sasa, pca, dccm, secondary_structure, hbonds]
statistics: [ks, welch_t, cohen_d, ensemble_similarity]
report_formats: [md, html]
```

### 3.2 What gets compared, and how (the science)

| Metric | Overlay | Difference / insight |
|---|---|---|
| RMSD | all systems vs. their own ref, one axis | stability shift on binding |
| **RMSF (per residue)** | overlaid lines by residue | **ΔRMSF map** → where ligand rigidifies/loosens the protein (binding-site footprint) |
| Rg / SASA / end-to-end | overlaid time-series + violin/box | compaction / burial changes |
| Secondary structure | %helix/sheet/coil bars per system | folding stabilization |
| **PCA** | project **all** ensembles onto the **common (reference) subspace** | conformational-space shift; overlap in essential subspace |
| **FEL** | side-by-side + **ΔG difference map** on shared PC axes | basins gained/lost on binding |
| **DCCM** | per-system + **difference matrix** | coupling induced/broken by the partner |
| **Contacts** | occupancy difference map | contacts formed/lost |
| H-bonds / salt bridges | overlaid counts + shared occupancy table | network remodeling |

### 3.3 Design & pseudocode

```python
# moldynx/compare/config.py
@dataclass
class CompareConfig:
    systems: dict[str, RunConfig]        # label -> per-system config
    common_selection: str = "protein and name CA"
    reference: str | None = None         # label used as the common frame of reference
    analyses: list[str] | None = None
    statistics: list[str] = ("ks", "welch_t", "cohen_d", "ensemble_similarity")
    output_dir: Path = Path("compare_results")
    stride: int = 1
    report_formats: list[str] = ("md", "html")
```

```python
# moldynx/compare/pipeline.py
def run_comparison(cfg: CompareConfig):
    # 1) Run each system through the normal pipeline, but RESTRICTED to the common
    #    selection so metrics are apples-to-apples. Reuse existing analyses.
    per_system = {}
    for label, sys_cfg in cfg.systems.items():
        sys_cfg.analyses = cfg.analyses
        sys_cfg.params.setdefault("_common", {})["selection"] = cfg.common_selection
        per_system[label] = run_pipeline(sys_cfg)      # existing engine, unchanged

    ref = cfg.reference or next(iter(cfg.systems))

    # 2) Build a COMMON PCA subspace on the reference, project every system into it.
    common_space = build_common_subspace(per_system[ref], cfg.common_selection)
    projections = {lab: common_space.project(res) for lab, res in per_system.items()}

    # 3) Overlay every comparable metric + compute difference maps.
    overlays = OverlayBuilder(per_system, reference=ref)
    overlays.timeseries("rmsd.csv", "rmsd_backbone_nm", ylabel="RMSD (nm)")
    overlays.per_residue("rmsf.csv", "rmsf_nm", diff=True)        # ΔRMSF map vs control
    overlays.distribution("rog.csv", "rg_nm", kind="violin")
    overlays.matrix_difference("dccm.csv")                        # DCCM(sys) - DCCM(ref)
    overlays.matrix_difference("contact_occupancy.csv")
    overlays.fel_difference(projections, reference=ref)          # ΔG landscape shift
    overlays.secondary_structure_bars(per_system)

    # 4) Statistics: per-metric hypothesis tests + ensemble similarity.
    stats = compare_statistics(per_system, cfg.statistics)       # KS, Welch t, Cohen d
    similarity = ensemble_similarity(projections)                # CES/DRES/HES (mdaencore)

    # 5) Comparative report (reuses the report renderer with a comparison layout).
    generate_comparison_report(cfg, per_system, overlays, stats, similarity)
```

```python
# moldynx/compare/overlay.py  — one figure, many systems
PALETTE_BY_LABEL = cycle_palette()   # deterministic, colour-blind-safe per label

def timeseries(self, csv, col, ylabel):
    fig, ax = new_axes()
    for label, res in self.per_system.items():
        df = read_csv(res.output_dir / "results" / csv)
        ax.plot(df.time_ns, df[col], label=label, color=PALETTE_BY_LABEL[label], alpha=.85)
    ax.legend(title="system"); ax.set_ylabel(ylabel); ax.set_xlabel("Time (ns)")
    save_figure(fig, self.figdir / f"compare_{stem(csv)}")

def per_residue(self, csv, col, diff=True):
    ref = read_csv(self.ref_dir / csv).set_index("resid")[col]
    fig, ax = new_axes()
    for label, res in self.per_system.items():
        s = read_csv(res.output_dir / "results" / csv).set_index("resid")[col]
        y = (s - ref) if (diff and label != self.reference) else s
        ax.plot(s.index, y, label=label)
    ax.axhline(0, ls="--", color="grey")           # ΔRMSF baseline
    ax.set_ylabel(f"Δ{col} vs {self.reference}" if diff else col)
    save_figure(fig, self.figdir / f"compare_{stem(csv)}")
```

```python
# moldynx/compare/subspace.py  — the key to a MEANINGFUL PCA comparison
class CommonSubspace:
    """PCA eigenvectors from the reference ensemble; project any system into them."""
    def __init__(self, ref_ca_coords):                 # (n, 3N), aligned
        Xc = ref_ca_coords - ref_ca_coords.mean(0)
        evals, evecs = np.linalg.eigh(Xc.T @ Xc / (len(Xc) - 1))
        self.mean = ref_ca_coords.mean(0)
        self.evecs = evecs[:, ::-1]                    # ref principal axes
    def project(self, ca_coords, k=2):
        return (ca_coords - self.mean) @ self.evecs[:, :k]   # same axes for all systems
```

```python
# moldynx/compare/statistics.py
def compare_statistics(per_system, tests):
    rows = []
    for metric in ("rmsd_backbone_nm", "rg_nm", "sasa_total_nm2"):
        series = {lab: load_metric(res, metric) for lab, res in per_system.items()}
        ref = list(series)[0]
        for lab, x in series.items():
            if lab == ref: continue
            rows.append({
                "metric": metric, "system": lab, "vs": ref,
                "delta_mean": x.mean() - series[ref].mean(),
                "ks_p":    ks_2samp(series[ref], x).pvalue,     # distribution shift
                "welch_p": ttest_ind(series[ref], x, equal_var=False).pvalue,
                "cohen_d": cohens_d(series[ref], x),            # effect size
            })
    return pd.DataFrame(rows)

def ensemble_similarity(projections):
    # Influence: MDAnalysis encore / mdaencore  (CES/DRES/HES, JS-divergence based).
    # Quantifies how different each ensemble is from the reference in [0,1].
    import mdaencore as encore     # optional dependency
    return encore.ces([u_for(lab) for lab in projections])      # similarity matrix
```

### 3.4 Tasks (v0.3 comparison)
- [ ] `moldynx/compare/` package: `config.py`, `pipeline.py`, `overlay.py`, `subspace.py`, `statistics.py`, `report.py`.
- [ ] `moldynx compare` CLI subcommand + YAML schema.
- [ ] Restrict per-system analyses to a shared `common_selection` (add a "common" param honored by RMSD/RMSF/Rg/PCA/…).
- [ ] Common-subspace PCA projection + FEL difference.
- [ ] ΔRMSF / ΔDCCM / Δcontact difference maps.
- [ ] Per-metric KS / Welch-t / Cohen's d table + violin/box overlays.
- [ ] Optional `mdaencore` ensemble similarity (CES/DRES/HES) with graceful skip.
- [ ] Comparative Markdown/HTML report with a "what changed" auto-summary.
- [ ] Tests on two tiny fixture trajectories (same protein, ± a dummy ligand).

---

## 4. Architecture upgrades

### 4.1 `FrameAnalysis` base with parallel execution (influence: MDAnalysis `AnalysisBase`)

```python
# moldynx/core/frame_analysis.py
class FrameAnalysis(BaseAnalysis):
    """Opt-in lifecycle base for per-frame analyses that can run in parallel."""
    parallelizable = True

    def _prepare(self, ctx): ...                       # allocate self.results
    def _single_frame(self, ts, ag): ...               # per-frame; independent
    def _conclude(self, ctx): ...                       # finalize/normalize

    def run(self, ctx):
        backend  = ctx.config.backend                   # "serial"|"multiprocessing"|"dask"
        n_workers = ctx.config.threads
        frames = list(sliced(ctx))                      # (start, stop, step)
        if backend == "serial" or not self.parallelizable or n_workers == 1:
            return self._run_serial(ctx, frames)
        parts  = split(frames, n_parts=n_workers)       # balanced groups
        chunks = executor(backend).map(self._compute_part, parts)   # workers
        self.results = self._aggregate(chunks)          # ResultsGroup-style combine
        self._conclude(ctx)
        return self.summary()

    def _aggregate(self, chunks):                       # e.g. vstack time-series
        return {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}
```

*Migration:* RMSD/RMSF/Rg/SASA/COM/contacts/DCCM become `FrameAnalysis` subclasses;
interface/report modules stay on the simpler `BaseAnalysis`. Backward compatible.

### 4.2 Functional API (influence: mdtraj `compute_*`)

```python
# moldynx/compute.py  — notebook-friendly, no CLI needed
import moldynx.compute as mfc
rmsd  = mfc.rmsd(u, select="backbone", ref_frame=0)          # -> np.ndarray (nm)
rg    = mfc.radius_of_gyration(u, select="protein")
q     = mfc.native_contacts(u, cutoff=8.0)
fel   = mfc.free_energy_landscape(pc1, pc2, T=310)
sim   = mfc.ensemble_similarity([u_apo, u_holo])
# every CLI analysis is a thin wrapper over these.
```

### 4.3 Chunked streaming reader (influence: mdtraj `iterload`)

```python
# moldynx/io/stream.py
def iter_chunks(topology, trajectory, selection="all", chunk=200, stride=1):
    u = mda.Universe(topology, trajectory)
    ag = u.select_atoms(selection)
    buf = []
    for ts in u.trajectory[::stride]:
        buf.append(ag.positions.copy())
        if len(buf) == chunk:
            yield np.asarray(buf); buf.clear()
    if buf: yield np.asarray(buf)
# Lets analyses process 100 GB trajectories in bounded memory without full extraction.
```

### 4.4 Ensemble similarity + convergence module (influence: encore/mdaencore)
`moldynx/analysis/ensemble_similarity.py` (`ensemble` group): CES/DRES/HES between
systems (comparison mode) and `*_convergence` within a run (a rigorous alternative
to the current RMSIP/block-averaging convergence check).

### 4.5 Kernel routing — prefer mdtraj C kernels for hot paths
Pairwise-RMSD matrix (clustering), contacts, SASA, DSSP, baker-hubbard H-bonds →
mdtraj; keep MDAnalysis for selections, PBC transforms, `.edr`, salt bridges.

### 4.6 Transform pipeline (influence: MDAnalysis `transformations`)

```yaml
transforms: [unwrap, center_in_box(protein), fit_rot_trans(protein and name CA)]
```

### 4.7 Citation duty (influence: MDAnalysis Duecredit)
Each analysis declares `citations = ["DOI/…"]`; the report gains a **References**
section and `manifest.json` records methods used → publishable provenance.

---

## 5. "More options" — user-facing knobs (rolling, from v0.3)

- **Reference choice:** `--reference first|average|<file.pdb>|<system-label>`.
- **Selections:** per-analysis `selection`, `align_on`, `common` overrides.
- **Units:** `--units nm|angstrom` everywhere.
- **PBC:** `--pbc none|whole|nojump|cluster` (transform pipeline).
- **Frames:** `--start/--end/--stride`, or `--frames "0:1000:2,1500,2000"`.
- **Output:** `--figure-format png,pdf,svg`, `--dpi`, `--palette`, `--theme light|dark`.
- **Analyses:** `--only-category structure`, `--exclude`, per-analysis params via YAML.
- **Report:** `--report md,html,pdf`, `--no-figures`, `moldynx report <dir>` to regen.

---

## 6. Versioned roadmap

### v0.3 — Comparison & parallelism  *(flagship)*
- ⭐ `moldynx compare` (§3) — overlays, difference maps, common-subspace PCA, stats, ensemble similarity.
- `FrameAnalysis` + parallel `--threads` backend (§4.1); wire `multiprocessing`.
- Functional `moldynx.compute.*` API (§4.2).
- Content-addressed **cache/resume** (skip unchanged analyses; checksum of inputs+params).
- `moldynx reproduce manifest.json` (re-run a past run bit-for-bit).
- Expanded options (§5): reference/units/frames/selections.

### v0.4 — Coverage: membrane + multi-engine + validated complexes
- **Membrane suite:** area-per-lipid, bilayer thickness, S_CD order parameters, lipid–protein contacts, insertion depth/tilt (`membrane` category, gated on `has_lipid`).
- **Multi-engine:** broaden discovery/detection to AMBER (`prmtop`/`nc`), NAMD (`psf`/`dcd`), CHARMM, OpenMM (`.h5`) — anything MDAnalysis reads.
- **Validate & ship** interface/ligand/nucleic with committed **fixture trajectories** + golden-file CI.
- Chunked streaming reader (§4.3) as default for solvent analyses.

### v0.5 — Advanced dynamics & interactive reports
- Dihedral-PCA (dPCA), **linear mutual-information DCCM** (nonlinear coupling), porcupine mode plots.
- Optional **Markov State Models** (deeptime): lag-time, implied timescales, macrostates.
- Water analyses: residence time, bridging waters, hydration density.
- **Interactive HTML** reports (plotly): zoomable time-series, hover tables, collapsible sections.

### v0.6 — Robustness & UX polish
- `rich` CLI (progress, tables, `--plan` pretty output).
- `moldynx config init` scaffolder; user-extensible detection dictionaries.
- Missing-hydrogen handling for H-bond analyses (optional protonation step).
- Benchmark suite + performance regression guard.

### v1.0 — Maturity
- Docs site (mkdocs-material + mkdocstrings API docs), tutorials/notebooks.
- **PyPI + conda-forge** packages; `mypy`/`ruff` pre-commit; codecov badge.
- Stable public API + semantic-versioning guarantees; Zenodo DOI per release.

---

## 7. Cross-cutting: testing, CI, docs, release

- **Fixtures:** commit tiny (<1 MB) trajectories (protein-only, protein+ligand, 2-chain, DNA) so CI runs *real* end-to-end analyses and comparison, not just unit tests.
- **Golden-file regression:** freeze key numeric outputs; fail CI on drift beyond tolerance.
- **Cross-validation:** spot-check RMSD/Rg/SASA against `gmx` native tools and mdtraj.
- **Docs:** every analysis's docstring → API site; a "compare apo vs holo" tutorial as the headline example.
- **Release:** CHANGELOG, semver tags, auto-archive to Zenodo (already wired).

---

## 8. Suggested execution order (next 3 PRs)

1. **PR: comparison MVP** — `moldynx compare` with RMSD/RMSF/Rg/SASA overlays + ΔRMSF map + KS/Welch/Cohen table + comparative report. (Delivers the headline value fastest.)
2. **PR: `FrameAnalysis` + `--threads`** — parallel execution for the per-frame analyses; benchmark on the 4.35 GB reference trajectory.
3. **PR: common-subspace PCA + FEL/DCCM difference + ensemble similarity** — completes the "what changed conformationally" story.

*Contributions welcome — each analysis is a drop-in `BaseAnalysis`/`FrameAnalysis`;
see [`docs/QUICKSTART.md`](docs/QUICKSTART.md) and the plugin template.*
