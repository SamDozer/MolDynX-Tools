"""
End-to-end pipeline: discover -> validate -> detect -> select -> run -> report.

Each analysis is isolated: a failure is recorded in the manifest and the pipeline
continues, so one broken module never aborts the whole run.
"""

from __future__ import annotations

import time
import traceback
from pathlib import Path

from moldynx.core.config import RunConfig
from moldynx.core.context import AnalysisContext
from moldynx.core.provenance import RunManifest
from moldynx.core.registry import registry
from moldynx.core.system import detect_system, SystemType
from moldynx.io.discovery import discover_files, pick_trajectory
from moldynx.io.validation import validate_fileset


def build_plan(config: RunConfig, ask=None):
    """Resolve files + system + selected analyses without running anything."""
    fs = discover_files(config.input_dir)
    fs.trajectory = pick_trajectory(fs, interactive=config.interactive, ask=ask)
    if config.trajectory:
        fs.trajectory = Path(config.trajectory)
    if config.topology:
        fs.topology = Path(config.topology)
    validation = validate_fileset(fs)

    system = None
    selected, skipped = [], []
    if validation.ok:
        topo = str(config.topology or fs.topology or fs.structure)
        system = detect_system(topo)
        if config.system_type:
            system.system_type = SystemType(config.system_type)
        registry.ensure_builtins_loaded()
        registry.discover_plugins(config.plugin_dirs)
        requested = config.analyses
        selected, skipped = registry.select(
            system, fs.available_keys(), requested=requested, run_all=config.run_all)
        if config.exclude:
            selected = [c for c in selected if c.name not in config.exclude]
    return fs, validation, system, selected, skipped


def run_pipeline(config: RunConfig | None = None, *, ask=None,
                 plan_only: bool = False, log=print, **kwargs):
    """
    Run the full pipeline.

    Parameters
    ----------
    config : a RunConfig, or None to build one from kwargs (e.g. input_dir=...).
    ask : optional callable(question, options)->index for interactive prompts.
    plan_only : if True, resolve and return the plan without executing.
    """
    if config is None:
        config = RunConfig(**kwargs)

    fs, validation, system, selected, skipped = build_plan(config, ask=ask)

    log("=" * 74)
    log(" MolDynX Tools -- MD analysis pipeline")
    log("=" * 74)
    log(f" input : {config.input_dir}")
    log(f" output: {config.output_dir}")
    if not validation.ok:
        log("\n[input validation FAILED]\n" + validation.report())
        raise SystemExit(2)
    for w in validation.warnings:
        log(f" [warn] {w}")
    log("\n" + system.summary())
    log(f"\n Selected {len(selected)} analyses: "
        f"{', '.join(c.name for c in selected)}")
    if skipped:
        log(" Skipped:")
        for name, reason in skipped:
            log(f"   - {name}: {reason}")

    if plan_only:
        return {"fileset": fs, "system": system,
                "selected": [c.name for c in selected], "skipped": skipped}

    config.ensure_dirs()
    manifest = RunManifest()
    manifest.set_config(config)
    manifest.set_system(system)
    manifest.set_inputs(fs)

    ctx = AnalysisContext(config, system, fs, provenance=manifest)
    results: dict[str, dict] = {}
    for cls in selected:
        log("-" * 74)
        log(f">>> {cls.name}: {cls.label}")
        t0 = time.time()
        try:
            summary = cls().run(ctx)
            dt = time.time() - t0
            results[cls.name] = summary
            manifest.record_analysis(cls.name, "ok", dt, summary=summary)
            log(f"[ok] {cls.name} in {dt:.1f}s")
        except Exception:
            dt = time.time() - t0
            err = traceback.format_exc()
            manifest.record_analysis(cls.name, "error", dt,
                                     error=err.splitlines()[-1])
            results[cls.name] = {"error": err.splitlines()[-1]}
            log(f"[ERROR] {cls.name}:\n{err}")

    manifest_paths = manifest.write(config.output_dir)
    log(f"\n[manifest] {manifest_paths[0]}")

    # -- report ----------------------------------------------------------- #
    try:
        from moldynx.report.generator import generate_report
        report_paths = generate_report(ctx, results, manifest,
                                       formats=config.report_formats)
        for p in report_paths:
            log(f"[report] {p}")
    except Exception:
        log("[report] generation failed:\n" + traceback.format_exc())

    log("=" * 74)
    log(" Done.")
    return {"system": system, "results": results, "manifest": manifest.data,
            "output_dir": str(config.output_dir)}
