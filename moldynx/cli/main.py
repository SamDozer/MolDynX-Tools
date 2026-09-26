"""
MolDynX Tools (moldynx) command-line entry point.

Subcommands
-----------
    moldynx intake   --input DIR [--output DIR]             which files, which run, what is missing
    moldynx analyze  --input DIR [--output DIR] [options]   run the pipeline
    moldynx detect   --input DIR                            detect system only
    moldynx list-analyses [--system TYPE]                   list available analyses
    moldynx version

MolDynX never writes into the simulation folder: outputs default to
./moldynx_results/<folder name>.
"""

from __future__ import annotations

import argparse
import sys

from moldynx import __version__


def _add_common(p):
    p.add_argument("--input", "-i", help="Simulation directory (searched recursively).")
    p.add_argument("--output", "-o", help="Output directory.")
    p.add_argument("--config", help="YAML config file (CLI flags override it).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="moldynx", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"MolDynX Tools (moldynx) {__version__}")
    sub = parser.add_subparsers(dest="command")

    a = sub.add_parser("analyze", help="Run the analysis pipeline.")
    _add_common(a)
    a.add_argument("--traj", help="Override the discovered trajectory.")
    a.add_argument("--top", help="Override the discovered topology.")
    a.add_argument("--edr", help="Override the discovered energy file.")
    a.add_argument("--ndx", help="Override the discovered index file.")
    a.add_argument("--system-type", dest="system_type", help="Override auto-detected system type.")
    a.add_argument("--chains", help="Chain selection, e.g. 'A' or 'A+B'.")
    a.add_argument("--ligand", help="Ligand resname override.")
    a.add_argument("--analyses", help="Comma-separated analyses to run (default: auto).")
    a.add_argument("--exclude", help="Comma-separated analyses to skip.")
    a.add_argument("--all", action="store_true", help="Run every applicable analysis.")
    a.add_argument("--stride", type=int, help="Use every Nth frame.")
    a.add_argument("--start", type=int, help="First frame.")
    a.add_argument("--end", type=int, help="Last frame.")
    a.add_argument("--threads", type=int, help="Worker threads (where supported).")
    a.add_argument("--seed", type=int, help="Random seed.")
    a.add_argument("--report", help="Report formats, comma-separated (md,html,pdf).")
    a.add_argument("--plugin-dir", help="Extra directory of plugin modules.")
    a.add_argument("--interactive", action="store_true", help="Ask questions when ambiguous.")
    a.add_argument("--plan", action="store_true", help="Dry-run: show what would run and why.")
    a.add_argument("--allow-ambiguous", dest="allow_ambiguous", action="store_true", default=None,
                   help="Proceed (with warnings) when the production run is ambiguous.")
    a.add_argument("--include-dir", dest="include_dir", action="append",
                   help="Folder never to be treated as derived output (repeatable).")
    a.add_argument("--pbc", choices=["auto", "none", "whole", "nojump"],
                   help="Periodic-boundary treatment of the solute trajectory (default auto).")

    it = sub.add_parser("intake", help="Report which files form the run and what is missing.")
    _add_common(it)
    it.add_argument("--no-deep", dest="deep", action="store_false",
                    help="Names/stages only; do not read trajectory, run-input or log headers.")
    it.add_argument("--detect", action="store_true",
                    help="Also load the run input and describe chains, components and ions.")
    it.add_argument("--allow-ambiguous", dest="allow_ambiguous", action="store_true")
    it.add_argument("--include-dir", dest="include_dir", action="append")

    be = sub.add_parser("binding-energy",
                        help="MM-GBSA/PBSA with gmx_MMPBSA: prepare, optionally run, analyse.")
    _add_common(be)
    be.add_argument("--execute", action="store_true",
                    help="Run the prepared gmx_MMPBSA script now (hours; Linux or WSL).")

    ds = sub.add_parser("dataset", help="Assemble a shareable, self-verifying dataset from a run.")
    ds.add_argument("--run", required=True, help="Output folder of 'moldynx analyze'.")
    ds.add_argument("--out", required=True, help="Dataset folder to create.")
    ds.add_argument("--include-raw", action="store_true",
                    help="Copy and SHA-256 the canonical raw files (can be many GB).")
    ds.add_argument("--title", help="Dataset title for the README.")
    vf = sub.add_parser("verify", help="Run a dataset's verify_dataset.py.")
    vf.add_argument("dataset")
    vf.add_argument("--full", action="store_true", help="Also hash the raw files.")
    pk = sub.add_parser("package", help="Zip a dataset, then extract and verify the zip.")
    pk.add_argument("dataset")
    pk.add_argument("--zip", help="Zip path (default: <dataset>.zip).")

    d = sub.add_parser("detect", help="Detect and print the system composition.")
    _add_common(d)

    ls = sub.add_parser("list-analyses", help="List registered analyses.")
    ls.add_argument("--system", help="Filter to a system type.")

    sub.add_parser("version", help="Print version.")
    return parser


def _cmd_analyze(args) -> int:
    from moldynx.core.config import RunConfig
    from moldynx.core.pipeline import run_pipeline
    from moldynx.cli.interactive import ask_choice
    cfg = RunConfig.from_args(args, yaml_path=args.config)
    if cfg.input_dir is None:
        print("error: --input (or 'input_dir' in --config) is required.")
        return 2
    if not args.output and "output_dir" not in _yaml_keys(args.config):
        cfg.output_dir = default_output_dir(cfg.input_dir)
    if args.include_dir:
        cfg.include_dirs = list(cfg.include_dirs) + list(args.include_dir)
    ask = ask_choice if args.interactive else None
    run_pipeline(cfg, ask=ask, plan_only=args.plan)
    return 0


def default_output_dir(input_dir) -> "Path":
    """``./moldynx_results/<input folder name>`` -- never inside the simulation folder."""
    from pathlib import Path
    return Path.cwd() / "moldynx_results" / Path(input_dir).resolve().name


def _yaml_keys(path) -> set:
    if not path:
        return set()
    import yaml
    from pathlib import Path
    return set((yaml.safe_load(Path(path).read_text()) or {}).keys())


def _cmd_binding_energy(args) -> int:
    """Prepare (or analyse) the gmx_MMPBSA calculation; --execute runs it in between."""
    import json
    import subprocess
    from moldynx.analysis.mmpbsa import MMPBSA
    from moldynx.core.config import RunConfig
    from moldynx.core.context import AnalysisContext
    from moldynx.core.pipeline import build_plan
    from moldynx.io.gromacs import windows_to_wsl
    cfg = RunConfig.from_args(args, yaml_path=args.config)
    if cfg.input_dir is None:
        print("error: --input (or 'input_dir' in --config) is required.")
        return 2
    if not args.output and "output_dir" not in _yaml_keys(args.config):
        cfg.output_dir = default_output_dir(cfg.input_dir)
    fs, val, system, _sel, _skip = build_plan(cfg)
    if not val.ok:
        print(val.report())
        return 2
    ctx = AnalysisContext(cfg, system, fs)
    out = MMPBSA().run(ctx)
    print(json.dumps({k: v for k, v in out.items() if k not in ("headline",)}, indent=1,
                     default=str)[:3000])
    if args.execute and out.get("status") == "prepared":
        d = out["directory"]
        cmd = (["wsl.exe", "--", "bash", "-c", f"cd '{windows_to_wsl(d)}' && bash run_mmpbsa.sh"]
               if sys.platform.startswith("win") else ["bash", "run_mmpbsa.sh"])
        print(f"[binding-energy] running gmx_MMPBSA in {d} (this takes hours) ...")
        rc = subprocess.call(cmd, cwd=None if sys.platform.startswith("win") else d)
        if rc != 0:
            print(f"[binding-energy] run script exited with {rc}; see {d}/run.log")
            return rc
        out = MMPBSA().run(ctx)
        print(json.dumps(out.get("headline"), indent=1, default=float))
    return 0


def _cmd_intake(args) -> int:
    from moldynx.io.intake import run_intake, write_intake
    if not args.input:
        print("error: --input is required.")
        return 2
    res = run_intake(args.input, deep=args.deep, detect=args.detect,
                     allow_ambiguous=args.allow_ambiguous, include_dirs=args.include_dir)
    out = args.output or default_output_dir(args.input) / "intake"
    paths = write_intake(res, out)
    fs, val = res.fileset, res.validation
    print(f"trajectory : {fs.trajectory}")
    print(f"run input  : {fs.topology}")
    print(f"evidence   : {fs.evidence.get('trajectory_choice', '—')}")
    print("\n" + val.capability_table())
    if val.errors or val.warnings:
        print("\n" + val.report())
    for p in paths:
        print(f"[intake] {p}")
    return 0 if val.ok else 2


def _cmd_detect(args) -> int:
    from moldynx.io.discovery import discover_files
    from moldynx.io.validation import validate_fileset
    from moldynx.core.system import detect_system
    if not args.input:
        print("error: --input is required.")
        return 2
    fs = discover_files(args.input)
    val = validate_fileset(fs)
    print("Discovered files:")
    for k, v in fs.to_dict().items():
        if v:
            print(f"  {k:11s}: {v}")
    print("\nValidation:\n" + val.report())
    if val.ok:
        topo = str(fs.topology or fs.structure)
        system = detect_system(topo)
        print("\n" + system.summary())
    return 0 if val.ok else 2


def _cmd_list(args) -> int:
    from moldynx.core.registry import registry
    from moldynx.core.system import SystemType
    registry.ensure_builtins_loaded()
    filt = SystemType(args.system) if args.system else None
    print(f"{'name':22s} {'category':13s} systems")
    print("-" * 70)
    for name in registry.names():
        cls = registry.get(name)
        if filt and "*" not in cls.supported_systems and filt not in cls.supported_systems:
            continue
        systems = "all" if "*" in cls.supported_systems else \
            ", ".join(sorted(s.value for s in cls.supported_systems))
        print(f"{cls.name:22s} {cls.category:13s} {systems}")
    return 0


def _force_utf8() -> None:
    """Make console output robust to non-ASCII (e.g. 'Cα') on Windows cp1252."""
    import sys
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


def main(argv=None) -> int:
    _force_utf8()
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command or args.command == "version":
        if not args.command:
            parser.print_help()
            return 1
        print(f"MolDynX Tools (moldynx) {__version__}")
        return 0
    if args.command == "dataset":
        from moldynx.dataset import assemble
        print(f"[dataset] {assemble(args.run, args.out, include_raw=args.include_raw, title=args.title)}")
        return 0
    if args.command == "verify":
        from moldynx.dataset import verify
        return verify(args.dataset, full=args.full)
    if args.command == "package":
        from pathlib import Path
        from moldynx.dataset import package
        r = package(args.dataset, args.zip or str(Path(args.dataset).resolve()) + ".zip")
        print(f"[package] {r['zip']}: {r['files']} files, {r['bytes'] / 1e6:.1f} MB, "
              f"CRC {'ok' if r['crc_ok'] else 'FAILED'}, verification exit {r['verify_exit']}")
        print(r["verify_output"])
        return 0 if r["crc_ok"] and r["verify_exit"] == 0 else 1
    return {"analyze": _cmd_analyze, "detect": _cmd_detect, "intake": _cmd_intake,
            "binding-energy": _cmd_binding_energy,
            "list-analyses": _cmd_list}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
