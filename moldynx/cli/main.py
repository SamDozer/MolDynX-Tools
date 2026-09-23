"""
MolDynX Tools (moldynx) command-line entry point.

Subcommands
-----------
    moldynx analyze  --input DIR [--output DIR] [options]   run the pipeline
    moldynx detect   --input DIR                            detect system only
    moldynx list-analyses [--system TYPE]                   list available analyses
    moldynx version
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
    if cfg.output_dir is None:
        cfg.output_dir = cfg.input_dir / "moldynx_results"
    ask = ask_choice if args.interactive else None
    run_pipeline(cfg, ask=ask, plan_only=args.plan)
    return 0


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
    return {"analyze": _cmd_analyze, "detect": _cmd_detect,
            "list-analyses": _cmd_list}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
