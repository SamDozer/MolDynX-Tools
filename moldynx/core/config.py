"""
Run configuration: merges command-line arguments with an optional YAML file.

Precedence (highest first): explicit CLI flags > YAML file > built-in defaults.
A single ``config.yaml`` can therefore describe an entire reproducible workflow
(inputs, selections, per-analysis parameters, report formats) and be re-run with
``moldynx analyze --config config.yaml``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path


@dataclass
class RunConfig:
    # -- inputs / outputs ------------------------------------------------- #
    input_dir: Path | None = None
    output_dir: Path = Path("moldynx_results")
    # explicit overrides for auto-discovered files (optional)
    trajectory: Path | None = None
    topology: Path | None = None
    energy: Path | None = None
    index: Path | None = None
    gmx_top: Path | None = None

    # -- system / selections --------------------------------------------- #
    system_type: str | None = None            # override auto-detection
    chains: str | None = None                 # e.g. "A", "A+B"
    ligand: str | None = None                 # ligand resname override

    # -- trajectory slicing / performance -------------------------------- #
    start: int = 0
    end: int | None = None
    stride: int = 1
    threads: int = 1
    core_extract: bool = True                 # cache a solute-only trajectory

    # -- analysis selection ---------------------------------------------- #
    analyses: list[str] | None = None         # None => auto-select for system
    run_all: bool = False
    exclude: list[str] = field(default_factory=list)
    plugin_dirs: list[str] = field(default_factory=list)
    params: dict = field(default_factory=dict)  # {analysis_name: {param: value}}

    # -- reporting / repro ----------------------------------------------- #
    report_formats: list[str] = field(default_factory=lambda: ["md", "html"])
    dpi: int = 300
    seed: int = 0
    interactive: bool = False

    def __post_init__(self):
        for k in ("input_dir", "output_dir", "trajectory", "topology",
                  "energy", "index", "gmx_top"):
            v = getattr(self, k)
            if v is not None:
                setattr(self, k, Path(v))

    # -- output sub-directories ------------------------------------------ #
    @property
    def data_dir(self) -> Path: return self.output_dir / "data"
    @property
    def results_dir(self) -> Path: return self.output_dir / "results"
    @property
    def figures_dir(self) -> Path: return self.output_dir / "figures"
    @property
    def tables_dir(self) -> Path: return self.output_dir / "tables"
    @property
    def report_dir(self) -> Path: return self.output_dir / "report"

    def ensure_dirs(self) -> None:
        for d in (self.output_dir, self.data_dir, self.results_dir,
                  self.figures_dir, self.tables_dir, self.report_dir):
            d.mkdir(parents=True, exist_ok=True)

    def params_for(self, name: str) -> dict:
        return dict(self.params.get(name, {}))

    def to_dict(self) -> dict:
        d = asdict(self)
        for k, v in d.items():
            if isinstance(v, Path):
                d[k] = str(v)
        return d

    # -- construction ---------------------------------------------------- #
    @classmethod
    def from_yaml(cls, path: str | Path) -> "RunConfig":
        import yaml
        data = yaml.safe_load(Path(path).read_text()) or {}
        known = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in data.items() if k in known})

    @classmethod
    def from_args(cls, args, yaml_path: str | Path | None = None) -> "RunConfig":
        """Build from an argparse namespace, layering YAML underneath if given."""
        base = cls.from_yaml(yaml_path) if yaml_path else cls()
        mapping = {
            "input_dir": "input", "output_dir": "output", "trajectory": "traj",
            "topology": "top", "energy": "edr", "index": "ndx",
            "system_type": "system_type", "chains": "chains", "ligand": "ligand",
            "start": "start", "end": "end", "stride": "stride", "threads": "threads",
            "run_all": "all", "interactive": "interactive", "seed": "seed",
        }
        for attr, argname in mapping.items():
            val = getattr(args, argname, None)
            if val is not None:
                setattr(base, attr, val)
        if getattr(args, "analyses", None):
            base.analyses = [a.strip() for a in args.analyses.split(",")]
        if getattr(args, "exclude", None):
            base.exclude = [a.strip() for a in args.exclude.split(",")]
        if getattr(args, "plugin_dir", None):
            base.plugin_dirs = list(base.plugin_dirs) + [args.plugin_dir]
        if getattr(args, "report", None):
            base.report_formats = [f.strip() for f in args.report.split(",")]
        base.__post_init__()
        return base
