"""
Reproducibility manifest.

Captures everything another researcher needs to reproduce a run: software and
library versions, OS/CPU, the exact command, git commit, input-file fingerprints,
random seeds, all parameters, per-analysis status/runtime, and total wall-time.
Written as both JSON and YAML.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def _lib_version(mod: str) -> str | None:
    try:
        m = __import__(mod)
        return getattr(m, "__version__", None)
    except Exception:
        return None


def _git_commit() -> dict:
    try:
        root = Path(__file__).resolve().parents[2]
        commit = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
        dirty = bool(subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain"],
            stderr=subprocess.DEVNULL).decode().strip())
        return {"commit": commit, "dirty": dirty}
    except Exception:
        return {"commit": None, "dirty": None}


def file_fingerprint(path: str | Path, chunk: int = 1 << 20) -> dict:
    """
    Cheap but robust fingerprint for possibly-huge trajectories: size + SHA-256
    of the first and last ``chunk`` bytes (avoids hashing 100 GB in full).
    """
    path = Path(path)
    if not path.exists():
        return {"path": str(path), "exists": False}
    size = path.stat().st_size
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read(chunk))
        if size > 2 * chunk:
            fh.seek(-chunk, 2)
            h.update(fh.read(chunk))
    return {"path": str(path), "exists": True, "size_bytes": size,
            "sha256_ends": h.hexdigest()}


class RunManifest:
    """Collects provenance during a pipeline run and writes it out."""

    def __init__(self) -> None:
        self._t0 = time.time()
        self.data: dict = {
            "moldynx_version": _lib_version("moldynx"),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "command": " ".join(sys.argv),
            "python": sys.version.split()[0],
            "platform": {
                "system": platform.system(), "release": platform.release(),
                "machine": platform.machine(), "processor": platform.processor(),
            },
            "libraries": {m: _lib_version(m) for m in (
                "numpy", "scipy", "pandas", "matplotlib", "MDAnalysis",
                "mdtraj", "sklearn", "networkx", "prolif", "yaml")},
            "git": _git_commit(),
            "inputs": {},
            "system": {},
            "config": {},
            "seed": None,
            "analyses": {},
        }

    # -- population ------------------------------------------------------- #
    def set_config(self, cfg) -> None:
        self.data["config"] = cfg.to_dict()
        self.data["seed"] = cfg.seed

    def set_system(self, system) -> None:
        self.data["system"] = system.to_dict()

    def set_inputs(self, fileset) -> None:
        d = fileset.to_dict()
        fp = {}
        for role in ("trajectory", "topology", "energy"):
            p = getattr(fileset, role)
            if p is not None:
                fp[role] = file_fingerprint(p)
        self.data["inputs"] = {"resolved": d, "fingerprints": fp}

    def record_analysis(self, name: str, status: str, runtime_s: float,
                        summary: dict | None = None, error: str | None = None) -> None:
        self.data["analyses"][name] = {
            "status": status, "runtime_s": round(runtime_s, 2),
            "summary": summary, "error": error,
        }

    # -- output ----------------------------------------------------------- #
    def write(self, outdir: Path) -> list[Path]:
        self.data["total_runtime_s"] = round(time.time() - self._t0, 2)
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        written = []
        jp = outdir / "manifest.json"
        jp.write_text(json.dumps(self.data, indent=2, default=str), encoding="utf-8")
        written.append(jp)
        try:
            import yaml
            yp = outdir / "manifest.yaml"
            yp.write_text(yaml.safe_dump(self.data, sort_keys=False, allow_unicode=True),
                          encoding="utf-8")
            written.append(yp)
        except Exception:
            pass
        return written
