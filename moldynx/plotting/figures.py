"""Figure saving: always PNG + vector PDF at 300 dpi."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt


def save_figure(fig, path_no_ext: Path | str, dpi: int = 300,
                formats: tuple[str, ...] = ("png", "pdf")) -> list[Path]:
    """Save ``fig`` to every format at high resolution; returns written paths."""
    path_no_ext = Path(path_no_ext)
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    written = []
    for ext in formats:
        out = path_no_ext.with_suffix(f".{ext}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        written.append(out)
    plt.close(fig)
    return written
