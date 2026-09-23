"""Minimal interactive prompt helpers (only used with --interactive)."""

from __future__ import annotations


def ask_choice(question: str, options: list[str]) -> int:
    """Ask the user to pick one option; returns its index (0-based)."""
    print(f"\n{question}")
    for i, opt in enumerate(options, 1):
        print(f"  [{i}] {opt}")
    while True:
        raw = input("Enter choice number: ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return int(raw) - 1
        print("  invalid choice, try again.")


def confirm(question: str, default: bool = True) -> bool:
    suffix = "[Y/n]" if default else "[y/N]"
    raw = input(f"{question} {suffix} ").strip().lower()
    if not raw:
        return default
    return raw.startswith("y")
