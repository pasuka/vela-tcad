#!/usr/bin/env python
"""Emit a PN2D BV Jacobian block-audit CSV.

When a probe executable is available, this script delegates to it and records
finite analytic-vs-finite-difference Jacobian block norms. Without a probe it
falls back to the CSV contract with ``nan`` norms so downstream tooling can keep
the same shape while the C++ probe is unavailable.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


FIELDS = [
    "bias_V",
    "block",
    "analytic_norm",
    "fd_norm",
    "diff_norm",
    "rel_diff",
]

BLOCKS = (
    "srh_auger",
    "sg_avalanche",
    "transport",
    "poisson",
    "dirichlet_or_gauge",
)


def write_placeholder_report(output: Path, bias_values: list[float]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for bias in bias_values:
            for block in BLOCKS:
                writer.writerow(
                    {
                        "bias_V": f"{bias:.6g}",
                        "block": block,
                        "analytic_norm": "nan",
                        "fd_norm": "nan",
                        "diff_norm": "nan",
                        "rel_diff": "nan",
                    }
                )


def default_probe_candidates() -> list[Path]:
    repo = Path(__file__).resolve().parents[1]
    exe = ".exe" if sys.platform.startswith("win") else ""
    return [
        repo / "build-release" / f"pn2d_jacobian_block_audit{exe}",
        repo / "build" / f"pn2d_jacobian_block_audit{exe}",
    ]


def find_default_probe() -> Path | None:
    for candidate in default_probe_candidates():
        if candidate.is_file():
            return candidate
    return None


def write_probe_report(output: Path,
                       bias_values: list[float],
                       probe_command: list[str | Path]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [str(item) for item in probe_command]
    command.extend(["--output", str(output)])
    for bias in bias_values:
        command.extend(["--bias", f"{bias:.17g}"])
    subprocess.run(command, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write the PN2D BV Jacobian block-audit CSV."
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bias", action="append", type=float, required=True)
    parser.add_argument(
        "--probe-exe",
        type=Path,
        default=None,
        help="Optional pn2d_jacobian_block_audit executable. If omitted, known build paths are probed.",
    )
    parser.add_argument(
        "--allow-placeholder",
        action="store_true",
        help="Write nan placeholder rows when no probe executable is found.",
    )
    args = parser.parse_args()

    probe = args.probe_exe if args.probe_exe is not None else find_default_probe()
    if probe is not None:
        write_probe_report(args.output, args.bias, [probe])
    elif args.allow_placeholder:
        write_placeholder_report(args.output, args.bias)
    else:
        raise SystemExit(
            "pn2d_jacobian_block_audit probe executable not found; build it or pass "
            "--allow-placeholder for contract-only output."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
