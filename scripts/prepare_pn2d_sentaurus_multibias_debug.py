#!/usr/bin/env python3
"""Prepare a Sentaurus pn2d IV deck that writes multi-bias field snapshots."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


DEFAULT_INPUT = Path("reference_tcad/pn2d_sentaurus2018/source/pn2d_iv_sdevice.cmd")
DEFAULT_OUTPUT = Path(
    "build/reference_tcad/pn2d_sentaurus2018/sentaurus_debug/"
    "pn2d_iv_multibias_debug_sdevice.cmd"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Source Sentaurus IV .cmd deck.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Generated debug .cmd deck.")
    parser.add_argument(
        "--file-prefix",
        default="pn2d_iv_multibias",
        help="Sentaurus Plot(FilePrefix=...) value for multi-bias snapshots.",
    )
    parser.add_argument(
        "--time-points",
        default="0,0.25,0.3,0.5,0.8,1.0",
        help="Comma-separated quasistationary Time points in normalized sweep coordinates.",
    )
    return parser.parse_args()


def parse_time_points(raw: str) -> list[str]:
    points: list[str] = []
    for part in raw.split(","):
        value = part.strip()
        if not value:
            continue
        float(value)
        points.append(value)
    if not points:
        raise ValueError("--time-points must contain at least one numeric value")
    return points


def ensure_plot_variables(text: str, variables: list[str]) -> str:
    match = re.search(r"\bPlot\s*\{", text)
    if match is None:
        raise ValueError("Sentaurus deck does not contain a Plot block")

    block_start = match.end()
    depth = 1
    pos = block_start
    while pos < len(text) and depth > 0:
        if text[pos] == "{":
            depth += 1
        elif text[pos] == "}":
            depth -= 1
        pos += 1
    if depth != 0:
        raise ValueError("Sentaurus Plot block braces are unbalanced")

    block_end = pos - 1
    block = text[block_start:block_end]
    missing = [var for var in variables if not re.search(rf"\b{re.escape(var)}\b", block)]
    if not missing:
        return text

    insertion = "\n  # Mobility diagnostics\n\n" + "\n".join(f"  {var}" for var in missing) + "\n"
    return text[:block_end] + insertion + text[block_end:]


def insert_multibias_plot(text: str, file_prefix: str, time_points: list[str]) -> str:
    plot_line = (
        f'    Plot(FilePrefix="{file_prefix}" '
        f'Time=({" ; ".join(time_points).replace(" ; ", ";")}) NoOverWrite)'
    )
    if f'FilePrefix="{file_prefix}"' in text:
        return text

    pattern = re.compile(
        r"(?P<coupled>\n(?P<indent>\s*)Coupled\s*\{\s*\n\s*Poisson\s+Electron\s+Hole\s*\n\s*\})"
        r"(?P<trailer>\s*\n\s*\}\s*\n\s*\}\s*)$",
        re.MULTILINE,
    )
    match = pattern.search(text)
    if match is None:
        raise ValueError("Could not locate the final Quasistationary Coupled block")
    return (
        text[: match.end("coupled")]
        + "\n"
        + plot_line
        + text[match.start("trailer") :]
    )


def prepare_deck(text: str, file_prefix: str, time_points: list[str]) -> str:
    text = ensure_plot_variables(text, ["eMobility", "hMobility"])
    return insert_multibias_plot(text, file_prefix, time_points)


def main() -> int:
    args = parse_args()
    time_points = parse_time_points(args.time_points)
    generated = prepare_deck(args.input.read_text(), args.file_prefix, time_points)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(generated)
    summary = {
        "input": str(args.input),
        "output": str(args.output),
        "file_prefix": args.file_prefix,
        "time_points": time_points,
        "added_plot_variables": ["eMobility", "hMobility"],
        "note": (
            "Run this deck from a Sentaurus environment in the source directory, then import each "
            "generated TDR snapshot with build/sentaurus_import."
        ),
    }
    (args.output.parent / "pn2d_iv_multibias_debug_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
