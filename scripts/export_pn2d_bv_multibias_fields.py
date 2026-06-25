#!/usr/bin/env python3
"""Export target PN2D BV Sentaurus multibias TDR snapshots to field CSV directories."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_DIR = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "sentaurus_vm_runs" / "pn2d_bv_validation_refresh_0p05v" / "source"
DEFAULT_OUT_ROOT = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "sentaurus_multibias_0p25v"
DEFAULT_IMPORTER = REPO / "build-release" / "sentaurus_import.exe"


def target_biases(start: float = -16.0, stop: float = -20.0, step: float = -0.25) -> list[float]:
    if step == 0.0:
        raise ValueError("step must be non-zero")
    values: list[float] = []
    value = start
    if step < 0.0:
        while value >= stop - 1.0e-12:
            values.append(round(value, 12))
            value += step
    else:
        while value <= stop + 1.0e-12:
            values.append(round(value, 12))
            value += step
    if not values:
        raise ValueError("empty target bias list")
    return values


def tdr_index_for_bias(bias: float, final_bias: float = -20.0, intervals: int = 400) -> int:
    if intervals <= 0:
        raise ValueError("intervals must be positive")
    if final_bias == 0.0:
        raise ValueError("final_bias must be non-zero")
    index = int(round((bias / final_bias) * intervals))
    if index < 0 or index > intervals:
        raise ValueError(f"bias {bias:g} V maps outside 0..{intervals}")
    reconstructed = final_bias * index / intervals
    if abs(reconstructed - bias) > 1.0e-8:
        raise ValueError(
            f"bias {bias:g} V is not represented by {intervals} intervals to {final_bias:g} V; "
            f"nearest index {index} is {reconstructed:g} V"
        )
    return index


def signed_bias_token(value: float) -> str:
    text = f"{value:g}".replace("+", "").replace(".", "p")
    return text


def export_dir_name(bias: float) -> str:
    return f"sentaurus_{signed_bias_token(bias)}v"


def run_checked(argv: Sequence[str], dry_run: bool) -> None:
    if dry_run:
        return
    subprocess.run(list(argv), check=True)


def export_targets(
    source_dir: Path,
    out_root: Path,
    importer: Path,
    biases: list[float],
    final_bias: float,
    intervals: int,
    clean: bool,
    dry_run: bool,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    out_root.mkdir(parents=True, exist_ok=True)
    for bias in biases:
        index = tdr_index_for_bias(bias, final_bias, intervals)
        tdr = source_dir / f"pn2d_bv_multibias_{index:04d}_des.tdr"
        export_dir = out_root / export_dir_name(bias)
        if not tdr.is_file():
            raise FileNotFoundError(f"missing multibias TDR for {bias:g} V: {tdr}")
        if clean and export_dir.exists() and not dry_run:
            shutil.rmtree(export_dir)
        command = [
            str(importer),
            "--tdr",
            str(tdr),
            "--export-dir",
            str(export_dir),
            "--compensated-doping-policy",
            "reported",
        ]
        run_checked(command, dry_run)
        rows.append({
            "bias_V": bias,
            "tdr_index": index,
            "tdr": str(tdr),
            "export_dir": str(export_dir),
            "command": command,
            "dry_run": dry_run,
        })
    return rows


def write_summary(out_root: Path, rows: list[dict[str, Any]], intervals: int, final_bias: float) -> Path:
    summary = {
        "schema": "vela.pn2d_bv_multibias_field_export.v1",
        "intervals": intervals,
        "final_bias_V": final_bias,
        "count": len(rows),
        "biases_V": [row["bias_V"] for row in rows],
        "rows": rows,
    }
    path = out_root / "export_summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    md = out_root / "README.md"
    lines = [
        "# PN2D BV 0.25V Sentaurus Field Exports",
        "",
        f"Source sweep: 0 to {final_bias:g} V with {intervals} intervals.",
        f"Exported target biases: {len(rows)}.",
        "",
        "| bias V | TDR index | export dir |",
        "| ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(f"| {row['bias_V']} | {row['tdr_index']} | `{row['export_dir']}` |")
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--tdr-importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--start", type=float, default=-16.0)
    parser.add_argument("--stop", type=float, default=-20.0)
    parser.add_argument("--step", type=float, default=-0.25)
    parser.add_argument("--final-bias", type=float, default=-20.0)
    parser.add_argument("--intervals", type=int, default=400)
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    biases = target_biases(args.start, args.stop, args.step)
    rows = export_targets(
        args.source_dir.resolve(),
        args.out_root.resolve(),
        args.tdr_importer.resolve(),
        biases,
        args.final_bias,
        args.intervals,
        args.clean,
        args.dry_run,
    )
    summary_path = write_summary(args.out_root.resolve(), rows, args.intervals, args.final_bias)
    print(json.dumps({"summary": str(summary_path), "count": len(rows)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())