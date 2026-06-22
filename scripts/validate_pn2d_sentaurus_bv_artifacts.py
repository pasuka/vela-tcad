#!/usr/bin/env python3
"""Validate fetched PN2D Sentaurus BV artifacts.

This validator is intentionally lightweight: it checks the artifact contract
needed by the BV validation workflow without requiring Sentaurus, TDR parsing,
or a license. PLT bias parsing is best-effort so minimal synthetic fixtures can
exercise the structural checks.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Sequence


MULTIBIAS_RE = re.compile(r"^pn2d_bv_multibias_(?P<index>\d{4})_des\.tdr$")
LOG_PROBLEM_RE = re.compile(
    r"\b(fatal|aborted|segmentation fault|uncaught exception|failed|failure)\b",
    re.IGNORECASE,
)
LICENSE_PROBLEM_RE = re.compile(
    r"\b(no\s+license|license\s+(error|failed|failure|denied|unavailable)|checkout\s+failed)\b",
    re.IGNORECASE,
)
LOG_CANDIDATES = ("pn2d_bv.log", "pn2d_bv.log_des.log", "run_pn2d_bv.out")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--require-final-bias", type=float, default=-20.0)
    parser.add_argument("--expected-multibias-count", type=int, default=201)
    return parser.parse_args(argv)


def fail(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


def check_required_file(source_dir: Path, name: str) -> Path:
    path = source_dir / name
    if not path.is_file():
        raise FileNotFoundError(f"missing required BV artifact: {path}")
    return path


def check_required_log(source_dir: Path) -> Path:
    for name in LOG_CANDIDATES:
        path = source_dir / name
        if path.is_file():
            return path
    raise FileNotFoundError(
        "missing required BV log artifact: expected one of " + ", ".join(LOG_CANDIDATES)
    )


def check_log_clean(log_path: Path) -> list[str]:
    warnings: list[str] = []
    for line_no, line in enumerate(log_path.read_text(errors="replace").splitlines(), 1):
        lowered = line.lower()
        if LICENSE_PROBLEM_RE.search(line):
            raise ValueError(
                f"{log_path} contains license problem marker on line {line_no}: {line.strip()}"
            )
        if LOG_PROBLEM_RE.search(line):
            raise ValueError(
                f"{log_path} contains log problem marker on line {line_no}: {line.strip()}"
            )
    if log_path.stat().st_size == 0:
        warnings.append(f"{log_path.name} is empty")
    return warnings


def discover_multibias_tdrs(source_dir: Path) -> tuple[list[Path], list[int]]:
    indexed: list[tuple[int, Path]] = []
    for path in source_dir.glob("pn2d_bv_multibias_*_des.tdr"):
        match = MULTIBIAS_RE.match(path.name)
        if match is None:
            continue
        indexed.append((int(match.group("index")), path))
    indexed.sort()
    return [path for _, path in indexed], [index for index, _ in indexed]


def check_multibias_set(source_dir: Path, expected_count: int) -> tuple[list[Path], list[int]]:
    files, indexes = discover_multibias_tdrs(source_dir)
    if len(files) != expected_count:
        raise ValueError(
            f"expected {expected_count} pn2d_bv_multibias TDR files, found {len(files)}"
        )
    if not indexes:
        raise ValueError("no pn2d_bv_multibias TDR files found")
    expected_last = expected_count - 1
    if indexes[0] != 0 or indexes[-1] != expected_last:
        raise ValueError(
            f"expected pn2d_bv_multibias indexes 0000..{expected_last:04d}, "
            f"found {indexes[0]:04d}..{indexes[-1]:04d}"
        )
    missing = sorted(set(range(expected_count)) - set(indexes))
    if missing:
        preview = ", ".join(f"{index:04d}" for index in missing[:8])
        raise ValueError(f"missing pn2d_bv_multibias index(es): {preview}")
    return files, indexes


def parse_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def parse_csv_biases(plt_path: Path) -> list[float] | None:
    text = plt_path.read_text(errors="replace")
    if "," not in text:
        return None
    with plt_path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return None
        bias_column = next((name for name in reader.fieldnames if name.strip() == "Anode OuterVoltage"), None)
        if bias_column is None:
            return None
        biases: list[float] = []
        for row in reader:
            value = parse_float((row.get(bias_column) or "").strip())
            if value is not None:
                biases.append(value)
    return biases


def parse_ise_text_biases(plt_path: Path) -> list[float] | None:
    text = plt_path.read_text(errors="replace")
    datasets_match = re.search(r"datasets\s*=\s*\[(?P<body>.*?)\]", text, re.DOTALL)
    data_match = re.search(r"Data\s*\{(?P<body>.*?)\}", text, re.DOTALL)
    if datasets_match is None or data_match is None:
        return None
    datasets = re.findall(r'"([^"]+)"', datasets_match.group("body"))
    try:
        bias_index = datasets.index("Anode OuterVoltage")
    except ValueError:
        return None
    values = [
        float(match.group(0))
        for match in re.finditer(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?", data_match.group("body"))
    ]
    width = len(datasets)
    if width <= 0 or len(values) < width:
        return None
    return [values[offset + bias_index] for offset in range(0, len(values) - width + 1, width)]


def check_plt_bias(plt_path: Path, required_bias: float) -> list[str]:
    warnings: list[str] = []
    biases = parse_csv_biases(plt_path)
    if biases is None:
        biases = parse_ise_text_biases(plt_path)
    if not biases:
        return [f"could not parse bias values from {plt_path.name}; skipped final-bias check"]

    tolerance = 1.0e-6
    if min(abs(value - required_bias) for value in biases) > tolerance:
        raise ValueError(
            f"{plt_path.name} does not reach required bias {required_bias:g} V "
            f"within {tolerance:g} V; parsed range {min(biases):g}..{max(biases):g} V"
        )
    return warnings


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    source_dir = args.source_dir.resolve()
    try:
        if not source_dir.is_dir():
            return fail(f"source directory does not exist: {source_dir}")
        plt_path = check_required_file(source_dir, "pn2d_bv.plt")
        log_path = check_required_log(source_dir)
        warnings = []
        warnings.extend(check_log_clean(log_path))
        files, indexes = check_multibias_set(source_dir, args.expected_multibias_count)
        warnings.extend(check_plt_bias(plt_path, args.require_final_bias))
    except Exception as exc:  # noqa: BLE001 - CLI should report clean failures.
        return fail(str(exc))

    print(f"source_dir: {source_dir}")
    print(f"plt: {plt_path.name}")
    print(f"log: {log_path.name}")
    print(
        f"multibias_tdr_count: {len(files)} "
        f"({files[0].name}..{files[-1].name}; indexes {indexes[0]:04d}..{indexes[-1]:04d})"
    )
    print(f"required_final_bias_V: {args.require_final_bias:g}")
    for warning in warnings:
        print(f"warning: {warning}")
    print("status: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
