#!/usr/bin/env python3
"""Import Sentaurus text artifacts into Vela-friendly reference files.

TDR/HDF5 import is handled by the C++ sentaurus_import tool. This script covers
the text pieces around it: DF-ISE .plt curves and SDevice .cmd summaries.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any


UNSUPPORTED_PHYSICS = [
    "Thermodynamic",
    "IALMob",
    "Trap",
    "Avalanche",
    "eTemperature",
    "hTemperature",
]


def parse_quoted_list(text: str, key: str) -> list[str]:
    match = re.search(rf"{re.escape(key)}\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        return []
    return re.findall(r'"([^"]+)"', match.group(1))


def parse_values_block(text: str, column_count: int) -> list[list[float]]:
    match = re.search(r"Values\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        match = re.search(r"Data\s*\{(.*?)\}", text, re.DOTALL | re.IGNORECASE)
    if not match:
        raise ValueError("Sentaurus .plt file does not contain a Data block")

    if "Values" not in match.group(0):
        numbers = [float(value) for value in re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", match.group(1))]
        if column_count <= 0 or len(numbers) % column_count != 0:
            raise ValueError(
                f"Sentaurus .plt Data block has {len(numbers)} values, "
                f"not a multiple of {column_count} datasets"
            )
        return [
            numbers[index:index + column_count]
            for index in range(0, len(numbers), column_count)
        ]

    rows: list[list[float]] = []
    for raw in match.group(1).splitlines():
        line = raw.strip().rstrip(",")
        if not line:
            continue
        rows.append([float(part) for part in line.split()])
    return rows


def format_float(value: float) -> str:
    return f"{value:.12g}"


def export_plt_reference(input_path: Path,
                         output_path: Path,
                         bias_column: str,
                         current_column: str) -> None:
    text = input_path.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    if not datasets:
        raise ValueError("Sentaurus .plt file does not declare datasets")
    try:
        bias_index = datasets.index(bias_column)
        current_index = datasets.index(current_column)
    except ValueError as exc:
        raise ValueError(f"requested column is missing from .plt datasets: {exc}") from exc

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["bias_V", "current_total"])
        for row in parse_values_block(text, len(datasets)):
            writer.writerow([
                format_float(row[bias_index]),
                format_float(row[current_index]),
            ])


def remove_comments(text: str) -> str:
    cleaned_lines = []
    for raw in text.splitlines():
        line = raw.split("#", maxsplit=1)[0]
        line = line.split("//", maxsplit=1)[0]
        cleaned_lines.append(line)
    return "\n".join(cleaned_lines)


def parse_file_block(text: str) -> dict[str, str]:
    match = re.search(r"File\s*\{(.*?)\}", text, re.DOTALL | re.IGNORECASE)
    if not match:
        return {}
    result: dict[str, str] = {}
    for key, value in re.findall(r"(\w+)\s*=\s*\"([^\"]+)\"", match.group(1)):
        result[key] = value
    return result


def parse_electrodes(text: str) -> list[dict[str, float | str]]:
    match = re.search(r"Electrode\s*\{(.*?)\n\}", text, re.DOTALL | re.IGNORECASE)
    if not match:
        return []
    electrodes = []
    for entry in re.findall(r"\{([^{}]*Name\s*=\s*\"[^\"]+\"[^{}]*)\}", match.group(1), re.DOTALL):
        name_match = re.search(r"Name\s*=\s*\"([^\"]+)\"", entry)
        voltage_match = re.search(r"Voltage\s*=\s*([-+0-9.eE]+)", entry)
        if name_match:
            electrodes.append({
                "name": name_match.group(1),
                "voltage": float(voltage_match.group(1)) if voltage_match else 0.0,
            })
    return electrodes


def parse_sweeps(text: str) -> list[dict[str, Any]]:
    sweeps = []
    pattern = re.compile(
        r"Quasistationary\s*\((.*?)\)\s*\{",
        re.DOTALL | re.IGNORECASE,
    )
    for match in pattern.finditer(text):
        body = match.group(1)
        goal = re.search(
            r"Goal\s*\{[^{}]*Name\s*=\s*\"([^\"]+)\"[^{}]*Voltage\s*=\s*([-+0-9.eE]+)",
            body,
            re.DOTALL | re.IGNORECASE,
        )
        if goal:
            sweeps.append({
                "mode": "quasistationary",
                "contact": goal.group(1),
                "stop": float(goal.group(2)),
            })
    return sweeps


def parse_unsupported_physics(text: str) -> list[str]:
    found = []
    for token in UNSUPPORTED_PHYSICS:
        if re.search(rf"\b{re.escape(token)}\b", text):
            found.append(token)
    return found


def parse_cmd(input_path: Path) -> dict[str, Any]:
    text = remove_comments(input_path.read_text(errors="ignore"))
    return {
        "files": parse_file_block(text),
        "electrodes": parse_electrodes(text),
        "sweeps": parse_sweeps(text),
        "unsupported_physics": parse_unsupported_physics(text),
    }


def write_python_runner(summary: dict[str, Any],
                        output_path: Path,
                        mesh_json: str,
                        output_csv: str) -> None:
    electrodes = {str(item["name"]): float(item["voltage"]) for item in summary.get("electrodes", [])}
    sweeps = summary.get("sweeps", [])
    final_sweep = sweeps[-1] if sweeps else {"contact": "drain", "stop": 0.0}
    contacts = [
        {"name": name, "bias": voltage}
        for name, voltage in electrodes.items()
    ]
    deck = {
        "simulation_type": "dc_sweep",
        "mesh_file": mesh_json,
        "output_csv": output_csv,
        "scaling": {"mode": "unit_scaling"},
        "contacts": contacts,
        "doping": [],
        "solver": {
            "method": "gummel",
            "max_iter": 100,
            "reltol": 1.0e-5,
        },
        "sweep": {
            "mode": "iv",
            "contact": final_sweep["contact"],
            "current_contact": final_sweep["contact"],
            "start": electrodes.get(str(final_sweep["contact"]), 0.0),
            "stop": final_sweep["stop"],
            "step": final_sweep["stop"] / 80.0 if final_sweep["stop"] != 0.0 else 1.0,
            "write_vtk": False,
        },
        "sentaurus_import": {
            "unsupported_physics": summary.get("unsupported_physics", []),
            "sweeps": sweeps,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "from vela.curves import run_iv_curve\n\n"
        "CONFIG = "
        + json.dumps(deck, indent=2)
        + "\n\n"
        "if __name__ == \"__main__\":\n"
        "    run_iv_curve(CONFIG)\n"
    )


def cmd_command(args: argparse.Namespace) -> None:
    summary = parse_cmd(args.input)
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    if args.python_runner:
        write_python_runner(summary, args.python_runner, args.mesh_json, args.output_csv)


def plt_command(args: argparse.Namespace) -> None:
    export_plt_reference(args.input, args.output, args.bias_column, args.current_column)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    plt = sub.add_parser("plt", help="Convert a Sentaurus DF-ISE .plt curve to reference CSV")
    plt.add_argument("--input", type=Path, required=True)
    plt.add_argument("--output", type=Path, required=True)
    plt.add_argument("--bias-column", required=True)
    plt.add_argument("--current-column", required=True)
    plt.set_defaults(func=plt_command)

    cmd = sub.add_parser("cmd", help="Summarize an SDevice .cmd file and generate a Python runner")
    cmd.add_argument("--input", type=Path, required=True)
    cmd.add_argument("--summary-json", type=Path)
    cmd.add_argument("--python-runner", type=Path)
    cmd.add_argument("--mesh-json", default="mesh.json")
    cmd.add_argument("--output-csv", default="outputs/sentaurus_import.csv")
    cmd.set_defaults(func=cmd_command)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
