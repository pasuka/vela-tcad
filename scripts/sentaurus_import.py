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
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any


UNSUPPORTED_PHYSICS = [
    "Thermodynamic",
    "IALMob",
    "Trap",
    "Traps",
    "Avalanche",
    "eTemperature",
    "hTemperature",
]

REPO = Path(__file__).resolve().parents[1]
NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")
PLACEHOLDER_RE = re.compile(r"@([^@\s]+)@")


def coerce_scalar(value: str) -> float | str:
    stripped = value.strip().strip('"')
    if NUMBER_RE.fullmatch(stripped):
        return float(stripped)
    return stripped


def numeric_or_none(value: object) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str) and NUMBER_RE.fullmatch(value):
        return float(value)
    return None


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


def parse_template_vars(items: list[str] | None) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"template variable must be KEY=VALUE: {item}")
        key, value = item.split("=", maxsplit=1)
        if not key:
            raise ValueError(f"template variable key is empty: {item}")
        result[key] = value
    return result


def apply_template_variables(text: str, variables: dict[str, str]) -> str:
    def replace(match: re.Match[str]) -> str:
        key = match.group(1)
        return variables.get(key, match.group(0))

    return PLACEHOLDER_RE.sub(replace, text)


def find_matching(text: str, open_index: int, open_char: str, close_char: str) -> int:
    depth = 0
    for index in range(open_index, len(text)):
        char = text[index]
        if char == open_char:
            depth += 1
        elif char == close_char:
            depth -= 1
            if depth == 0:
                return index
    raise ValueError(f"unbalanced {open_char}{close_char} block")


def iter_named_blocks(text: str, name: str) -> list[tuple[str, str, int, int]]:
    blocks: list[tuple[str, str, int, int]] = []
    pattern = re.compile(rf"\b{re.escape(name)}\b\s*(\([^{{]*?\))?\s*\{{", re.DOTALL | re.IGNORECASE)
    for match in pattern.finditer(text):
        open_index = match.end() - 1
        close_index = find_matching(text, open_index, "{", "}")
        qualifier = (match.group(1) or "").strip()
        blocks.append((qualifier, text[open_index + 1:close_index], match.start(), close_index + 1))
    return blocks


def first_block(text: str, name: str) -> tuple[str, str, int, int] | None:
    blocks = iter_named_blocks(text, name)
    return blocks[0] if blocks else None


def parse_assignments(text: str) -> dict[str, float | str]:
    result: dict[str, float | str] = {}
    pattern = re.compile(
        r"([A-Za-z_][\w]*(?:\([^)]*\))?)\s*=\s*(\"[^\"]*\"|@\S+?@|[-+.\w|]+)",
        re.IGNORECASE,
    )
    for key, value in pattern.findall(text):
        result[key] = coerce_scalar(value)
    return result


def assignment_spans(text: str) -> list[tuple[int, int]]:
    pattern = re.compile(
        r"[A-Za-z_][\w]*(?:\([^)]*\))?\s*=\s*(?:\"[^\"]*\"|@\S+?@|[-+.\w|]+)",
        re.IGNORECASE,
    )
    return [match.span() for match in pattern.finditer(text)]


def remove_spans(text: str, spans: list[tuple[int, int]]) -> str:
    chars = list(text)
    for start, end in spans:
        for index in range(start, end):
            chars[index] = " "
    return "".join(chars)


def parse_file_block(text: str) -> dict[str, str]:
    block = first_block(text, "File")
    if not block:
        return {}
    result: dict[str, str] = {}
    for key, value in re.findall(r"(\w+)\s*=\s*\"([^\"]+)\"", block[1]):
        result[key] = value
    return result


def parse_electrodes(text: str) -> list[dict[str, float | str]]:
    block = first_block(text, "Electrode")
    if not block:
        return []
    electrodes = []
    for entry in re.findall(r"\{([^{}]*Name\s*=\s*\"[^\"]+\"[^{}]*)\}", block[1], re.DOTALL):
        name_match = re.search(r"Name\s*=\s*\"([^\"]+)\"", entry)
        voltage_match = re.search(r"Voltage\s*=\s*([-+0-9.eE]+)", entry)
        if name_match:
            electrodes.append({
                "name": name_match.group(1),
                "voltage": float(voltage_match.group(1)) if voltage_match else 0.0,
            })
    return electrodes


def parse_thermodes(text: str) -> list[dict[str, float | str]]:
    block = first_block(text, "Thermode")
    if not block:
        return []
    thermodes = []
    for entry in re.findall(r"\{([^{}]*Name\s*=\s*\"[^\"]+\"[^{}]*)\}", block[1], re.DOTALL):
        values = parse_assignments(entry)
        name = values.pop("Name", None)
        if name is None:
            continue
        item: dict[str, float | str] = {"name": str(name)}
        if "Temperature" in values:
            item["temperature"] = values["Temperature"]
        if "SurfaceResistance" in values:
            item["surface_resistance"] = values["SurfaceResistance"]
        thermodes.append(item)
    return thermodes


def parse_scope(qualifier: str) -> dict[str, Any]:
    if not qualifier:
        return {"kind": "global"}
    inner = qualifier.strip()[1:-1] if qualifier.startswith("(") and qualifier.endswith(")") else qualifier
    values = parse_assignments(inner)
    for key, kind in [("MaterialInterface", "material_interface"), ("Material", "material")]:
        if key in values:
            return {"kind": kind, "name": values[key]}
    return {"kind": "qualified", "raw": inner.strip(), "parameters": values}


def parse_models(body: str) -> list[str]:
    without_assignments = remove_spans(body, assignment_spans(body))
    tokens = re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", without_assignments)
    ignored = {"Conc", "EnergyMid", "BondConc", "ActEnergy"}
    return sorted(dict.fromkeys(token for token in tokens if token not in ignored))


def parse_physics_blocks(text: str) -> list[dict[str, Any]]:
    result = []
    for qualifier, body, _start, _end in iter_named_blocks(text, "Physics"):
        result.append({
            "scope": parse_scope(qualifier),
            "parameters": parse_assignments(body),
            "models": parse_models(body),
        })
    return result


def parse_plot_fields(text: str) -> list[dict[str, Any]]:
    block = first_block(text, "Plot")
    if not block:
        return []
    fields = []
    for token in re.findall(r"[A-Za-z_][A-Za-z0-9_]*(?:/Vector)?", block[1]):
        name, sep, suffix = token.partition("/")
        fields.append({"name": name, "vector": bool(sep and suffix.lower() == "vector")})
    return fields


def parse_math_block(text: str) -> dict[str, Any]:
    block = first_block(text, "Math")
    if not block:
        return {"parameters": {}, "flags": []}
    body = block[1]
    parameters = parse_assignments(body)
    remainder = remove_spans(body, assignment_spans(body))
    flags = sorted(dict.fromkeys(re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", remainder)))
    return {"parameters": parameters, "flags": flags}


def parse_coupled(text: str) -> list[dict[str, Any]]:
    result = []
    pattern = re.compile(r"Coupled\s*(\([^)]*\))?\s*\{([^{}]*)\}", re.DOTALL | re.IGNORECASE)
    for match in pattern.finditer(text):
        equations = re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", match.group(2))
        result.append({
            "type": "Coupled",
            "parameters": parse_assignments(match.group(1) or ""),
            "equations": equations,
            "start": match.start(),
            "end": match.end(),
        })
    return result


def parse_current_plot(text: str) -> dict[str, Any] | None:
    match = re.search(r"\bCurrentPlot\b\s*\(", text, re.IGNORECASE)
    if not match:
        return None
    open_index = match.end() - 1
    close_index = find_matching(text, open_index, "(", ")")
    body = text[open_index + 1:close_index]
    result: dict[str, Any] = {"raw": " ".join(body.split())}
    intervals = re.search(r"\bIntervals\s*=\s*(\d+)", body, re.IGNORECASE)
    if intervals:
        result["intervals"] = int(intervals.group(1))
    range_match = re.search(r"Range\s*=\s*\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)", body, re.IGNORECASE)
    if range_match:
        result["range"] = [float(range_match.group(1)), float(range_match.group(2))]
    return result


def iter_quasistationary(solve_body: str) -> list[dict[str, Any]]:
    sweeps = []
    pattern = re.compile(r"\bQuasistationary\b\s*\(", re.IGNORECASE)
    for match in pattern.finditer(solve_body):
        params_open = match.end() - 1
        params_close = find_matching(solve_body, params_open, "(", ")")
        brace_match = re.search(r"\{", solve_body[params_close + 1:])
        if not brace_match:
            continue
        body_open = params_close + 1 + brace_match.start()
        body_close = find_matching(solve_body, body_open, "{", "}")
        sweeps.append({
            "params": solve_body[params_open + 1:params_close],
            "body": solve_body[body_open + 1:body_close],
            "start": match.start(),
            "end": body_close + 1,
        })
    return sweeps


def parse_initial_steps(solve_body: str, sweep_ranges: list[tuple[int, int]]) -> list[dict[str, Any]]:
    chars = list(solve_body)
    for start, end in sweep_ranges:
        for index in range(start, end):
            chars[index] = " "
    text = "".join(chars)
    coupled_steps = parse_coupled(text)
    occupied = [(step["start"], step["end"]) for step in coupled_steps]
    events: list[tuple[int, dict[str, Any]]] = []
    for step in coupled_steps:
        events.append((step["start"], {
            "type": "Coupled",
            "parameters": step["parameters"],
            "equations": step["equations"],
        }))
    for match in re.finditer(r"\bPoisson\b", text, re.IGNORECASE):
        if any(start <= match.start() < end for start, end in occupied):
            continue
        events.append((match.start(), {"type": "Poisson", "equations": ["Poisson"]}))
    return [item for _position, item in sorted(events, key=lambda pair: pair[0])]


def parse_solve_block(text: str) -> dict[str, Any]:
    block = first_block(text, "Solve")
    if not block:
        return {"initial_steps": [], "new_current_files": []}
    body = block[1]
    sweep_blocks = iter_quasistationary(body)
    new_current_files = [
        match.group(1)
        for match in re.finditer(r"NewCurrentFile\s*=\s*\"([^\"]+)\"", body, re.IGNORECASE)
    ]
    sweeps = []
    for index, sweep in enumerate(sweep_blocks):
        params = sweep["params"]
        action = sweep["body"]
        goal = re.search(
            r"Goal\s*\{[^{}]*Name\s*=\s*\"([^\"]+)\"[^{}]*Voltage\s*=\s*(\"[^\"]*\"|@\S+?@|[-+0-9.eE]+)",
            params,
            re.DOTALL | re.IGNORECASE,
        )
        if not goal:
            continue
        goal_block = goal.group(0)
        step_control = parse_assignments(params.replace(goal_block, " "))
        stop_value = coerce_scalar(goal.group(2))
        equations = parse_coupled(action)
        current_file = None
        preceding = body[:sweep["start"]]
        current_matches = list(re.finditer(r"NewCurrentFile\s*=\s*\"([^\"]+)\"", preceding, re.IGNORECASE))
        if current_matches:
            current_file = current_matches[-1].group(1)
        entry: dict[str, Any] = {
            "mode": "quasistationary",
            "contact": goal.group(1),
            "stop": stop_value,
            "step_control": step_control,
            "equations": equations[0]["equations"] if equations else [],
        }
        current_plot = parse_current_plot(action)
        if current_plot:
            entry["current_plot"] = current_plot
        if current_file:
            entry["current_file_prefix"] = current_file
        sweeps.append(entry)

    return {
        "initial_steps": parse_initial_steps(body, [(item["start"], item["end"]) for item in sweep_blocks]),
        "new_current_files": new_current_files,
        "sweeps": sweeps,
    }


def parse_sweeps(text: str) -> list[dict[str, Any]]:
    return parse_solve_block(text).get("sweeps", [])


def parse_unsupported_physics(text: str) -> list[str]:
    found = []
    for token in UNSUPPORTED_PHYSICS:
        if re.search(rf"\b{re.escape(token)}\b", text):
            found.append(token)
    return found


def unsupported_report(tokens: list[str]) -> list[dict[str, str]]:
    reasons = {
        "Thermodynamic": "Vela does not yet solve the full self-heating temperature equation from SDevice.",
        "IALMob": "IALMob surface-orientation mobility is not represented by the current mobility model.",
        "Trap": "Interface and bulk trap kinetics are imported as metadata only.",
        "Traps": "Interface and bulk trap kinetics are imported as metadata only.",
        "Avalanche": "Avalanche generation is not calibrated to Sentaurus model semantics yet.",
        "eTemperature": "Carrier temperature transport is not supported in the Vela runner.",
        "hTemperature": "Carrier temperature transport is not supported in the Vela runner.",
    }
    return [{"feature": token, "reason": reasons.get(token, "unsupported Sentaurus physics feature")} for token in tokens]


def parse_cmd(input_path: Path, template_vars: dict[str, str] | None = None) -> dict[str, Any]:
    variables = template_vars or {}
    raw_text = remove_comments(input_path.read_text(errors="ignore"))
    text = apply_template_variables(raw_text, variables)
    solve = parse_solve_block(text)
    unsupported = parse_unsupported_physics(text)
    return {
        "template_variables": variables,
        "unresolved_placeholders": sorted(dict.fromkeys(PLACEHOLDER_RE.findall(text))),
        "files": parse_file_block(text),
        "electrodes": parse_electrodes(text),
        "thermodes": parse_thermodes(text),
        "physics": parse_physics_blocks(text),
        "plot_fields": parse_plot_fields(text),
        "math": parse_math_block(text),
        "solve": solve,
        "sweeps": solve.get("sweeps", []),
        "unsupported_physics": unsupported,
        "unsupported_report": unsupported_report(unsupported),
    }


def build_runner_deck(summary: dict[str, Any], mesh_json: str, output_csv: str) -> dict[str, Any]:
    electrodes = {str(item["name"]): float(item["voltage"]) for item in summary.get("electrodes", [])}
    sweeps = summary.get("sweeps", [])
    final_sweep = sweeps[-1] if sweeps else {"contact": "drain", "stop": 0.0}
    final_stop = numeric_or_none(final_sweep.get("stop"))
    if final_stop is None:
        final_stop = 0.0
    start = electrodes.get(str(final_sweep["contact"]), 0.0)
    intervals = final_sweep.get("current_plot", {}).get("intervals") if isinstance(final_sweep.get("current_plot"), dict) else None
    step = (final_stop - start) / float(intervals) if intervals else final_stop / 80.0 if final_stop != 0.0 else 1.0
    contacts = [
        {"name": name, "bias": voltage}
        for name, voltage in electrodes.items()
    ]
    return {
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
            "start": start,
            "stop": final_stop,
            "step": step,
            "write_vtk": False,
        },
        "sentaurus_import": {
            "unsupported_physics": summary.get("unsupported_physics", []),
            "unsupported_report": summary.get("unsupported_report", []),
            "sweeps": sweeps,
            "physics": summary.get("physics", []),
            "plot_fields": summary.get("plot_fields", []),
            "math": summary.get("math", {}),
            "solve": summary.get("solve", {}),
        },
    }


def write_deck_json(deck: dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(deck, indent=2) + "\n")


def write_python_runner(summary: dict[str, Any],
                        output_path: Path,
                        mesh_json: str,
                        output_csv: str,
                        deck_json: Path | None = None) -> None:
    deck = build_runner_deck(summary, mesh_json, output_csv)
    if deck_json is not None:
        write_deck_json(deck, deck_json)
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
    summary = parse_cmd(args.input, parse_template_vars(args.template_var))
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    if args.python_runner:
        write_python_runner(summary, args.python_runner, args.mesh_json, args.output_csv, args.deck_json)
    elif args.deck_json:
        write_deck_json(build_runner_deck(summary, args.mesh_json, args.output_csv), args.deck_json)


def plt_command(args: argparse.Namespace) -> None:
    export_plt_reference(args.input, args.output, args.bias_column, args.current_column)


def default_tdr_importer() -> str:
    build_dir = Path(__file__).resolve().parents[1] / "build"
    exe_name = "sentaurus_import.exe" if sys.platform.startswith("win") else "sentaurus_import"
    return str(build_dir / exe_name)


def run_command(command: list[str], cwd: Path) -> None:
    proc = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {proc.returncode}: {' '.join(command)}\n"
            f"stdout={proc.stdout.strip()}\n"
            f"stderr={proc.stderr.strip()}"
        )


def relative_generated(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def project_artifacts(project_dir: Path, node: int) -> dict[str, Path]:
    return {
        "tdr": project_dir / f"n{node}_des.tdr",
        "plt": project_dir / f"IdVd_n{node}_des.plt",
        "cmd": project_dir / f"pp{node}_des.cmd",
    }


def write_reference_manifest(output_dir: Path,
                             node: int,
                             device: str,
                             artifacts: dict[str, Path],
                             generated: list[str],
                             warnings: list[str]) -> None:
    reference_curves = sorted(item for item in generated if item.startswith("reference_curves/"))
    vela_decks = sorted(
        item for item in generated
        if item.startswith("vela/") and item.endswith(".json")
    )
    manifest = {
        "schema": "vela.reference_tcad.sentaurus_project.v1",
        "device": device,
        "node": node,
        "source_artifacts": {
            key: {
                "path": str(path),
                "filename": path.name,
            }
            for key, path in artifacts.items()
        },
        "generated": sorted(dict.fromkeys(generated)),
        "reference_curves": reference_curves,
        "vela_decks": vela_decks,
        "warnings": warnings,
        "commit_policy": {
            "raw_sentaurus_artifacts": False,
            "generated_reference_tree": "local-or-explicit-review",
        },
    }
    (output_dir / "reference_tcad_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def project_command(args: argparse.Namespace) -> None:
    project_dir = args.project_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    artifacts = project_artifacts(project_dir, args.node)
    for label, path in artifacts.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing Sentaurus {label} artifact: {path}")

    generated: list[str] = []
    tdr_importer = shlex.split(args.tdr_importer or default_tdr_importer(),
                               posix=not sys.platform.startswith("win"))
    if not tdr_importer:
        raise ValueError("tdr importer command is empty")
    run_command(
        [
            *tdr_importer,
            "--tdr",
            str(artifacts["tdr"]),
            "--inventory-json",
            str(output_dir / "tdr_inventory.json"),
            "--export-dir",
            str(output_dir),
        ],
        cwd=REPO,
    )
    for name in [
        "nodes.csv",
        "elements.csv",
        "contacts.csv",
        "doping.csv",
        "metadata.json",
        "field_manifest.json",
        "tdr_inventory.json",
    ]:
        if (output_dir / name).exists():
            generated.append(name)

    reference_curve = output_dir / "reference_curves" / f"{args.device}_idvd_reference.csv"
    export_plt_reference(
        artifacts["plt"],
        reference_curve,
        args.bias_column,
        args.current_column,
    )
    generated.append(relative_generated(reference_curve, output_dir))

    cmd_summary = output_dir / "cmd_summary.json"
    runner = output_dir / "run_idvd.py"
    runner_deck = output_dir / "run_idvd_deck.json"
    summary = parse_cmd(artifacts["cmd"])
    cmd_summary.write_text(json.dumps(summary, indent=2) + "\n")
    write_python_runner(
        summary,
        runner,
        args.mesh_json,
        args.output_csv,
        runner_deck,
    )
    generated.extend(["cmd_summary.json", "run_idvd.py", "run_idvd_deck.json"])

    vela_dir = output_dir / "vela"
    run_command(
        [
            sys.executable,
            str(REPO / "scripts" / "convert_tcad_export.py"),
            "--input-dir",
            str(output_dir),
            "--output-dir",
            str(vela_dir),
            "--device",
            args.device,
            "--simulation-types",
            "iv",
        ],
        cwd=REPO,
    )
    for name in ["mesh.json", "simulation_iv.json"]:
        if (vela_dir / name).exists():
            generated.append(f"vela/{name}")

    warnings = [
        "unsupported physics: " + ", ".join(summary["unsupported_physics"])
    ] if summary.get("unsupported_physics") else []
    generated.append("reference_tcad_manifest.json")
    write_reference_manifest(
        output_dir,
        args.node,
        args.device,
        artifacts,
        generated,
        warnings,
    )
    import_summary = {
        "node": args.node,
        "device": args.device,
        "sources": {key: str(value) for key, value in artifacts.items()},
        "generated": sorted(dict.fromkeys(generated)),
        "warnings": warnings,
    }
    (output_dir / "import_summary.json").write_text(json.dumps(import_summary, indent=2) + "\n")
    print(json.dumps(import_summary, indent=2))


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
    cmd.add_argument("--deck-json", type=Path)
    cmd.add_argument("--mesh-json", default="mesh.json")
    cmd.add_argument("--output-csv", default="outputs/sentaurus_import.csv")
    cmd.add_argument(
        "--template-var",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Expand a Sentaurus Workbench placeholder such as @Vg@ before parsing",
    )
    cmd.set_defaults(func=cmd_command)

    project = sub.add_parser("project", help="Import a Sentaurus project directory into a neutral reference tree")
    project.add_argument("--project-dir", type=Path, required=True)
    project.add_argument("--node", type=int, required=True)
    project.add_argument("--device", choices=["pn_diode", "nmos2d", "pmos2d", "ldmos2d", "igbt2d"], required=True)
    project.add_argument("--output-dir", type=Path, required=True)
    project.add_argument("--tdr-importer", default=None, help="Command used to run the C++ TDR importer")
    project.add_argument("--bias-column", default="drain OuterVoltage")
    project.add_argument("--current-column", default="drain TotalCurrent")
    project.add_argument("--mesh-json", default="vela/mesh.json")
    project.add_argument("--output-csv", default="outputs/ldmos_idvd.csv")
    project.set_defaults(func=project_command)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
