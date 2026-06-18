#!/usr/bin/env python3
"""Evaluate Vela Newton residuals for PN2D BV Vela/Sentaurus external states."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


STATE_FIELD_MAP = {
    "ElectrostaticPotential": "Potential",
    "eQuasiFermiPotential": "ElectronQuasiFermi",
    "hQuasiFermiPotential": "HoleQuasiFermi",
}

NODE_FIELDS = [
    "state",
    "source",
    "bias_V",
    "node_id",
    "x",
    "y",
    "psi",
    "phin",
    "phip",
    "psi_residual",
    "phin_residual",
    "phip_residual",
    "abs_psi_residual",
    "abs_phin_residual",
    "abs_phip_residual",
    "donors_m3",
    "acceptors_m3",
    "net_doping_m3",
    "ni_eff_m3",
    "block_psi",
    "block_phin",
    "block_phip",
    "block_combined",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--states",
        default="vela:-10,vela:-13.2,sentaurus:-13.2",
        help="Comma-separated source:bias list, e.g. vela:-10,sentaurus:-13.2.",
    )
    parser.add_argument("--nodes", default="351,986")
    parser.add_argument("--bias-contact", default="Anode")
    parser.add_argument("--ground-contact", default="Cathode")
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args()


def parse_states(raw: str) -> list[tuple[str, float]]:
    result = []
    valid_sources = {
        "vela",
        "sentaurus",
        "hybrid_vpsi_sqf",
        "hybrid_spsi_vqf",
        "hybrid_spsi_shift_vqf",
    }
    for item in raw.split(","):
        if not item.strip():
            continue
        source, value = item.split(":", 1)
        source = source.strip().lower()
        if source not in valid_sources:
            raise ValueError(f"unknown state source: {source}")
        result.append((source, float(value.strip())))
    return result


def parse_nodes(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def state_name(source: str, bias: float) -> str:
    token = f"{bias:g}".replace("-", "m").replace(".", "p")
    return f"{source}_{token}v"


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    result: dict[float, Path] = {}
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if match:
            result[bias_key(float(match.group("bias")))] = path
    return result


def sentaurus_dir(root: Path, bias: float) -> Path:
    candidates = [
        root / f"sentaurus_{signed_bias_token(bias)}v",
        root / f"sentaurus_{bias:g}v",
    ]
    for candidate in candidates:
        if (candidate / "fields").exists():
            return candidate
    raise FileNotFoundError(f"missing Sentaurus export for {bias:g} V under {root}")


def write_scalar_csv(path: Path, values: list[float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "component0"])
        for node_id, value in enumerate(values):
            writer.writerow([node_id, f"{value:.17g}"])


def read_scalar_csv(path: Path) -> list[float]:
    values_by_id: dict[int, float] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            values_by_id[int(row["node_id"])] = float(row["component0"])
    if not values_by_id:
        return []
    return [values_by_id[node_id] for node_id in range(max(values_by_id) + 1)]


def load_vela_state(vtk: Path, node_count: int) -> dict[str, list[float]]:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    values_by_field: dict[str, list[float]] = {}
    for output_name, vtk_name in STATE_FIELD_MAP.items():
        if vtk_name not in scalars:
            raise RuntimeError(f"{vtk} is missing VTK scalar {vtk_name}")
        values = scalars[vtk_name]
        if len(values) != node_count:
            raise RuntimeError(
                f"{vtk_name} has {len(values)} values, expected {node_count}")
        values_by_field[output_name] = values
    return values_by_field


def load_sentaurus_state(export_dir: Path, node_count: int) -> dict[str, list[float]]:
    values_by_field: dict[str, list[float]] = {}
    for output_name in STATE_FIELD_MAP:
        source = export_dir / "fields" / f"{output_name}_region0.csv"
        if not source.exists():
            raise FileNotFoundError(source)
        values = read_scalar_csv(source)
        if len(values) != node_count:
            raise RuntimeError(
                f"{source} has {len(values)} values, expected {node_count}")
        values_by_field[output_name] = values
    return values_by_field


def write_state(fields_dir: Path, values_by_field: dict[str, list[float]]) -> None:
    for output_name in STATE_FIELD_MAP:
        write_scalar_csv(fields_dir / f"{output_name}_region0.csv", values_by_field[output_name])


def shifted(values: list[float], delta: float) -> list[float]:
    return [value + delta for value in values]


def psi_shift_for_nodes(reference: list[float], candidate: list[float], nodes: list[int]) -> float:
    if not nodes:
        return 0.0
    return sum(reference[node_id] - candidate[node_id] for node_id in nodes) / len(nodes)


def prepare_state(
    source: str,
    bias: float,
    fields_dir: Path,
    node_count: int,
    selected_nodes: list[int],
    vtks: dict[float, Path],
    sentaurus_root: Path,
) -> tuple[str, dict[str, Any]]:
    vtk = vtks.get(bias_key(bias))
    export = sentaurus_dir(sentaurus_root, bias)
    metadata: dict[str, Any] = {}
    if source == "vela":
        if vtk is None:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        write_state(fields_dir, load_vela_state(vtk, node_count))
        return str(vtk), metadata
    if source == "sentaurus":
        write_state(fields_dir, load_sentaurus_state(export, node_count))
        return str(export), metadata

    if vtk is None:
        raise FileNotFoundError(f"missing Vela VTK for hybrid {bias:g} V")
    vela = load_vela_state(vtk, node_count)
    sentaurus = load_sentaurus_state(export, node_count)
    hybrid: dict[str, list[float]]
    if source == "hybrid_vpsi_sqf":
        hybrid = {
            "ElectrostaticPotential": vela["ElectrostaticPotential"],
            "eQuasiFermiPotential": sentaurus["eQuasiFermiPotential"],
            "hQuasiFermiPotential": sentaurus["hQuasiFermiPotential"],
        }
        metadata["components"] = "Vela psi + Sentaurus phin/phip"
    elif source == "hybrid_spsi_vqf":
        hybrid = {
            "ElectrostaticPotential": sentaurus["ElectrostaticPotential"],
            "eQuasiFermiPotential": vela["eQuasiFermiPotential"],
            "hQuasiFermiPotential": vela["hQuasiFermiPotential"],
        }
        metadata["components"] = "Sentaurus psi + Vela phin/phip"
    elif source == "hybrid_spsi_shift_vqf":
        delta = psi_shift_for_nodes(
            vela["ElectrostaticPotential"],
            sentaurus["ElectrostaticPotential"],
            selected_nodes,
        )
        hybrid = {
            "ElectrostaticPotential": shifted(sentaurus["ElectrostaticPotential"], delta),
            "eQuasiFermiPotential": vela["eQuasiFermiPotential"],
            "hQuasiFermiPotential": vela["hQuasiFermiPotential"],
        }
        metadata["components"] = "Sentaurus psi shifted to selected-node Vela mean + Vela phin/phip"
        metadata["psi_shift_V"] = delta
        metadata["psi_shift_nodes"] = selected_nodes
    else:
        raise ValueError(source)
    write_state(fields_dir, hybrid)
    return f"{vtk}; {export}", metadata


def build_probe_config(
    base_config: Path,
    output_csv: Path,
    fields_dir: Path,
    bias: float,
    bias_contact: str,
    ground_contact: str,
) -> dict[str, Any]:
    cfg = json.loads(base_config.read_text())
    base_dir = base_config.parent
    cfg["simulation_type"] = "newton_residual_probe"
    cfg["output_csv"] = str(output_csv.resolve())
    cfg["state_fields_dir"] = str(fields_dir.resolve())
    for key in ["mesh_file", "node_doping_file", "materials_file"]:
        if key in cfg:
            path = Path(cfg[key])
            cfg[key] = str((path if path.is_absolute() else base_dir / path).resolve())
    cfg.pop("sweep", None)
    for contact in cfg.get("contacts", []):
        name = contact.get("name")
        if name == bias_contact:
            contact["bias"] = bias
        elif name == ground_contact:
            contact["bias"] = 0.0
    return cfg


def run_runner(runner: Path, config: Path) -> dict[str, Any]:
    result = subprocess.run(
        [str(runner.resolve()), "--config", str(config.resolve())],
        cwd=config.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"runner failed for {config}: {result.stderr.strip() or result.stdout.strip()}")
    return json.loads(result.stdout)


def read_residual_rows(path: Path) -> dict[int, dict[str, str]]:
    with path.open(newline="") as handle:
        return {int(row["node_id"]): row for row in csv.DictReader(handle)}


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    return value


def main() -> int:
    args = parse_args()
    if args.runner is None and not args.prepare_only:
        raise RuntimeError("--runner is required unless --prepare-only is set")
    mesh = json.loads(args.base_config.read_text())
    mesh_path = args.base_config.parent / mesh["mesh_file"]
    mesh_data = json.loads(mesh_path.read_text())
    node_count = len(mesh_data["nodes"])
    vtks = discover_vela_vtks(args.vela_vtk_root)
    selected_nodes = parse_nodes(args.nodes)
    node_rows: list[dict[str, Any]] = []
    states_summary: list[dict[str, Any]] = []

    for source, bias in parse_states(args.states):
        name = state_name(source, bias)
        fields_dir = args.out_dir / "states" / name / "fields"
        source_path, metadata = prepare_state(
            source,
            bias,
            fields_dir,
            node_count,
            selected_nodes,
            vtks,
            args.sentaurus_root,
        )

        output_csv = args.out_dir / "residuals" / f"{name}_residual.csv"
        config_path = args.out_dir / "configs" / f"{name}_newton_residual_probe.json"
        config_path.parent.mkdir(parents=True, exist_ok=True)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        cfg = build_probe_config(
            args.base_config,
            output_csv,
            fields_dir,
            bias,
            args.bias_contact,
            args.ground_contact,
        )
        config_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
        state_summary: dict[str, Any] = {
            "state": name,
            "source": source,
            "bias_V": bias,
            "source_path": source_path,
            "fields_dir": str(fields_dir),
            "config": str(config_path),
            "output_csv": str(output_csv),
            **metadata,
        }
        if not args.prepare_only:
            assert args.runner is not None
            status = run_runner(args.runner, config_path)
            state_summary["status"] = status
            residual_rows = read_residual_rows(output_csv)
            blocks = status.get("block_residuals", {})
            for node_id in selected_nodes:
                row = dict(residual_rows[node_id])
                row.update({
                    "state": name,
                    "source": source,
                    "bias_V": bias,
                    "block_psi": blocks.get("psi"),
                    "block_phin": blocks.get("phin"),
                    "block_phip": blocks.get("phip"),
                    "block_combined": blocks.get("combined"),
                })
                node_rows.append(row)
        states_summary.append(state_summary)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if node_rows:
        write_csv(args.out_dir / "newton_residual_state_nodes.csv", node_rows, NODE_FIELDS)
    (args.out_dir / "newton_residual_state_summary.json").write_text(
        json.dumps(clean_json({
            "states": states_summary,
            "selected_nodes": selected_nodes,
            "node_rows": len(node_rows),
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
