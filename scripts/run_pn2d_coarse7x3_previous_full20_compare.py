#!/usr/bin/env python3
"""Run the PN2D coarse 7x3 Sentaurus/Vela previous-full20 diagnostic.

The script intentionally keeps the coarse Sentaurus case separate from the
fine PN2D reference. It can reuse already-fetched Sentaurus artifacts, import
them into neutral CSV/Vela decks, patch the generated Vela BV deck to the
previous-full20 physics path, and write diagnostic comparison artifacts.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
COARSE_CASE = "pn2d_sentaurus2018_coarse7x3"
DEFAULT_CASE_ROOT = REPO / "reference_tcad" / COARSE_CASE
DEFAULT_BUILD_ROOT = REPO / "build-release" / "reference_tcad" / COARSE_CASE
DEFAULT_REFERENCE_CONFIG = DEFAULT_CASE_ROOT / f"{COARSE_CASE}_reference.json"
DEFAULT_RUNNER = REPO / "build-release" / ("vela_example_runner.exe" if os.name == "nt" else "vela_example_runner")
DEFAULT_IMPORTER = REPO / "build-release" / ("sentaurus_import.exe" if os.name == "nt" else "sentaurus_import")

PREVIOUS_FULL20_BIAS_POINTS = [
    0.0,
    -0.002,
    -0.00438666666667,
    -0.00725066666667,
    -0.01067792,
    -0.0147791998222,
    -0.0196870646761,
    -0.0255601429513,
    -0.032588259954,
    -0.0409985733005,
    -0.0510629149385,
    -0.0631065770986,
    -0.0775188261502,
    -0.094765484182,
    -0.11540398496,
    -0.140101390891,
    -0.169655953322,
    -0.205022913031,
    -0.247345374816,
    -0.297991254085,
    -0.358428670013,
    -0.43055065302,
    -0.516616219409,
    -0.619321128633,
    -0.741882320307,
    -0.888138675705,
    -1.06267125981,
    -1.27094681018,
    -1.51948896696,
    -1.81608260737,
    -2.17001768494,
    -2.59120042724,
    -3.09240789057,
    -3.68884477194,
    -4.39860466078,
    -5.24321892849,
    -6.24321892849,
    -7.24321892849,
    -8.24321892849,
    -9.24321892849,
    -10.2432189285,
    -11.2432189285,
    -12.2432189285,
    -13.2432189285,
    -14.2432189285,
    -15.2432189285,
    -16.2432189285,
    -17.2432189285,
    -18.2432189285,
    -19.2432189285,
    -20.0,
]

REQUIRED_SENT_FIELDS = [
    "eAlphaAvalanche",
    "hAlphaAvalanche",
    "eVelocity",
    "hVelocity",
    "eIonIntegral",
    "hIonIntegral",
    "MeanIonIntegral",
]

ELEMENTARY_CHARGE_C = 1.602176634e-19

NODE_FIELD_MAP: dict[str, dict[str, Any]] = {
    "potential": {"sentaurus": ["ElectrostaticPotential", "Potential"], "vela": ["Potential"], "scale": 1.0},
    "electric_field": {"sentaurus": ["ElectricField"], "vela": ["ElectricField"], "scale": 1.0},
    "electron_density": {"sentaurus": ["eDensity"], "vela": ["Electrons"], "scale": 1.0e-6},
    "hole_density": {"sentaurus": ["hDensity"], "vela": ["Holes"], "scale": 1.0e-6},
    "electron_qf": {"sentaurus": ["eQuasiFermiPotential", "eQuasiFermi"], "vela": ["ElectronQuasiFermi"], "scale": 1.0},
    "hole_qf": {"sentaurus": ["hQuasiFermiPotential", "hQuasiFermi"], "vela": ["HoleQuasiFermi"], "scale": 1.0},
    "electron_current": {"sentaurus": ["eCurrentDensity", "eCurrent"], "vela": ["ElectronCurrentDensity"], "scale": 1.0},
    "hole_current": {"sentaurus": ["hCurrentDensity", "hCurrent"], "vela": ["HoleCurrentDensity"], "scale": 1.0},
    "total_current": {"sentaurus": ["TotalCurrentDensity", "TotalCurrent"], "vela": ["TotalCurrentDensity"], "scale": 1.0},
    "srh": {"sentaurus": ["SRHRecombination", "srhRecombination"], "vela": ["SRHRecombination"], "scale": 1.0e-6},
    "avalanche": {"sentaurus": ["AvalancheGeneration", "ImpactIonization"], "vela": ["AvalancheGeneration"], "scale": 1.0e-6},
    "electron_mobility": {"sentaurus": ["eMobility"], "vela": ["ElectronMobility"], "scale": 1.0e4},
    "hole_mobility": {"sentaurus": ["hMobility"], "vela": ["HoleMobility"], "scale": 1.0e4},
    "electron_velocity": {"sentaurus": ["eVelocity"], "vela": ["ElectronVelocity"], "scale": 1.0},
    "hole_velocity": {"sentaurus": ["hVelocity"], "vela": ["HoleVelocity"], "scale": 1.0},
    "electron_alpha_avalanche": {"sentaurus": ["eAlphaAvalanche"], "vela": ["ElectronAlphaAvalanche"], "sentaurus_scale": 100.0, "scale": 1.0},
    "hole_alpha_avalanche": {"sentaurus": ["hAlphaAvalanche"], "vela": ["HoleAlphaAvalanche"], "sentaurus_scale": 100.0, "scale": 1.0},
    "electron_ion_integral": {"sentaurus": ["eIonIntegral"], "vela": ["ElectronIonIntegral"], "scale": 1.0},
    "hole_ion_integral": {"sentaurus": ["hIonIntegral"], "vela": ["HoleIonIntegral"], "scale": 1.0},
    "mean_ion_integral": {"sentaurus": ["MeanIonIntegral"], "vela": ["MeanIonIntegral"], "scale": 1.0},
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def run_command(cmd: list[str], cwd: Path) -> None:
    print(":: running", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def resolve_config_path(base_config: Path, raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw)
    if not path.is_absolute():
        path = base_config.parent / path
    return path.resolve()


def write_previous_full20_config(
    *,
    base_config: Path,
    out_dir: Path,
    output_csv_name: str = "coarse_previous_full20.csv",
    bias_points: list[float] | None = None,
    config_name: str = "simulation_coarse_previous_full20.json",
    vtk_subdir: str = "vtk",
    diagnostics_suffix: str = "",
) -> Path:
    config = read_json(base_config)
    out_dir.mkdir(parents=True, exist_ok=True)

    mesh = resolve_config_path(base_config, str(config.get("mesh_file", "mesh.json")))
    doping = resolve_config_path(base_config, str(config.get("node_doping_file", "doping.csv")))
    materials = resolve_config_path(base_config, config.get("materials_file"))
    if materials is None:
        candidate = base_config.parent / "pn2d_sentaurus2018_iv_materials.json"
        if candidate.exists():
            materials = candidate.resolve()
    if materials is None or not materials.exists():
        raise FileNotFoundError(
            f"{base_config} does not provide an existing materials_file; coarse PN2D diagnostics must not use default ni"
        )
    if mesh is None or not mesh.exists():
        raise FileNotFoundError(f"mesh_file does not exist for {base_config}: {mesh}")
    if doping is None or not doping.exists():
        raise FileNotFoundError(f"node_doping_file does not exist for {base_config}: {doping}")

    config["mesh_file"] = str(mesh)
    config["node_doping_file"] = str(doping)
    config["materials_file"] = str(materials)
    config["output_csv"] = output_csv_name

    solver = config.setdefault("solver", {})
    solver["method"] = "gummel_newton"
    solver["max_iter"] = 40
    solver["reltol"] = 1.0e-8
    solver["abstol"] = 1.0e-9
    solver["damping_psi"] = 0.2
    solver["damping_factor"] = 1.0
    solver["max_update"] = 0
    solver["line_search"] = True
    solver["warm_start"] = True
    solver["contact_boundary_reconstruction"] = "dominant_signed_contact_mean"
    solver["contact_boundary_minority_electron_relaxation"] = False
    solver["quasi_fermi_update_limit_V"] = 0.1
    solver["mobility"] = {
        "model": "masetti_field",
        "high_field_driving_force": "quasi_fermi_gradient",
        "jacobian_field_derivatives": False,
    }
    solver["recombination"] = ["srh"]
    solver["bandgap_narrowing"] = "old_slotboom"
    solver["impact_ionization"] = {
        "model": "van_overstraeten",
        "driving_force": "quasi_fermi_gradient",
        "generation": "current_density",
        "current_approximation": "grad_qf",
        "quasi_fermi_gradient_discretization": "cell_gradient",
    }
    solver["diagnostics"] = True

    sweep = config.setdefault("sweep", {})
    sweep["mode"] = "bv_reverse"
    sweep["start"] = 0.0
    sweep["stop"] = -20.0
    sweep["step"] = -0.05
    sweep["bias_points"] = list(PREVIOUS_FULL20_BIAS_POINTS if bias_points is None else bias_points)
    sweep["contact"] = "Anode"
    sweep["current_contact"] = "Anode"
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = f"{vtk_subdir}/dc_sweep"
    sweep["write_state_file"] = f"coarse_previous_full20{diagnostics_suffix}_last_state.csv"
    sweep["diagnostics"] = {
        "sg_avalanche_edges": {"enabled": True, "csv_file": str(out_dir / f"sg_avalanche_edges{diagnostics_suffix}.csv")},
        "continuity_balance": {
            "enabled": True,
            "contacts": ["Anode", "Cathode"],
            "csv_file": str(out_dir / f"continuity_balance{diagnostics_suffix}.csv"),
        },
        "newton_history": {"enabled": True, "csv_file": str(out_dir / f"newton_history{diagnostics_suffix}.csv")},
    }

    config_path = out_dir / config_name
    write_json(config_path, config)
    return config_path


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def normalize_sentaurus_source_dir(path: Path) -> Path:
    if (path / "sentaurus_vm_run_manifest.json").exists() and (path / "source").is_dir():
        return path / "source"
    return path


def ensure_sentaurus_source_inputs(source_dir: Path) -> None:
    template_dir = DEFAULT_CASE_ROOT / "source"
    for name in ["pn2d_sde.cmd", "pn2d_bv_sdevice.cmd", "models.par", "pn2d_sentaurus2018_iv_materials.json"]:
        target = source_dir / name
        if target.exists():
            continue
        source = template_dir / name
        if source.exists():
            shutil.copyfile(source, target)

def multibias_index_for_bias(bias: float) -> int:
    return int(round(abs(bias) / 20.0 * 400.0))


def export_selected_multibias(
    *,
    source_dir: Path,
    out_root: Path,
    tdr_importer: Path,
    biases: list[float],
) -> dict[float, Path]:
    exports: dict[float, Path] = {}
    for bias in biases:
        index = multibias_index_for_bias(bias)
        tdr = source_dir / f"pn2d_bv_multibias_{index:04d}_des.tdr"
        if not tdr.exists():
            continue
        out_dir = out_root / f"sentaurus_{signed_bias_token(bias)}v"
        out_dir.mkdir(parents=True, exist_ok=True)
        run_command([
            str(tdr_importer),
            "--tdr",
            str(tdr),
            "--inventory-json",
            str(out_dir / "tdr_inventory.json"),
            "--export-dir",
            str(out_dir),
        ], REPO)
        exports[bias_key(bias)] = out_dir
    return exports


def manifest_field_names(export_dir: Path) -> set[str]:
    manifest = export_dir / "field_manifest.json"
    if not manifest.exists():
        return set()
    data = read_json(manifest)
    return {str(item.get("name", "")) for item in data.get("fields", [])}


def validate_required_sentaurus_fields(export_dirs: list[Path]) -> dict[str, Any]:
    if not export_dirs:
        return {"status": "missing_exports", "missing_fields": REQUIRED_SENT_FIELDS}
    field_sets = [manifest_field_names(path) for path in export_dirs if (path / "field_manifest.json").exists()]
    if not field_sets:
        return {"status": "missing_field_manifest", "missing_fields": REQUIRED_SENT_FIELDS}
    common = set.intersection(*field_sets) if field_sets else set()
    missing = sorted(field for field in REQUIRED_SENT_FIELDS if field not in common)
    return {
        "status": "pass" if not missing else "missing_required_fields",
        "missing_fields": missing,
        "checked_export_count": len(field_sets),
    }


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def load_nodes(export_dir: Path) -> list[dict[str, Any]]:
    rows = read_csv_rows(export_dir / "nodes.csv")
    return [
        {"id": int(row["id"]), "x_um": float(row["x_um"]), "y_um": float(row["y_um"])}
        for row in rows
    ]


def load_sentaurus_scalar(export_dir: Path, aliases: list[str]) -> tuple[str, dict[int, float]] | None:
    for alias in aliases:
        for path in [export_dir / "fields" / f"{alias}_region0.csv", export_dir / f"{alias}_region0.csv"]:
            if not path.exists():
                continue
            values: dict[int, float] = {}
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                component_cols = [name for name in reader.fieldnames or [] if name != "node_id"]
                for row in reader:
                    comps = [float(row[name]) for name in component_cols if row.get(name, "") != ""]
                    values[int(row["node_id"])] = comps[0] if len(comps) == 1 else math.sqrt(sum(c * c for c in comps))
            return alias, values
    return None


def parse_vtk(path: Path) -> dict[str, list[float]]:
    lines = path.read_text().splitlines()
    scalars: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            values: list[float] = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {"SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA"}:
                    break
                values.extend(float(value) for value in next_parts)
                i += 1
            scalars[name] = values
            continue
        i += 1
    return scalars


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    result: dict[float, Path] = {}
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if match:
            result[bias_key(float(match.group("bias")))] = path
    return result


def node_field_compare(
    *,
    sentaurus_exports: dict[float, Path],
    vela_vtks: dict[float, Path],
    biases: list[float],
    exact_vela_bias: bool = False,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for bias in biases:
        sent_dir = sentaurus_exports.get(bias_key(bias))
        if sent_dir is None:
            continue
        nodes = load_nodes(sent_dir)
        exact_bias_key = bias_key(bias)
        if exact_vela_bias:
            vela_bias = exact_bias_key if exact_bias_key in vela_vtks else None
        else:
            vela_bias = min(vela_vtks, key=lambda item: abs(item - bias)) if vela_vtks else None
        vtk_scalars = parse_vtk(vela_vtks[vela_bias]) if vela_bias is not None else {}
        for quantity, spec in NODE_FIELD_MAP.items():
            loaded = load_sentaurus_scalar(sent_dir, spec["sentaurus"])
            sent_field = loaded[0] if loaded else ""
            sent_values = loaded[1] if loaded else {}
            vela_field = next((name for name in spec["vela"] if name in vtk_scalars), "")
            vela_values = vtk_scalars.get(vela_field, [])
            for node in nodes:
                node_id = int(node["id"])
                raw_sent_value = sent_values.get(node_id)
                sent_scale = float(spec.get("sentaurus_scale", 1.0))
                vela_scale = float(spec.get("vela_scale", spec.get("scale", 1.0)))
                sent_value = None if raw_sent_value is None else raw_sent_value * sent_scale
                vela_value = None
                if vela_values and node_id < len(vela_values):
                    vela_value = float(vela_values[node_id]) * vela_scale
                ratio = None
                if sent_value not in (None, 0.0) and vela_value is not None:
                    ratio = vela_value / sent_value
                rows.append({
                    "bias_V": bias,
                    "vela_bias_V": vela_bias,
                    "node_id": node_id,
                    "x_um": node["x_um"],
                    "y_um": node["y_um"],
                    "quantity": quantity,
                    "sentaurus_field": sent_field,
                    "sentaurus_value": sent_value,
                    "vela_field": vela_field,
                    "vela_value_scaled_to_sentaurus_units": vela_value,
                    "vela_over_sentaurus": ratio,
                    "abs_diff": None if sent_value is None or vela_value is None else vela_value - sent_value,
                })
    return rows


def nearest_export(sentaurus_exports: dict[float, Path], bias: float) -> tuple[float, Path] | None:
    if not sentaurus_exports:
        return None
    key = min(sentaurus_exports, key=lambda item: abs(item - bias))
    return key, sentaurus_exports[key]


def current_support_compare(
    *,
    sg_edges_csv: Path,
    sentaurus_exports: dict[float, Path],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for edge in read_csv_rows(sg_edges_csv):
        try:
            bias = float(edge.get("bias_V", "nan"))
            node0 = int(edge["node0"])
            node1 = int(edge["node1"])
        except (KeyError, ValueError):
            continue
        matched = nearest_export(sentaurus_exports, bias)
        sent_values: dict[str, dict[int, float]] = {}
        matched_bias = None
        if matched is not None:
            matched_bias, export_dir = matched
            for name, spec in {
                "sent_e_velocity": {"aliases": ["eVelocity"]},
                "sent_h_velocity": {"aliases": ["hVelocity"]},
                "sent_e_alpha": {"aliases": ["eAlphaAvalanche"], "scale": 100.0},
                "sent_h_alpha": {"aliases": ["hAlphaAvalanche"], "scale": 100.0},
                "sent_e_ion_integral": {"aliases": ["eIonIntegral"]},
                "sent_h_ion_integral": {"aliases": ["hIonIntegral"]},
                "sent_mean_ion_integral": {"aliases": ["MeanIonIntegral"]},
                "sent_e_current": {"aliases": ["eCurrentDensity", "eCurrent"]},
                "sent_h_current": {"aliases": ["hCurrentDensity", "hCurrent"]},
            }.items():
                loaded = load_sentaurus_scalar(export_dir, spec["aliases"])
                scale = float(spec.get("scale", 1.0))
                sent_values[name] = {node: value * scale for node, value in loaded[1].items()} if loaded else {}

        def endpoint_avg(values: dict[int, float]) -> float | None:
            vals = [values[node] for node in (node0, node1) if node in values]
            return sum(vals) / len(vals) if vals else None

        row: dict[str, Any] = {
            "bias_V": bias,
            "nearest_sentaurus_bias_V": matched_bias,
            "edge_id": edge.get("edge_id", ""),
            "node0": node0,
            "node1": node1,
            "edge_class": edge.get("edge_class", ""),
            "source_integral": edge.get("source_integral", ""),
            "electron_source_integral": edge.get("electron_source_integral", ""),
            "hole_source_integral": edge.get("hole_source_integral", ""),
            "electron_flux_abs": edge.get("electron_flux_abs", edge.get("electron_flux_proxy", "")),
            "hole_flux_abs": edge.get("hole_flux_abs", edge.get("hole_flux_proxy", "")),
            "electron_alpha_m_inv": edge.get("electron_alpha_m_inv", ""),
            "hole_alpha_m_inv": edge.get("hole_alpha_m_inv", ""),
        }
        for key, values in sent_values.items():
            row[key] = endpoint_avg(values)
        rows.append(row)
    return rows


def optional_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def build_active_edge_top_rows(
    rows: list[dict[str, Any]],
    *,
    target_bias: float,
    nearest_vela_bias: float | None = None,
    limit: int = 30,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for row in rows:
        sent_bias = optional_float(row.get("nearest_sentaurus_bias_V"))
        vela_bias = optional_float(row.get("bias_V"))
        if sent_bias is None or abs(sent_bias - target_bias) > 1.0e-9:
            continue
        if nearest_vela_bias is not None and (vela_bias is None or abs(vela_bias - nearest_vela_bias) > 1.0e-9):
            continue
        electron_flux = optional_float(row.get("electron_flux_abs")) or 0.0
        hole_flux = optional_float(row.get("hole_flux_abs")) or 0.0
        electron_alpha = optional_float(row.get("electron_alpha_m_inv")) or 0.0
        hole_alpha = optional_float(row.get("hole_alpha_m_inv")) or 0.0
        electron_alpha_flux = electron_flux * electron_alpha
        hole_alpha_flux = hole_flux * hole_alpha
        combined_alpha_flux = electron_alpha_flux + hole_alpha_flux
        source_integral = optional_float(row.get("source_integral"))
        if source_integral is None:
            source_integral = (optional_float(row.get("electron_source_integral")) or 0.0) + \
                (optional_float(row.get("hole_source_integral")) or 0.0)
        sent_e_current = optional_float(row.get("sent_e_current"))
        sent_h_current = optional_float(row.get("sent_h_current"))
        sent_e_flux = abs(sent_e_current) * 1.0e4 / ELEMENTARY_CHARGE_C if sent_e_current is not None else None
        sent_h_flux = abs(sent_h_current) * 1.0e4 / ELEMENTARY_CHARGE_C if sent_h_current is not None else None
        sent_e_alpha = optional_float(row.get("sent_e_alpha"))
        sent_h_alpha = optional_float(row.get("sent_h_alpha"))
        item = {
            "rank": 0,
            "target_sentaurus_bias_V": target_bias,
            "bias_V": vela_bias,
            "nearest_sentaurus_bias_V": sent_bias,
            "edge_id": str(row.get("edge_id", "")),
            "node0": row.get("node0", ""),
            "node1": row.get("node1", ""),
            "edge_class": row.get("edge_class", ""),
            "electron_flux_abs": electron_flux,
            "hole_flux_abs": hole_flux,
            "electron_alpha_m_inv": electron_alpha,
            "hole_alpha_m_inv": hole_alpha,
            "electron_alpha_flux": electron_alpha_flux,
            "hole_alpha_flux": hole_alpha_flux,
            "combined_alpha_flux": combined_alpha_flux,
            "source_integral": source_integral,
            "electron_source_integral": optional_float(row.get("electron_source_integral")),
            "hole_source_integral": optional_float(row.get("hole_source_integral")),
            "sent_e_current": sent_e_current,
            "sent_h_current": sent_h_current,
            "sent_e_flux_from_J_over_q": sent_e_flux,
            "sent_h_flux_from_J_over_q": sent_h_flux,
            "sent_e_alpha": sent_e_alpha,
            "sent_h_alpha": sent_h_alpha,
            "electron_flux_over_sentaurus": electron_flux / sent_e_flux if sent_e_flux not in (None, 0.0) else None,
            "hole_flux_over_sentaurus": hole_flux / sent_h_flux if sent_h_flux not in (None, 0.0) else None,
            "electron_alpha_over_sentaurus": electron_alpha / sent_e_alpha if sent_e_alpha not in (None, 0.0) else None,
            "hole_alpha_over_sentaurus": hole_alpha / sent_h_alpha if sent_h_alpha not in (None, 0.0) else None,
        }
        selected.append(item)
    selected.sort(key=lambda item: (abs(float(item["combined_alpha_flux"])), abs(float(item["source_integral"]))), reverse=True)
    for rank, item in enumerate(selected[:limit], start=1):
        item["rank"] = rank
    return selected[:limit]


def build_active_edge_top_tables(
    support_rows: list[dict[str, Any]],
    *,
    biases: list[float],
    per_bias_limit: int = 30,
) -> list[dict[str, Any]]:
    top_rows: list[dict[str, Any]] = []
    for target_bias in biases:
        candidate_biases = sorted({
            value for value in (
                optional_float(row.get("bias_V"))
                for row in support_rows
                if optional_float(row.get("nearest_sentaurus_bias_V")) is not None and
                abs((optional_float(row.get("nearest_sentaurus_bias_V")) or 0.0) - target_bias) <= 1.0e-9
            ) if value is not None
        }, key=lambda value: abs(value - target_bias))
        if not candidate_biases:
            continue
        top_rows.extend(build_active_edge_top_rows(
            support_rows,
            target_bias=target_bias,
            nearest_vela_bias=candidate_biases[0],
            limit=per_bias_limit,
        ))
    return top_rows

def load_curve(path: Path, current_cols: list[str]) -> list[tuple[float, float]]:
    rows = read_csv_rows(path)
    if not rows:
        return []
    fields = list(rows[0])
    bias_col = next((name for name in fields if name in {"bias_V", "voltage", "V"}), fields[0])
    current_col = next((name for name in current_cols if name in fields), fields[-1])
    result = []
    for row in rows:
        try:
            result.append((float(row[bias_col]), float(row[current_col])))
        except ValueError:
            continue
    return sorted(result)


def write_curve_plot(
    reference_csv: Path,
    candidate_csv: Path,
    out_dir: Path,
    *,
    stem: str = "coarse_bv_iv_sentaurus_vs_vela_previous_full20",
    reference_alias_name: str = "sentaurus_coarse_bv_reference.csv",
) -> dict[str, Any]:
    mpl_config = out_dir / ".matplotlib"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))
    import matplotlib.pyplot as plt

    ref = load_curve(reference_csv, ["current_total", "current"])
    cand = load_curve(candidate_csv, ["current_total_A_per_um", "current_total", "current"])
    reference_alias = out_dir / reference_alias_name
    png = out_dir / f"{stem}.png"
    svg = out_dir / f"{stem}.svg"
    points_csv = out_dir / f"{stem}_points.csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(reference_csv, reference_alias)

    with points_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bias_V", "sentaurus_current_A", "vela_current_A_per_um"])
        writer.writeheader()
        cand_by_bias = {bias_key(b): v for b, v in cand}
        for b, value in ref:
            writer.writerow({
                "bias_V": b,
                "sentaurus_current_A": value,
                "vela_current_A_per_um": cand_by_bias.get(bias_key(b), ""),
            })

    if ref and cand:
        fig, ax = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
        ax.semilogy([abs(b) for b, _ in ref], [max(abs(i), 1.0e-300) for _, i in ref], label="Sentaurus coarse |I| (A)")
        ax.semilogy([abs(b) for b, _ in cand], [max(abs(i), 1.0e-300) for _, i in cand], marker="o", label="Vela previous full20 |I| (A/um)")
        ax.set_xlabel("Reverse bias |V| (V)")
        ax.set_ylabel("Absolute terminal current")
        ax.set_title("PN2D coarse BV diagnostic current comparison")
        ax.grid(True, which="both", linewidth=0.45, alpha=0.35)
        ax.legend()
        fig.savefig(png, dpi=220)
        fig.savefig(svg)
        plt.close(fig)
    return {
        "sentaurus_coarse_bv_reference_csv": str(reference_alias),
        "png": str(png),
        "svg": str(svg),
        "points_csv": str(points_csv),
    }


def write_rootcause_summary(path: Path, summary: dict[str, Any]) -> None:
    def existing_path(raw: Any) -> Path | None:
        if not raw:
            return None
        candidate = Path(str(raw))
        return candidate if candidate.exists() else None

    def nearest(points: list[tuple[float, float]], target: float) -> tuple[float, float] | None:
        if not points:
            return None
        return min(points, key=lambda item: abs(item[0] - target))

    curve = summary.get("curve_plot", {}) if isinstance(summary.get("curve_plot"), dict) else {}
    ref_path = existing_path(curve.get("sentaurus_coarse_bv_reference_csv"))
    cand_path = existing_path(summary.get("coarse_previous_full20_csv"))
    ref_curve = load_curve(ref_path, ["current_total", "current"]) if ref_path else []
    cand_curve = load_curve(cand_path, ["current_total_A_per_um", "current_total", "current"]) if cand_path else []

    current_rows: list[str] = []
    for target in summary.get("target_biases", []):
        try:
            bias = float(target)
        except (TypeError, ValueError):
            continue
        ref = nearest(ref_curve, bias)
        cand = nearest(cand_curve, bias)
        if ref is None or cand is None:
            continue
        ref_bias, ref_current = ref
        cand_bias, cand_current = cand
        ratio = ""
        if ref_current != 0.0:
            ratio = f"{abs(cand_current) / abs(ref_current):.6g}"
        current_rows.append(
            f"| {bias:g} | {ref_bias:.6g} | {ref_current:.6e} | "
            f"{cand_bias:.6g} | {cand_current:.6e} | {ratio} | {cand_bias - ref_bias:.6g} |"
        )

    node_path = existing_path(summary.get("node_field_compare_csv"))
    node_rows = read_csv_rows(node_path) if node_path else []
    field_rows: list[str] = []
    for quantity in [
        "electron_velocity",
        "hole_velocity",
        "electron_alpha_avalanche",
        "hole_alpha_avalanche",
        "electron_ion_integral",
        "hole_ion_integral",
        "mean_ion_integral",
    ]:
        selected = [row for row in node_rows if row.get("quantity") == quantity]
        sent_count = sum(1 for row in selected if row.get("sentaurus_value") not in {"", None})
        vela_count = sum(1 for row in selected if row.get("vela_value_scaled_to_sentaurus_units") not in {"", None})
        field_rows.append(f"| {quantity} | {len(selected)} | {sent_count} | {vela_count} |")

    support_path = existing_path(summary.get("current_support_compare_csv"))
    support_rows = read_csv_rows(support_path) if support_path else []
    support_fields = [
        "sent_e_velocity",
        "sent_h_velocity",
        "sent_e_alpha",
        "sent_h_alpha",
        "sent_e_ion_integral",
        "sent_h_ion_integral",
        "sent_mean_ion_integral",
        "sent_e_current",
        "sent_h_current",
        "electron_alpha_m_inv",
        "hole_alpha_m_inv",
        "electron_flux_abs",
        "hole_flux_abs",
    ]
    support_lines = [
        f"| {field} | {sum(1 for row in support_rows if row.get(field) not in {'', None})} |"
        for field in support_fields
    ]

    active_path = existing_path(summary.get("active_edge_top_csv"))
    active_rows = read_csv_rows(active_path) if active_path else []
    active_lines: list[str] = []
    for row in active_rows:
        try:
            target_bias = float(row.get("target_sentaurus_bias_V", "nan"))
            rank = int(float(row.get("rank", "nan")))
        except (TypeError, ValueError):
            continue
        if abs(target_bias + 10.0) > 1.0e-9 or rank > 6:
            continue
        active_lines.append(
            f"| {rank} | {float(row.get('bias_V', 'nan')):.6g} | {row.get('edge_id', '')} | "
            f"{row.get('node0', '')}-{row.get('node1', '')} | {row.get('edge_class', '')} | "
            f"{float(row.get('combined_alpha_flux', 'nan')):.6e} | "
            f"{float(row.get('source_integral', 'nan')):.6e} | "
            f"{float(row.get('electron_flux_over_sentaurus', 'nan')):.6g} | "
            f"{float(row.get('electron_alpha_over_sentaurus', 'nan')):.6g} | "
            f"{float(row.get('hole_flux_over_sentaurus', 'nan')):.6g} | "
            f"{float(row.get('hole_alpha_over_sentaurus', 'nan')):.6g} |"
        )

    lines = [
        "# PN2D coarse 7x3 previous-full20 diagnostic",
        "",
        "This report is a diagnostic coarse-mesh comparison. Do not use it as final BV curve parity evidence.",
        "",
        "## Gates",
        "",
        f"- Sentaurus field gate: `{summary.get('sentaurus_field_gate', {}).get('status')}`",
        f"- Missing Sentaurus fields: `{summary.get('sentaurus_field_gate', {}).get('missing_fields')}`",
        f"- Vela run status: `{summary.get('vela_run_status')}`",
        "",
        "## Nearest-Bias Current Comparison",
        "",
        "| target V | Sentaurus V | Sentaurus I | Vela V | Vela I | abs(Vela/Sentaurus) | Vela-Sentaurus dV |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        *current_rows,
        "",
        "## Added Sentaurus Field Coverage",
        "",
        "The Sentaurus exports include velocity, avalanche alpha, and ion-integral fields. Vela VTK now exports matching direct node diagnostics for velocity, alpha, and local alpha-length ion-integral support; remaining gaps should be interpreted as nearest-bias or field-availability differences.",
        "",
        "| quantity | node rows | Sentaurus values | Vela direct node values |",
        "| --- | ---: | ---: | ---: |",
        *field_rows,
        "",
        "## Active Edge Top (-10 V target)",
        "",
        "Top Vela active edges ranked by combined `alpha*flux`, using nearest Vela bias to Sentaurus -10 V.",
        "",
        "| rank | Vela V | edge | nodes | class | combined alpha*flux | source integral | e flux/Sent | e alpha/Sent | h flux/Sent | h alpha/Sent |",
        "| ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        *active_lines,
        "",
        "## Current-Support Coverage",
        "",
        "The edge table joins Vela SG flux/source diagnostics with endpoint-averaged Sentaurus current, velocity, alpha, and ion-integral support on the same imported mesh node ids.",
        "",
        "| field | non-empty rows |",
        "| --- | ---: |",
        *support_lines,
        "",
        "## Artifacts",
        "",
    ]
    for key in [
        "coarse_previous_full20_csv",
        "curve_plot",
        "node_field_compare_csv",
        "current_support_compare_csv",
        "active_edge_top_csv",
    ]:
        if key in summary:
            lines.append(f"- {key}: `{summary[key]}`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_CASE_ROOT / "source")
    parser.add_argument("--reference-config", type=Path, default=DEFAULT_REFERENCE_CONFIG)
    parser.add_argument("--import-output-dir", type=Path, default=DEFAULT_BUILD_ROOT / "imported_reference")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_BUILD_ROOT / "reports" / "coarse_previous_full20")
    parser.add_argument("--sentaurus-export-root", type=Path, default=None)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--tdr-importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--biases", default="0,-1,-5,-10,-16,-18,-20")
    parser.add_argument(
        "--aligned-fixed-biases",
        default="",
        help="Run Vela only at these exact bias points and write *_aligned outputs.",
    )
    parser.add_argument("--skip-import", action="store_true")
    parser.add_argument("--skip-vela-run", action="store_true")
    parser.add_argument("--skip-multibias-export", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    aligned_bias_points = parse_biases(args.aligned_fixed_biases) if args.aligned_fixed_biases else None
    biases = aligned_bias_points if aligned_bias_points is not None else parse_biases(args.biases)
    run_suffix = "_aligned" if aligned_bias_points is not None else ""
    run_stem = f"coarse_previous_full20{run_suffix}"
    vtk_subdir = f"vtk{run_suffix}"
    args.source_dir = normalize_sentaurus_source_dir(args.source_dir).resolve()
    ensure_sentaurus_source_inputs(args.source_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    vtk_dir = args.out_dir / vtk_subdir
    vtk_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_import:
        run_command([
            sys.executable,
            str(REPO / "scripts" / "sentaurus_import.py"),
            "reference",
            "--config",
            str(args.reference_config),
            "--source-dir",
            str(args.source_dir),
            "--output-dir",
            str(args.import_output_dir),
            "--tdr-importer",
            str(args.tdr_importer),
            "--runner",
            str(args.runner),
            "--skip-vela-run",
        ], REPO)

    sentaurus_export_root = args.sentaurus_export_root or (args.out_dir / "sentaurus_multibias")
    sentaurus_exports: dict[float, Path] = {}
    if not args.skip_multibias_export and args.tdr_importer.exists():
        sentaurus_exports = export_selected_multibias(
            source_dir=args.source_dir,
            out_root=sentaurus_export_root,
            tdr_importer=args.tdr_importer,
            biases=biases,
        )
    elif sentaurus_export_root.exists():
        for bias in biases:
            candidate = sentaurus_export_root / f"sentaurus_{signed_bias_token(bias)}v"
            if (candidate / "field_manifest.json").exists():
                sentaurus_exports[bias_key(bias)] = candidate

    base_config = args.import_output_dir / "vela" / "simulation_bv.json"
    config_path = write_previous_full20_config(
        base_config=base_config,
        out_dir=args.out_dir,
        output_csv_name=f"{run_stem}.csv",
        bias_points=aligned_bias_points,
        config_name=f"simulation_{run_stem}.json",
        vtk_subdir=vtk_subdir,
        diagnostics_suffix=run_suffix,
    )
    csv_path = args.out_dir / f"{run_stem}.csv"

    vela_run_status = "skipped"
    if not args.skip_vela_run:
        run_command([str(args.runner), "--config", str(config_path)], args.out_dir)
        vela_run_status = "attempted"

    vela_vtks = discover_vela_vtks(vtk_dir)
    node_rows = node_field_compare(
        sentaurus_exports=sentaurus_exports,
        vela_vtks=vela_vtks,
        biases=biases,
        exact_vela_bias=aligned_bias_points is not None,
    )
    node_csv = args.out_dir / f"coarse_node_field_compare{run_suffix}.csv"
    write_csv_rows(node_csv, node_rows, [
        "bias_V", "vela_bias_V", "node_id", "x_um", "y_um", "quantity", "sentaurus_field",
        "sentaurus_value", "vela_field", "vela_value_scaled_to_sentaurus_units",
        "vela_over_sentaurus", "abs_diff",
    ])
    write_json(args.out_dir / f"coarse_node_field_compare{run_suffix}.json", {"rows": node_rows[:1000], "row_count": len(node_rows)})

    support_rows = current_support_compare(
        sg_edges_csv=args.out_dir / f"sg_avalanche_edges{run_suffix}.csv",
        sentaurus_exports=sentaurus_exports,
    )
    support_csv = args.out_dir / f"coarse_current_support_compare{run_suffix}.csv"
    write_csv_rows(support_csv, support_rows, [
        "bias_V", "nearest_sentaurus_bias_V", "edge_id", "node0", "node1", "edge_class",
        "source_integral", "electron_source_integral", "hole_source_integral",
        "electron_flux_abs", "hole_flux_abs", "electron_alpha_m_inv", "hole_alpha_m_inv",
        "sent_e_velocity", "sent_h_velocity", "sent_e_alpha", "sent_h_alpha",
        "sent_e_ion_integral", "sent_h_ion_integral", "sent_mean_ion_integral",
        "sent_e_current", "sent_h_current",
    ])
    write_json(args.out_dir / f"coarse_current_support_compare{run_suffix}.json", {"rows": support_rows[:1000], "row_count": len(support_rows)})

    active_top_rows = build_active_edge_top_tables(support_rows, biases=biases, per_bias_limit=30)
    active_top_csv = args.out_dir / f"coarse_active_edge_top{run_suffix}.csv"
    write_csv_rows(active_top_csv, active_top_rows, [
        "rank", "target_sentaurus_bias_V", "bias_V", "nearest_sentaurus_bias_V",
        "edge_id", "node0", "node1", "edge_class",
        "electron_flux_abs", "hole_flux_abs", "electron_alpha_m_inv", "hole_alpha_m_inv",
        "electron_alpha_flux", "hole_alpha_flux", "combined_alpha_flux", "source_integral",
        "electron_source_integral", "hole_source_integral",
        "sent_e_current", "sent_h_current", "sent_e_flux_from_J_over_q", "sent_h_flux_from_J_over_q",
        "sent_e_alpha", "sent_h_alpha", "electron_flux_over_sentaurus", "hole_flux_over_sentaurus",
        "electron_alpha_over_sentaurus", "hole_alpha_over_sentaurus",
    ])
    write_json(args.out_dir / f"coarse_active_edge_top{run_suffix}.json", {"rows": active_top_rows, "row_count": len(active_top_rows)})

    curve_reference = args.import_output_dir / "reference_curves" / f"{COARSE_CASE}_bv_reference.csv"
    curve = write_curve_plot(
        curve_reference,
        csv_path,
        args.out_dir,
        stem=f"coarse_bv_iv_sentaurus_vs_vela_previous_full20{run_suffix}",
        reference_alias_name=f"sentaurus_coarse_bv_reference{run_suffix}.csv",
    ) if curve_reference.exists() and csv_path.exists() else {}
    field_gate = validate_required_sentaurus_fields(list(sentaurus_exports.values()))
    summary = {
        "schema": "vela.pn2d_coarse7x3_previous_full20_diagnostic.v1",
        "config": str(config_path),
        "coarse_previous_full20_csv": str(csv_path),
        "sentaurus_field_gate": field_gate,
        "vela_run_status": vela_run_status,
        "sentaurus_exports": {str(bias): str(path) for bias, path in sentaurus_exports.items()},
        "target_biases": biases,
        "curve_plot": curve,
        "node_field_compare_csv": str(node_csv),
        "current_support_compare_csv": str(support_csv),
        "active_edge_top_csv": str(active_top_csv),
    }
    write_json(args.out_dir / f"coarse_diagnostic_summary{run_suffix}.json", summary)
    write_rootcause_summary(args.out_dir / f"coarse_rootcause_summary{run_suffix}.md", summary)
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
