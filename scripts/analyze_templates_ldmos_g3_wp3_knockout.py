#!/usr/bin/env python3
"""Import and qualify one completed Templates/LDMOS WP3 T4 knockout run.

This driver is deliberately endpoint-only.  It imports the two independently
converged Sentaurus TDRs, replays each state through the matching Vela physics,
and writes the exact raw document consumed by the frozen R4 qualifier.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable, Sequence

try:
    from scripts.audit_templates_ldmos_g3_idvg_shift_kcl import (
        drain_nodes,
        drain_sg_cut,
        read_csv,
        write_sentaurus_state_csv,
    )
    from scripts.audit_templates_ldmos_g3_residual_scaling import (
        current_scale_from_sg,
    )
    from scripts.audit_templates_ldmos_g3_hfs_gradient_ab import mesh_node_sets
    from scripts.qualify_templates_ldmos_g3_wp3_knockout import (
        FIXED_HOTSPOTS,
        VALID_VARIANTS,
        qualify,
    )
    from scripts.sentaurus_import import parse_quoted_list, parse_values_block
except ModuleNotFoundError:
    from audit_templates_ldmos_g3_idvg_shift_kcl import (  # type: ignore
        drain_nodes,
        drain_sg_cut,
        read_csv,
        write_sentaurus_state_csv,
    )
    from audit_templates_ldmos_g3_residual_scaling import (  # type: ignore
        current_scale_from_sg,
    )
    from audit_templates_ldmos_g3_hfs_gradient_ab import mesh_node_sets  # type: ignore
    from qualify_templates_ldmos_g3_wp3_knockout import (  # type: ignore
        FIXED_HOTSPOTS,
        VALID_VARIANTS,
        qualify,
    )
    from sentaurus_import import parse_quoted_list, parse_values_block  # type: ignore


Q_C = 1.602176634e-19
ENDPOINTS = ((1.0, "vg1p0"), (1.166666666666667, "vg1p166667"))
REPO = Path(__file__).resolve().parents[1]
SUMMARY_SCHEMA = REPO / "schemas/vela.templates_ldmos.g3_wp3_knockout_summary.v1.schema.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def assert_summary_schema_contract(summary: dict[str, Any], schema: dict[str, Any]) -> None:
    """Validate the closed, pre-registered R4 surface without a new dependency."""
    required = set(schema["required"])
    properties = schema["properties"]
    if schema.get("additionalProperties") is not False or set(summary) != required:
        raise ValueError(
            f"summary keys differ from closed R4 schema: {sorted(set(summary) ^ required)}"
        )
    if required != set(properties):
        raise ValueError("R4 schema required/property sets differ")
    for name in ("schema", "benchmark", "incident_edge_count",
                 "denominator_port_current_source", "residual_and_flux_unit"):
        if summary[name] != properties[name]["const"]:
            raise ValueError(f"summary field {name} violates schema const")
    for name in ("variant", "bias_V", "verdict"):
        if summary[name] not in properties[name]["enum"]:
            raise ValueError(f"summary field {name} violates schema enum")
    sha = summary["artifact_sha256"]
    if set(sha) != set(properties["artifact_sha256"]["required"]):
        raise ValueError("artifact SHA-256 keys violate schema")
    if any(
        not isinstance(value, str) or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
        for value in sha.values()
    ):
        raise ValueError("artifact SHA-256 value violates schema")


def canonical_hash(rows: Iterable[Iterable[Any]]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update((",".join(str(value) for value in row) + "\n").encode("utf-8"))
    return digest.hexdigest()


def coordinate_signature_from_mesh(mesh: dict[str, Any]) -> tuple[int, str]:
    rows = [
        (int(node["id"]), format(float(node["x"]), ".17g"),
         format(float(node["y"]), ".17g"))
        for node in sorted(mesh["nodes"], key=lambda item: int(item["id"]))
    ]
    return len(rows), canonical_hash(rows)


def coordinate_signature_from_import(path: Path) -> tuple[int, str]:
    rows = sorted(read_csv(path), key=lambda item: int(item["id"]))
    canonical = [
        (int(row["id"]), format(float(row["x_um"]), ".17g"),
         format(float(row["y_um"]), ".17g"))
        for row in rows
    ]
    return len(canonical), canonical_hash(canonical)


def edge_signature(rows: list[dict[str, str]]) -> str:
    canonical = sorted(
        (int(row["edge_id"]), min(int(row["node0"]), int(row["node1"])),
         max(int(row["node0"]), int(row["node1"])))
        for row in rows
    )
    if len(canonical) != len(set(canonical)):
        raise ValueError("transport edge probe contains duplicate edge records")
    return canonical_hash(canonical)


def matching_physics(base: dict[str, Any], variant: str) -> dict[str, Any]:
    if variant not in VALID_VARIANTS:
        raise ValueError(f"unregistered variant: {variant}")
    config = deepcopy(base)
    config.pop("sweep", None)
    solver = config["solver"]
    if variant == "k1_hfs_off":
        solver["mobility"] = {"model": "constant"}
    elif variant == "k2_boltzmann":
        solver["carrier_statistics"] = {"model": "boltzmann"}
    elif variant == "k3_bgn_off":
        solver["bandgap_narrowing"] = {"model": "none"}
    mobility_text = json.dumps(solver.get("mobility", {})).lower()
    if "ialmob" in mobility_text or config.get("sweep", {}).get("predictor") not in (
        None, False, "off",
    ):
        raise ValueError("T4 replay must keep IALMob and predictor disabled")
    return config


def endpoint_tdrs(raw: Path, variant: str) -> dict[float, Path]:
    result: dict[float, Path] = {}
    for bias, token in ENDPOINTS:
        matches = sorted(raw.glob(f"{variant}_{token}*.tdr"))
        if len(matches) != 1:
            raise ValueError(
                f"{variant} Vg={bias}: expected exactly one endpoint TDR, "
                f"found {[path.name for path in matches]}"
            )
        result[bias] = matches[0]
    return result


def sentaurus_endpoint_currents(raw: Path) -> dict[float, float]:
    candidates: dict[float, list[float]] = {bias: [] for bias, _ in ENDPOINTS}
    for path in sorted(raw.glob("*.plt")):
        text = path.read_text(encoding="utf-8", errors="ignore")
        datasets = parse_quoted_list(text, "datasets")
        if not datasets or "drain TotalCurrent" not in datasets:
            continue
        bias_name = next(
            (name for name in ("gate OuterVoltage", "gate Voltage") if name in datasets),
            None,
        )
        if bias_name is None:
            continue
        bias_index = datasets.index(bias_name)
        current_index = datasets.index("drain TotalCurrent")
        for row in parse_values_block(text, len(datasets)):
            value = float(row[bias_index])
            for endpoint, _ in ENDPOINTS:
                if abs(value - endpoint) <= 5.0e-8:
                    candidates[endpoint].append(float(row[current_index]))
    result: dict[float, float] = {}
    for bias, values in candidates.items():
        if not values:
            raise ValueError(f"Sentaurus PLT contains no exact Vg={bias} endpoint current")
        reference = values[-1]
        spread = max(abs(value - reference) for value in values)
        if spread > 1.0e-10 * max(abs(reference), 1.0e-30):
            raise ValueError(f"conflicting Sentaurus currents at Vg={bias}: {values}")
        result[bias] = reference
    return result


def runner_environment() -> dict[str, str]:
    environment = dict(os.environ)
    if os.name == "nt":
        environment["PATH"] = os.pathsep.join([
            r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin",
            environment.get("PATH", ""),
        ])
    environment["VELA_LINEAR_SOLVER"] = "sparselu"
    return environment


def run_command(argv: Sequence[str], *, cwd: Path, stdout: Path) -> None:
    completed = subprocess.run(
        list(argv), cwd=cwd, env=runner_environment(), text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    stdout.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(argv)}\n{completed.stdout}"
        )


def import_tdr(importer: Path, tdr: Path, output: Path) -> None:
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty import: {output}")
    output.mkdir(parents=True, exist_ok=True)
    run_command([
        str(importer.resolve()), "--tdr", str(tdr.resolve()),
        "--export-dir", str(output.resolve()), "--coordinate-unit", "um",
        "--compensated-doping-policy", "reported",
    ], cwd=output, stdout=output / "sentaurus_import.stdout.txt")


def set_gate_bias(config: dict[str, Any], bias: float) -> None:
    gates = [contact for contact in config["contacts"] if contact["name"].lower() == "gate"]
    if len(gates) != 1:
        raise ValueError("replay config must contain exactly one gate contact")
    gates[0]["bias"] = bias


def run_probe(
    runner: Path, base: dict[str, Any], state: Path, bias: float,
    simulation_type: str, output: Path,
) -> tuple[Path, Path]:
    output.mkdir(parents=True, exist_ok=True)
    config = deepcopy(base)
    config["simulation_type"] = simulation_type
    config["state_file"] = str(state.resolve())
    csv_path = output / ("carrier_terms.csv" if simulation_type == "newton_carrier_term_probe" else "sg_edges.csv")
    config["output_csv"] = str(csv_path.resolve())
    set_gate_bias(config, bias)
    config_path = output / ("carrier_probe.json" if simulation_type == "newton_carrier_term_probe" else "sg_probe.json")
    write_json(config_path, config)
    run_command(
        [str(runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((output / f"{simulation_type}.log").resolve())],
        cwd=output, stdout=output / f"{simulation_type}.stdout.txt",
    )
    if not csv_path.is_file():
        raise FileNotFoundError(csv_path)
    return config_path, csv_path


def state_map(path: Path) -> dict[int, dict[str, float]]:
    return {
        int(row["node_id"]): {
            key: float(row[key])
            for key in ("psi", "phin", "phip", "electrons_m3", "holes_m3")
        }
        for row in read_csv(path)
    }


def reconstructed_densities(rows: list[dict[str, str]]) -> dict[int, tuple[float, float]]:
    samples: dict[int, list[tuple[float, float]]] = {}
    for row in rows:
        for suffix in ("0", "1"):
            node = int(row[f"node{suffix}"])
            samples.setdefault(node, []).append((
                float(row[f"electron_density{suffix}_m3"]),
                float(row[f"hole_density{suffix}_m3"]),
            ))
    result: dict[int, tuple[float, float]] = {}
    for node, values in samples.items():
        first = values[0]
        for value in values[1:]:
            for candidate, reference in zip(value, first, strict=True):
                if abs(candidate - reference) / max(
                    abs(candidate), abs(reference), 1.0
                ) > 1.0e-11:
                    raise ValueError(f"inconsistent reconstructed density on node {node}")
        result[node] = first
    return result


def log_density_errors(
    state: dict[int, dict[str, float]], reconstructed: dict[int, tuple[float, float]],
    nodes: set[int],
) -> dict[str, list[float]]:
    missing = sorted(nodes - set(reconstructed))
    if missing:
        raise ValueError(f"density reconstruction misses free silicon nodes: {missing[:12]}")
    return {
        "electron": [
            math.log10(max(reconstructed[node][0], 1.0))
            - math.log10(max(state[node]["electrons_m3"], 1.0))
            for node in sorted(nodes)
        ],
        "hole": [
            math.log10(max(reconstructed[node][1], 1.0))
            - math.log10(max(state[node]["holes_m3"], 1.0))
            for node in sorted(nodes)
        ],
    }


def doping_map(path: Path) -> dict[int, float]:
    return {
        int(row["node_id"]): (
            float(row["donors_cm3"]) - float(row["acceptors_cm3"])
        ) * 1.0e6
        for row in read_csv(path)
    }


def contact_gate_inputs(
    mesh: dict[str, Any], config: dict[str, Any],
    state: dict[int, dict[str, float]], reconstructed: dict[int, tuple[float, float]],
    doping: dict[int, float],
) -> dict[str, list[float]]:
    mesh_contacts = {str(item["name"]).lower(): item for item in mesh["contacts"]}
    phin: list[float] = []
    phip: list[float] = []
    neutrality: list[float] = []
    for contact in config["contacts"]:
        if str(contact.get("type", "")).lower() != "ohmic":
            continue
        name = str(contact["name"]).lower()
        voltage = float(contact["bias"])
        for node in mesh_contacts[name]["node_ids"]:
            node = int(node)
            phin.append(state[node]["phin"] - voltage)
            phip.append(state[node]["phip"] - voltage)
            # The R4 contact-row gate qualifies the independently converged
            # endpoint boundary state.  SG edge densities are bulk-kernel
            # reconstructions and intentionally do not apply Vela's Ohmic
            # Dirichlet neutrality reconstruction; using them here would mix
            # gate (a) into gate (b).
            n = state[node]["electrons_m3"]
            p = state[node]["holes_m3"]
            net = doping[node]
            neutrality.append((n - p - net) / max(n, p, abs(net), 1.0))
    if not phin:
        raise ValueError("no ohmic contact nodes were selected for the R4 gate")
    return {
        "phin_minus_contact_V": phin,
        "phip_minus_contact_V": phip,
        "neutrality_relative_residual": neutrality,
    }


def baseline_normalized(path: Path) -> dict[float, float]:
    result = {}
    for row in read_csv(path):
        bias = float(row["bias_V"])
        for endpoint, _ in ENDPOINTS:
            if abs(bias - endpoint) <= 5.0e-8:
                result[endpoint] = float(row["seven_node_residual_over_incident_abs_flux"])
    if set(result) != {item[0] for item in ENDPOINTS}:
        raise ValueError("baseline table does not contain both T4 endpoints")
    return result


def variant_top7(residuals: dict[int, float], free_silicon: set[int]) -> list[int]:
    return sorted(free_silicon, key=lambda node: (-abs(residuals[node]), node))[:7]


def analyze_endpoint(
    *, variant: str, bias: float, tdr: Path, export: Path, output: Path,
    runner: Path, replay: dict[str, Any], mesh_path: Path, mesh: dict[str, Any],
    importer: Path, sentaurus_current: float, baseline: float,
    prepared_variant: Path, grid: Path, baseline_edges: Path,
) -> dict[str, Any]:
    import_tdr(importer, tdr, export)
    state = output / "sentaurus_state.csv"
    write_sentaurus_state_csv(export, mesh_path, state)
    state_values = state_map(state)
    carrier_config, carrier_csv = run_probe(
        runner, replay, state, bias, "newton_carrier_term_probe", output / "replay",
    )
    sg_config, sg_csv = run_probe(
        runner, replay, state, bias, "sg_edge_flux_probe", output / "replay",
    )
    carrier_rows = read_csv(carrier_csv)
    sg_rows = read_csv(sg_csv)
    scale = current_scale_from_sg(sg_rows)["A_per_um_per_scaled"]
    residuals = {
        int(row["node_id"]): float(row["electron_residual"]) * scale
        for row in carrier_rows
    }
    sets = mesh_node_sets(mesh_path)
    top = variant_top7(residuals, sets["free_silicon"])
    incident = [
        row for row in sg_rows
        if int(row["node0"]) in FIXED_HOTSPOTS or int(row["node1"]) in FIXED_HOTSPOTS
    ]
    incident_ids = [int(row["edge_id"]) for row in incident]
    if len(incident_ids) != 40 or len(set(incident_ids)) != 40:
        raise ValueError(f"{variant} Vg={bias}: fixed hotspot set does not have 40 edges")
    reconstructed = reconstructed_densities(sg_rows)
    density_errors = log_density_errors(
        state_values, reconstructed, sets["free_silicon"],
    )
    doping = doping_map(Path(replay["node_doping_file"]))
    contact = contact_gate_inputs(mesh, replay, state_values, reconstructed, doping)
    vela_current = drain_sg_cut(sg_rows, drain_nodes(mesh_path))["total_A_per_um"]
    mesh_count, mesh_coordinate_hash = coordinate_signature_from_mesh(mesh)
    imported_count, imported_coordinate_hash = coordinate_signature_from_import(
        export / "nodes.csv"
    )
    variant_edge_hash = edge_signature(sg_rows)
    baseline_edge_hash = edge_signature(read_csv(baseline_edges))
    raw = {
        "variant": variant,
        "bias_V": bias,
        "np_log10_errors_dex": density_errors,
        "contact": contact,
        "port_current_A_per_um": {
            "sentaurus": sentaurus_current,
            "vela": vela_current,
        },
        "mesh": {
            "vertex_count_equal": mesh_count == imported_count,
            "coordinate_sha256_equal": mesh_coordinate_hash == imported_coordinate_hash,
            "transport_edge_sha256_equal": baseline_edge_hash == variant_edge_hash,
        },
        "residual": {
            "fixed_seven_rows_A_per_um": [residuals[node] for node in FIXED_HOTSPOTS],
            "variant_top7_node_ids": top,
            "variant_top7_rows_A_per_um": [residuals[node] for node in top],
            "incident_silicon_edge_ids": incident_ids,
            "incident_silicon_edge_flux_A_per_um": [
                float(row["electron_particle_line_flux_per_m_s"]) * Q_C * 1.0e-6
                for row in incident
            ],
        },
        "baseline_normalized_residual": baseline,
        "artifact_sha256": {
            "deck": sha256_file(prepared_variant / "IdVg.cmd"),
            "grid": sha256_file(grid),
            "parameter": sha256_file(prepared_variant / "sdevice.par"),
            "sentaurus_tdr": sha256_file(tdr),
            "imported_state": sha256_file(state),
            "vela_replay_config": sha256_file(carrier_config),
            "transport_edges": sha256_file(sg_csv),
        },
    }
    raw_path = output / "analysis" / "raw.json"
    summary_path = output / "analysis" / "summary.json"
    write_json(raw_path, raw)
    summary = qualify(raw)
    assert_summary_schema_contract(
        summary, json.loads(SUMMARY_SCHEMA.read_text(encoding="utf-8")),
    )
    write_json(summary_path, summary)
    write_json(output / "analysis" / "provenance.json", {
        "schema": "vela.templates_ldmos.g3_wp3_knockout_endpoint_provenance.v1",
        "variant": variant,
        "bias_V": bias,
        "tdr": str(tdr.resolve()),
        "export": str(export.resolve()),
        "state": str(state.resolve()),
        "carrier_config": str(carrier_config.resolve()),
        "carrier_terms": str(carrier_csv.resolve()),
        "sg_config": str(sg_config.resolve()),
        "sg_edges": str(sg_csv.resolve()),
        "mesh_coordinate_sha256": mesh_coordinate_hash,
        "imported_coordinate_sha256": imported_coordinate_hash,
        "baseline_transport_edge_sha256": baseline_edge_hash,
        "variant_transport_edge_sha256": variant_edge_hash,
    })
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=sorted(VALID_VARIANTS), required=True)
    parser.add_argument("--capture-dir", type=Path, required=True)
    parser.add_argument("--prepared-variant", type=Path, required=True)
    parser.add_argument("--prepared-manifest", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--baseline-table", type=Path, required=True)
    parser.add_argument("--baseline-edge-vg1", type=Path, required=True)
    parser.add_argument("--baseline-edge-vg1p166667", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--importer", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    capture = args.capture_dir.resolve()
    manifest_path = capture / "state_capture_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if any(item["status"] != "pass" for item in manifest["stage_results"]):
        raise ValueError("Sentaurus state capture did not pass")
    prepared_variant = args.prepared_variant.resolve()
    if manifest["state_deck_manifest_sha256"] != sha256_file(
        prepared_variant / "state_deck_manifest.json"
    ):
        raise ValueError("capture state-deck manifest hash does not match prepared variant")
    prepared = json.loads(args.prepared_manifest.resolve().read_text(encoding="utf-8"))
    if prepared["frozen_plan_commit"] != "01f20ace88dc2b68a7bc3788534244dadf0d18f2":
        raise ValueError("T4 prepared manifest is not frozen plan 01f20ac")
    if prepared["status"] != "prepared_not_run":
        raise ValueError("unexpected prepared-manifest status")
    grid = Path(prepared["grid_file"]).resolve()
    if sha256_file(grid) != prepared["grid_file_sha256"]:
        raise ValueError("local sealed grid no longer matches prepared manifest")

    raw = capture / "raw"
    tdrs = endpoint_tdrs(raw, args.variant)
    currents = sentaurus_endpoint_currents(raw)
    baseline = baseline_normalized(args.baseline_table.resolve())
    base = json.loads(args.baseline_config.resolve().read_text(encoding="utf-8"))
    replay = matching_physics(base, args.variant)
    mesh_path = args.mesh.resolve()
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty T4 analysis: {output}")
    output.mkdir(parents=True, exist_ok=True)
    edge_paths = {
        1.0: args.baseline_edge_vg1.resolve(),
        1.166666666666667: args.baseline_edge_vg1p166667.resolve(),
    }
    summaries = []
    for bias, token in ENDPOINTS:
        point = output / token
        summaries.append(analyze_endpoint(
            variant=args.variant, bias=bias, tdr=tdrs[bias],
            export=point / "sentaurus_import", output=point,
            runner=args.runner.resolve(), replay=replay, mesh_path=mesh_path,
            mesh=mesh, importer=args.importer.resolve(),
            sentaurus_current=currents[bias], baseline=baseline[bias],
            prepared_variant=prepared_variant, grid=grid,
            baseline_edges=edge_paths[bias],
        ))
    aggregate = {
        "schema": "vela.templates_ldmos.g3_wp3_knockout_variant_analysis.v1",
        "variant": args.variant,
        "capture_manifest": str(manifest_path.resolve()),
        "endpoints": summaries,
        "all_contract_equivalent": all(
            item["verdict"] != "contract_not_equivalent" for item in summaries
        ),
        "verdicts": [item["verdict"] for item in summaries],
        "prohibited_actions": {
            "ialmob": False,
            "predictor": False,
            "idvg_31_point": False,
            "production_default_change": False,
            "reclose": False,
        },
    }
    write_json(output / "analysis" / "variant_summary.json", aggregate)
    print(json.dumps(aggregate, indent=2))
    return 0 if aggregate["all_contract_equivalent"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
