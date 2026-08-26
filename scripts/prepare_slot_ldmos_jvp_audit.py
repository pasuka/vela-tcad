#!/usr/bin/env python3
"""Prepare a localized SLOT-LDMOS Newton JVP audit deck."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any


DEFAULT_AMPLITUDES_V = (1.0e-4, 1.0e-6, 1.0e-8)


def build_document(
    base: dict[str, Any],
    *,
    state_file: str,
    output_csv: str,
    hotspot_node: int,
    drain_bias_V: float,
    amplitudes_V: tuple[float, ...] = DEFAULT_AMPLITUDES_V,
    source_jacobian: str = "local_ad",
    mobility_model: str | None = None,
    disable_impact: bool = False,
    disable_recombination: bool = False,
    basis_nodes: list[int] | None = None,
) -> dict[str, Any]:
    document = copy.deepcopy(base)
    document["simulation_type"] = "newton_jvp_probe"
    document["state_file"] = state_file
    document["output_csv"] = output_csv
    document.pop("sweep", None)
    for contact in document["contacts"]:
        if contact["name"] == "drain":
            contact["bias"] = drain_bias_V
    impact = document["solver"].get("impact_ionization")
    if impact is not None:
        impact["source_jacobian"] = source_jacobian
    if mobility_model is not None:
        document["solver"]["mobility"]["model"] = mobility_model
    if disable_impact:
        if impact is not None:
            impact["coupling_mode"] = "postprocess_only"
    if disable_recombination:
        document["solver"]["recombination"] = []
    directions: list[dict[str, Any]] = []
    active_nodes = basis_nodes if basis_nodes is not None else [hotspot_node]
    for amplitude in amplitudes_V:
        for node in active_nodes:
          for mode in ("psi", "phin", "phip"):
            directions.append(
                {
                    "name": f"node{node}_{mode}_h{amplitude:.0e}",
                    "mode": mode,
                    "amplitude_V": amplitude,
                    "exclude_contacts": False,
                    "node_ids": [node],
                    "node_index_base": 0,
                    "adjacent_cell_rings": 0,
                }
            )
    document["directions"] = directions
    document["_jvp_audit"] = {
        "hotspot_node_zero_based": hotspot_node,
        "drain_bias_V": drain_bias_V,
        "amplitudes_V": list(amplitudes_V),
        "source_jacobian": source_jacobian,
        "response_scope": "global Jv; perturbation localized to hotspot node",
        "basis_nodes_zero_based": active_nodes,
    }
    return document


def adjacent_basis_nodes(
    mesh_document: dict[str, Any], hotspot_node: int, rings: int
) -> list[int]:
    selected = {hotspot_node}
    cells = mesh_document.get("cells", mesh_document.get("triangles", []))
    for _ in range(rings):
        expanded = set(selected)
        for cell in cells:
            nodes = set(cell["node_ids"])
            if nodes & selected:
                expanded.update(nodes)
        selected = expanded
    return sorted(selected)


def prepare(
    bundle: Path,
    base_config: str,
    state_file: str,
    output_config: str,
    output_csv: str,
    hotspot_node: int,
    drain_bias_V: float,
    source_jacobian: str,
    mobility_model: str | None,
    disable_impact: bool,
    disable_recombination: bool,
    basis_adjacent_rings: int,
) -> Path:
    bundle = bundle.resolve()
    with (bundle / base_config).open(encoding="utf-8") as handle:
        base = json.load(handle)
    basis_nodes = None
    amplitudes = DEFAULT_AMPLITUDES_V
    if basis_adjacent_rings > 0:
        mesh_path = bundle / base["mesh_file"]
        mesh_document = json.loads(mesh_path.read_text(encoding="utf-8"))
        basis_nodes = adjacent_basis_nodes(
            mesh_document, hotspot_node, basis_adjacent_rings
        )
        amplitudes = (1.0e-6,)
    document = build_document(
        base,
        state_file=state_file,
        output_csv=output_csv,
        hotspot_node=hotspot_node,
        drain_bias_V=drain_bias_V,
        source_jacobian=source_jacobian,
        mobility_model=mobility_model,
        disable_impact=disable_impact,
        disable_recombination=disable_recombination,
        amplitudes_V=amplitudes,
        basis_nodes=basis_nodes,
    )
    destination = bundle / output_config
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument(
        "--base-config",
        default="simulation_ialmob_external_bv_ialmob_off.json",
    )
    parser.add_argument("--state-file", required=True)
    parser.add_argument(
        "--output-config",
        default="diagnostics/jvp/node10236_local_ad.json",
    )
    parser.add_argument(
        "--output-csv",
        default="diagnostics/jvp/node10236_local_ad.csv",
    )
    parser.add_argument("--hotspot-node", type=int, default=10236)
    parser.add_argument("--drain-bias-V", type=float, required=True)
    parser.add_argument(
        "--source-jacobian",
        choices=("local_ad", "finite_difference", "frozen"),
        default="local_ad",
    )
    parser.add_argument("--mobility-model")
    parser.add_argument("--disable-impact", action="store_true")
    parser.add_argument("--disable-recombination", action="store_true")
    parser.add_argument("--basis-adjacent-rings", type=int, default=0)
    args = parser.parse_args()
    destination = prepare(
        args.bundle,
        args.base_config,
        args.state_file,
        args.output_config,
        args.output_csv,
        args.hotspot_node,
        args.drain_bias_V,
        args.source_jacobian,
        args.mobility_model,
        args.disable_impact,
        args.disable_recombination,
        args.basis_adjacent_rings,
    )
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
