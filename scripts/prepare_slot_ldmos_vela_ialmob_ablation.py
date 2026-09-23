#!/usr/bin/env python3
"""Prepare strict Vela Slot-LDMOS Masetti/Enhanced-Lombardi A/B decks."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


DEFAULT_BASE_CONFIG = Path("simulation_05_avalanche_on_60v.json")
DEFAULT_STAGE05_IV = Path("outputs/stages/05_avalanche_on_60v/iv.csv")
SILICON_REGION = "Silicon_1"
OXIDE_REGION = "Oxide_1"


class IALMobPreparationError(ValueError):
    """Raised when a strict IALMob A/B contract cannot be established."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_interface_edges(
    mesh: dict[str, Any], region_a: str, region_b: str
) -> int:
    regions = {int(row["id"]): row["name"] for row in mesh["regions"]}
    if region_a not in regions.values() or region_b not in regions.values():
        raise IALMobPreparationError(
            f"mesh must contain regions {region_a!r} and {region_b!r}"
        )
    owners: dict[tuple[int, int], list[str]] = {}
    for cell in mesh["triangles"]:
        nodes = [int(value) for value in cell["node_ids"]]
        region = regions[int(cell["region_id"])]
        for index in range(3):
            edge = tuple(sorted((nodes[index], nodes[(index + 1) % 3])))
            owners.setdefault(edge, []).append(region)
    wanted = {region_a, region_b}
    return sum(1 for adjacent in owners.values() if set(adjacent) == wanted)


def replace_output_prefix(value: Any, source: str, target: str) -> Any:
    if isinstance(value, dict):
        return {
            key: replace_output_prefix(item, source, target)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [replace_output_prefix(item, source, target) for item in value]
    if isinstance(value, str):
        return value.replace(source, target)
    return value


def build_probe_case(
    base: dict[str, Any],
    case: str,
    initial_inner_voltage_V: float,
    resume_boundary_control: bool = False,
) -> dict[str, Any]:
    if case not in {"ialmob_off", "ialmob_on"}:
        raise IALMobPreparationError(f"unknown case {case!r}")
    source_output = "outputs/stages/05_avalanche_on_60v"
    target_output = f"outputs/ialmob_ablation/probe_60v/{case}"
    document = replace_output_prefix(copy.deepcopy(base), source_output, target_output)
    document["_comment"] = (
        "Strict 60 V outer-load-line IALMob A/B probe. Both cases share the "
        "same mesh, Stage 05 initial state, avalanche model, resistor, and solver."
    )
    document["_ialmob_ablation"] = {
        "case": case,
        "controlled_delta": (
            "masetti_field versus masetti_field_lombardi at "
            f"{SILICON_REGION}/{OXIDE_REGION}"
        ),
        "handoff": "direct_newton_from_shared_converged_stage05_state",
    }

    mobility = document["solver"]["mobility"]
    mobility.pop("surface", None)
    if case == "ialmob_off":
        mobility["model"] = "masetti_field"
    else:
        mobility["model"] = "masetti_field_lombardi"
        mobility["surface"] = {
            "surface_interface": [SILICON_REGION, OXIDE_REGION]
        }

    sweep = document["sweep"]
    document["solver"]["handoff"]["gummel_max_iter"] = 0
    sweep["bias_points"] = [60.0]
    sweep["start"] = 60.0
    sweep["stop"] = 60.0
    sweep["initial_state_file"] = (
        "outputs/stages/05_avalanche_on_60v/final_state.h5"
    )
    sweep["external_circuit"]["initial_inner_voltage_V"] = (
        initial_inner_voltage_V
    )
    sweep["boundary_control"]["resume"] = resume_boundary_control
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = f"{target_output}/vtk/state"
    return document


def normalized_physics(document: dict[str, Any]) -> dict[str, Any]:
    """Return a copy suitable for proving that only mobility differs."""
    result = copy.deepcopy(document)
    result.pop("_comment", None)
    result.pop("_ialmob_ablation", None)
    result["solver"].pop("mobility", None)

    def normalize(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: normalize(item) for key, item in value.items()}
        if isinstance(value, list):
            return [normalize(item) for item in value]
        if isinstance(value, str):
            return value.replace("ialmob_off", "CASE").replace("ialmob_on", "CASE")
        return value

    return normalize(result)


def read_stage05_terminal(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    converged = [row for row in rows if row.get("converged") == "1"]
    if not converged:
        raise IALMobPreparationError("Stage 05 IV contains no converged row")
    row = converged[-1]
    if abs(float(row["outer_voltage_V"]) - 60.0) > 1.0e-9:
        raise IALMobPreparationError("Stage 05 terminal row must be the 60 V point")
    return row


def prepare(bundle: Path, resume_boundary_control: bool = False) -> dict[str, Any]:
    bundle = bundle.resolve()
    with (bundle / "mesh.json").open(encoding="utf-8") as handle:
        mesh = json.load(handle)
    interface_edges = count_interface_edges(mesh, SILICON_REGION, OXIDE_REGION)
    if interface_edges <= 0:
        raise IALMobPreparationError("selected Si/oxide regions share no mesh edge")

    with (bundle / DEFAULT_BASE_CONFIG).open(encoding="utf-8") as handle:
        base = json.load(handle)
    terminal = read_stage05_terminal(bundle / DEFAULT_STAGE05_IV)
    initial_inner = float(terminal["inner_voltage_V"])

    documents = {
        case: build_probe_case(
            base, case, initial_inner, resume_boundary_control
        )
        for case in ("ialmob_off", "ialmob_on")
    }
    if normalized_physics(documents["ialmob_off"]) != normalized_physics(
        documents["ialmob_on"]
    ):
        raise IALMobPreparationError(
            "IALMob A/B documents differ outside mobility and isolated outputs"
        )

    cases: list[dict[str, Any]] = []
    for case, document in documents.items():
        case_output = bundle / "outputs" / "ialmob_ablation" / "probe_60v" / case
        for directory in (
            case_output,
            case_output / "boundary_control_checkpoints",
            case_output / "states",
            case_output / "vtk",
        ):
            directory.mkdir(parents=True, exist_ok=True)
        filename = f"simulation_ialmob_probe_60v_{case}.json"
        path = bundle / filename
        path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        cases.append(
            {
                "case": case,
                "config": filename,
                "sha256": sha256(path),
                "mobility": document["solver"]["mobility"],
            }
        )

    manifest = {
        "schema": "vela.slot_ldmos.ialmob_ablation.v1",
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "shared_initial_state": (
            "outputs/stages/05_avalanche_on_60v/final_state.h5"
        ),
        "shared_initial_inner_voltage_V": initial_inner,
        "surface_interface": [SILICON_REGION, OXIDE_REGION],
        "surface_interface_edge_count": interface_edges,
        "resume_boundary_control": resume_boundary_control,
        "cases": cases,
    }
    manifest_path = bundle / "ialmob_probe_60v_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse matching boundary-control checkpoints from an interrupted probe.",
    )
    args = parser.parse_args()
    print(json.dumps(prepare(args.bundle, args.resume), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
