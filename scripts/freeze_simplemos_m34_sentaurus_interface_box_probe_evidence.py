#!/usr/bin/env python3
"""Freeze qualified M34 Sentaurus interface box-probe evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil


REPO = Path(__file__).resolve().parents[1]
BUILD = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m34_interface_box_probe"
)
PORTABLE = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "sentaurus_interface_box_probe"
)
EVIDENCE = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m34_sentaurus_interface_box_probe_evidence.json"
)
CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m34_sentaurus_interface_box_probe_contract_v1.json"
)
RUNNER = REPO / "scripts/run_simplemos_m34_sentaurus_interface_box_probe.py"
FREEZER = REPO / "scripts/freeze_simplemos_m34_sentaurus_interface_box_probe_evidence.py"
TEST = REPO / "tests/regression/test_simplemos_m34_sentaurus_interface_box_probe.py"
FILES = (
    "m34_sentaurus_interface_box_probe_report.json",
    "m34_interface_edge_coefficients.csv",
    "m34_interface_node_measures.csv",
    "m34_vertex_mapping_summary.csv",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO).as_posix()


def main() -> int:
    source = BUILD / "analysis"
    report = json.loads(
        (source / "m34_sentaurus_interface_box_probe_report.json").read_text(
            encoding="utf-8"
        )
    )
    if report["status"] != "qualified":
        raise ValueError("M34 report is not qualified")
    probe = report["probe"]
    mapping = report["vertex_mapping"]
    interface = report["si_oxide_interface"]
    conclusions = report["conclusions"]
    if probe["triangle_record_count"] != 2746:
        raise ValueError("M34 triangle coverage changed")
    if probe["input_to_debug_local_permutation"] != [1, 0, 2]:
        raise ValueError("M34 local permutation changed")
    if probe["input_to_debug_measure_local_permutation"] != [0, 2, 1]:
        raise ValueError("M34 Measure local permutation changed")
    if not mapping["contact_nodes_explain_delta"] or mapping["double_edges"] != 0:
        raise ValueError("M34 vertex expansion is no longer contact-closed")
    if interface["edge_count"] != 20 or interface["node_count"] != 21:
        raise ValueError("M34 Si/Oxide interface support changed")
    if abs(interface["vela_region_over_sentaurus_region_poisson_min"] - 1.0) > 1e-12:
        raise ValueError("M34 region-local Poisson coefficient no longer matches")
    if abs(interface["vela_region_over_sentaurus_region_poisson_max"] - 1.0) > 1e-12:
        raise ValueError("M34 region-local Poisson coefficient no longer matches")
    if conclusions["explicit_si_oxide_double_nodes_supported"]:
        raise ValueError("M34 unexpectedly claims explicit Si/Oxide double nodes")
    if conclusions["production_default_changed"]:
        raise ValueError("M34 must not change the production default")

    PORTABLE.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        shutil.copy2(source / name, PORTABLE / name)
    raw = BUILD / "sentaurus_raw/sentaurus_bundle/n23"
    sources = [CONTRACT, RUNNER, FREEZER, TEST]
    evidence = {
        "schema": "vela.simplemos.sdevice.m34_sentaurus_interface_box_probe_evidence.v1",
        "status": "passed",
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "artifact_hashes": {
            portable(PORTABLE / name): sha256(PORTABLE / name) for name in FILES
        },
        "raw_oracle_hashes": {
            "MeasureCoefficients.debug": sha256(raw / "MeasureCoefficients.debug"),
            "m34_box_probe.console.log": sha256(raw / "m34_box_probe.console.log"),
            "input_fps.tdr": sha256(raw / "input_fps.tdr"),
        },
        "summary": {
            "triangle_records": probe["triangle_record_count"],
            "input_to_debug_local_permutation": probe[
                "input_to_debug_local_permutation"
            ],
            "input_to_debug_measure_local_permutation": probe[
                "input_to_debug_measure_local_permutation"
            ],
            "input_grid_vertices": mapping["input_grid_vertices"],
            "sdevice_internal_vertices": mapping["sdevice_internal_vertices"],
            "internal_vertex_delta": mapping["internal_vertex_delta"],
            "unique_contact_nodes": mapping["unique_contact_nodes"],
            "double_edges": mapping["double_edges"],
            "si_oxide_interface_edges": interface["edge_count"],
            "si_oxide_interface_nodes": interface["node_count"],
            "region_poisson_ratio_min": interface[
                "vela_region_over_sentaurus_region_poisson_min"
            ],
            "region_poisson_ratio_max": interface[
                "vela_region_over_sentaurus_region_poisson_max"
            ],
            "legacy_poisson_ratio_median": interface[
                "vela_legacy_over_sentaurus_region_poisson_median"
            ],
            "sentaurus_total_over_si_measure_median": interface[
                "sentaurus_total_over_si_measure_median"
            ],
            "barycentric_total_over_si_measure_median": interface[
                "barycentric_total_over_si_measure_median"
            ],
            "explicit_si_oxide_double_nodes_supported": conclusions[
                "explicit_si_oxide_double_nodes_supported"
            ],
            "production_default_changed": conclusions[
                "production_default_changed"
            ],
        },
    }
    EVIDENCE.write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(evidence["summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
