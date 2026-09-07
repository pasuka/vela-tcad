#!/usr/bin/env python3
"""Replay general-Tri3 avalanche source with staged alpha/current replacement."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_general_tri3_element_edge_avalanche import (
    geometry_rows,
    parse_log,
)
from scripts.diagnose_pn2d_general_tri3_imported_state import (
    abs_dex,
    error_summary,
)
from scripts.sentaurus_avalanche_replay import (
    currentplot_targets,
)
from scripts.pn2d_general_tri3_contract import EXACT_BIASES_V, SCHEMA_ID


CURRENT_CLOSURE_SCHEMA = (
    "pn2d_general_tri3_element_edge_current_closure/v1"
)
IMPORTED_STATE_SCHEMA = "pn2d_general_tri3_imported_state/v1"
OUTPUT_SCHEMA = "pn2d_general_tri3_element_edge_source/v1"
Q_SENT_C = 1.6021918e-19
MEASURED_SOURCE_TO_A_UM = 1.0e-12
CARRIERS = ("electron", "hole")
CURRENT_METHODS = (
    "gss_laux_truncated_support",
    "charon_whitney_hcurl_cell_average",
    "genius_least_squares_tangent",
    "box_active_edge_exact",
)
CANDIDATES = (
    "vela_contact_fallback_alpha_vela_current",
    "vela_global_qf_alpha_vela_current",
    "vela_global_electric_alpha_vela_current",
    "native_vertex_alpha_vela_current",
    "vela_contact_fallback_alpha_sentaurus_box_current",
    "vela_global_qf_alpha_sentaurus_box_current",
    "vela_global_electric_alpha_sentaurus_box_current",
    "native_vertex_alpha_sentaurus_box_current",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--imported-state-root", type=Path, required=True)
    parser.add_argument("--current-closure-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--current-method",
        choices=CURRENT_METHODS,
        default="gss_laux_truncated_support",
    )
    return parser.parse_args()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0e-300)


def integrated_qg_A_um(
    alpha_cm_inv: float,
    current_density_A_cm2: float,
    measure_um2: float,
) -> float:
    return (
        alpha_cm_inv
        * current_density_A_cm2
        * measure_um2
        * MEASURED_SOURCE_TO_A_UM
    )


def matching_driver(element_class: str) -> str:
    return (
        "electric_field"
        if element_class == "contact"
        else "quasi_fermi_gradient"
    )


def main() -> int:
    args = parse_args()
    raw_root = args.raw_root.resolve()
    imported_root = args.imported_state_root.resolve()
    current_root = args.current_closure_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    raw_manifest_path = raw_root / "manifest.json"
    raw_manifest = json.loads(raw_manifest_path.read_text(encoding="ascii"))
    if raw_manifest.get("schema") != SCHEMA_ID:
        raise ValueError("raw schema mismatch")
    if raw_manifest.get("status") != "passed":
        raise ValueError("raw root is incomplete")
    imported_manifest_path = imported_root / "analysis_manifest.json"
    imported_manifest = json.loads(
        imported_manifest_path.read_text(encoding="ascii")
    )
    if imported_manifest.get("schema") != IMPORTED_STATE_SCHEMA:
        raise ValueError("imported-state schema mismatch")
    if imported_manifest.get("source_manifest_sha256") != digest(
        raw_manifest_path
    ):
        raise ValueError("imported-state/raw manifest mismatch")
    current_manifest_path = current_root / "analysis_manifest.json"
    current_manifest = json.loads(
        current_manifest_path.read_text(encoding="ascii")
    )
    if current_manifest.get("schema") != CURRENT_CLOSURE_SCHEMA:
        raise ValueError("current-closure schema mismatch")
    if current_manifest.get("source_raw_manifest_sha256") != digest(
        raw_manifest_path
    ):
        raise ValueError("current-closure/raw manifest mismatch")

    case_names = tuple(raw_manifest["cases"])
    if len(case_names) != 1:
        raise ValueError(f"expected one case, got {case_names}")
    case_name = case_names[0]
    variant_root = raw_root / case_name / "implicit_default"
    log_path = variant_root / "fetched" / "run_implicit_default.out"
    plt_path = (
        variant_root
        / "fetched"
        / "runtime_general_tri3_avalanche_probe_implicit_default.plt"
    )
    groups = parse_log(log_path)
    geometry = geometry_rows(groups, EXACT_BIASES_V[0])
    geometry_by_element = {
        int(row["element"]): row for row in geometry
    }
    vertices = {
        (float(row["bias_V"]), int(row["vertex"])): row
        for row in groups["vertices"]
    }
    measures_by_key: dict[
        tuple[float, int],
        list[dict[str, Any]],
    ] = defaultdict(list)
    for row in groups["measures"]:
        measures_by_key[
            (float(row["bias_V"]), int(row["element"]))
        ].append(row)

    alpha_rows = read_csv(imported_root / "element_alpha_replay.csv")
    vela_alpha = {
        (
            float(row["bias_V"]),
            int(row["element"]),
            row["carrier"],
            row["forced_driver"],
        ): float(row["vela_alpha_cm_inv"])
        for row in alpha_rows
    }
    vector_rows = read_csv(
        current_root / "matching_support_cell_vectors.csv"
    )
    selected_vectors = {
        (
            float(row["bias_V"]),
            int(row["element"]),
            row["carrier"],
        ): row
        for row in vector_rows
        if row["method"] == args.current_method
    }

    element_vertex_rows: list[dict[str, Any]] = []
    cell_rows: list[dict[str, Any]] = []
    state_integrals: dict[
        tuple[float, str, str],
        float,
    ] = defaultdict(float)
    node_buckets: dict[
        tuple[float, str, str, int],
        dict[str, float],
    ] = defaultdict(
        lambda: {
            "predicted_qg_A_um": 0.0,
            "measure_um2": 0.0,
            "native_qg_A_um": 0.0,
        }
    )

    for key in sorted(measures_by_key):
        bias, element_id = key
        element_class = (
            "contact"
            if geometry_by_element[element_id]["contact_adjacent"]
            else "interior"
        )
        driver = matching_driver(element_class)
        measures = sorted(
            measures_by_key[key],
            key=lambda row: int(row["local_vertex"]),
        )
        for carrier in CARRIERS:
            vector = selected_vectors[(bias, element_id, carrier)]
            currents = {
                "vela": float(vector["vela_magnitude_A_cm2"]),
                "sentaurus_box": float(
                    vector["sentaurus_box_magnitude_A_cm2"]
                ),
            }
            element_alpha = vela_alpha[
                (bias, element_id, carrier, driver)
            ]
            global_qf_alpha = vela_alpha[
                (bias, element_id, carrier, "quasi_fermi_gradient")
            ]
            global_electric_alpha = vela_alpha[
                (bias, element_id, carrier, "electric_field")
            ]
            cell_accumulator = {
                candidate: 0.0 for candidate in CANDIDATES
            }
            for measure in measures:
                vertex_id = int(measure["vertex"])
                vertex = vertices[(bias, vertex_id)]
                measure_um2 = float(measure["measure_um2"])
                native_alpha = float(
                    vertex[
                        "alpha_n_cm_inv"
                        if carrier == "electron"
                        else "alpha_p_cm_inv"
                    ]
                )
                native_qg = float(
                    measure[
                        "qg_n_A_um"
                        if carrier == "electron"
                        else "qg_p_A_um"
                    ]
                )
                alpha_by_candidate = {
                    "vela_contact_fallback_alpha_vela_current": element_alpha,
                    "vela_global_qf_alpha_vela_current": global_qf_alpha,
                    "vela_global_electric_alpha_vela_current": global_electric_alpha,
                    "native_vertex_alpha_vela_current": native_alpha,
                    "vela_contact_fallback_alpha_sentaurus_box_current": element_alpha,
                    "vela_global_qf_alpha_sentaurus_box_current": global_qf_alpha,
                    "vela_global_electric_alpha_sentaurus_box_current": global_electric_alpha,
                    "native_vertex_alpha_sentaurus_box_current": native_alpha,
                }
                current_by_candidate = {
                    candidate: (
                        currents["sentaurus_box"]
                        if "sentaurus_box" in candidate
                        else currents["vela"]
                    )
                    for candidate in CANDIDATES
                }
                candidate_values = {
                    candidate: integrated_qg_A_um(
                        alpha_by_candidate[candidate],
                        current_by_candidate[candidate],
                        measure_um2,
                    )
                    for candidate in CANDIDATES
                }
                for candidate, predicted in candidate_values.items():
                    error = abs_dex(abs(predicted), abs(native_qg))
                    element_vertex_rows.append(
                        {
                            "case": case_name,
                            "bias_V": bias,
                            "element": element_id,
                            "local_vertex": int(measure["local_vertex"]),
                            "vertex": vertex_id,
                            "element_class": element_class,
                            "angle_class": geometry_by_element[element_id][
                                "angle_class"
                            ],
                            "carrier": carrier,
                            "candidate": candidate,
                            "driver": driver,
                            "measure_um2": measure_um2,
                            "alpha_cm_inv": alpha_by_candidate[candidate],
                            "current_magnitude_A_cm2": current_by_candidate[
                                candidate
                            ],
                            "predicted_qg_A_um": predicted,
                            "native_qg_A_um": native_qg,
                            "absolute_error_A_um": abs(
                                predicted - native_qg
                            ),
                            "absolute_error_dex": (
                                "" if error is None else error
                            ),
                            "observation_label": (
                                "box_operator_reconstruction"
                                if "sentaurus_box" in candidate
                                else "vela_recomputation"
                            ),
                        }
                    )
                    cell_accumulator[candidate] += predicted
                    state_integrals[(bias, carrier, candidate)] += predicted
                    bucket = node_buckets[
                        (bias, carrier, candidate, vertex_id)
                    ]
                    bucket["predicted_qg_A_um"] += predicted
                    bucket["measure_um2"] += measure_um2
                    bucket["native_qg_A_um"] += native_qg
            native_cell_qg = sum(
                float(
                    measure[
                        "qg_n_A_um"
                        if carrier == "electron"
                        else "qg_p_A_um"
                    ]
                )
                for measure in measures
            )
            for candidate, predicted in cell_accumulator.items():
                error = abs_dex(abs(predicted), abs(native_cell_qg))
                cell_rows.append(
                    {
                        "case": case_name,
                        "bias_V": bias,
                        "element": element_id,
                        "element_class": element_class,
                        "angle_class": geometry_by_element[element_id][
                            "angle_class"
                        ],
                        "carrier": carrier,
                        "candidate": candidate,
                        "predicted_qg_A_um": predicted,
                        "native_readmeasure_qg_A_um": native_cell_qg,
                        "absolute_error_A_um": abs(
                            predicted - native_cell_qg
                        ),
                        "absolute_error_dex": (
                            "" if error is None else error
                        ),
                    }
                )

    node_rows: list[dict[str, Any]] = []
    for (bias, carrier, candidate, vertex_id), bucket in sorted(
        node_buckets.items()
    ):
        measure_um2 = bucket["measure_um2"]
        predicted_generation = (
            bucket["predicted_qg_A_um"]
            / (Q_SENT_C * measure_um2 * MEASURED_SOURCE_TO_A_UM)
        )
        native_generation = (
            bucket["native_qg_A_um"]
            / (Q_SENT_C * measure_um2 * MEASURED_SOURCE_TO_A_UM)
        )
        error = abs_dex(
            abs(predicted_generation),
            abs(native_generation),
        )
        node_rows.append(
            {
                "case": case_name,
                "bias_V": bias,
                "vertex": vertex_id,
                "carrier": carrier,
                "candidate": candidate,
                "predicted_generation_cm3_s": predicted_generation,
                "native_generation_cm3_s": native_generation,
                "absolute_error_dex": "" if error is None else error,
            }
        )

    currentplot = {
        float(row["bias_V"]): row
        for row in currentplot_targets(
            plt_path,
            tuple(float(value) for value in EXACT_BIASES_V),
        )
    }
    raw_integrals = {
        (float(bias), carrier): sum(
            float(
                row[
                    "qg_n_A_um"
                    if carrier == "electron"
                    else "qg_p_A_um"
                ]
            )
            for row in groups["measures"]
            if float(row["bias_V"]) == float(bias)
        )
        for bias in EXACT_BIASES_V
        for carrier in CARRIERS
    }
    currentplot_names = {
        "electron": "IntegreAvalancheIntegral eAvalancheGeneration",
        "hole": "IntegrhAvalancheIntegral hAvalancheGeneration",
    }
    reference_maximum = {
        carrier: max(
            abs(raw_integrals[(float(bias), carrier)])
            for bias in EXACT_BIASES_V
        )
        for carrier in CARRIERS
    }
    state_rows: list[dict[str, Any]] = []
    identity_rows: list[dict[str, Any]] = []
    for bias in (float(value) for value in EXACT_BIASES_V):
        for carrier in CARRIERS:
            readmeasure = raw_integrals[(bias, carrier)]
            currentplot_qg = (
                float(currentplot[bias][currentplot_names[carrier]])
                * Q_SENT_C
                * MEASURED_SOURCE_TO_A_UM
            )
            identity_rows.append(
                {
                    "case": case_name,
                    "bias_V": bias,
                    "carrier": carrier,
                    "readmeasure_qg_A_um": readmeasure,
                    "currentplot_qg_A_um": currentplot_qg,
                    "absolute_error_A_um": abs(
                        readmeasure - currentplot_qg
                    ),
                    "relative_error": relative_error(
                        readmeasure,
                        currentplot_qg,
                    ),
                }
            )
            status = (
                "below_state_relative_floor"
                if abs(readmeasure)
                <= max(reference_maximum[carrier] * 1.0e-12, 1.0e-300)
                else "active"
            )
            for candidate in CANDIDATES:
                predicted = state_integrals[(bias, carrier, candidate)]
                error = abs_dex(abs(predicted), abs(readmeasure))
                state_rows.append(
                    {
                        "case": case_name,
                        "bias_V": bias,
                        "carrier": carrier,
                        "candidate": candidate,
                        "predicted_qg_A_um": predicted,
                        "native_readmeasure_qg_A_um": readmeasure,
                        "absolute_error_A_um": abs(
                            predicted - readmeasure
                        ),
                        "absolute_error_dex": (
                            "" if error is None else error
                        ),
                        "signal_status": status,
                    }
                )

    outputs = {
        "element_vertex_source.csv": element_vertex_rows,
        "cell_source_integrals.csv": cell_rows,
        "node_generation.csv": node_rows,
        "state_source_integrals.csv": state_rows,
        "source_identity_closure.csv": identity_rows,
    }
    for name, rows in outputs.items():
        write_csv(output_root / name, rows)

    candidate_summary: dict[str, Any] = {}
    for candidate in CANDIDATES:
        candidate_summary[candidate] = {}
        for carrier in CARRIERS:
            selected = [
                row for row in state_rows
                if row["candidate"] == candidate
                and row["carrier"] == carrier
                and row["signal_status"] == "active"
                and row["absolute_error_dex"] != ""
            ]
            candidate_summary[candidate][carrier] = {
                "active_integral_error_dex": error_summary(
                    [float(row["absolute_error_dex"]) for row in selected]
                ),
                "active_state_count": len(selected),
            }
    manifest = {
        "schema": OUTPUT_SCHEMA,
        "status": "valid",
        "case_name": case_name,
        "source_raw_manifest_sha256": digest(raw_manifest_path),
        "source_imported_state_manifest_sha256": digest(
            imported_manifest_path
        ),
        "source_current_closure_manifest_sha256": digest(
            current_manifest_path
        ),
        "source_log_sha256": digest(log_path),
        "source_currentplot_sha256": digest(plt_path),
        "exact_biases_V": list(EXACT_BIASES_V),
        "current_vector_method": args.current_method,
        "driver_contract": {
            "interior": "quasi_fermi_gradient",
            "contact": "electric_field_fallback_candidate",
        },
        "candidate_summary": candidate_summary,
        "source_identity": {
            "maximum_readmeasure_currentplot_relative_error": max(
                float(row["relative_error"]) for row in identity_rows
            )
        },
        "outputs": {
            name: digest(output_root / name) for name in outputs
        },
    }
    (output_root / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
