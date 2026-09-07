#!/usr/bin/env python3
"""Summarize the exact PN2D high-bias Sentaurus process-variable oracle."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.sentaurus_avalanche_replay import (
    currentplot_targets,
)
from scripts.pn2d_high_bias_process_contract import EXACT_HIGH_BIAS_V
from scripts.pn2d_sentaurus_process_run_contract import validate_case
from scripts.run_pn2d_high_bias_oracle_variant_vm import (
    ORACLE_VARIANTS,
    oracle_deck,
)
from scripts.run_pn2d_high_bias_process_probe_vm import PROCESS_FIELDS


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
PAIR = re.compile(rf"([A-Za-z0-9_]+)=({NUMBER})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-a", type=Path, required=True)
    parser.add_argument("--root-b", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def records(path: Path) -> list[dict[str, Any]]:
    result = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith("AVAL_PROBE_"):
            continue
        kind = line.split(maxsplit=1)[0].removeprefix("AVAL_PROBE_").lower()
        row: dict[str, Any] = {"kind": kind}
        for name, value in PAIR.findall(line):
            row[name] = float(value)
        result.append(row)
    return result


def vector_norm(row: dict[str, Any], x: str, y: str) -> float:
    return math.hypot(float(row[x]), float(row[y]))


def maximum(rows: list[dict[str, Any]], field: str) -> float:
    return max(abs(float(row[field])) for row in rows)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def variant_files(root: Path, variant: str) -> tuple[Path, Path, Path]:
    case = root / variant
    manifest = case / "manifest.json"
    fetched = case / "fetched"
    run = fetched / f"run_{variant}.out"
    plt = fetched / f"runtime_general_tri3_avalanche_probe_{variant}.plt"
    return manifest, run, plt


def variant_deck(root: Path, variant: str) -> Path:
    return (
        root
        / variant
        / "bundle"
        / f"runtime_general_tri3_avalanche_probe_{variant}.cmd"
    )


def implicit_template(root: Path) -> str:
    deck = variant_deck(root, "implicit_default").read_text(encoding="ascii")
    stem = "runtime_general_tri3_avalanche_probe_implicit_default"
    if deck.count(stem) != 3:
        raise ValueError(f"{root}: implicit output stem count mismatch")
    deck = deck.replace(stem, "runtime_element_avalanche_probe_default")
    if deck.count("pn2d_msh.tdr") != 2:
        raise ValueError(f"{root}: implicit mesh reference count mismatch")
    deck = deck.replace("pn2d_msh.tdr", "pn2d_minimal6.tdr")
    expanded = "\n".join(
        ("  hAlphaAvalanche", *("  " + name for name in PROCESS_FIELDS), "}")
    )
    if deck.count(expanded) != 1:
        raise ValueError(f"{root}: process Plot field block mismatch")
    return deck.replace(expanded, "  hAlphaAvalanche\n}", 1)


def validate_variant_bundle_contract(
    roots: tuple[Path, Path],
    manifests: dict[tuple[Path, str], dict[str, Any]],
) -> None:
    common_files = (
        "pn2d_msh.tdr",
        "models.par",
        "runtime_element_avalanche_probe.tcl",
    )
    for root in roots:
        baseline_hashes = manifests[(root, "implicit_default")]["bundle_sha256"]
        template = implicit_template(root)
        for variant in ORACLE_VARIANTS:
            hashes = manifests[(root, variant)]["bundle_sha256"]
            for name in common_files:
                if (
                    name not in hashes
                    or name not in baseline_hashes
                    or hashes[name] != baseline_hashes[name]
                ):
                    raise ValueError(
                        f"{root}/{variant}: uncontrolled common input difference: {name}"
                    )
            actual = variant_deck(root, variant).read_text(encoding="ascii")
            expected = oracle_deck(template, variant, EXACT_HIGH_BIAS_V)
            if actual != expected:
                raise ValueError(
                    f"{root}/{variant}: deck differs beyond declared variant controls"
                )


def reproducibility_row(
    variant: str,
    bundle_a: dict[str, str],
    bundle_b: dict[str, str],
    records_a: list[dict[str, Any]],
    records_b: list[dict[str, Any]],
    currents_a: list[dict[str, float]],
    currents_b: list[dict[str, float]],
) -> dict[str, Any]:
    if bundle_a != bundle_b:
        raise ValueError(f"{variant}: paired bundle hashes differ")
    normalized_a = [row for row in records_a if row["kind"] != "begin"]
    normalized_b = [row for row in records_b if row["kind"] != "begin"]
    if normalized_a != normalized_b:
        raise ValueError(f"{variant}: paired runtime records differ")
    if currents_a != currents_b:
        raise ValueError(f"{variant}: paired CurrentPlot rows differ")
    return {
        "variant": variant,
        "runtime_records_equal": 1,
        "currentplot_rows_equal": 1,
        "bundle_hashes_equal": 1,
    }


def process_summary(
    variant: str,
    parsed: list[dict[str, Any]],
    current_rows: list[dict[str, float]],
) -> list[dict[str, Any]]:
    groups: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    for row in parsed:
        if "bias_V" in row:
            groups[(str(row["kind"]), float(row["bias_V"]))].append(row)
    current = {float(row["bias_V"]): row for row in current_rows}
    result: list[dict[str, Any]] = []
    for bias in EXACT_HIGH_BIAS_V:
        vertices = groups[("vertex", bias)]
        elements = groups[("element", bias)]
        process = groups[("process", bias)]
        integrals = groups[("integral", bias)]
        if len(vertices) != 33 or len(elements) != 32 or len(process) != 33:
            raise ValueError(
                f"{variant}/{bias:g}: record counts "
                f"{len(vertices)}/{len(elements)}/{len(process)}"
            )
        if len(integrals) != 1:
            raise ValueError(f"{variant}/{bias:g}: integral count mismatch")
        curve = current[bias]
        result.append(
            {
                "variant": variant,
                "bias_V": bias,
                "max_abs_psi_V": maximum(vertices, "psi_V"),
                "max_abs_eQFP_V": maximum(vertices, "eQFP_V"),
                "max_abs_hQFP_V": maximum(vertices, "hQFP_V"),
                "max_n_cm3": maximum(vertices, "n_cm3"),
                "max_p_cm3": maximum(vertices, "p_cm3"),
                "max_efield_V_cm": max(
                    vector_norm(row, "efield_x_V_cm", "efield_y_V_cm")
                    for row in elements
                ),
                "max_grad_qf_n_V_cm": max(
                    vector_norm(
                        row, "grad_qf_n_x_V_cm", "grad_qf_n_y_V_cm"
                    )
                    for row in elements
                ),
                "max_grad_qf_p_V_cm": max(
                    vector_norm(
                        row, "grad_qf_p_x_V_cm", "grad_qf_p_y_V_cm"
                    )
                    for row in elements
                ),
                "max_mu_n_cm2_Vs": maximum(elements, "mu_n_cm2_Vs"),
                "max_mu_p_cm2_Vs": maximum(elements, "mu_p_cm2_Vs"),
                "max_velocity_n_cm_s": maximum(process, "velocity_n_cm_s"),
                "max_velocity_p_cm_s": maximum(process, "velocity_p_cm_s"),
                "max_current_n_A_cm2": max(
                    vector_norm(
                        row, "current_n_x_A_cm2", "current_n_y_A_cm2"
                    )
                    for row in elements
                ),
                "max_current_p_A_cm2": max(
                    vector_norm(
                        row, "current_p_x_A_cm2", "current_p_y_A_cm2"
                    )
                    for row in elements
                ),
                "max_alpha_n_cm_inv": maximum(vertices, "alpha_n_cm_inv"),
                "max_alpha_p_cm_inv": maximum(vertices, "alpha_p_cm_inv"),
                "max_generation_total_cm3_s": maximum(
                    vertices, "generation_total_cm3_s"
                ),
                "max_ion_n": maximum(process, "ion_n"),
                "max_ion_p": maximum(process, "ion_p"),
                "max_ion_mean": maximum(process, "ion_mean"),
                "source_integral_A_um": abs(
                    float(integrals[0]["qg_total_A_um"])
                ),
                "anode_total_current_A": abs(
                    float(curve["Anode TotalCurrent"])
                ),
            }
        )
    return result


def adjacent(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    metrics = [
        name
        for name in rows[0]
        if name not in {"variant", "bias_V"}
    ]
    result = []
    for before, after in zip(rows, rows[1:]):
        delta_v = abs(float(after["bias_V"]) - float(before["bias_V"]))
        for metric in metrics:
            first = abs(float(before[metric]))
            second = abs(float(after[metric]))
            ratio = second / first if first > 0.0 else None
            slope = (
                math.log(second / first) / delta_v
                if first > 0.0 and second > 0.0
                else None
            )
            result.append(
                {
                    "variant": before["variant"],
                    "from_bias_V": before["bias_V"],
                    "to_bias_V": after["bias_V"],
                    "metric": metric,
                    "ratio": ratio,
                    "log_slope_per_V": slope,
                }
            )
    return result


def main() -> int:
    args = parse_args()
    roots = (args.root_a.resolve(), args.root_b.resolve())
    output = args.output_root.resolve()
    all_summary: list[dict[str, Any]] = []
    all_adjacent: list[dict[str, Any]] = []
    reproducibility: list[dict[str, Any]] = []
    manifests: dict[tuple[Path, str], dict[str, Any]] = {}
    for variant in ORACLE_VARIANTS:
        root_records = []
        root_currents = []
        for root in roots:
            case = root / variant
            manifest = validate_case(
                case,
                experiment="pn2d_exact_high_bias_oracle_variant",
                variant=variant,
                exact_biases=EXACT_HIGH_BIAS_V,
            )
            manifests[(root, variant)] = manifest
            _, run_path, plt_path = variant_files(root, variant)
            root_records.append(records(run_path))
            root_currents.append(currentplot_targets(plt_path, EXACT_HIGH_BIAS_V))
        reproducibility.append(
            reproducibility_row(
                variant,
                manifests[(roots[0], variant)]["bundle_sha256"],
                manifests[(roots[1], variant)]["bundle_sha256"],
                root_records[0],
                root_records[1],
                root_currents[0],
                root_currents[1],
            )
        )
        summary = process_summary(variant, root_records[0], root_currents[0])
        all_summary.extend(summary)
        all_adjacent.extend(adjacent(summary))
    validate_variant_bundle_contract(roots, manifests)
    write_csv(output / "process_summary.csv", all_summary)
    write_csv(output / "adjacent_growth.csv", all_adjacent)
    write_csv(output / "reproducibility.csv", reproducibility)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
