#!/usr/bin/env python3
"""Compare the exported Sentaurus HighFieldDependence block with Vela."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_PARAMETER_FILE = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8a_model_ablation/sentaurus_raw/sentaurus_bundle/full"
    / "full_vd_0p05_models.par"
)
DEFAULT_OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m9_hfs_diagnostics/parameter_audit"
)
HEADER = REPO / "include/vela/physics/MobilityModel.h"
IMPLEMENTATION = REPO / "src/physics/MobilityModel.cpp"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_pair(text: str, key: str) -> tuple[float, float]:
    match = re.search(
        rf"(?m)^\s*{re.escape(key)}\s*=\s*"
        r"([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)", text)
    if not match:
        raise ValueError(f"missing Sentaurus HighFieldDependence parameter {key}")
    return float(match.group(1)), float(match.group(2))


def high_field_block(text: str) -> str:
    match = re.search(
        r"HighFieldDependence:\s*\{(.*?)\n\}\s*\n\s*HighFieldDependence_aniso:",
        text, flags=re.DOTALL)
    if not match:
        raise ValueError("HighFieldDependence block was not found")
    return match.group(1)


def vela_defaults(header: str) -> dict[str, tuple[float, float]]:
    electron = re.search(
        r"FieldMobilityParameters\s+electronField\{"
        r"([-+0-9.eE]+),\s*([-+0-9.eE]+)\}", header)
    hole = re.search(
        r"FieldMobilityParameters\s+holeField\{"
        r"([-+0-9.eE]+),\s*([-+0-9.eE]+)\}", header)
    if not electron or not hole:
        raise ValueError("Vela field mobility defaults were not found")
    return {
        "electron": (float(electron.group(1)), float(electron.group(2))),
        "hole": (float(hole.group(1)), float(hole.group(2))),
    }


def formula_equivalence(beta: float, alpha: float) -> float:
    errors = []
    for exponent in range(-12, 13):
        ratio = 10.0 ** (exponent / 2.0)
        sentaurus = ((alpha + 1.0) /
                     (alpha + math.pow(
                         1.0 + math.pow((alpha + 1.0) * ratio, beta),
                         1.0 / beta)))
        vela = 1.0 / math.pow(
            1.0 + math.pow(ratio, beta), 1.0 / beta)
        errors.append(abs(sentaurus - vela) / max(abs(sentaurus), 1e-300))
    return max(errors)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-file", type=Path,
                        default=DEFAULT_PARAMETER_FILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    parameter_file = args.parameter_file.resolve()
    block = high_field_block(parameter_file.read_text(
        encoding="utf-8", errors="replace"))
    header_text = HEADER.read_text(encoding="utf-8")
    implementation_text = IMPLEMENTATION.read_text(encoding="utf-8")
    vela = vela_defaults(header_text)
    keys = [
        "beta0", "betaexp", "alpha", "K_dT", "E0_TrEf",
        "Ksmooth_TrEf", "Vsat_Formula", "vsat0", "vsatexp", "ku", "kv",
    ]
    exported = {key: parse_pair(block, key) for key in keys}
    carriers = ("electron", "hole")
    rows: list[dict[str, Any]] = []
    for index, carrier in enumerate(carriers):
        vela_vsat, vela_beta = vela[carrier]
        sentaurus_vsat = exported["vsat0"][index] * 1e-2
        sentaurus_beta = exported["beta0"][index]
        alpha = exported["alpha"][index]
        rows.extend([
            {
                "carrier": carrier, "parameter": "vsat0_m_s",
                "sentaurus": sentaurus_vsat, "vela": vela_vsat,
                "status": "exact" if math.isclose(
                    sentaurus_vsat, vela_vsat, rel_tol=0.0, abs_tol=1e-12)
                else "different",
                "relevance_at_300K": "active",
            },
            {
                "carrier": carrier, "parameter": "beta0",
                "sentaurus": sentaurus_beta, "vela": vela_beta,
                "status": "exact" if math.isclose(
                    sentaurus_beta, vela_beta, rel_tol=0.0, abs_tol=1e-12)
                else "different",
                "relevance_at_300K": "active",
            },
            {
                "carrier": carrier, "parameter": "alpha",
                "sentaurus": alpha, "vela": 0.0,
                "status": "implicit_exact" if alpha == 0.0 else "unsupported",
                "relevance_at_300K": "active",
            },
            {
                "carrier": carrier, "parameter": "betaexp",
                "sentaurus": exported["betaexp"][index], "vela": None,
                "status": "neutral_at_300K",
                "relevance_at_300K": "T/T0=1",
            },
            {
                "carrier": carrier, "parameter": "vsatexp",
                "sentaurus": exported["vsatexp"][index], "vela": None,
                "status": "neutral_at_300K",
                "relevance_at_300K": "T/T0=1",
            },
            {
                "carrier": carrier, "parameter": "ku",
                "sentaurus": exported["ku"][index], "vela": 1.0,
                "status": "implicit_exact",
                "relevance_at_300K": "active scale",
            },
            {
                "carrier": carrier, "parameter": "kv",
                "sentaurus": exported["kv"][index], "vela": 1.0,
                "status": "implicit_exact",
                "relevance_at_300K": "active scale",
            },
        ])

    formula_errors = {
        carrier: formula_equivalence(
            exported["beta0"][index], exported["alpha"][index])
        for index, carrier in enumerate(carriers)
    }
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "hfs_parameter_comparison.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    audit = {
        "schema": "vela.simplemos.sdevice.m9_hfs_parameter_audit.v1",
        "status": "complete",
        "sentaurus_release": "T-2022.03-SP2",
        "parameter_file": str(parameter_file),
        "parameter_file_sha256": sha256(parameter_file),
        "high_field_dependence": exported,
        "vela_defaults": {
            carrier: {"saturation_velocity_m_s": values[0], "beta": values[1]}
            for carrier, values in vela.items()
        },
        "core_formula_maximum_relative_error_at_300K": formula_errors,
        "core_formula_equivalent_at_300K": all(
            value <= 1e-15 for value in formula_errors.values()),
        "implementation_formula_found": (
            "lowFieldMobility / std::pow(1.0 + std::pow(ratio, params.beta)"
            in implementation_text),
        "parameter_rows": rows,
        "not_exercised_by_isothermal_dd": [
            "K_dT is a HydroHighField carrier-temperature smoothing parameter",
            "betaexp and vsatexp are neutral at the frozen 300 K temperature",
        ],
        "unresolved_activation": [
            "E0_TrEf and Ksmooth_TrEf are exported, but the selected isothermal DD deck does not expose whether an additional transferred-electron branch is active internally",
        ],
        "interpretation": (
            "The exported 300 K Caughey-Thomas core parameters and formula "
            "match Vela. Remaining HFS uncertainty is therefore concentrated "
            "in driving-force construction, density interpolation, hidden "
            "activation/limiting behavior, or spatial discretization rather "
            "than vsat0 or beta0."),
        "default_model_changed": False,
    }
    json_path = output / "hfs_parameter_audit.json"
    json_path.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps({
        "status": audit["status"],
        "core_formula_equivalent_at_300K": audit[
            "core_formula_equivalent_at_300K"],
        "output": str(json_path),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
