#!/usr/bin/env python3
"""Gate and localize Genius NPN BJT spatial-state differences."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


DENSITY_FLOOR_CM3 = 1.0e-300


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_scalar(path: Path, column: str) -> list[float]:
    rows = read_csv(path)
    values = [0.0] * len(rows)
    for row in rows:
        node_id = int(row["node_id"])
        if not 0 <= node_id < len(rows):
            raise ValueError(f"non-contiguous node id {node_id} in {path}")
        values[node_id] = float(row[column])
    return values


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("cannot calculate a percentile of an empty population")
    position = (len(ordered) - 1) * fraction
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def error_statistics(errors: list[float], mask: list[bool]) -> dict[str, float | int]:
    selected = [abs(value) for value, keep in zip(errors, mask, strict=True) if keep]
    if not selected:
        raise ValueError("spatial-field mask selected no nodes")
    return {
        "node_count": len(selected),
        "rmse": math.sqrt(sum(value * value for value in selected) / len(selected)),
        "p95_absolute_error": percentile(selected, 0.95),
        "maximum_absolute_error": max(selected),
    }


def field_mask(reference: list[float], contract: dict[str, object]) -> list[bool]:
    mask = contract["mask"]
    mask_type = str(mask["type"])
    if mask_type == "all_common_nodes":
        return [True] * len(reference)
    if mask_type == "sentaurus_reference_min_cm3":
        minimum = float(mask["minimum_cm3"])
        if not math.isfinite(minimum) or minimum <= 0.0:
            raise ValueError("density-mask minimum must be finite and positive")
        return [value >= minimum for value in reference]
    raise ValueError(f"unsupported spatial-field mask: {mask_type}")


def evaluate_gate(
    errors: list[float], reference: list[float], contract: dict[str, object]
) -> dict[str, object]:
    mask = field_mask(reference, contract)
    observed = error_statistics(errors, mask)
    thresholds = {
        "rmse": float(contract["maximum_rmse"]),
        "p95_absolute_error": float(contract["maximum_p95_absolute_error"]),
        "maximum_absolute_error": float(contract["maximum_absolute_error"]),
    }
    if any(not math.isfinite(value) or value < 0.0 for value in thresholds.values()):
        raise ValueError("spatial-field thresholds must be finite and non-negative")
    checks = {
        name: float(observed[name]) <= threshold
        for name, threshold in thresholds.items()
    }
    return {
        "units": contract["units"],
        "error": contract["error"],
        "mask": contract["mask"],
        "selected_node_count": sum(mask),
        "excluded_node_count": len(mask) - sum(mask),
        "thresholds": thresholds,
        "observed": observed,
        "checks": checks,
        "pass": all(checks.values()),
    }


def density_decade_errors(actual: list[float], reference: list[float]) -> list[float]:
    return [
        math.log10(max(value, DENSITY_FLOOR_CM3))
        - math.log10(max(target, DENSITY_FLOOR_CM3))
        for value, target in zip(actual, reference, strict=True)
    ]


def density_bins(
    reference: list[float], errors: list[float]
) -> list[dict[str, float | int | str | None]]:
    edges = (-math.inf, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, math.inf)
    output = []
    logs = [math.log10(max(value, DENSITY_FLOOR_CM3)) for value in reference]
    for lower, upper in zip(edges[:-1], edges[1:], strict=True):
        selected = [
            abs(error)
            for value, error in zip(logs, errors, strict=True)
            if lower <= value < upper
        ]
        if not selected:
            continue
        lower_label = "-inf" if not math.isfinite(lower) else f"{lower:g}"
        upper_label = "+inf" if not math.isfinite(upper) else f"{upper:g}"
        output.append(
            {
                "reference_log10_cm3_range": f"[{lower_label}, {upper_label})",
                "node_count": len(selected),
                "rmse_decade": math.sqrt(
                    sum(value * value for value in selected) / len(selected)
                ),
                "p95_absolute_error_decade": percentile(selected, 0.95),
                "maximum_absolute_error_decade": max(selected),
            }
        )
    return output


def region_statistics(
    coordinates: list[tuple[float, float]], errors: list[float]
) -> dict[str, dict[str, float | int]]:
    regions = {
        "top_device_y_le_0p5um": lambda x, y: y <= 0.5,
        "bulk_y_gt_0p5um": lambda x, y: y > 0.5,
        "emitter_window_x2p75_4p25_y_le_0p5um":
            lambda x, y: 2.75 <= x <= 4.25 and y <= 0.5,
        "outside_emitter_window_y_le_0p5um":
            lambda x, y: y <= 0.5 and not 2.75 <= x <= 4.25,
    }
    return {
        name: error_statistics(
            errors,
            [selector(x, y) for x, y in coordinates],
        )
        for name, selector in regions.items()
    }


def top_error_nodes(
    coordinates: list[tuple[float, float]],
    reference: list[float],
    actual: list[float],
    errors: list[float],
    count: int = 20,
) -> list[dict[str, float | int]]:
    ranked = sorted(range(len(errors)), key=lambda index: abs(errors[index]), reverse=True)
    return [
        {
            "node_id": index,
            "x_um": coordinates[index][0],
            "y_um": coordinates[index][1],
            "sentaurus_holes_cm3": reference[index],
            "vela_holes_cm3": actual[index],
            "log10_error_decade": errors[index],
        }
        for index in ranked[:count]
    ]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, report: dict[str, object]) -> None:
    bias = report["bias"]
    status = report["comparison_status"]
    lines = [
        "# Genius NPN BJT spatial-state acceptance",
        "",
        f"This {status.replace('_', ' ')} comparison uses the exact common mesh at "
        f"VBE={bias['VBE_V']:.2f} V and VCE={bias['VCE_V']:.2f} V.",
        "Potential is checked on all nodes. Carrier-density decade errors are gated only where the Sentaurus reference density is at least 1e10 cm^-3.",
        "Full-domain density statistics remain available as characterization and cannot override the registered gate.",
        "",
        "| Field | Masked nodes | RMSE | P95 | Maximum | Pass |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    labels = {
        "electrostatic_potential": "Electrostatic potential [V]",
        "electron_density": "Electron density [decade]",
        "hole_density": "Hole density [decade]",
    }
    for name in labels:
        item = report["fields"][name]
        observed = item["gate"]["observed"]
        lines.append(
            f"| {labels[name]} | {item['gate']['selected_node_count']} | "
            f"{observed['rmse']:.6g} | {observed['p95_absolute_error']:.6g} | "
            f"{observed['maximum_absolute_error']:.6g} | {item['gate']['pass']} |"
        )
    localization = report["hole_density_localization"]
    lines.extend(
        [
            "",
            f"Overall spatial-state pass: **{report['overall_pass']}**",
            "",
            "## Hole-density localization",
            "",
            f"The largest full-domain error is {localization['full_domain']['maximum_absolute_error']:.6g} decade; "
            f"the registered reference-significant mask reduces it to "
            f"{report['fields']['hole_density']['gate']['observed']['maximum_absolute_error']:.6g} decade.",
            "",
            "The highest errors are concentrated in reference-low-density nodes; see `hole_density_top_errors.csv` and the density-bin table in the JSON report.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--sentaurus-fields-root", type=Path, required=True)
    parser.add_argument("--vela-state", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--vce", type=float, default=3.0)
    parser.add_argument(
        "--characterization-only",
        action="store_true",
        help="Report the registered metrics without treating threshold misses as process failure",
    )
    args = parser.parse_args()

    threshold_path = args.reference_root / "contracts" / "comparison_thresholds.json"
    contract_document = json.loads(threshold_path.read_text(encoding="utf-8"))
    contract = contract_document["wp3_wp5_vela_comparison"]["spatial_state_gate"]
    if contract.get("status") != "asserted":
        raise ValueError("spatial-state gate must be asserted")

    nodes = read_csv(args.sentaurus_fields_root / "nodes.csv")
    coordinates = [(float(row["x_um"]), float(row["y_um"])) for row in nodes]
    field_root = args.sentaurus_fields_root / "fields"
    sentaurus_potential = read_scalar(
        field_root / "ElectrostaticPotential_region0.csv", "component0"
    )
    sentaurus_electrons = read_scalar(field_root / "eDensity_region0.csv", "component0")
    sentaurus_holes = read_scalar(field_root / "hDensity_region0.csv", "component0")
    vela_potential = read_scalar(args.vela_state, "psi")
    vela_electrons = [
        value / 1.0e6 for value in read_scalar(args.vela_state, "electrons_m3")
    ]
    vela_holes = [value / 1.0e6 for value in read_scalar(args.vela_state, "holes_m3")]
    lengths = {
        len(coordinates), len(sentaurus_potential), len(sentaurus_electrons),
        len(sentaurus_holes), len(vela_potential), len(vela_electrons), len(vela_holes)
    }
    if len(lengths) != 1:
        raise ValueError(f"spatial input length mismatch: {sorted(lengths)}")

    errors = {
        "electrostatic_potential": [
            actual - reference
            for actual, reference in zip(vela_potential, sentaurus_potential, strict=True)
        ],
        "electron_density": density_decade_errors(vela_electrons, sentaurus_electrons),
        "hole_density": density_decade_errors(vela_holes, sentaurus_holes),
    }
    references = {
        "electrostatic_potential": sentaurus_potential,
        "electron_density": sentaurus_electrons,
        "hole_density": sentaurus_holes,
    }
    fields: dict[str, object] = {}
    for name in errors:
        gate = evaluate_gate(errors[name], references[name], contract["fields"][name])
        fields[name] = {
            "gate": gate,
            "full_domain_characterization": error_statistics(
                errors[name], [True] * len(errors[name])
            ),
        }

    hole_localization = {
        "full_domain": fields["hole_density"]["full_domain_characterization"],
        "by_sentaurus_reference_density": density_bins(sentaurus_holes, errors["hole_density"]),
        "by_geometric_region": region_statistics(coordinates, errors["hole_density"]),
    }
    report = {
        "schema_version": 1,
        "device": "Genius NPN BJT",
        "bias": {"VBE_V": 0.7, "VCE_V": args.vce},
        "comparison_status": (
            "characterization_only" if args.characterization_only else "asserted"
        ),
        "sampling": contract["sampling"],
        "contract_reason": contract["reason"],
        "common_node_count": len(coordinates),
        "fields": fields,
        "hole_density_localization": hole_localization,
        "overall_pass": all(fields[name]["gate"]["pass"] for name in fields),
        "source_paths": {
            "threshold_contract": str(threshold_path),
            "sentaurus_field_manifest": str(args.sentaurus_fields_root / "field_manifest.json"),
            "sentaurus_nodes": str(args.sentaurus_fields_root / "nodes.csv"),
            "vela_state": str(args.vela_state),
        },
        "source_sha256": {
            "threshold_contract": sha256(threshold_path),
            "sentaurus_field_manifest": sha256(args.sentaurus_fields_root / "field_manifest.json"),
            "sentaurus_nodes": sha256(args.sentaurus_fields_root / "nodes.csv"),
            "vela_state": sha256(args.vela_state),
        },
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / "spatial_comparison_summary.json"
    report_path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    write_markdown(args.output_dir / "spatial_comparison_summary.md", report)
    top_rows = top_error_nodes(
        coordinates, sentaurus_holes, vela_holes, errors["hole_density"]
    )
    write_csv(args.output_dir / "hole_density_top_errors.csv", top_rows)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["overall_pass"] or args.characterization_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
