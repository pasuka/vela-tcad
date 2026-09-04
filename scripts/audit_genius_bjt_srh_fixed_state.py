#!/usr/bin/env python3
"""Decompose Genius BJT SRH differences on frozen SDevice carrier states."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
COARSE_INDICES = (0, 10, 20, 30)
RATE_COLUMNS = (
    "qf_reconstructed_production_ni_generalized_srh_cm3_s",
    "exact_np_production_ni_generalized_srh_cm3_s",
    "exact_np_production_ni_classical_srh_cm3_s",
    "exact_np_sdevice_ni_generalized_srh_cm3_s",
    "exact_np_sdevice_ni_classical_srh_cm3_s",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = (len(ordered) - 1) * fraction
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def region_name(y_um: float) -> str:
    if y_um <= 0.42:
        return "emitter_base"
    if y_um <= 0.88:
        return "base_collector"
    return "collector_bulk"


def run_probe(
    *, runner: Path, base_config: dict[str, object], fields: Path,
    mesh: Path, doping: Path, output_root: Path, case: str, vce: float,
    variant: str,
) -> tuple[dict[str, object], list[dict[str, str]]]:
    config = copy.deepcopy(base_config)
    config.pop("sweep", None)
    config.update(
        {
            "simulation_type": "srh_fixed_state_probe",
            "mesh_file": str(mesh.resolve()),
            "node_doping_file": str(doping.resolve()),
            "materials_file": str(
                (FIXTURE / "vela" / "materials_sentaurus2022.json").resolve()
            ),
            "fixed_state_fields_dir": str(fields.resolve()),
            "output_csv": str((output_root / f"{case}_{variant}_nodes.csv").resolve()),
            "output_summary": str(
                (output_root / f"{case}_{variant}_summary.json").resolve()
            ),
        }
    )
    for contact in config["contacts"]:
        if contact["name"] == "collector":
            contact["bias"] = vce
    config["solver"]["recombination"] = ["srh"]
    if variant == "bgn_no_fermi_correction":
        config["solver"]["bandgap_narrowing"]["fermi_statistics_correction"] = False
    elif variant == "bgn_none":
        config["solver"]["bandgap_narrowing"] = {
            "model": "none", "fermi_statistics_correction": False
        }
    elif variant != "production":
        raise ValueError(f"unknown probe variant: {variant}")

    config_path = output_root / f"{case}_{variant}_probe.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        cwd=REPO,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    (output_root / f"{case}_{variant}.stdout.log").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output_root / f"{case}_{variant}.stderr.log").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr or completed.stdout)
    summary = json.loads(Path(config["output_summary"]).read_text(encoding="utf-8"))
    rows = read_csv(Path(config["output_csv"]))
    return summary, rows


def metric(rows: list[dict[str, str]], candidate: str) -> dict[str, float]:
    ref = [float(row["sdevice_srh_cm3_s"]) for row in rows]
    got = [float(row[candidate]) for row in rows]
    area = [float(row["lumped_area_um2"]) for row in rows]
    signed_ref = sum(a * value for a, value in zip(area, ref, strict=True))
    signed_got = sum(a * value for a, value in zip(area, got, strict=True))
    abs_ref = sum(a * abs(value) for a, value in zip(area, ref, strict=True))
    abs_got = sum(a * abs(value) for a, value in zip(area, got, strict=True))
    l1 = sum(
        a * abs(value - reference)
        for a, value, reference in zip(area, got, ref, strict=True)
    )
    return {
        "signed_integral_ratio": signed_got / signed_ref,
        "absolute_integral_ratio": abs_got / abs_ref,
        "normalized_l1_error": l1 / abs_ref,
    }


def regional_metrics(rows: list[dict[str, str]], candidate: str) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for region in ("emitter_base", "base_collector", "collector_bulk"):
        selected = [row for row in rows if region_name(float(row["y"])) == region]
        values = metric(selected, candidate)
        output.append(
            {"region": region, "candidate": candidate, "node_count": len(selected), **values}
        )
    return output


def field_diagnostics(rows: list[dict[str, str]]) -> dict[str, object]:
    ni_errors = [
        abs(math.log10(float(row["production_ni_eff_cm3"]) /
                       float(row["sdevice_ni_eff_cm3"])))
        for row in rows
        if float(row["production_ni_eff_cm3"]) > 0.0
        and float(row["sdevice_ni_eff_cm3"]) > 0.0
    ]
    closure_errors = [
        abs(math.log10(float(row["generalized_np_closure_ratio"])))
        for row in rows
        if float(row["generalized_np_closure_ratio"]) > 0.0
    ]
    electron_errors = [
        abs(math.log10(float(row["qf_reconstructed_electron_density_cm3"]) /
                       float(row["electron_density_cm3"])))
        for row in rows
        if float(row["qf_reconstructed_electron_density_cm3"]) > 0.0
        and float(row["electron_density_cm3"]) > 0.0
    ]
    hole_errors = [
        abs(math.log10(float(row["qf_reconstructed_hole_density_cm3"]) /
                       float(row["hole_density_cm3"])))
        for row in rows
        if float(row["qf_reconstructed_hole_density_cm3"]) > 0.0
        and float(row["hole_density_cm3"]) > 0.0
    ]
    return {
        "ni_eff_absolute_log10_error": {
            "median": percentile(ni_errors, 0.5),
            "p95": percentile(ni_errors, 0.95),
            "maximum": max(ni_errors),
        },
        "generalized_np_closure_absolute_log10_error": {
            "median": percentile(closure_errors, 0.5),
            "p95": percentile(closure_errors, 0.95),
            "maximum": max(closure_errors),
        },
        "qf_reconstructed_electron_absolute_log10_error": {
            "median": percentile(electron_errors, 0.5),
            "p95": percentile(electron_errors, 0.95),
            "maximum": max(electron_errors),
        },
        "qf_reconstructed_hole_absolute_log10_error": {
            "median": percentile(hole_errors, 0.5),
            "p95": percentile(hole_errors, 0.95),
            "maximum": max(hole_errors),
        },
    }


def markdown(report: dict[str, object]) -> str:
    lines = [
        "# Genius NPN BJT SRH fixed-state decomposition",
        "",
        "SDevice n, p, effective intrinsic density, and quasi-Fermi potentials are held fixed. "
        "Only the Vela SRH formula inputs are changed, so the table separates state, "
        "BGN-Fermi, and generalized-SRH effects without re-solving either device.",
        "",
        "| Grid | VCE (V) | Variant | Integral ratio | Normalized L1 |",
        "|---|---:|---|---:|---:|",
    ]
    for point in report["points"]:
        for candidate, values in point["metrics"].items():
            lines.append(
                f"| {point['grid']} | {point['VCE_V']:.1f} | {candidate} | "
                f"{values['absolute_integral_ratio']:.6f} | "
                f"{values['normalized_l1_error']:.6f} |"
            )
    attribution = report["vce3_attribution"]
    lines.extend(
        [
            "",
            "## VCE=3 V attribution",
            "",
            f"- Potential/QF-only production operator ratio: {attribution['qf_state_production_operator_ratio']:.6f}",
            f"- Legacy rescaled carrier-term ratio (not an equivalent field replay): {attribution['legacy_rescaled_carrier_term_ratio']:.6f}",
            f"- Exact n,p with production ni_eff ratio: {attribution['exact_np_production_ni_generalized_ratio']:.6f}",
            f"- Exact n,p,ni_eff generalized-SRH ratio: {attribution['exact_np_sdevice_ni_generalized_ratio']:.6f}",
            f"- Exact n,p,ni_eff classical-SRH ratio: {attribution['exact_np_sdevice_ni_classical_ratio']:.6f}",
            "",
            "## Interpretation",
            "",
            *[f"- {item}" for item in report["interpretation"]],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runner", type=Path, default=REPO / "build-release" / "vela_example_runner.exe"
    )
    parser.add_argument(
        "--coarse-sentaurus-root", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "sentaurus",
    )
    parser.add_argument(
        "--previous-srh-audit", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "srh_alignment" / "srh_alignment_audit.json",
    )
    parser.add_argument(
        "--refined-root", type=Path,
        default=BUILD_ROOT / "mesh_sensitivity" / "local_refined",
    )
    parser.add_argument(
        "--output-root", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "srh_fixed_state",
    )
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    coarse_config = json.loads(
        (FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json").read_text(
            encoding="utf-8"
        )
    )
    previous = json.loads(args.previous_srh_audit.read_text(encoding="utf-8"))
    previous_qf_ratios = {
        int(point["index"]): float(
            point["frozen_sdevice_state_operator_ratio_vela_over_sdevice"]
        )
        for point in previous["points"]
    }
    points: list[dict[str, object]] = []
    regional: list[dict[str, object]] = []
    production_rows_by_case: dict[str, list[dict[str, str]]] = {}
    for index in COARSE_INDICES:
        case = f"coarse_vce_{index:03d}"
        fields = args.coarse_sentaurus_root / f"vce_{index:03d}" / "fields"
        variant_results: dict[str, object] = {}
        for variant in ("production", "bgn_no_fermi_correction", "bgn_none"):
            _, rows = run_probe(
                runner=args.runner,
                base_config=coarse_config,
                fields=fields,
                mesh=FIXTURE / "vela" / "input" / "mesh.json",
                doping=FIXTURE / "vela" / "input" / "doping.csv",
                output_root=args.output_root,
                case=case,
                vce=index / 10.0,
                variant=variant,
            )
            variant_results[variant] = {
                candidate: metric(rows, candidate) for candidate in RATE_COLUMNS
            }
            if variant == "production":
                production_rows_by_case[case] = rows
                for candidate in RATE_COLUMNS:
                    regional.extend(
                        {"grid": "coarse", "VCE_V": index / 10.0, **item}
                        for item in regional_metrics(rows, candidate)
                    )
        points.append(
            {
                "grid": "coarse",
                "index": index,
                "VCE_V": index / 10.0,
                "qf_state_production_operator_ratio": variant_results["production"][
                    "qf_reconstructed_production_ni_generalized_srh_cm3_s"
                ]["absolute_integral_ratio"],
                "legacy_rescaled_carrier_term_ratio": previous_qf_ratios[index],
                "legacy_method_absolute_disagreement": abs(
                    variant_results["production"][
                        "qf_reconstructed_production_ni_generalized_srh_cm3_s"
                    ]["absolute_integral_ratio"] - previous_qf_ratios[index]
                ),
                "metrics": variant_results["production"],
                "bgn_ab": {
                    variant: result[
                        "exact_np_production_ni_generalized_srh_cm3_s"
                    ]
                    for variant, result in variant_results.items()
                },
                "field_diagnostics": field_diagnostics(
                    production_rows_by_case[case]
                ),
            }
        )

    refined_config = json.loads(
        (args.refined_root / "vela" / "configs" / "m1_spatial_vce3.json").read_text(
            encoding="utf-8"
        )
    )
    refined_variants: dict[str, object] = {}
    refined_rows: list[dict[str, str]] = []
    for variant in ("production", "bgn_no_fermi_correction", "bgn_none"):
        _, rows = run_probe(
            runner=args.runner,
            base_config=refined_config,
            fields=args.refined_root / "sentaurus_exports" / "vce3" / "fields",
            mesh=args.refined_root / "vela" / "input" / "mesh.json",
            doping=args.refined_root / "vela" / "input" / "doping.csv",
            output_root=args.output_root,
            case="refined_vce_030",
            vce=3.0,
            variant=variant,
        )
        refined_variants[variant] = {
            candidate: metric(rows, candidate) for candidate in RATE_COLUMNS
        }
        if variant == "production":
            refined_rows = rows
            for candidate in RATE_COLUMNS:
                regional.extend(
                    {"grid": "refined", "VCE_V": 3.0, **item}
                    for item in regional_metrics(rows, candidate)
                )
    points.append(
        {
            "grid": "refined",
            "index": 30,
            "VCE_V": 3.0,
            "qf_state_production_operator_ratio": refined_variants["production"][
                "qf_reconstructed_production_ni_generalized_srh_cm3_s"
            ]["absolute_integral_ratio"],
            "metrics": refined_variants["production"],
            "bgn_ab": {
                variant: result[
                    "exact_np_production_ni_generalized_srh_cm3_s"
                ]
                for variant, result in refined_variants.items()
            },
            "field_diagnostics": field_diagnostics(refined_rows),
        }
    )

    coarse_vce3 = next(
        point for point in points if point["grid"] == "coarse" and point["index"] == 30
    )
    metrics = coarse_vce3["metrics"]
    report: dict[str, object] = {
        "schema_version": 1,
        "scope": "SDevice exact n/p/ni_eff fixed-state SRH decomposition, generalized/classical and BGN-Fermi A/B, with coarse/refined mesh comparison",
        "points": points,
        "regional_breakdown": regional,
        "vce3_attribution": {
            "qf_state_production_operator_ratio": coarse_vce3[
                "qf_state_production_operator_ratio"
            ],
            "legacy_rescaled_carrier_term_ratio": coarse_vce3[
                "legacy_rescaled_carrier_term_ratio"
            ],
            "exact_np_production_ni_generalized_ratio": metrics[
                "exact_np_production_ni_generalized_srh_cm3_s"
            ]["absolute_integral_ratio"],
            "exact_np_sdevice_ni_generalized_ratio": metrics[
                "exact_np_sdevice_ni_generalized_srh_cm3_s"
            ]["absolute_integral_ratio"],
            "exact_np_sdevice_ni_classical_ratio": metrics[
                "exact_np_sdevice_ni_classical_srh_cm3_s"
            ]["absolute_integral_ratio"],
        },
    }
    attribution = report["vce3_attribution"]
    report["interpretation"] = [
        "The potential/QF-only direct production result is the valid operator baseline. The earlier rescaled carrier-term audit is retained only as provenance and is not an equivalent SRH field replay.",
        "Replacing production ni_eff by the exported SDevice ni_eff isolates OldSlotboom plus Fermi-correction differences.",
        "The generalized-versus-classical result at identical n, p, and ni_eff isolates the SRH carrier-statistics formula.",
        "Regional metrics identify whether the residual is concentrated in the emitter-base junction, base-collector junction, or collector bulk.",
        "The refined-grid 3 V replay separates mesh sensitivity from the frozen local constitutive formula.",
    ]
    report["attribution_deltas"] = {
        "exact_density_ratio_change":
            attribution["exact_np_production_ni_generalized_ratio"]
            - attribution["qf_state_production_operator_ratio"],
        "sdevice_ni_ratio_change":
            attribution["exact_np_sdevice_ni_generalized_ratio"]
            - attribution["exact_np_production_ni_generalized_ratio"],
        "classical_formula_ratio_change":
            attribution["exact_np_sdevice_ni_classical_ratio"]
            - attribution["exact_np_sdevice_ni_generalized_ratio"],
    }
    (args.output_root / "srh_fixed_state_decomposition.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    with (args.output_root / "srh_fixed_state_regional.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(regional[0]))
        writer.writeheader()
        writer.writerows(regional)
    (args.output_root / "srh_fixed_state_decomposition.md").write_text(
        markdown(report), encoding="utf-8"
    )
    print(json.dumps(report["vce3_attribution"], allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
