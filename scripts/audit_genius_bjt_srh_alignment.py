#!/usr/bin/env python3
"""Separate Genius BJT SRH parameter, state, and integration-weight effects."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path


Q_C = 1.602176634e-19
REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_INDICES = (0, 10, 20, 30)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_scharfetter(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Scharfetter\s*\*.*?\n\{(.*?)\n\}", text, re.DOTALL)
    if not match:
        raise ValueError(f"Scharfetter section not found in {path}")
    section = match.group(1)

    def pair(name: str) -> list[float]:
        found = re.search(
            rf"^\s*{re.escape(name)}\s*=\s*([^,\s]+)\s*,\s*([^\s#]+)",
            section,
            re.MULTILINE,
        )
        if not found:
            raise ValueError(f"{name} not found in Scharfetter section")
        return [float(found.group(1)), float(found.group(2))]

    etrap = re.search(r"^\s*Etrap\s*=\s*([^\s#]+)", section, re.MULTILINE)
    if not etrap:
        raise ValueError("Etrap not found in Scharfetter section")
    return {
        "tau_min_s": pair("taumin"),
        "tau_max_s": pair("taumax"),
        "reference_doping_cm3": pair("Nref"),
        "gamma": pair("gamma"),
        "temperature_exponent": pair("Talpha"),
        "trap_energy_eV": float(etrap.group(1)),
    }


def independent_lumped_areas(mesh_path: Path) -> list[float]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = {int(node["id"]): (float(node["x"]), float(node["y"])) for node in mesh["nodes"]}
    areas = [0.0] * len(nodes)
    for triangle in mesh["triangles"]:
        ids = [int(value) for value in triangle["node_ids"]]
        a, b, c = (nodes[node_id] for node_id in ids)
        area = 0.5 * abs(
            (b[0] - a[0]) * (c[1] - a[1])
            - (c[0] - a[0]) * (b[1] - a[1])
        )
        for node_id in ids:
            areas[node_id] += area / 3.0
    return areas


def run_srh_probe(
    *, runner: Path, base_config: dict[str, object], state: Path, output_root: Path,
    token: str, label: str, vce: float, state_fields: bool,
) -> tuple[Path, float]:
    config = json.loads(json.dumps(base_config))
    config.pop("sweep", None)
    config.update(
        {
            "simulation_type": "newton_carrier_term_probe",
            "mesh_file": str((FIXTURE / "vela" / "input" / "mesh.json").resolve()),
            "node_doping_file": str((FIXTURE / "vela" / "input" / "doping.csv").resolve()),
            "materials_file": str((FIXTURE / "vela" / "materials_sentaurus2022.json").resolve()),
            "output_csv": str((output_root / f"{token}_{label}_terms.csv").resolve()),
            "carrier_term_probe": {"solved_equation_terms": True},
        }
    )
    if state_fields:
        config["state_fields_dir"] = str(state.resolve())
        config.pop("state_file", None)
    else:
        config["state_file"] = str(state.resolve())
        config.pop("state_fields_dir", None)
    for contact in config["contacts"]:
        if contact["name"] == "collector":
            contact["bias"] = vce
    config["solver"]["recombination"] = ["srh"]
    config_path = output_root / f"{token}_{label}_probe.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    process = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        cwd=REPO,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    (output_root / f"{token}_{label}.stdout.log").write_text(process.stdout, encoding="utf-8")
    (output_root / f"{token}_{label}.stderr.log").write_text(process.stderr, encoding="utf-8")
    if process.returncode:
        raise RuntimeError(process.stderr or process.stdout)
    terms_path = Path(config["output_csv"])
    term_sum = sum(float(row["electron_recombination"]) for row in read_csv(terms_path))
    return terms_path, term_sum


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=REPO / "build-release" / "vela_example_runner.exe")
    parser.add_argument(
        "--models-par", type=Path,
        default=BUILD_ROOT / "m1_current_diagnosis" / "Silicon_T2022.03-SP2_models.par",
    )
    parser.add_argument(
        "--accepted-root", type=Path, default=BUILD_ROOT / "m1_p0_acceptance"
    )
    parser.add_argument(
        "--sentaurus-root", type=Path, default=BUILD_ROOT / "m1_p0_multibias" / "sentaurus"
    )
    parser.add_argument(
        "--source-root", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "transport_sources" / "sources",
    )
    parser.add_argument(
        "--output-root", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "srh_alignment",
    )
    parser.add_argument("--indices", type=int, nargs="+", default=list(DEFAULT_INDICES))
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    base_config = json.loads(
        (FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json").read_text(encoding="utf-8")
    )
    configured = base_config["solver"]["srh_doping_dependence"]
    printed = parse_scharfetter(args.models_par)
    expected = {
        "tau_min_s": [configured["electron"]["tau_min_s"], configured["hole"]["tau_min_s"]],
        "tau_max_s": [configured["electron"]["tau_max_s"], configured["hole"]["tau_max_s"]],
        "reference_doping_cm3": [
            configured["electron"]["reference_doping_m3"],
            configured["hole"]["reference_doping_m3"],
        ],
        "gamma": [configured["electron"]["gamma"], configured["hole"]["gamma"]],
        "temperature_exponent": [
            configured["electron_temperature_exponent"], configured["hole_temperature_exponent"]
        ],
        "trap_energy_eV": 0.0,
    }
    parameter_checks = {
        key: (
            all(math.isclose(a, b, rel_tol=1.0e-14, abs_tol=0.0) for a, b in zip(printed[key], value, strict=True))
            if isinstance(value, list)
            else math.isclose(float(printed[key]), float(value), rel_tol=0.0, abs_tol=1.0e-15)
        )
        for key, value in expected.items()
    }

    mesh_path = FIXTURE / "vela" / "input" / "mesh.json"
    independent_areas = independent_lumped_areas(mesh_path)
    points: list[dict[str, object]] = []
    for index in args.indices:
        token = f"vce_{index:03d}"
        vce = index / 10.0
        source_rows = read_csv(args.source_root / f"{token}.csv")
        source_areas = [float(row["lumped_area_um2"]) for row in source_rows]
        area_differences = [
            abs(a - b) for a, b in zip(source_areas, independent_areas, strict=True)
        ]
        vela_srh = sum(float(row["srh_integrated_A_per_um"]) for row in source_rows)
        _, vela_term_sum = run_srh_probe(
            runner=args.runner,
            base_config=base_config,
            state=args.accepted_root / "states" / f"{token}.csv",
            output_root=args.output_root,
            token=token,
            label="vela_state",
            vce=vce,
            state_fields=False,
        )
        sentaurus_fields = args.sentaurus_root / token / "fields"
        _, sentaurus_state_term_sum = run_srh_probe(
            runner=args.runner,
            base_config=base_config,
            state=sentaurus_fields,
            output_root=args.output_root,
            token=token,
            label="sdevice_state",
            vce=vce,
            state_fields=True,
        )
        scale = vela_srh / vela_term_sum
        vela_operator_on_sdevice = sentaurus_state_term_sum * scale
        sentaurus_rates = [
            float(row["component0"])
            for row in read_csv(sentaurus_fields / "srhRecombination_region0.csv")
        ]
        sentaurus_integral = Q_C * sum(
            rate * area * 1.0e-12
            for rate, area in zip(sentaurus_rates, source_areas, strict=True)
        )
        points.append(
            {
                "index": index,
                "VCE_V": vce,
                "vela_srh_integral_A_per_um": vela_srh,
                "sdevice_exported_srh_integral_with_vela_weights_A_per_um": sentaurus_integral,
                "vela_operator_on_sdevice_qf_state_A_per_um": vela_operator_on_sdevice,
                "self_consistent_integral_ratio_vela_over_sdevice": vela_srh / sentaurus_integral,
                "frozen_sdevice_state_operator_ratio_vela_over_sdevice": vela_operator_on_sdevice / sentaurus_integral,
                "state_feedback_ratio_sdevice_qf_over_vela_state": vela_operator_on_sdevice / vela_srh,
                "integration_weight_audit": {
                    "total_area_um2": sum(source_areas),
                    "maximum_absolute_difference_um2": max(area_differences),
                    "maximum_relative_difference": max(
                        difference / max(abs(reference), 1.0e-300)
                        for difference, reference in zip(area_differences, independent_areas, strict=True)
                    ),
                },
            }
        )

    frozen_ratios = [float(point["frozen_sdevice_state_operator_ratio_vela_over_sdevice"]) for point in points]
    weight_errors = [float(point["integration_weight_audit"]["maximum_relative_difference"]) for point in points]
    report = {
        "schema_version": 1,
        "models_par": {"path": str(args.models_par.resolve()), "sha256": sha256(args.models_par)},
        "parameter_contract": {"sdevice": printed, "vela": expected, "checks": parameter_checks, "pass": all(parameter_checks.values())},
        "integration_contract": {
            "method": "compare production-exported node areas with independent triangle area/3 lumping",
            "maximum_relative_weight_difference": max(weight_errors),
            "pass": max(weight_errors) <= 1.0e-12,
        },
        "frozen_state_result": {
            "minimum_operator_ratio_vela_over_sdevice": min(frozen_ratios),
            "maximum_operator_ratio_vela_over_sdevice": max(frozen_ratios),
            "interpretation": "Residual gap after substituting SDevice psi/eQF/hQF into the Vela production SRH operator; includes carrier-statistics/BGN formula semantics but excludes Vela nonlinear-state error.",
        },
        "points": points,
    }
    report["verified"] = {
        "parameter_values_aligned": report["parameter_contract"]["pass"],
        "comparison_weights_match_production_lumping": report["integration_contract"]["pass"],
    }
    report["unresolved"] = [
        "SDevice internal element/source quadrature is not exposed by the plotted nodal SRH field.",
        "The frozen-state operator gap still combines generalized Fermi-SRH, effective-intrinsic-density/BGN, and proprietary SDevice source evaluation semantics.",
    ]
    (args.output_root / "srh_alignment_audit.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    chart_rows: list[dict[str, object]] = []
    for point in points:
        common = {
            "VCE_V": point["VCE_V"],
            "state_feedback_ratio": point["state_feedback_ratio_sdevice_qf_over_vela_state"],
            "parameter_pass": report["parameter_contract"]["pass"],
            "weight_pass": report["integration_contract"]["pass"],
        }
        chart_rows.extend(
            [
                {**common, "series": "Vela self-consistent", "ratio": point["self_consistent_integral_ratio_vela_over_sdevice"]},
                {**common, "series": "SDevice potential/QF frozen", "ratio": point["frozen_sdevice_state_operator_ratio_vela_over_sdevice"]},
            ]
        )
    with (args.output_root / "srh_alignment_points.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(chart_rows[0]))
        writer.writeheader()
        writer.writerows(chart_rows)
    print(json.dumps({
        "parameter_pass": report["parameter_contract"]["pass"],
        "integration_weight_pass": report["integration_contract"]["pass"],
        "frozen_ratio_min": min(frozen_ratios),
        "frozen_ratio_max": max(frozen_ratios),
    }))
    return 0 if all(parameter_checks.values()) and max(weight_errors) <= 1.0e-12 else 1


if __name__ == "__main__":
    raise SystemExit(main())
