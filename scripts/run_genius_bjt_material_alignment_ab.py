#!/usr/bin/env python3
"""Quantify Nc/Nv and OldSlotboom Fermi-correction effects for the BJT case."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
SDEVICE_NC_CM3 = 2.8567e19
SDEVICE_NV_CM3 = 3.1046e19


def absolute(path: Path) -> str:
    return str(path.resolve())


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def run(command: list[str], cwd: Path, allow_gate_failure: bool = False) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    if result.returncode != 0 and not allow_gate_failure:
        raise RuntimeError(result.stderr or result.stdout)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=REPO / "build-release" / "vela_example_runner.exe")
    parser.add_argument("--output-root", type=Path, default=BUILD_ROOT / "m1_material_alignment_ab")
    parser.add_argument(
        "--sentaurus-fields-root",
        type=Path,
        default=BUILD_ROOT / "m1_current_diagnosis" / "sentaurus_vce3",
    )
    args = parser.parse_args()
    output = args.output_root.resolve()
    base_config = json.loads(
        (FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json").read_text(encoding="utf-8")
    )
    base_material = json.loads(
        (FIXTURE / "vela" / "materials_sentaurus2022.json").read_text(encoding="utf-8")
    )
    base_config.pop("sweep", None)
    base_config.update(
        {
            "simulation_type": "newton_solve_from_state",
            "mesh_file": absolute(FIXTURE / "vela" / "input" / "mesh.json"),
            "node_doping_file": absolute(FIXTURE / "vela" / "input" / "doping.csv"),
            "state_file": absolute(BUILD_ROOT / "m1_accepted_states" / "states" / "vce_030.csv"),
        }
    )
    variants = (
        ("legacy_nc_nv", 2.8e19, 1.04e19, True),
        ("nc_nv_aligned", SDEVICE_NC_CM3, SDEVICE_NV_CM3, True),
        (
            "nc_nv_aligned_bgn_fermi_off",
            SDEVICE_NC_CM3,
            SDEVICE_NV_CM3,
            False,
        ),
    )
    records = []
    for name, nc_cm3, nv_cm3, bgn_fermi in variants:
        material = json.loads(json.dumps(base_material))
        material["materials"][0]["Nc_m3"] = nc_cm3
        material["materials"][0]["Nv_m3"] = nv_cm3
        material_path = output / name / "materials.json"
        write_json(material_path, material)

        config = json.loads(json.dumps(base_config))
        config["materials_file"] = absolute(material_path)
        config["solver"]["bandgap_narrowing"]["fermi_statistics_correction"] = bgn_fermi
        config["output_state_file"] = absolute(output / name / "state.csv")
        config["output_vtk"] = absolute(output / name / "state.vtk")
        config_path = output / name / "config.json"
        write_json(config_path, config)
        solve = run([absolute(args.runner), "--config", absolute(config_path)], REPO)
        status = json.loads(solve.stdout.splitlines()[-1])
        (output / name / "runner.stdout.log").write_text(solve.stdout, encoding="utf-8")
        (output / name / "runner.stderr.log").write_text(solve.stderr, encoding="utf-8")

        compare_dir = output / name / "comparison"
        compare = run(
            [
                sys.executable,
                absolute(REPO / "scripts" / "compare_genius_bjt_transport_fields.py"),
                "--reference-root",
                absolute(FIXTURE),
                "--sentaurus-fields-root",
                absolute(args.sentaurus_fields_root),
                "--vela-vtk",
                absolute(output / name / "state.vtk"),
                "--output-dir",
                absolute(compare_dir),
            ],
            REPO,
            allow_gate_failure=True,
        )
        if not (compare_dir / "transport_source_comparison.json").is_file():
            raise RuntimeError(compare.stderr or compare.stdout)
        comparison = json.loads(
            (compare_dir / "transport_source_comparison.json").read_text(encoding="utf-8")
        )
        records.append(
            {
                "variant": name,
                "Nc_cm3": nc_cm3,
                "Nv_cm3": nv_cm3,
                "old_slotboom_fermi_statistics_correction": bgn_fermi,
                "solver": {
                    "converged": status["converged"],
                    "iterations": status["iterations"],
                    "contact_currents_A_per_um": status["contact_currents_A_per_um"],
                },
                "transport_source_overall_pass": comparison["overall_pass"],
                "electron_current_density": comparison["current_density"]["electron"]["gate"]["observed"],
                "hole_current_density": comparison["current_density"]["hole"]["gate"]["observed"],
                "srh": comparison["recombination"]["srh"]["gate"]["observed"],
                "auger": comparison["recombination"]["auger"]["gate"]["observed"],
            }
        )
    baseline_comparison = json.loads(
        (FIXTURE / "comparison" / "transport_source_comparison.json").read_text(
            encoding="utf-8"
        )
    )
    accepted_currents = json.loads(
        (BUILD_ROOT / "m1_accepted_states" / "configs" / "vce_030_solve.json").read_text(
            encoding="utf-8"
        )
    )
    baseline_terminal_rows = (BUILD_ROOT / "m1_accepted_states" / "terminal_currents.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    terminal_header = baseline_terminal_rows[0].split(",")
    terminal_values = baseline_terminal_rows[-1].split(",")
    baseline_terminals = dict(zip(terminal_header, terminal_values, strict=True))
    baseline = {
        "variant": "current_contract",
        "Nc_cm3": base_material["materials"][0]["Nc_m3"],
        "Nv_cm3": base_material["materials"][0]["Nv_m3"],
        "old_slotboom_fermi_statistics_correction": accepted_currents["solver"][
            "bandgap_narrowing"
        ]["fermi_statistics_correction"],
        "contact_currents_A_per_um": {
            name: float(baseline_terminals[f"{name}_A_per_um"])
            for name in ("collector", "base", "emitter")
        },
        "transport_source_overall_pass": baseline_comparison["overall_pass"],
        "electron_current_density": baseline_comparison["current_density"]["electron"]["gate"]["observed"],
        "hole_current_density": baseline_comparison["current_density"]["hole"]["gate"]["observed"],
        "srh": baseline_comparison["recombination"]["srh"]["gate"]["observed"],
        "auger": baseline_comparison["recombination"]["auger"]["gate"]["observed"],
    }
    report = {
        "schema_version": 1,
        "bias": {"VBE_V": 0.7, "VCE_V": 3.0},
        "purpose": "separate Nc/Nv alignment from the OldSlotboom Fermi-statistics correction",
        "sentaurus_parameter_source": "T-2022.03-SP2 Silicon parameter dump retained by the Genius BJT investigation",
        "baseline": baseline,
        "variants": records,
    }
    write_json(output / "summary.json", report)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
