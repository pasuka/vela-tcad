#!/usr/bin/env python3
"""Prepare, run, and qualify the SimpleMOS M4 controlled mobility matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m3_workflow as m3  # noqa: E402
from sentaurus_import import parse_quoted_list, parse_values_block  # noqa: E402


DEFAULT_CONTRACT = (
    REPO / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m4_controlled_mobility_contract_v1.json"
)
DEFAULT_OUTPUT = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m4_controlled_mobility"
)
DEFAULT_REFERENCE_DIR = (
    REPO / "reference_tcad" / "simplemos_sentaurus2022"
    / "controlled_mobility"
)
DEFAULT_NEUTRAL_DIR = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m2_gate"
)
DEFAULT_MATERIALS = (
    REPO / "reference_tcad" / "transportmodels_sentaurus2022" / "vela"
    / "materials_sentaurus2022.json"
)
DEFAULT_REMOTE_TDR = (
    "~/sentaurus_runs/vela_oracle/simplemos_m2_nominal_20260826/source/"
    "n17_fps.tdr"
)
DEFAULT_REMOTE_ROOT = (
    "~/sentaurus_runs/vela_oracle/simplemos_m4_controlled_mobility_20260826"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = (Path(os.environ.get("SystemRoot", r"C:\Windows"))
                     / "System32" / "OpenSSH" / f"{name}.exe")
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
    )
    return completed.stdout or ""


def voltage_tag(value: float) -> str:
    return format(value, ".12g").replace("-", "m").replace(".", "p")


def validate_contract(contract: dict[str, Any]) -> None:
    variants = contract["variants"]
    if [item["id"] for item in variants] != ["A0", "A1", "A2", "A3"]:
        raise ValueError("M4 variants must be ordered A0, A1, A2, A3")
    expected_prefixes = [[], ["DopingDependence"],
                         ["DopingDependence", "HighFieldSaturation"],
                         ["DopingDependence", "HighFieldSaturation", "Enormal"]]
    for item, expected in zip(variants, expected_prefixes):
        if item["sentaurus_mobility_models"] != expected:
            raise ValueError(f"{item['id']}: controlled model ladder mismatch")
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    if len(lattice) != 51 or len(set(lattice)) != len(lattice):
        raise ValueError("M4 gate lattice must contain 51 unique points")
    for index, value in enumerate(lattice):
        if not math.isclose(value, index * 0.05, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError("M4 gate lattice must be exactly 0:0.05:2.5 V")
    protocol = contract["vela_execution_protocol"]
    if not protocol["gate_acceptance_lattice_unchanged"]:
        raise ValueError("M4 Vela protocol must preserve the gate lattice")
    drain = protocol["auxiliary_drain_ramp"]
    if drain["acceptance_role"] != "none":
        raise ValueError("auxiliary drain points cannot be acceptance points")


def auxiliary_drain_schedule(contract: dict[str, Any], variant: str,
                             drain_voltage_V: float) -> dict[str, Any]:
    base = dict(contract["vela_execution_protocol"]["auxiliary_drain_ramp"])
    overrides = base.pop("overrides", [])
    for item in overrides:
        if (item["variant"] == variant and math.isclose(
                float(item["drain_voltage_V"]), drain_voltage_V,
                rel_tol=0.0, abs_tol=1e-12)):
            base.update({key: value for key, value in item.items()
                         if key not in {"variant", "drain_voltage_V"}})
    return base


def sentaurus_deck(case: str, vd: float, models: list[str], intervals: int) -> str:
    mobility = ""
    if models:
        mobility = "   Mobility( " + " ".join(models) + " )\n"
    return f'''File {{
   Grid="n17_fps.tdr"
   Plot="{case}_des.tdr"
   Current="{case}"
   Output="{case}.log"
}}

Electrode {{
   {{ Name="source" Voltage=0.0 }}
   {{ Name="drain" Voltage=0.0 }}
   {{ Name="gate" Voltage=0.0 }}
   {{ Name="substrate" Voltage=0.0 }}
}}

Physics {{ EffectiveIntrinsicDensity(OldSlotboom NoFermi) }}
Physics(Material="Silicon") {{
{mobility}   Recombination(SRH(DopingDependence))
}}

Math {{ Extrapolate Iterations=20 ExitOnFailure }}

Solve {{
   Coupled(Iterations=100) {{ Poisson }}
   Coupled {{ Poisson Electron Hole }}
   Quasistationary(
      InitialStep=0.1 Increment=1.5 MinStep=1e-5 MaxStep=1
      Goal {{ Name="drain" Voltage={vd:.17g} }}
   ) {{ Coupled {{ Poisson Electron Hole }} }}
   NewCurrentPrefix="IdVg_"
   Quasistationary(
      DoZero InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05
      Goal {{ Name="gate" Voltage=2.5 }}
   ) {{
      Coupled {{ Poisson Electron Hole }}
      CurrentPlot(Time=(Range=(0 1) Intervals={intervals}))
   }}
}}
'''


def prepare_sentaurus_bundle(contract: dict[str, Any], contract_path: Path,
                              output_dir: Path) -> dict[str, Any]:
    bundle = output_dir / "sentaurus_bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    intervals = contract["bias_matrix"]["gate_lattice"]["point_count"] - 1
    cases: list[dict[str, Any]] = []
    for variant in contract["variants"]:
        for vd in contract["bias_matrix"]["drain_voltages_V"]:
            case = f"{str(variant['id']).lower()}_vd_{voltage_tag(float(vd))}"
            case_dir = bundle / case
            case_dir.mkdir(parents=True, exist_ok=True)
            deck_path = case_dir / f"{case}_des.cmd"
            deck_path.write_text(sentaurus_deck(
                case, float(vd), list(variant["sentaurus_mobility_models"]),
                intervals), encoding="utf-8", newline="\n")
            cases.append({
                "case": case,
                "variant": variant["id"],
                "drain_voltage_V": float(vd),
                "deck": str(deck_path.resolve()),
                "deck_sha256": sha256(deck_path),
                # SDevice appends the input deck's ``_des`` stem to a
                # NewCurrentPrefix-generated current file.
                "expected_plot": f"IdVg_{case}_des.plt",
            })
    manifest = {
        "schema": "vela.simplemos.sdevice.m4_sentaurus_matrix.v1",
        "status": "prepared",
        "contract_sha256": sha256(contract_path),
        "cases": cases,
    }
    write_json(output_dir / "sentaurus_matrix_manifest.json", manifest)
    return manifest


def run_sentaurus_vm(manifest: dict[str, Any], output_dir: Path,
                     ssh_target: str, ssh_bin: str, scp_bin: str,
                     remote_root: str, remote_tdr: str) -> str:
    banner = run(
        [ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
        capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, ssh_target, f"mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(output_dir / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])
    for index, case in enumerate(manifest["cases"], start=1):
        name = case["case"]
        print(f"[{index}/{len(manifest['cases'])}] running {name}", flush=True)
        command = (
            f"set -eu; cd {remote_root}/sentaurus_bundle/{name}; "
            f"cp {remote_tdr} n17_fps.tdr; "
            f"sdevice {name}_des.cmd > {name}.console.log 2>&1"
        )
        run([ssh_bin, ssh_target, command])
    archive_name = "simplemos_m4_results.tgz"
    run([ssh_bin, ssh_target,
         f"cd {remote_root} && tar -czf {archive_name} sentaurus_bundle"])
    raw = output_dir / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / archive_name
    run([scp_bin, f"{ssh_target}:{remote_root}/{archive_name}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (output_dir / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n")
    return banner


def exact_plt_curve(plot: Path, lattice: list[float], tolerance: float
                    ) -> list[dict[str, float]]:
    text = plot.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    gate_index = datasets.index("gate OuterVoltage")
    drain_index = datasets.index("drain TotalCurrent")
    source = parse_values_block(text, len(datasets))
    result: list[dict[str, float]] = []
    for target in lattice:
        matches = [row for row in source
                   if abs(float(row[gate_index]) - target) <= tolerance]
        if len(matches) != 1:
            raise ValueError(
                f"{plot}: expected exactly one row at Vg={target:g} V, "
                f"observed {len(matches)}")
        result.append({
            "gate_voltage_V": target,
            "drain_total_current_A_per_um": float(matches[0][drain_index]),
        })
    if not all(math.isfinite(row["drain_total_current_A_per_um"])
               for row in result):
        raise ValueError(f"{plot}: non-finite drain current")
    return result


def extract_references(contract: dict[str, Any], manifest: dict[str, Any],
                       contract_path: Path, output_dir: Path, reference_dir: Path,
                       banner: str) -> dict[str, Any]:
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    tolerance = float(contract["comparison"]["exact_bias_tolerance_V"])
    raw_bundle = output_dir / "sentaurus_raw" / "sentaurus_bundle"
    reference_dir.mkdir(parents=True, exist_ok=True)
    artifacts: list[dict[str, Any]] = []
    for case in manifest["cases"]:
        plot = raw_bundle / case["case"] / case["expected_plot"]
        rows = exact_plt_curve(plot, lattice, tolerance)
        output = reference_dir / f"{case['case']}_reference.csv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        artifacts.append({
            "case": case["case"],
            "variant": case["variant"],
            "drain_voltage_V": case["drain_voltage_V"],
            "path": output.name,
            "sha256": sha256(output),
            "point_count": len(rows),
        })
    reference_manifest = {
        "schema": "vela.simplemos.sdevice.m4_reference_curves.v1",
        "status": "qualified",
        "sentaurus_banner": banner,
        "source_tdr_sha256": contract["upstream_tdr_sha256"],
        "contract": Path(os.path.relpath(
            contract_path, reference_dir)).as_posix(),
        "contract_sha256": sha256(contract_path),
        "interpolation": "forbidden",
        "artifacts": artifacts,
        "raw_binary_commit_policy": "excluded",
    }
    write_json(reference_dir / "reference_manifest.json", reference_manifest)
    return reference_manifest


def neutral_mesh(neutral_dir: Path) -> dict[str, Any]:
    with (neutral_dir / "nodes.csv").open(newline="", encoding="utf-8-sig") as handle:
        nodes = [{"id": int(row["id"]), "x": float(row["x_um"]),
                  "y": float(row["y_um"])} for row in csv.DictReader(handle)]
    with (neutral_dir / "elements.csv").open(newline="", encoding="utf-8-sig") as handle:
        element_rows = list(csv.DictReader(handle))
    region_names = list(dict.fromkeys(row["region"] for row in element_rows))
    region_id = {name: index for index, name in enumerate(region_names)}
    material_map = {"Si": "Si", "SiO2": "SiO2", "Nitride": "Nitride"}
    triangles = [{
        "id": int(row["id"]),
        "region_id": region_id[row["region"]],
        "node_ids": [int(row["node0"]), int(row["node1"]), int(row["node2"])],
    } for row in element_rows]
    regions = []
    for name in region_names:
        rows = [row for row in element_rows if row["region"] == name]
        material = rows[0]["material"]
        if material not in material_map:
            raise ValueError(f"unsupported neutral material {material!r}")
        regions.append({
            "id": region_id[name], "name": name,
            "material": material_map[material],
            "cell_ids": [int(row["id"]) for row in rows],
        })
    with (neutral_dir / "contacts.csv").open(newline="", encoding="utf-8-sig") as handle:
        contact_rows = list(csv.DictReader(handle))
    contacts = [{
        "id": index,
        "name": row["name"],
        "region_id": region_id[row["region"]],
        "node_ids": [int(value) for value in row["node_ids"].split(";") if value],
    } for index, row in enumerate(contact_rows)]
    return {"nodes": nodes, "triangles": triangles,
            "regions": regions, "contacts": contacts}


def base_config(mesh: Path, doping: Path, materials: Path,
                mobility: dict[str, Any]) -> dict[str, Any]:
    return {
        "simulation_type": "dc_sweep",
        "mesh_file": str(mesh.resolve()),
        "node_doping_file": str(doping.resolve()),
        "materials_file": str(materials.resolve()),
        "output_csv": "unused.csv",
        "scaling": {"mode": "unit_scaling"},
        "doping": [],
        "contacts": [{"name": name, "bias": 0.0, "type": "ohmic"}
                     for name in ("source", "drain", "gate", "substrate")],
        "solver": {
            "method": "newton", "max_iter": 100, "reltol": 1e-7,
            "abstol": 1e-12, "damping_psi": 0.35, "warm_start": True,
            "quasi_fermi_update_limit_V": 0.025,
            "srh_density_coupling": "sentaurus_default",
            "quasi_fermi_reference": "contact_basin",
            "continuity_row_scaling": {
                "enabled": True, "flux_fraction": 1e-3,
                "scale_floor": 1e-30, "min_source_scale": 1e-18,
                "min_weight": 1e-12, "max_weight": 1e12,
            },
            "mobility": mobility,
            "recombination": ["srh"],
            "srh_doping_dependence": {
                "enabled": True, "concentration_basis": "total_impurity",
                "temperature_dependence": False,
                "reference_temperature_K": 300.0,
                "electron_temperature_exponent": -1.5,
                "hole_temperature_exponent": -1.5,
                "electron": {"tau_min_s": 0.0, "tau_max_s": 3e-8,
                             "reference_doping_m3": 1e16, "gamma": 1.0},
                "hole": {"tau_min_s": 0.0, "tau_max_s": 3e-6,
                         "reference_doping_m3": 1e16, "gamma": 1.0},
            },
            "bandgap_narrowing": {
                "model": "old_slotboom", "fermi_statistics_correction": False},
            "impact_ionization": {"model": "none"},
        },
        "sweep": {"mode": "iv"},
    }


def prepare_vela_workflows(contract: dict[str, Any], neutral_dir: Path,
                           materials: Path, output_dir: Path,
                           contract_path: Path = DEFAULT_CONTRACT
                           ) -> dict[str, Any]:
    vela_dir = output_dir / "vela"
    vela_dir.mkdir(parents=True, exist_ok=True)
    mesh_path = vela_dir / "mesh.json"
    write_json(mesh_path, neutral_mesh(neutral_dir))
    doping_path = neutral_dir / "doping.csv"
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    stages: list[dict[str, Any]] = []
    for variant in contract["variants"]:
        mobility: dict[str, Any] = {
            "model": variant["vela_mobility_model"],
            "doping_concentration_basis": "total_impurity",
        }
        for source, target in (
            ("vela_high_field_driving_force", "high_field_driving_force"),
            ("vela_high_field_gradient_discretization",
             "high_field_gradient_discretization"),
        ):
            if source in variant:
                mobility[target] = variant[source]
        base_path = vela_dir / f"{variant['id'].lower()}_base.json"
        write_json(base_path, base_config(
            mesh_path, doping_path, materials, mobility))
        workflow_dir = vela_dir / str(variant["id"]).lower()
        workflow = m3.materialize(
            base_path, REPO / "reference_tcad" / "simplemos_sentaurus2022"
            / "simplemos_sdevice_validation_contract_v1.json",
            workflow_dir,
            [float(value) for value in contract["bias_matrix"]["drain_voltages_V"]],
        )
        for stage in workflow["stages"]:
            config_path = Path(stage["config"])
            config = read_json(config_path)
            if stage["phase"] == "equilibrium":
                config["solver"]["quasi_fermi_update_limit_V"] = 0.1
            if stage["phase"] == "drain_ramp":
                config["solver"].update({
                    "line_search_mode": "block_filter",
                    "residual_filter_gamma": 1e-4,
                    "residual_filter_envelope_factor": 2.0,
                })
            if stage["phase"] == "drain_ramp":
                drain_voltage = float(
                    stage["vela_physical_step_control"]["voltage_span_V"])
                schedule = auxiliary_drain_schedule(
                    contract, str(variant["id"]), drain_voltage)
                config["sweep"]["step"] = math.copysign(
                    min(abs(float(config["sweep"]["step"])),
                        float(schedule["initial_step_V"])),
                    float(config["sweep"]["stop"])
                    - float(config["sweep"]["start"]),
                )
                config["sweep"]["max_step"] = min(
                    float(config["sweep"]["max_step"]),
                    float(schedule["max_step_V"]))
                config["simplemos_m4"] = {
                    "purpose": "continuation_only",
                    "does_not_change_acceptance_lattice": True,
                    "auxiliary_drain_ramp": schedule,
                }
            if stage["phase"] != "gate_sweep":
                write_json(config_path, config)
                stage["config_sha256"] = sha256(config_path)
                continue
            config["sweep"]["bias_points"] = lattice
            config["sweep"]["step"] = 0.05
            config["sweep"]["min_step"] = min(
                float(config["sweep"]["min_step"]), 0.05)
            config["sweep"]["max_step"] = max(
                float(config["sweep"]["max_step"]), 0.05)
            write_json(config_path, config)
            stage["config_sha256"] = sha256(config_path)
            stage["exact_gate_lattice_V"] = lattice
        write_json(workflow_dir / "workflow_manifest.json", workflow)
        stages.extend(workflow["stages"])
    manifest = {
        "schema": "vela.simplemos.sdevice.m4_vela_matrix.v1",
        "status": "materialized",
        "mesh": str(mesh_path.resolve()),
        "mesh_sha256": sha256(mesh_path),
        "doping": str(doping_path.resolve()),
        "doping_sha256": sha256(doping_path),
        "materials": str(materials.resolve()),
        "materials_sha256": sha256(materials),
        "contract_sha256": sha256(contract_path),
        "vela_execution_protocol": contract["vela_execution_protocol"],
        "stages": stages,
    }
    write_json(vela_dir / "vela_matrix_manifest.json", manifest)
    return manifest


def execute_vela_workflows(contract: dict[str, Any], output_dir: Path,
                           runner: Path,
                           selected_variants: set[str] | None = None
                           ) -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    for variant in contract["variants"]:
        variant_id = str(variant["id"])
        if selected_variants is not None and variant_id not in selected_variants:
            continue
        workflow_dir = output_dir / "vela" / variant_id.lower()
        manifest_path = workflow_dir / "workflow_manifest.json"
        workflow = read_json(manifest_path)
        workflow = m3.execute(workflow, runner.resolve(), workflow_dir)
        results.append({
            "variant": variant_id,
            "status": workflow["status"],
            "manifest": str(manifest_path.resolve()),
            "manifest_sha256": sha256(manifest_path),
        })
    report = {
        "schema": "vela.simplemos.sdevice.m4_vela_execution.v1",
        "status": ("accepted" if all(item["status"] == "accepted"
                                     for item in results) else "fail"),
        "runner": str(runner.resolve()),
        "runner_sha256": sha256(runner),
        "variants": results,
    }
    write_json(output_dir / "vela" / "vela_execution_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    parser.add_argument("--neutral-dir", type=Path, default=DEFAULT_NEUTRAL_DIR)
    parser.add_argument("--materials", type=Path, default=DEFAULT_MATERIALS)
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--extract-existing", action="store_true")
    parser.add_argument("--execute-vela", action="store_true")
    parser.add_argument(
        "--execute-vela-variant", action="append", default=[],
        choices=["A0", "A1", "A2", "A3"],
        help="execute only the selected Vela variant; may be repeated")
    parser.add_argument("--runner", type=Path,
                        default=REPO / "build-release" / "vela_example_runner.exe")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--remote-tdr", default=DEFAULT_REMOTE_TDR)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if args.clean and output_dir.exists():
        shutil.rmtree(output_dir)
    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    manifest = prepare_sentaurus_bundle(contract, contract_path, output_dir)
    prepare_vela_workflows(
        contract, args.neutral_dir.resolve(), args.materials.resolve(), output_dir,
        contract_path)
    vela_report = None
    if args.execute_vela:
        if not args.runner.is_file():
            raise FileNotFoundError(args.runner)
        vela_report = execute_vela_workflows(
            contract, output_dir, args.runner.resolve(),
            set(args.execute_vela_variant) or None)
    banner = ""
    if args.live_sentaurus:
        banner = run_sentaurus_vm(
            manifest, output_dir, args.ssh_target, args.ssh_bin, args.scp_bin,
            args.remote_root, args.remote_tdr)
    elif args.extract_existing:
        banner_path = output_dir / "sentaurus_banner.txt"
        if not banner_path.is_file():
            raise FileNotFoundError(banner_path)
        banner = banner_path.read_text(encoding="utf-8").strip()
    if args.live_sentaurus or args.extract_existing:
        references = extract_references(
            contract, manifest, contract_path, output_dir,
            args.reference_dir.resolve(), banner)
        print(json.dumps({"status": references["status"],
                          "reference_curves": len(references["artifacts"]),
                          "vela_status": (vela_report or {}).get("status")}))
    else:
        print(json.dumps({"status": "prepared", "cases": len(manifest["cases"]),
                          "vela_status": (vela_report or {}).get("status")}))
    return 1 if vela_report is not None and vela_report["status"] == "fail" else 0


if __name__ == "__main__":
    raise SystemExit(main())
