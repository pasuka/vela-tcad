#!/usr/bin/env python3
"""Bind P1 generated grids to the real sd-sim check-only frontend (P2)."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from sentaurus_sde_cases import REPO, DEFAULT_MANIFEST, PRIMARY, digest, require, load_manifest, check_contract, compare_native
from jsonschema import Draft7Validator
from sentaurus_gmsh_pipeline import _strict_json
from vela_sdevice.cmd_parser import parse_cmd_text

CONFIG = REPO / "configs/reference_sde_checks.json"
SEMANTIC_KEYS = ("electrodes", "physics", "plot", "current_plot", "math", "solve", "stages")


def fingerprint(path: Path) -> dict:
    return {"path": path.relative_to(REPO).as_posix(), "sha256": digest(path),
            "sha256_lf": hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()}


def load_cases(path=CONFIG):
    data = json.loads(path.read_text())
    Draft7Validator(json.loads((REPO / "schemas/vela.reference_sde_checks.v1.schema.json").read_text())).validate(data)
    require(len({c["id"] for c in data["cases"]}) == len(data["cases"]), "duplicate check IDs")
    for case in data["cases"]:
        require(re.fullmatch(r"[a-z0-9_]+", case["id"]) is not None, "invalid check ID")
        for key in ("command", "parameter", "sidecar"):
            spec = case.get(key)
            if spec:
                p = (REPO / spec["path"]).resolve()
                require(p.is_relative_to(REPO.resolve()) and p.is_file(), "input path outside repository or missing")
                require(fingerprint(p)["sha256_lf"] == spec["sha256_lf"], f"{case['id']}: {key} input identity changed")
    return data


def effective_text(text: str) -> tuple[str, list[str]]:
    changed, n = re.subn(r'\bGrid\s*=\s*"[^"]+"', 'Grid="device.grd"', text)
    require(n == 1, "expected exactly one explicit Grid declaration")
    edits = ["Grid=device.grd"]
    if re.search(r"\bParameter\s*=", changed):
        changed, n = re.subn(r'\bParameter\s*=\s*"[^"]+"', 'Parameter="models.par"', changed)
        require(n == 1, "expected one explicit Parameter declaration")
        edits.append("Parameter=models.par (explicitly staged copy)")
    else:
        changed = changed.replace('Grid="device.grd"', 'Grid="device.grd" Parameter="models.par"', 1)
        edits.append("Parameter=models.par (added; original E301 boundary retained separately)")
    return changed, edits


def prepare_case(case: dict, meshes: Path, directory: Path) -> dict:
    require(not directory.exists(), "case directory must be fresh")
    source = REPO / case["command"]["path"]
    original = source.read_text(encoding="utf-8")
    effective, edits = effective_text(original)
    before = parse_cmd_text(original, path_name="device.cmd")
    after = parse_cmd_text(effective, path_name="device.cmd")
    for key in SEMANTIC_KEYS:
        require(before[key] == after[key], f"input staging altered {key}")
    directory.mkdir(parents=True)
    for ext in ("grd", "dat"):
        shutil.copy2(meshes / case["mesh_case"] / ("device." + ext), directory / ("device." + ext))
    shutil.copy2(REPO / case["parameter"]["path"], directory / "models.par")
    command = directory / "device.cmd"
    command.write_text(effective, encoding="utf-8", newline="\n")
    if case.get("sidecar"):
        shutil.copy2(REPO / case["sidecar"]["path"], command.with_suffix(".vela.json"))
    return {"source": fingerprint(source), "parameter": fingerprint(REPO / case["parameter"]["path"]),
            "parameter_policy": case["parameter"]["policy"], "edits": edits,
            "effective_command_sha256": digest(command),
            "sidecar_sha256": digest(command.with_suffix(".vela.json")) if case.get("sidecar") else None,
            "original_verdict": before["verdict"], "original_ir": before}


def invoke_check(launcher: Path, directory: Path) -> subprocess.CompletedProcess:
    # Both external tools are observable sentinels. --check must call neither.
    directory = directory.resolve()
    sentinel = directory / "forbidden_tool.py"
    sentinel.write_text("from pathlib import Path\nPath(__file__).with_suffix('.called').write_text('called')\nraise SystemExit(91)\n")
    env = dict(os.environ, VELA_PYTHON=sys.executable,
               VELA_EXAMPLE_RUNNER=str(sentinel), VELA_SENTAURUS_IMPORT=str(sentinel))
    result = subprocess.run([str(launcher.resolve()), "--check", str((directory / "device.cmd").resolve())],
                            cwd=directory, env=env, capture_output=True, text=True)
    require(not sentinel.with_suffix(".called").exists(), "check invoked solver or TDR importer")
    return result


def verify_check(case: dict, packet: dict, meshes: Path, directory: Path) -> dict:
    work = directory / "device_vela"
    run = json.loads((work / "run_manifest.json").read_text())
    require(run["exit_code"] == 0 and run["check_only"] and run["stages"] == [], "check invoked numerical stages")
    plan = json.loads((work / "stage_plan.json").read_text())
    ir = json.loads((work / "execution_ir.json").read_text())
    require(plan["cmd"]["sha256"] == packet["effective_command_sha256"] == ir["source"]["sha256"], "cmd provenance mismatch")
    for key in SEMANTIC_KEYS:
        require(packet["original_ir"][key] == ir[key], f"frontend changed command semantics: {key}")
    report = json.loads((work / "materials_report.json").read_text())
    require(report["parameter_file"]["sha256"] == packet["parameter"]["sha256"], "Parameter provenance mismatch")
    require(not report["blocking"] and not report["lossy"], "unexpected material approximation")
    require(report["default_library"]["sha256"] == digest(REPO/'share/vela/materials/models.par'),
            "runtime default material library differs from source")
    materials = json.loads((work / "materials.json").read_text())
    require(materials["unit_system"]["concentration"] == "cm^-3" and
            materials["unit_system"]["energy"] == "eV", "material units changed")
    qualification = json.loads((work / "mesh_qualification.json").read_text())
    require(qualification["source"]["format"] == "dfise" and
            qualification["source"]["sha256"] == digest(directory/'device.grd') and
            qualification["source"]["companion_sha256"] == digest(directory/'device.dat'), "grid provenance mismatch")
    require(qualification["coordinate_unit"]["applied_scale_to_um"] == 1.0 and
            qualification["coordinate_unit"]["declared"] == "um", "mesh units changed")
    require(set(qualification["materials"]) == set(case["expected_materials"]), "volume material identity mismatch")
    require(all(row["matched"] for row in qualification["qualification"]["electrodes_checked"]), "electrode gate failed")
    require({e["name"] for e in plan["electrodes"]} == set(case["contacts"]), "electrode identity changed")
    compare_native(meshes / case["mesh_case"], work)
    require([s["kind"] for s in plan["stages"]] == case["stage_kinds"], "stage sequence changed")
    source_goals = [s["goal"] for s in packet["original_ir"]["stages"] if s["kind"] == "sweep"]
    planned_sweeps = [s for s in plan["stages"] if s["kind"] == "sweep"]
    require(len(source_goals) == len(planned_sweeps), "sweep goals dropped")
    for goal, stage in zip(source_goals, planned_sweeps):
        require(stage["electrode_voltages_goal_V"][goal["contact"]] == goal["voltage_V"], "sweep voltage changed")
    decks = []
    voltages = {e["name"]: e["voltage_V"] for e in packet["original_ir"]["electrodes"]}
    goal_index = 0
    for i, stage in enumerate(plan["stages"]):
        require(stage["index"] == i and stage["input_state"] ==
                (plan["stages"][i-1]["output_state"] if i else None), "stage dependency changed")
        require(stage["electrode_voltages_start_V"] == voltages, "stage start voltages changed")
        goals = dict(voltages)
        if stage["kind"] == "sweep":
            goal = source_goals[goal_index]
            goal_index += 1
            goals[goal["contact"]] = goal["voltage_V"]
        require(stage["electrode_voltages_goal_V"] == goals, "stage end voltages changed")
        deck = json.loads((work / stage["deck"]).read_text())
        if stage["kind"] == "sweep":
            require(deck["sweep"]["contact"] == goal["contact"] and
                    deck["sweep"]["start"] == voltages[goal["contact"]] and
                    deck["sweep"]["stop"] == goal["voltage_V"], "emitted sweep goal changed")
        require(deck["format_version"] == 2 and deck["sweep"]["write_state_file"] == stage["output_state"], "state output mismatch")
        if i:
            require(deck["sweep"]["initial_state_file"] == stage["input_state"], "state input mismatch")
        require(all(deck[k] == v for k,v in {"mesh_file":"mesh.json","node_doping_file":"doping.csv","materials_file":"materials.json"}.items()), "stage inputs changed")
        byname = {c["name"]: c for c in deck["contacts"]}
        for name, expected in case["contacts"].items():
            require(byname[name]["bias"] == voltages[name], "contact start bias changed")
            for key, value in expected.items():
                require(byname[name].get(key) == value, f"contact boundary changed: {name}/{key}")
        if case.get("sidecar"):
            profile = deck["solver"]["impact_ionization"]
            require(profile["current_approximation"] == "element_edge_sg_gss_laux" and
                profile["source_mapping_mode"] == "element_vertex_box_measure" and
                profile["cell_reconstructed_midpoint_density"] == "bernoulli" and
                deck["mesh_geometry"]["node_volume_policy"] == "mixed_voronoi" and
                deck["mesh_geometry"]["require_non_obtuse"] and
                deck["solver"]["mobility"]["high_field_gradient_discretization"] == "transport_cell_vector",
                "atomic sidecar profile changed")
        voltages = goals
        decks.append({"index": i, "kind": stage["kind"], "sha256": digest(work/stage["deck"])})
    require(plan["inputs"].get("sidecar") == ("device.vela.json" if packet["sidecar_sha256"] else None), "sidecar identity mismatch")
    if packet["sidecar_sha256"]:
        require(digest(directory/'device.vela.json') == packet["sidecar_sha256"], "sidecar modified by frontend")
    require(not list(work.glob("stage*.csv")) and not list(work.glob("*.h5")), "check produced numerical outputs")
    return {"status":"checked", "stage_decks":decks, "grid_sha256":digest(directory/'device.grd'),
            "dat_sha256":digest(directory/'device.dat'), "materials_sha256":digest(work/'materials.json'),
            "node_count":qualification["nodes"]["count"], "cell_count":qualification["elements"]["count"],
            "runtime":run["version"], "numerical_qualification":"not_run"}


def run_suite(meshes: Path, output: Path, launcher: Path) -> dict:
    require(launcher.is_file(), "build sd-sim first or specify --launcher")
    require(not output.exists() or not any(output.iterdir()), "use a fresh empty output directory")
    require(not meshes.resolve().is_relative_to(output.resolve()), "output contains input grids")
    require(not REPO.resolve().is_relative_to(output.resolve()), "output contains repository")
    cases = load_cases()["cases"]
    original_hashes = {str(REPO/c[key]["path"]): digest(REPO/c[key]["path"])
        for c in cases for key in ("command", "parameter", "sidecar") if c.get(key)}
    p1 = load_manifest(DEFAULT_MANIFEST)
    summary = json.loads((meshes / 'case_summary.json').read_text())
    require(summary["status"] == "passed" and summary["manifest_sha256"] == digest(DEFAULT_MANIFEST), "P1 identity mismatch")
    require({r["id"] for r in summary["cases"] if r["status"] == "generated"} ==
            {c["mesh_case"] for c in cases}, "missing P1 generated case")
    for row in summary["cases"]:
        if row["status"] == "generated":
            for name in PRIMARY:
                require(digest(meshes/row["id"]/name) == row["runs"][-1]["sha256"][name], "P1 generated file modified")
            contract = next(c["contract"] for c in p1["cases"] if c["id"] == row["id"])
            check_contract(meshes/row["id"],contract)
    output.mkdir(parents=True,exist_ok=True)
    results=[]
    for case in cases:
        row={"id":case["id"], "mesh_case":case["mesh_case"]}
        try:
            directory=output/case["id"]
            packet=prepare_case(case,meshes,directory)
            row["inputs"]={k:v for k,v in packet.items() if k!='original_ir'}
            result=invoke_check(launcher,directory)
            require(result.returncode == 0, result.stderr.strip())
            row.update(verify_check(case,packet,meshes,directory))
        except Exception as error:
            row.update(status="failed",error=str(error))
        finally:
            for key in ("command","parameter","sidecar"):
                if case.get(key):
                    p = REPO/case[key]["path"]
                    require(digest(p) == original_hashes[str(p)], "repository input bytes changed")
        results.append(row)
    package = next((p for p in (launcher.parent/'bin/vela_sdevice.pyz', launcher.parent/'vela_sdevice.pyz') if p.is_file()), None)
    require(package is not None, "cannot fingerprint launcher zipapp")
    report={"schema":"vela.reference_sde_checks_run.v1", "status":"passed" if all(r['status']=='checked' for r in results) else "failed",
            "check_contract_sha256":digest(CONFIG), "p1_manifest_sha256":digest(DEFAULT_MANIFEST),
            "launcher_sha256":digest(launcher), "package_sha256":digest(package),
            "p1_summary_sha256":digest(meshes/'case_summary.json'), "numerical_qualification":"not_run", "cases":results}
    _strict_json(output/'check_summary.json',report)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--meshes-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--launcher',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        result=run_suite(args.meshes_dir,args.output_dir,args.launcher)
    except Exception as error:
        print(str(error),file=sys.stderr);return 2
    print(json.dumps({'status':result['status'],'checked':sum(r['status']=='checked' for r in result['cases']),
                      'numerical_qualification':'not_run'},indent=2))
    return 0 if result['status']=='passed' else 2


if __name__=='__main__':
    raise SystemExit(main())
