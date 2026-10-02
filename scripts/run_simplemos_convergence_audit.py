"""Frozen read-only audit of M65/M77 gate endpoints and stalled Vg=0 states."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

from build_simplemos_convergence_audit_runner import REPO, LOCAL, RUNNER

ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUT = ROOT / "convergence_acceptance_audit"
CONTRACT = ROOT / "simplemos_convergence_audit_contract_v1.json"
FREEZE = ROOT / "simplemos_convergence_audit_freeze_v1.json"
PROFILE = {"mode": "report", "eps_row": 1e-3, "scale_floor": 1e-30,
           "min_source_scale": 0., "min_source_scale_fraction": 1e-3,
           "min_source_global_fraction": 1e-6, "min_carrier_density_m3": 0.,
           "min_flux_scale_fraction": 1e-8, "min_flux_scale": 0.}
GLOBAL = {"mode": "report", "tolerance": 1e-3, "source_floor": 1e-12}


def read(path):
    return json.loads(path.read_text())


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path):
    return Path(path).relative_to(REPO).as_posix()


def write(path, value):
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def write_csv(path, values):
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(values[0]))
        writer.writeheader(); writer.writerows(values)


def vtk_state(path):
    text = path.read_text()
    point_data = text.split("POINT_DATA ", 1)[1]
    count = int(point_data.splitlines()[0])
    result = {}
    for name in ("Potential", "ElectronQuasiFermi", "HoleQuasiFermi", "Electrons", "Holes"):
        block = point_data.split("SCALARS " + name + " ", 1)[1].split("LOOKUP_TABLE default", 1)[1]
        data = [float(v) for v in block.split()[:count]]
        if len(data) != count or not all(math.isfinite(v) for v in data):
            raise ValueError("Invalid VTK state")
        result[name] = data
    # VTK uses the internal unit system (cm^-3 density for these decks), while
    # the restart CSV explicitly uses SI density.
    return [{"node_id": i, "psi": result["Potential"][i], "phin": result["ElectronQuasiFermi"][i],
             "phip": result["HoleQuasiFermi"][i], "electrons_m3": result["Electrons"][i]*1e6,
             "holes_m3": result["Holes"][i]*1e6} for i in range(count)]


def deck_for(source, state, vd, vg):
    deck = copy.deepcopy(read(source))
    for key in ("sweep", "output_csv", "log_file"):
        deck.pop(key, None)
    for c in deck["contacts"]:
        c["bias"] = vd if c["name"] == "drain" else vg if c["name"] == "gate" else 0.
    deck["state_file"] = str(state.resolve())
    return deck


def prepare():
    if CONTRACT.exists():
        raise FileExistsError(CONTRACT)
    cases = []
    inputs = {RUNNER, REPO / "build-release/vela_example_runner.exe", REPO / "build-release/libvela_core.a",
              REPO / "src/solver/NewtonSolver.cpp", REPO / "include/vela/solver/NewtonSolver.h",
              Path(__file__), REPO / "scripts/build_simplemos_convergence_audit_runner.py",
              LOCAL / "NewtonSolver.cpp", LOCAL / "runner.cpp", LOCAL / "include/vela/solver/NewtonSolver.h"}
    for label, folder in (("m65", "m65_ni/vela"), ("m77", "m77_combined_poisson_charge_volume_post_main")):
        for path in sorted((REPO / "build-release" / folder).glob("*/vd_*/20_gate_sweep/config.json")):
            config = read(path)
            curve = rows(path.parent / "curve.csv")
            selected = [r for r in curve if r["newton_convergence_reason"] == "stall_residual_floor"]
            if curve[-1] not in selected:
                selected.append(curve[-1])
            for r in selected:
                vg = float(r["bias_V"]); vd = next(c["bias"] for c in config["contacts"] if c["name"] == "drain")
                device = path.parts[-4]
                stage = "endpoint" if r is curve[-1] else "stall"
                case = f"{label}_{device}_{path.parts[-3]}_{stage}"
                target = LOCAL / "cases" / case
                if stage == "endpoint":
                    state = path.parent / "state.csv"
                else:
                    matches = [p for p in (path.parent / "vtk").glob("*.vtk")
                               if abs(float(p.stem.rsplit("_", 1)[1][:-1])-vg)<1e-10]
                    if len(matches) != 1:
                        raise ValueError("Ambiguous VTK")
                    state = target / "state.csv"
                    write_csv(state, vtk_state(matches[0])); inputs.add(matches[0])
                inputs.update((path, state, path.parent / "curve.csv"))
                inputs.update(Path(config[key]) for key in ("mesh_file", "node_doping_file", "materials_file"))
                cases.append({"case": case, "branch": label, "device": device, "vd": vd, "vg": vg,
                              "stage": stage, "source_config": rel(path), "state": rel(state),
                              "source_reason": r["newton_convergence_reason"],
                              "reference_current_A_per_um": float(r["current_total_A_per_um"])})
    write(CONTRACT, {"status": "frozen_before_read_only_audit", "cases": cases,
        "local_profile": PROFILE, "global_profile": GLOBAL,
        "rationale": "Audit both source-significant and flux-significant free carrier rows. Density screening disabled because diagnostic density naming does not establish SI conversion in unit_scaling decks. Profiles are diagnostics, not retroactive acceptance claims.",
        "current_identity": {"relative": 1e-8, "absolute_A_per_um": 1e-22},
        "counts": {"read_only_states": len(cases), "nonlinear_solves": 0, "sentaurus_runs": 0},
        "next": "Review residual/source magnitudes and qualification coverage before freezing four-workpoint strict convergence contrast; retain failed controls and historical inputs."})
    inputs.add(CONTRACT)
    write(FREEZE, {"input_hashes": {rel(p): sha(p) for p in sorted(inputs)}})
    print("Frozen read-only states:", len(cases))


def verify(path=FREEZE):
    for name, digest in read(path)["input_hashes"].items():
        if sha(REPO / name) != digest:
            raise ValueError("Frozen file changed: " + name)


def execute(deck, config, runner=RUNNER):
    write(config, deck)
    p = subprocess.run([str(runner), "--config", str(config), "--log", "off"], capture_output=True, text=True)
    config.with_suffix(".stdout.txt").write_text(p.stdout)
    config.with_suffix(".stderr.txt").write_text(p.stderr)
    status = json.loads(p.stdout.strip().splitlines()[-1]) if p.stdout.strip() else {}
    status["exit_code"] = p.returncode
    write(config.with_suffix(".status.json"), status)
    return status


def audit_case(case):
    target = LOCAL / "cases" / case["case"]
    source, state = REPO/case["source_config"], REPO/case["state"]
    deck = deck_for(source, state, case["vd"], case["vg"])
    deck["solver"]["carrier_row_convergence"] = PROFILE
    deck["solver"]["global_continuity_closure"] = GLOBAL
    term = execute(dict(deck, simulation_type="newton_carrier_term_probe", output_csv=str(target/"terms.csv"),
                        carrier_term_probe={"solved_equation_terms": True}), target/"audit.json")
    if term["exit_code"] or not term["read_only"]:
        raise ValueError("Read-only audit failed")
    contacts = {}
    for contact in ("source", "drain", "substrate", "gate"):
        status = execute(dict(deck, simulation_type="terminal_current_functional_probe", contact=contact), target/(contact+".json"))
        if status["exit_code"]:
            raise ValueError("Current audit failed")
        contacts[contact] = status
    current = contacts["drain"]["current_A_per_um"]
    identity = abs(current-case["reference_current_A_per_um"])
    local, glob = term["carrier_row_convergence"], term["global_continuity_closure"]
    result = {**case, "current_A_per_um": current,
        "identity_error_A_per_um": identity,
        "identity_pass": identity <= max(1e-22, 1e-8*abs(case["reference_current_A_per_um"])),
        "kcl_A_per_um": math.fsum(s["current_A_per_um"] for s in contacts.values()),
        "kcl_over_Id": abs(math.fsum(s["current_A_per_um"] for s in contacts.values()))/max(abs(current),1e-20),
        "local_pass": local["satisfied"], "local_max_ratio": local["max_ratio"],
        "local_violations": local["violation_count"], "qualified_rows": local["qualified_row_count"],
        "ignored_rows": local["ignored_row_count"], "global_pass": glob["satisfied"],
        "psi_residual": term["block_residuals"]["psi"], "electron_residual": term["block_residuals"]["phin"],
        "hole_residual": term["block_residuals"]["phip"]}
    for carrier in ("electron", "hole"):
        for key, value in glob[carrier].items():
            result[carrier+"_global_"+key] = value
        # Current scale is the operator normalization, independently replayed
        # against ContactCurrent; only magnitude is used for source relevance.
        result[carrier+"_source_magnitude_A_per_um"] = abs(glob[carrier]["integrated_source"]*contacts["drain"]["current_scale"])
        result[carrier+"_mismatch_magnitude_A_per_um"] = abs(glob[carrier]["mismatch"]*contacts["drain"]["current_scale"])
    print(case["case"], "local", result["local_violations"], "global", result["global_pass"], flush=True)
    return result


def run():
    verify()
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(audit_case, read(CONTRACT)["cases"]))
    write_csv(OUT/"state_audit.csv", results)
    write(OUT/"audit_summary.json", {"states": len(results), "identity_pass": sum(r["identity_pass"] for r in results),
        "local_pass": sum(r["local_pass"] for r in results), "global_pass": sum(r["global_pass"] for r in results),
        "max_kcl_over_Id": max(r["kcl_over_Id"] for r in results), "nonlinear_solves": 0,
        "strict_contrast_not_yet_run": True})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run", "verify"))
    a = p.parse_args()
    {"prepare": prepare, "run": run, "verify": verify}[a.action]()
