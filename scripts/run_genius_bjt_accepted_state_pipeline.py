#!/usr/bin/env python3
"""Build the unique 31-point accepted-state chain for the Genius BJT M1 case."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from decimal import Decimal
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_OUTPUT = BUILD_ROOT / "m1_accepted_states"
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
NODE_COUNT = 5611
MAX_POST_ACCEPT_UPDATE_V = 1.0e-6
MAX_KCL_A_PER_UM = 1.0e-10
ACTIVE_CARRIER_DENSITY_M3 = 1.0e16


def absolute(path: Path) -> str:
    return str(path.resolve())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def common_config() -> dict:
    cfg = load_json(FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json")
    cfg.pop("sweep", None)
    cfg["mesh_file"] = absolute(FIXTURE / "vela" / "input" / "mesh.json")
    cfg["node_doping_file"] = absolute(FIXTURE / "vela" / "input" / "doping.csv")
    cfg["materials_file"] = absolute(FIXTURE / "vela" / "materials_sentaurus2022.json")
    cfg["solver"]["carrier_row_qualified_stall_acceptance"] = True
    cfg["solver"]["carrier_row_convergence"] = {
        "mode": "enforce",
        "eps_row": 1.0e-3,
        "scale_floor": 1.0e-30,
        "min_source_scale": 1.0e-18,
        "min_source_scale_fraction": 1.0e-3,
        "min_source_global_fraction": 1.0e-6,
        "min_carrier_density_m3": ACTIVE_CARRIER_DENSITY_M3,
        "min_flux_scale_fraction": 1.0e-6,
        "min_flux_scale": 0.0,
        "min_newton_max_iter": 200,
    }
    cfg["solver"]["global_continuity_closure"] = {
        "mode": "enforce",
        "tolerance": 1.0e-6,
        "source_floor": 1.0e-18,
    }
    return cfg


def set_collector_bias(cfg: dict, bias: float) -> None:
    for contact in cfg["contacts"]:
        if contact["name"] == "collector":
            contact["bias"] = bias
            return
    raise KeyError("collector contact not found")


def run(runner: Path, cfg: dict, config_path: Path, log_prefix: Path) -> dict:
    write_json(config_path, cfg)
    process = subprocess.run(
        [absolute(runner), "--config", absolute(config_path)],
        cwd=REPO,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    log_prefix.with_suffix(".stdout.log").write_text(process.stdout, encoding="utf-8")
    log_prefix.with_suffix(".stderr.log").write_text(process.stderr, encoding="utf-8")
    if process.returncode != 0:
        raise RuntimeError(f"runner failed for {config_path}; see {log_prefix}.stderr.log")
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError(f"runner returned no JSON for {config_path}")
    return json.loads(lines[-1])


def step_summary(path: Path, state_rows: list[dict[str, str]]) -> dict:
    rows = read_csv(path)
    state = {int(row["node_id"]): row for row in state_rows}
    keys = ("delta_psi_V", "delta_phin_V", "delta_phip_V")
    maxima = {key: max(abs(float(row[key])) for row in rows) for key in keys}
    l2 = math.sqrt(sum(float(row[key]) ** 2 for row in rows for key in keys))
    active_phip = [
        abs(float(row["delta_phip_V"]) - float(row["delta_psi_V"])) for row in rows
        if float(state[int(row["node_id"])]["holes_m3"]) >= ACTIVE_CARRIER_DENSITY_M3
    ]
    active_phin = [
        abs(float(row["delta_phin_V"]) - float(row["delta_psi_V"])) for row in rows
        if float(state[int(row["node_id"])]["electrons_m3"]) >= ACTIVE_CARRIER_DENSITY_M3
    ]
    driving_electron = [
        abs(float(row["delta_phin_V"]) - float(row["delta_psi_V"])) for row in rows
    ]
    driving_hole = [
        abs(float(row["delta_phip_V"]) - float(row["delta_psi_V"])) for row in rows
    ]
    return {
        "row_count": len(rows),
        "max_abs_delta_psi_V": maxima["delta_psi_V"],
        "max_abs_delta_phin_V": maxima["delta_phin_V"],
        "max_abs_delta_phip_V": maxima["delta_phip_V"],
        "max_abs_delta_any_V": max(maxima.values()),
        "max_abs_delta_phin_minus_psi_V": max(driving_electron, default=0.0),
        "max_abs_delta_phip_minus_psi_V": max(driving_hole, default=0.0),
        "max_abs_delta_relevant_V": max(
            maxima["delta_psi_V"], max(active_phin, default=0.0), max(active_phip, default=0.0)
        ),
        "active_electron_node_count": len(active_phin),
        "active_hole_node_count": len(active_phip),
        "delta_l2_V": l2,
        "sha256": sha256(path),
    }


def recombination_summary(path: Path) -> dict:
    rows = read_csv(path)
    electron = [float(row["electron_recombination"]) for row in rows]
    hole = [float(row["hole_recombination"]) for row in rows]
    return {
        "row_count": len(rows),
        "electron_recombination_sum": sum(electron),
        "hole_recombination_sum": sum(hole),
        "electron_recombination_max_abs": max(map(abs, electron)),
        "hole_recombination_max_abs": max(map(abs, hole)),
        "sha256": sha256(path),
    }


def probe_point(runner: Path, output: Path, index: int, bias: float) -> tuple[dict, dict]:
    token = f"vce_{index:03d}"
    state = output / "states" / f"{token}.csv"
    state_rows = read_csv(state)
    step_results = {}
    for label, limit in (("raw", 0.0), ("capped", 0.1)):
        path = output / f"steps_{label}" / f"{token}.csv"
        probe = common_config()
        set_collector_bias(probe, bias)
        probe.update({"simulation_type": "newton_step_probe", "state_file": absolute(state), "output_csv": absolute(path)})
        probe["solver"]["quasi_fermi_update_limit_V"] = limit
        run(runner, probe, output / "configs" / f"{token}_step_{label}.json", output / "logs" / f"{token}_step_{label}")
        step_results[label] = step_summary(path, state_rows)
    terms_csv = output / "recombination" / f"{token}.csv"
    terms = common_config()
    set_collector_bias(terms, bias)
    terms.update({"simulation_type": "newton_carrier_term_probe", "state_file": absolute(state), "output_csv": absolute(terms_csv), "carrier_term_probe": {"solved_equation_terms": True}})
    run(runner, terms, output / "configs" / f"{token}_terms.json", output / "logs" / f"{token}_terms")
    updates = {"index": index, "vce_V": bias, **{f"raw_{key}": value for key, value in step_results["raw"].items()}, **{f"capped_{key}": value for key, value in step_results["capped"].items()}}
    return updates, {"index": index, "vce_V": bias, **recombination_summary(terms_csv)}


def main() -> int:
    args = parse_args()
    output = args.output_root.resolve()
    for name in ("states", "fields", "recombination", "steps_raw", "steps_capped", "configs", "logs"):
        (output / name).mkdir(parents=True, exist_ok=True)
    biases = [float(Decimal(index) / Decimal(10)) for index in range(31)]
    records, terminal_rows = [], []
    previous = BUILD_ROOT / "vela_wp3_wp5" / "m1_vbe070_state.csv"
    previous_bias = 0.0
    for index, bias in enumerate(biases):
        token = f"vce_{index:03d}"
        state, vtk = output / "states" / f"{token}.csv", output / "fields" / f"{token}.vtk"
        candidate = previous
        if index:
            candidate = output / "states" / f"{token}_continuation.csv"
            cfg = common_config(); cfg["simulation_type"] = "dc_sweep"; set_collector_bias(cfg, previous_bias)
            cfg["output_csv"] = absolute(output / "logs" / f"{token}_continuation.csv")
            cfg["sweep"] = {"mode": "iv", "contact": "collector", "current_contact": "collector", "start": previous_bias, "stop": bias, "step": 0.05, "initial_step": 0.05, "min_step": 1e-8, "max_step": 0.05, "growth_factor": 1.0, "shrink_factor": 0.5, "max_retries": 20, "write_vtk": False, "initial_state_file": absolute(previous), "write_state_file": absolute(candidate)}
            run(args.runner, cfg, output / "configs" / f"{token}_continuation.json", output / "logs" / f"{token}_continuation_run")
            if not candidate.is_file(): raise RuntimeError(f"continuation did not produce {candidate}")
        cfg = common_config(); set_collector_bias(cfg, bias)
        cfg.update({"simulation_type": "newton_solve_from_state", "state_file": absolute(candidate), "output_state_file": absolute(state), "output_vtk": absolute(vtk)})
        status = run(args.runner, cfg, output / "configs" / f"{token}_solve.json", output / "logs" / f"{token}_solve")
        if not status.get("converged", False): raise RuntimeError(f"fixed-bias solve failed at VCE={bias}: {status['failure_reason']}")
        rows, currents = read_csv(state), status["contact_currents_A_per_um"]
        finite = len(rows) == NODE_COUNT and all(math.isfinite(float(row[key])) for row in rows for key in ("psi", "phin", "phip", "electrons_m3", "holes_m3"))
        kcl = sum(float(value) for value in currents.values())
        carrier_rows = status["carrier_row_convergence"]
        closure = status["global_continuity_closure"]
        record = {"index": index, "vce_V": bias, "solver_converged": True, "iterations": status["iterations"], "convergence_reason": status["convergence_reason"], "final_residual": status["final_residual"], "carrier_rows_satisfied": carrier_rows["satisfied"], "carrier_rows_qualified": carrier_rows["qualified_row_count"], "carrier_rows_ignored": carrier_rows["ignored_row_count"], "carrier_row_max_ratio": carrier_rows["max_ratio"], "global_continuity_satisfied": closure["satisfied"], "global_electron_closure_ratio": closure["electron"]["ratio"], "global_hole_closure_ratio": closure["hole"]["ratio"], "state_row_count": len(rows), "finite_primary_fields": finite, "state_sha256": sha256(state), "parent_state_sha256": sha256(previous), "kcl_A_per_um": kcl, "state_file": state.relative_to(output).as_posix(), "field_file": vtk.relative_to(output).as_posix(), "recombination_file": f"recombination/{token}.csv"}
        records.append(record)
        terminal_rows.append({"index": index, "vce_V": bias, "collector_A_per_um": currents["collector"], "base_A_per_um": currents["base"], "emitter_A_per_um": currents["emitter"], "kcl_A_per_um": kcl, "state_sha256": record["state_sha256"]})
        write_csv(output / "state_solve_ledger.csv", records); write_csv(output / "terminal_currents.csv", terminal_rows)
        previous, previous_bias = state, bias
        print(f"state {token}: converged, KCL={kcl:.3e} A/um", flush=True)

    update_rows, recombination_rows = [], []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(probe_point, args.runner, output, index, bias): index for index, bias in enumerate(biases)}
        for future in as_completed(futures):
            updates, recombination = future.result(); update_rows.append(updates); recombination_rows.append(recombination)
            print(f"derived vce_{futures[future]:03d}", flush=True)
    update_rows.sort(key=lambda row: row["index"]); recombination_rows.sort(key=lambda row: row["index"])
    write_csv(output / "update_diagnostics.csv", update_rows); write_csv(output / "recombination_summary.csv", recombination_rows)
    by_index = {row["index"]: row for row in update_rows}
    for record in records:
        update = by_index[record["index"]]
        record["raw_max_update_V"] = update["raw_max_abs_delta_any_V"]
        record["raw_max_relevant_update_V"] = update["raw_max_abs_delta_relevant_V"]
        record["capped_max_update_V"] = update["capped_max_abs_delta_any_V"]
        record["accepted"] = record["finite_primary_fields"] and record["carrier_rows_satisfied"] and record["global_continuity_satisfied"] and abs(record["kcl_A_per_um"]) <= MAX_KCL_A_PER_UM and record["raw_max_relevant_update_V"] <= MAX_POST_ACCEPT_UPDATE_V
    write_csv(output / "accepted_states.csv", records)

    manifest = {
        "schema_version": 2,
        "contract": {
            "biases": "VBE=0.70 V; VCE=0.00..3.00 V in exact 0.10 V increments",
            "state_lineage": "each fixed-bias state starts from the immediately preceding accepted state",
            "solver_converged_required": True,
            "state_rows_required": NODE_COUNT,
            "finite_primary_fields_required": True,
            "max_post_accept_raw_update_V": MAX_POST_ACCEPT_UPDATE_V,
            "post_accept_update_definition": "max |delta_psi| and active-carrier |delta_phi_carrier-delta_psi|",
            "carrier_update_density_mask_m3": ACTIVE_CARRIER_DENSITY_M3,
            "carrier_row_convergence_required": True,
            "global_continuity_closure_required": True,
            "max_terminal_kcl_A_per_um": MAX_KCL_A_PER_UM,
        },
        "accepted_count": sum(bool(row["accepted"]) for row in records),
        "point_count": len(records),
        "overall_pass": len(records) == 31 and all(row["accepted"] for row in records),
        "runner": {"path": absolute(args.runner), "sha256": sha256(args.runner)},
        "inputs": {
            "base_config_sha256": sha256(FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json"),
            "mesh_sha256": sha256(FIXTURE / "vela" / "input" / "mesh.json"),
            "doping_sha256": sha256(FIXTURE / "vela" / "input" / "doping.csv"),
            "materials_sha256": sha256(FIXTURE / "vela" / "materials_sentaurus2022.json"),
        },
        "points": records,
    }
    write_json(output / "manifest.json", manifest)
    print(json.dumps({
        "overall_pass": manifest["overall_pass"],
        "accepted_count": manifest["accepted_count"],
        "point_count": manifest["point_count"],
        "output": absolute(output),
    }))
    return 0 if manifest["overall_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
