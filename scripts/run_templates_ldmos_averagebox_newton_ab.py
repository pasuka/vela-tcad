#!/usr/bin/env python3
"""Run the gated Templates/LDMOS external-AverageBox Newton qualification.

The candidate changes only carrier-transport edge couples.  Poisson geometry,
source volumes, physics, IALMob, and predictor remain identical to the supplied
baseline.  A full same-bias reclose is attempted only after the fixed-state and
accepted-one-step gates pass.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any


FROZEN_HOTSPOTS = (3747, 10233, 4492, 4538)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def vector_norm(rows: list[dict[str, str]], field: str) -> dict[str, float]:
    values = [float(row[field]) for row in rows]
    return {
        "l2": math.sqrt(sum(value * value for value in values)),
        "maximum_abs": max((abs(value) for value in values), default=0.0),
    }


def status_from_stdout(stdout: str) -> dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            return value
    raise RuntimeError(f"runner did not emit a JSON status:\n{stdout}")


def run_config(
    runner: Path, config: dict[str, Any], output: Path, name: str,
    *, require_success: bool,
) -> tuple[dict[str, Any], Path, int]:
    config_path = output / f"{name}.json"
    log_path = output / f"{name}.log"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    log_path.write_text(completed.stdout, encoding="utf-8")
    status = status_from_stdout(completed.stdout)
    if require_success and completed.returncode != 0:
        raise RuntimeError(json.dumps({
            "return_code": completed.returncode,
            "config": str(config_path),
            "log": str(log_path),
            "status": status,
        }, indent=2))
    return status, config_path, completed.returncode


def with_profile(
    base: dict[str, Any], profile_csv: Path, expected_edges: int,
) -> dict[str, Any]:
    candidate = json.loads(json.dumps(base))
    geometry = candidate.setdefault("mesh_geometry", {})
    geometry["carrier_transport_couple_profile"] = (
        "templates_ldmos_external_averagebox"
    )
    geometry["external_averagebox_couples_file"] = str(profile_csv.resolve())
    geometry["external_averagebox_expected_edges"] = expected_edges
    return candidate


def without_profile(base: dict[str, Any]) -> dict[str, Any]:
    baseline = json.loads(json.dumps(base))
    geometry = baseline.get("mesh_geometry")
    if isinstance(geometry, dict):
        for key in (
            "carrier_transport_couple_profile",
            "external_averagebox_couples_file",
            "external_averagebox_expected_edges",
        ):
            geometry.pop(key, None)
    return baseline


def with_contact_boundary_reconstruction(
    base: dict[str, Any], mode: str | None,
) -> dict[str, Any]:
    config = json.loads(json.dumps(base))
    if mode is not None:
        config.setdefault("solver", {})["contact_boundary_reconstruction"] = mode
    return config


def step_probe_config(
    base: dict[str, Any], output_csv: Path,
) -> dict[str, Any]:
    config = json.loads(json.dumps(base))
    config["simulation_type"] = "newton_step_probe"
    config["output_csv"] = str(output_csv.resolve())
    return config


def residual_probe_config(
    base: dict[str, Any], state: Path, output_csv: Path,
) -> dict[str, Any]:
    config = json.loads(json.dumps(base))
    config["simulation_type"] = "newton_residual_probe"
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str(output_csv.resolve())
    return config


def solve_config(
    base: dict[str, Any], output_state: Path, max_iter: int,
) -> dict[str, Any]:
    config = json.loads(json.dumps(base))
    config["simulation_type"] = "newton_solve_from_state"
    config.pop("output_csv", None)
    config["output_state_file"] = str(output_state.resolve())
    config["solver"]["max_iter"] = max_iter
    config["solver"]["diagnostics"] = True
    config["solver"]["verbose"] = False
    return config


def per_node_ratio(
    baseline: list[dict[str, str]], candidate: list[dict[str, str]],
    field: str,
) -> dict[int, float]:
    by_node = {int(row["node_id"]): float(row[field]) for row in baseline}
    result: dict[int, float] = {}
    for row in candidate:
        node = int(row["node_id"])
        result[node] = abs(float(row[field])) / max(abs(by_node[node]), 1.0e-300)
    return result


def compare_step_probes(
    baseline_rows: list[dict[str, str]], candidate_rows: list[dict[str, str]],
) -> dict[str, Any]:
    baseline_phin = vector_norm(baseline_rows, "phin_residual")
    candidate_phin = vector_norm(candidate_rows, "phin_residual")
    baseline_psi = {int(row["node_id"]): float(row["psi_residual"])
                    for row in baseline_rows}
    maximum_psi_difference = max(
        abs(float(row["psi_residual"]) - baseline_psi[int(row["node_id"])])
        for row in candidate_rows
    )
    ratios = per_node_ratio(
        baseline_rows, candidate_rows, "phin_residual"
    )
    hotspot_ratios = {str(node): ratios[node] for node in FROZEN_HOTSPOTS}
    l2_ratio = candidate_phin["l2"] / max(baseline_phin["l2"], 1.0e-300)
    maximum_ratio = (
        candidate_phin["maximum_abs"]
        / max(baseline_phin["maximum_abs"], 1.0e-300)
    )
    gate = (
        l2_ratio <= 0.5
        and maximum_ratio <= 0.5
        and all(value <= 0.5 for value in hotspot_ratios.values())
        and maximum_psi_difference <= 1.0e-20
    )
    return {
        "baseline_phin_residual": baseline_phin,
        "candidate_phin_residual": candidate_phin,
        "candidate_over_baseline": {
            "l2": l2_ratio,
            "maximum_abs": maximum_ratio,
        },
        "frozen_hotspot_abs_ratios": hotspot_ratios,
        "maximum_absolute_initial_psi_residual_difference":
            maximum_psi_difference,
        "gate": {
            "requires_l2_ratio_at_most": 0.5,
            "requires_maximum_ratio_at_most": 0.5,
            "requires_each_hotspot_ratio_at_most": 0.5,
            "requires_initial_psi_residual_difference_at_most": 1.0e-20,
            "passed": gate,
        },
    }


def max_state_delta(initial: Path, updated: Path) -> dict[str, float]:
    initial_rows = read_rows(initial)
    updated_rows = read_rows(updated)
    initial_by_node = {int(row["node_id"]): row for row in initial_rows}
    fields = ("psi", "phin", "phip")
    return {
        field: max(
            abs(float(row[field]) - float(initial_by_node[int(row["node_id"])][field]))
            for row in updated_rows
        )
        for field in fields
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--transport-couples", type=Path, required=True)
    parser.add_argument("--expected-edges", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sentaurus-current-A-per-um", type=float)
    parser.add_argument("--legacy-reclose-current-A-per-um", type=float)
    parser.add_argument(
        "--contact-boundary-reconstruction",
        choices=("dominant_signed_contact_mean", "legacy_node_local"),
        help=(
            "Optional explicit Ohmic-contact reconstruction override. "
            "Omission preserves the supplied deck/default."
        ),
    )
    parser.add_argument("--run-reclose", action="store_true")
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base = with_contact_boundary_reconstruction(
        json.loads(args.base_config.read_text(encoding="utf-8")),
        args.contact_boundary_reconstruction,
    )
    baseline = without_profile(base)
    candidate = with_profile(
        base, args.transport_couples.resolve(), args.expected_edges
    )

    step_status: dict[str, Any] = {}
    step_rows: dict[str, list[dict[str, str]]] = {}
    for name, config in (("baseline", baseline), ("candidate", candidate)):
        csv_path = output / f"{name}_step.csv"
        status, _, _ = run_config(
            args.runner.resolve(), step_probe_config(config, csv_path),
            output, f"{name}_step", require_success=True,
        )
        step_status[name] = status
        step_rows[name] = read_rows(csv_path)
    fixed_gate = compare_step_probes(step_rows["baseline"], step_rows["candidate"])

    summary: dict[str, Any] = {
        "schema": "vela.templates_ldmos.averagebox_newton_ab.v1",
        "contract": {
            "bias": "Vd=0.1 V, Vg=0.5 V",
            "candidate": "templates_ldmos_external_averagebox_carrier_transport_only",
            "poisson_geometry": "unchanged_mesh_default",
            "source_volume": "unchanged_barycentric",
            "ialmob": "off",
            "predictor": "off",
            "contact_boundary_reconstruction": base.get("solver", {}).get(
                "contact_boundary_reconstruction", "solver_default"
            ),
            "production_default_changed": False,
        },
        "profile": {
            "path": str(args.transport_couples.resolve()),
            "sha256": sha256(args.transport_couples.resolve()),
            "expected_edges": args.expected_edges,
        },
        "fixed_state_step_probe": {
            "statuses": step_status,
            **fixed_gate,
        },
    }

    if fixed_gate["gate"]["passed"]:
        one_step: dict[str, Any] = {}
        accepted = True
        for name, config in (("baseline", baseline), ("candidate", candidate)):
            state_path = output / f"{name}_one_newton_state.csv"
            status, config_path, return_code = run_config(
                args.runner.resolve(), solve_config(config, state_path, 1),
                output, f"{name}_one_newton", require_success=False,
            )
            if not state_path.is_file():
                accepted = False
                one_step[name] = {
                    "status": status, "return_code": return_code,
                    "config": str(config_path), "state_written": False,
                }
                continue
            residual_csv = output / f"{name}_one_newton_residual.csv"
            residual_status, _, _ = run_config(
                args.runner.resolve(),
                residual_probe_config(config, state_path, residual_csv),
                output, f"{name}_one_newton_residual", require_success=True,
            )
            initial_combined = float(
                step_status[name]["block_residuals"]["combined"]
            )
            final_combined = float(
                residual_status["block_residuals"]["combined"]
            )
            decreased = final_combined < initial_combined
            accepted = accepted and decreased and status.get("iterations") == 1
            one_step[name] = {
                "status": status,
                "return_code": return_code,
                "state_written": True,
                "maximum_state_delta_V": max_state_delta(
                    Path(config["state_file"]), state_path
                ),
                "initial_combined_residual": initial_combined,
                "post_one_step_combined_residual": final_combined,
                "post_over_initial": final_combined / max(initial_combined, 1.0e-300),
                "accepted_residual_decrease": decreased,
                "residual_csv": str(residual_csv),
            }
        summary["one_newton_ab"] = {
            "branches": one_step,
            "gate": {
                "requires_each_branch_one_iteration": True,
                "requires_each_branch_combined_residual_decrease": True,
                "passed": accepted,
            },
        }

        if args.run_reclose and accepted:
            reclose_state = output / "candidate_reclose_state.csv"
            status, config_path, return_code = run_config(
                args.runner.resolve(), solve_config(candidate, reclose_state, 400),
                output, "candidate_reclose", require_success=False,
            )
            contacts = status.get("contact_currents_A_per_um", {})
            contact_values = [float(value) for value in contacts.values()]
            current_scale = max((abs(value) for value in contact_values), default=1.0e-300)
            reclose: dict[str, Any] = {
                "status": status,
                "return_code": return_code,
                "config": str(config_path),
                "state_written": reclose_state.is_file(),
                "maximum_state_delta_V": (
                    max_state_delta(Path(candidate["state_file"]), reclose_state)
                    if reclose_state.is_file() else None
                ),
                "terminal_kcl_relative_error": (
                    abs(sum(contact_values)) / current_scale
                ),
                "gate": {
                    "requires_converged": True,
                    "passed": bool(status.get("converged", False)),
                },
            }
            drain_current = float(contacts.get("drain", math.nan))
            if args.sentaurus_current_A_per_um is not None:
                ratio_to_sentaurus = (
                    abs(drain_current) / abs(args.sentaurus_current_A_per_um)
                )
                reclose["sentaurus_current_A_per_um"] = (
                    args.sentaurus_current_A_per_um
                )
                reclose["drain_abs_over_sentaurus"] = ratio_to_sentaurus
                reclose["drain_log10_abs_error_dex"] = abs(
                    math.log10(ratio_to_sentaurus)
                )
                reclose["stage3_current_gate"] = {
                    "requires_log10_abs_error_dex_at_most": 0.20,
                    "passed": abs(math.log10(ratio_to_sentaurus)) <= 0.20,
                }
            if args.legacy_reclose_current_A_per_um is not None:
                reclose["legacy_reclose_current_A_per_um"] = (
                    args.legacy_reclose_current_A_per_um
                )
                reclose["candidate_abs_over_legacy"] = (
                    abs(drain_current)
                    / abs(args.legacy_reclose_current_A_per_um)
                )
            summary["candidate_same_bias_reclose"] = reclose
    else:
        summary["one_newton_ab"] = {"skipped": "fixed_state_gate_failed"}
        if args.run_reclose:
            summary["candidate_same_bias_reclose"] = {
                "skipped": "fixed_state_gate_failed"
            }

    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "summary": str(summary_path),
        "fixed_state_gate": fixed_gate["gate"],
        "one_newton_gate": summary.get("one_newton_ab", {}).get("gate"),
        "reclose_gate": summary.get("candidate_same_bias_reclose", {}).get("gate"),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
