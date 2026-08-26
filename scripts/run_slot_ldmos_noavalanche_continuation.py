#!/usr/bin/env python3
"""Run bounded, adaptive SLOT-LDMOS avalanche-off voltage continuation."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import time
from collections import deque
from pathlib import Path


DEFAULT_TARGETS = [0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0]


def tag(value: float) -> str:
    return f"{value:.12g}".replace("-", "m").replace(".", "p")


def load_last_row(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return rows[-1] if rows else {}


def write_secant_restart(previous_path: Path, current_path: Path,
                          output_path: Path, previous_bias: float,
                          current_bias: float, target_bias: float,
                          max_ratio: float = 2.0) -> float:
    """Apply the C++ sweep predictor semantics to restart CSV fields."""
    denominator = current_bias - previous_bias
    if denominator == 0.0:
        raise ValueError("secant restart requires distinct previous/current biases")
    ratio = max(-max_ratio, min(max_ratio, (target_bias - current_bias) / denominator))
    with previous_path.open(newline="", encoding="utf-8-sig") as stream:
        previous_rows = list(csv.DictReader(stream))
    with current_path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        current_rows = list(reader)
        fieldnames = reader.fieldnames
    if fieldnames is None or len(previous_rows) != len(current_rows):
        raise ValueError("secant restart state shapes do not match")
    extrapolated = {"psi", "phin", "phip"}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for previous, current in zip(previous_rows, current_rows, strict=True):
            if previous.get("node_id") != current.get("node_id"):
                raise ValueError("secant restart node ordering does not match")
            row = dict(current)
            for field in extrapolated:
                row[field] = format(
                    float(current[field]) + ratio * (
                        float(current[field]) - float(previous[field])), ".17g")
            if "electron_qf_reference_V" in row and "electron_qf_increment_V" in row:
                row["electron_qf_increment_V"] = format(
                    float(row["phin"]) - float(row["electron_qf_reference_V"]), ".17g")
            if "hole_qf_reference_V" in row and "hole_qf_increment_V" in row:
                row["hole_qf_increment_V"] = format(
                    float(row["phip"]) - float(row["hole_qf_reference_V"]), ".17g")
            writer.writerow(row)
    return ratio


def render_attempt(base: dict, bundle: Path, target: float,
                   initial_state: Path, attempt_index: int,
                   max_newton_iterations: int | None = None,
                   stall_residual_floor: float | None = None,
                   qf_update_limit_mode: str | None = None,
                   line_search_mode: str | None = None,
                   solver_method: str | None = None,
                   gummel_max_iterations: int | None = None) -> tuple[dict, Path, Path]:
    document = json.loads(json.dumps(base))
    run_rel = Path("outputs") / "stages" / "01_delaunay_noavalanche_adaptive" / (
        f"attempt_{attempt_index:03d}_v_{tag(target)}")
    state_rel = run_rel / "final_state.csv"
    document["output_csv"] = str(run_rel / "iv.csv")
    document["sweep"]["bias_points"] = [target]
    document["sweep"]["initial_state_file"] = str(initial_state.resolve())
    document["sweep"]["write_state_every_point_prefix"] = str(run_rel / "states" / "state")
    document["sweep"]["write_state_file"] = str(state_rel)
    if max_newton_iterations is not None:
        solver = document["solver"]
        solver["max_iter"] = max_newton_iterations
        if "handoff" in solver:
            solver["handoff"]["newton_max_iter"] = max_newton_iterations
        if "carrier_row_convergence" in solver:
            solver["carrier_row_convergence"]["min_newton_max_iter"] = (
                max_newton_iterations)
    if stall_residual_floor is not None:
        document["solver"]["stall_residual_floor"] = stall_residual_floor
    if qf_update_limit_mode is not None:
        document["solver"]["quasi_fermi_update_limit_mode"] = qf_update_limit_mode
    if line_search_mode is not None:
        document["solver"]["line_search_mode"] = line_search_mode
    if solver_method is not None:
        document["solver"]["method"] = solver_method
    if gummel_max_iterations is not None:
        document["solver"].setdefault("handoff", {})["gummel_max_iter"] = (
            gummel_max_iterations)
    if "carrier_row_convergence" in document["solver"]:
        document["solver"]["carrier_row_convergence"]["diagnostic_csv"] = str(
            (bundle / run_rel / "carrier_row_convergence.csv").resolve())
    config_path = bundle / f"adaptive_noavalanche_attempt_{attempt_index:03d}.json"
    return document, config_path, bundle / run_rel


def render_poisson_prewarm(base: dict, bundle: Path, target: float,
                           initial_state: Path, attempt_index: int,
                           max_iterations: int) -> tuple[dict, Path, Path]:
    """Render a Poisson-only contact-bias basin solve before coupled Newton."""
    document = json.loads(json.dumps(base))
    run_rel = Path("outputs") / "stages" / "01_delaunay_noavalanche_adaptive" / (
        f"prewarm_{attempt_index:03d}_v_{tag(target)}")
    document["output_csv"] = str(run_rel / "iv.csv")
    document["sweep"]["bias_points"] = [target]
    document["sweep"]["initial_state_file"] = str(initial_state.resolve())
    document["sweep"]["write_state_every_point_prefix"] = str(
        run_rel / "states" / "state")
    document["sweep"]["write_state_file"] = str(run_rel / "final_state.csv")
    solver = document["solver"]
    solver["method"] = "poisson_only"
    solver["poisson_only_quasi_fermi_initialization"] = "preserve_initial"
    solver["max_iter"] = max_iterations
    if "handoff" in solver:
        solver["handoff"]["newton_max_iter"] = max_iterations
    config_path = bundle / f"adaptive_noavalanche_prewarm_{attempt_index:03d}.json"
    return document, config_path, bundle / run_rel


def render_batch_predictor(base: dict, bundle: Path, biases: list[float],
                           initial_state: Path, max_newton_iterations: int,
                           stall_residual_floor: float | None,
                           qf_update_limit_mode: str | None,
                           line_search_mode: str | None) -> tuple[dict, Path, Path]:
    """Render one multi-point sweep so the C++ secant history stays alive."""
    if len(biases) < 2:
        raise ValueError("batch predictor requires at least two bias points")
    document, _, _ = render_attempt(
        base, bundle, biases[-1], initial_state, 1,
        max_newton_iterations, stall_residual_floor,
        qf_update_limit_mode, line_search_mode, "newton", None)
    run_rel = Path("outputs") / "stages" / "01_delaunay_noavalanche_adaptive" / (
        f"batch_predictor_v_{tag(biases[-1])}")
    document["output_csv"] = str(run_rel / "iv.csv")
    document["sweep"]["bias_points"] = biases
    document["sweep"]["write_state_every_point_prefix"] = str(
        run_rel / "states" / "state")
    document["sweep"]["write_state_file"] = str(run_rel / "final_state.csv")
    continuation = document["sweep"].setdefault("continuation", {})
    continuation["predictor"] = {
        "mode": "secant",
        "fields": ["psi", "phin", "phip"],
        "max_extrapolation_ratio": 2.0,
    }
    if "carrier_row_convergence" in document["solver"]:
        document["solver"]["carrier_row_convergence"]["diagnostic_csv"] = str(
            (bundle / run_rel / "carrier_row_convergence.csv").resolve())
    return (
        document,
        bundle / "adaptive_noavalanche_batch_predictor.json",
        bundle / run_rel,
    )


def execute_config(runner: Path, repo_root: Path, config_path: Path,
                   run_dir: Path, timeout_seconds: float) -> tuple[int, bool, float]:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "states").mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    timed_out = False
    exit_code = -1
    try:
        completed = subprocess.run(
            [str(runner), "--config", str(config_path)],
            cwd=str(repo_root.resolve()), timeout=timeout_seconds,
            check=False, capture_output=True, text=True)
        exit_code = completed.returncode
        (run_dir / "stdout.log").write_text(completed.stdout, encoding="utf-8")
        (run_dir / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        (run_dir / "stdout.log").write_text(exc.stdout or "", encoding="utf-8")
        (run_dir / "stderr.log").write_text(exc.stderr or "", encoding="utf-8")
    return exit_code, timed_out, time.monotonic() - started


def run(args: argparse.Namespace) -> int:
    base_config = args.base_config.resolve()
    bundle = base_config.parent
    base = json.loads(base_config.read_text(encoding="utf-8"))
    runner = args.runner.resolve()
    initial_state = args.initial_state.resolve()
    if not initial_state.exists():
        raise FileNotFoundError(initial_state)

    if args.batch_predictor:
        biases = [args.initial_bias, *args.targets]
        document, config_path, run_dir = render_batch_predictor(
            base, bundle, biases, initial_state, args.max_newton_iterations,
            args.stall_residual_floor, args.qf_update_limit_mode,
            args.line_search_mode)
        if args.closure_feasibility_restoration is not None:
            document["solver"].setdefault("global_continuity_closure", {})[
                "feasibility_restoration"] = args.closure_feasibility_restoration
        if args.closure_feasibility_restoration_factor is not None:
            document["solver"].setdefault("global_continuity_closure", {})[
                "feasibility_restoration_factor"] = (
                    args.closure_feasibility_restoration_factor)
        config_path.write_text(json.dumps(document, indent=2), encoding="utf-8")
        exit_code, timed_out, elapsed = execute_config(
            runner, args.repo_root, config_path, run_dir, args.timeout_seconds)
        rows = []
        iv_path = run_dir / "iv.csv"
        if iv_path.exists():
            with iv_path.open(newline="", encoding="utf-8-sig") as stream:
                rows = list(csv.DictReader(stream))
        summary = {
            "schema": "vela.slot_ldmos.noavalanche_batch_predictor.v1",
            "bias_points_V": biases,
            "exit_code": exit_code,
            "timed_out": timed_out,
            "elapsed_s": elapsed,
            "converged_biases_V": [
                float(row["bias_V"]) for row in rows if row.get("converged") == "1"],
            "rows": rows,
            "config_file": str(config_path.resolve()),
        }
        report_dir = bundle / "reports"
        report_dir.mkdir(parents=True, exist_ok=True)
        (report_dir / "batch_predictor_summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps({
            "exit_code": exit_code,
            "timed_out": timed_out,
            "converged_biases_V": summary["converged_biases_V"],
        }, sort_keys=True))
        return 0 if exit_code == 0 and not timed_out else 1

    requested = deque(args.targets)
    current_bias = args.initial_bias
    current_state = initial_state
    attempts: list[dict[str, object]] = []
    attempt_index = 0
    while requested:
        target = requested[0]
        attempt_index += 1
        if attempt_index > args.max_attempts:
            break
        attempt_initial_state = current_state
        if attempt_index == 1 and args.secant_previous_state is not None:
            previous_state = args.secant_previous_state.resolve()
            if not previous_state.exists():
                raise FileNotFoundError(previous_state)
            predicted_state = bundle / "outputs" / "stages" / (
                "01_delaunay_noavalanche_adaptive") / (
                f"secant_seed_v_{tag(target)}.csv")
            write_secant_restart(
                previous_state, current_state, predicted_state,
                args.secant_previous_bias, current_bias, target,
                args.secant_max_ratio)
            attempt_initial_state = predicted_state
        prewarm_converged = False
        prewarm_elapsed = 0.0
        prewarm_failure = ""
        coupled_initial_state = attempt_initial_state
        if args.poisson_prewarm:
            prewarm_document, prewarm_config, prewarm_dir = render_poisson_prewarm(
                base, bundle, target, current_state, attempt_index,
                args.poisson_prewarm_max_iterations)
            prewarm_config.write_text(
                json.dumps(prewarm_document, indent=2), encoding="utf-8")
            prewarm_exit, prewarm_timed_out, prewarm_elapsed = execute_config(
                runner, args.repo_root, prewarm_config, prewarm_dir,
                args.timeout_seconds)
            prewarm_row = load_last_row(prewarm_dir / "iv.csv")
            prewarm_state = prewarm_dir / "final_state.csv"
            prewarm_converged = (
                prewarm_exit == 0 and prewarm_row.get("converged") == "1"
                and prewarm_state.exists())
            prewarm_failure = (
                "timeout" if prewarm_timed_out
                else prewarm_row.get("failure_reason", ""))
            if prewarm_converged:
                coupled_initial_state = prewarm_state
        document, config_path, run_dir = render_attempt(
            base, bundle, target, coupled_initial_state, attempt_index,
            args.max_newton_iterations, args.stall_residual_floor,
            args.qf_update_limit_mode, args.line_search_mode,
            args.solver_method, args.gummel_max_iterations)
        if args.save_rejected_states:
            document["solver"]["diagnostic_history"] = True
            document["sweep"].setdefault("diagnostics", {})["newton_history"] = {
                "enabled": True,
                "attempts_csv_file": str(run_dir.relative_to(bundle) / "attempts.csv"),
                "iterations_csv_file": str(run_dir.relative_to(bundle) / "iterations.csv"),
                "rejected_state_directory": str(
                    run_dir.relative_to(bundle) / "rejected_states"),
            }
        if args.return_best_feasible_on_failure:
            document["solver"]["return_best_feasible_on_failure"] = True
        if args.closure_feasibility_restoration is not None:
            document["solver"].setdefault("global_continuity_closure", {})[
                "feasibility_restoration"] = args.closure_feasibility_restoration
        if args.closure_feasibility_restoration_factor is not None:
            document["solver"].setdefault("global_continuity_closure", {})[
                "feasibility_restoration_factor"] = (
                    args.closure_feasibility_restoration_factor)
        config_path.write_text(json.dumps(document, indent=2), encoding="utf-8")
        exit_code, timed_out, elapsed = execute_config(
            runner, args.repo_root, config_path, run_dir, args.timeout_seconds)
        iv_row = load_last_row(run_dir / "iv.csv")
        next_state = run_dir / "final_state.csv"
        converged = exit_code == 0 and iv_row.get("converged") == "1" and next_state.exists()
        attempts.append({
            "attempt": attempt_index, "start_bias_V": current_bias,
            "target_bias_V": target, "elapsed_s": elapsed,
            "poisson_prewarm_enabled": args.poisson_prewarm,
            "poisson_prewarm_converged": prewarm_converged,
            "poisson_prewarm_elapsed_s": prewarm_elapsed,
            "poisson_prewarm_failure_reason": prewarm_failure,
            "timed_out": timed_out, "exit_code": exit_code,
            "converged": converged, "failure_reason": iv_row.get("failure_reason", ""),
            "iterations": iv_row.get("iterations", ""),
            "current_A_per_um": iv_row.get("current_total_A_per_um", ""),
            "state_file": str(next_state.resolve()) if converged else "",
            "config_file": str(config_path.resolve()),
        })
        if converged:
            current_bias = target
            current_state = next_state
            requested.popleft()
            continue
        midpoint = 0.5 * (current_bias + target)
        if target - current_bias <= args.min_step_V or midpoint <= current_bias:
            break
        requested.appendleft(midpoint)

    report_dir = bundle / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    csv_path = report_dir / "adaptive_noavalanche_attempts.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(attempts[0]) if attempts else ["attempt"])
        writer.writeheader()
        writer.writerows(attempts)
    completed_all = not requested
    summary = {
        "schema": "vela.slot_ldmos.noavalanche_adaptive.v1",
        "completed": completed_all,
        "last_successful_bias_V": current_bias,
        "last_successful_state": str(current_state),
        "remaining_targets_V": list(requested),
        "timeout_seconds": args.timeout_seconds,
        "max_newton_iterations": args.max_newton_iterations,
        "stall_residual_floor": args.stall_residual_floor,
        "qf_update_limit_mode": args.qf_update_limit_mode,
        "line_search_mode": args.line_search_mode,
        "closure_feasibility_restoration": args.closure_feasibility_restoration,
        "closure_feasibility_restoration_factor": (
            args.closure_feasibility_restoration_factor),
        "solver_method": args.solver_method,
        "gummel_max_iterations": args.gummel_max_iterations,
        "poisson_prewarm": args.poisson_prewarm,
        "poisson_prewarm_max_iterations": args.poisson_prewarm_max_iterations,
        "min_step_V": args.min_step_V,
        "attempts": attempts,
    }
    (report_dir / "adaptive_noavalanche_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({key: summary[key] for key in (
        "completed", "last_successful_bias_V", "remaining_targets_V")}, sort_keys=True))
    return 0 if completed_all else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--initial-state", type=Path, required=True)
    parser.add_argument("--initial-bias", type=float, default=0.105)
    parser.add_argument("--targets", type=float, nargs="+", default=DEFAULT_TARGETS)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--max-newton-iterations", type=int, default=8)
    parser.add_argument("--stall-residual-floor", type=float)
    parser.add_argument(
        "--qf-update-limit-mode",
        choices=("componentwise_poisson_recorrect", "uniform_trust_region"))
    parser.add_argument(
        "--line-search-mode",
        choices=("merit", "block_filter", "block_closure_filter"))
    parser.add_argument(
        "--closure-feasibility-restoration",
        choices=("off", "linear_qf_projection"))
    parser.add_argument(
        "--closure-feasibility-restoration-factor", type=float)
    parser.add_argument(
        "--solver-method", choices=("newton", "gummel_newton"), default="newton")
    parser.add_argument("--gummel-max-iterations", type=int)
    parser.add_argument("--poisson-prewarm", action="store_true")
    parser.add_argument("--poisson-prewarm-max-iterations", type=int, default=20)
    parser.add_argument("--min-step-V", type=float, default=0.01)
    parser.add_argument("--max-attempts", type=int, default=100)
    parser.add_argument("--batch-predictor", action="store_true")
    parser.add_argument("--secant-previous-state", type=Path)
    parser.add_argument("--secant-previous-bias", type=float, default=0.0)
    parser.add_argument("--secant-max-ratio", type=float, default=2.0)
    parser.add_argument("--save-rejected-states", action="store_true")
    parser.add_argument("--return-best-feasible-on-failure", action="store_true")
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
