#!/usr/bin/env python3
"""Run bounded SLOT-LDMOS continuation, preserving the base deck's physics."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import hashlib
import subprocess
import time
from collections import deque
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slot_state_archive as states


DEFAULT_TARGETS = [0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0]


def physics_branch(base: dict) -> str:
    impact = base.get("solver", {}).get("impact_ionization", "none")
    model = impact.get("model", "none") if isinstance(impact, dict) else impact
    coupling = impact.get("coupling_mode", "self_consistent") if isinstance(impact, dict) else "self_consistent"
    return "avalanche" if model not in (None, "none", "off", "") and coupling == "self_consistent" else "noavalanche"


def stage_name(base: dict) -> str:
    branch = physics_branch(base)
    return f"{'02' if branch == 'avalanche' else '01'}_delaunay_{branch}_adaptive"


def require_direct_device_bias(base: dict) -> None:
    circuit = base.get("sweep", {}).get("external_circuit", {})
    if circuit and circuit.get("enabled", True) and circuit.get("mode", "none") != "none":
        raise ValueError("this continuation driver requires direct device bias; external-circuit audits need the solved inner voltage")


def tag(value: float) -> str:
    return f"{value:.12g}".replace("-", "m").replace(".", "p")


def load_last_row(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return rows[-1] if rows else {}


def load_rejected_final_state(path: Path) -> Path | None:
    """Return the solver's persisted final accepted iterate for a failed point."""
    value = load_last_row(path).get("rejected_final_state_file", "").strip()
    if not value:
        return None
    candidate = Path(value)
    return candidate.resolve() if candidate.exists() else None


def write_adaptive_checkpoint(bundle: Path, checkpoint: dict) -> None:
    """Persist one small recovery record per finished attempt, not per iteration.

    This records completed evidence only; it does not claim a process is alive.
    A restart must still validate the state/configuration before reusing it.
    """
    target = bundle / "reports" / "adaptive_checkpoint.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(checkpoint, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, target)


def redirect_diagnostics(document: dict, run_dir: Path) -> None:
    """Keep each attempt's diagnostics separate from the base deck evidence."""
    solver = document.get("solver", {})
    for section, fields in (
        ("global_continuity_closure", {"audit_jsonl_file": "closure_audit.jsonl"}),
        ("carrier_row_convergence", {"diagnostic_csv": "carrier_row_convergence.csv", "trace_csv": "carrier_trace.csv"}),
        ("local_update_diagnostics", {"csv_file": "local_updates.csv"}),
    ):
        value = solver.get(section)
        if isinstance(value, dict):
            for key, filename in fields.items():
                if value.get(key):
                    value[key] = str((run_dir / filename).resolve())
    history = document.get("sweep", {}).get("diagnostics", {}).get("newton_history")
    if isinstance(history, dict):
        for key, filename in (("attempts_csv_file", "attempts.csv"),
                              ("iterations_csv_file", "iterations.csv"),
                              ("rejected_state_directory", "rejected_states")):
            if history.get(key):
                history[key] = str((run_dir / filename).resolve())


def write_secant_restart(previous_path: Path, current_path: Path,
                          output_path: Path, previous_bias: float,
                          current_bias: float, target_bias: float,
                          max_ratio: float = 2.0, *, document: dict,
                          config_directory: Path, ratio_override: float | None = None) -> float:
    """Predict HDF5 unknowns without folding sub-ULP split QF increments."""
    if not all(math.isfinite(v) for v in (previous_bias, current_bias, target_bias, max_ratio)) or max_ratio <= 0:
        raise ValueError("invalid secant bias/ratio")
    if ratio_override is None:
        denominator = current_bias - previous_bias
        if denominator == 0:
            raise ValueError("secant restart requires distinct previous/current biases")
        ratio = max(-max_ratio, min(max_ratio, (target_bias-current_bias)/denominator))
    else:
        # Outer-voltage continuation remains defined at an inner-voltage fold
        # or plateau. The caller binds its own coordinate and bounded ratio.
        if not math.isfinite(ratio_override) or abs(ratio_override) > max_ratio:
            raise ValueError("invalid explicit secant ratio")
        ratio = ratio_override
    previous, _ = states.read(previous_path, document, config_directory)
    current, metadata = states.read(current_path, document, config_directory)
    if previous.keys() != current.keys():
        raise ValueError("secant restart state fields do not match")
    output = {key: value.copy() for key, value in current.items()}
    for field in ("psi", "phin", "phip"):
        output[field] = current[field] + ratio * (current[field] - previous[field])
    for carrier, field in (("electron", "phin"), ("hole", "phip")):
        ref, inc = f"{carrier}_qf_reference_V", f"{carrier}_qf_increment_V"
        if ref in output:
            for component in (ref, inc):
                output[component] = [math.fsum((float(now), ratio * (float(now)-float(old))))
                    for now, old in zip(current[component], previous[component], strict=True)]
            output[field] = [math.fsum((a,b)) for a,b in zip(output[ref], output[inc], strict=True)]
    states.write_candidate(output_path, output, metadata, document, [previous_path, current_path],
                           target_bias, "secant_candidate")
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
    document["state_format"] = "hdf5"
    run_rel = Path("outputs") / "stages" / stage_name(base) / (
        f"attempt_{attempt_index:03d}_v_{tag(target)}")
    state_rel = run_rel / "final_state.h5"
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
    config_path = bundle / f"adaptive_{physics_branch(base)}_attempt_{attempt_index:03d}.json"
    redirect_diagnostics(document, bundle / run_rel)
    return document, config_path, bundle / run_rel


def render_poisson_prewarm(base: dict, bundle: Path, target: float,
                           initial_state: Path, attempt_index: int,
                           max_iterations: int,
                           qf_initialization: str = "preserve_initial") -> tuple[dict, Path, Path]:
    """Render a Poisson-only contact-bias basin solve before coupled Newton."""
    document = json.loads(json.dumps(base))
    document["state_format"] = "hdf5"
    run_rel = Path("outputs") / "stages" / stage_name(base) / (
        f"prewarm_{attempt_index:03d}_v_{tag(target)}")
    document["output_csv"] = str(run_rel / "iv.csv")
    document["sweep"]["bias_points"] = [target]
    document["sweep"]["initial_state_file"] = str(initial_state.resolve())
    document["sweep"]["write_state_every_point_prefix"] = str(
        run_rel / "states" / "state")
    document["sweep"]["write_state_file"] = str(run_rel / "final_state.h5")
    solver = document["solver"]
    solver["method"] = "poisson_only"
    if qf_initialization not in ("preserve_initial", "contact_basin"):
        raise ValueError("Unknown Poisson prewarm quasi-Fermi initialization")
    solver["poisson_only_quasi_fermi_initialization"] = qf_initialization
    solver["max_iter"] = max_iterations
    if "handoff" in solver:
        solver["handoff"]["newton_max_iter"] = max_iterations
    config_path = bundle / f"adaptive_{physics_branch(base)}_prewarm_{attempt_index:03d}.json"
    redirect_diagnostics(document, bundle / run_rel)
    return document, config_path, bundle / run_rel


def render_cold_equilibrium(base: dict, bundle: Path) -> tuple[dict, Path, Path]:
    """Build a local zero-bias seed without consuming an imported solution."""
    require_direct_device_bias(base)
    # Keep the configured source and its derivatives. The local initialization
    # is a candidate; run() requalifies zero with the original coupled solver
    # and every active carrier row before allowing a positive-bias step.
    if any(float(c.get("bias", 0.0)) != 0.0 for c in base["contacts"]):
        raise ValueError("cold SLOT initialization requires zero terminal biases")
    document = json.loads(json.dumps(base))
    document["state_format"] = "hdf5"
    document.pop("state_file", None)
    document.pop("state_fields_dir", None)
    sweep = document["sweep"]
    sweep.pop("initial_state_file", None)
    run_dir = (bundle / "cold_equilibrium").resolve()
    document["output_csv"] = str(run_dir / "iv.csv")
    sweep.update(bias_points=[0.0], start=0.0, stop=0.0, max_retries=0,
                 write_state_file=str(run_dir / "final_state.h5"),
                 write_state_every_point_prefix=str(run_dir / "states/state"))
    document["solver"].update(method="poisson_only", max_iter=100,
        poisson_only_quasi_fermi_initialization="contact_basin")
    if "handoff" in document["solver"]:
        document["solver"]["handoff"]["newton_max_iter"] = 100
    redirect_diagnostics(document, run_dir)
    return document, bundle / "cold_equilibrium.json", run_dir


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
    run_rel = Path("outputs") / "stages" / stage_name(base) / (
        f"batch_predictor_v_{tag(biases[-1])}")
    document["output_csv"] = str(run_rel / "iv.csv")
    document["sweep"]["bias_points"] = biases
    document["sweep"]["write_state_every_point_prefix"] = str(
        run_rel / "states" / "state")
    document["sweep"]["write_state_file"] = str(run_rel / "final_state.h5")
    continuation = document["sweep"].setdefault("continuation", {})
    continuation["predictor"] = {
        "mode": "secant",
        "fields": ["psi", "phin", "phip"],
        "max_extrapolation_ratio": 2.0,
    }
    if "carrier_row_convergence" in document["solver"]:
        document["solver"]["carrier_row_convergence"]["diagnostic_csv"] = str(
            (bundle / run_rel / "carrier_row_convergence.csv").resolve())
    redirect_diagnostics(document, bundle / run_rel)
    return (
        document,
        bundle / f"adaptive_{physics_branch(base)}_batch_predictor.json",
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
        loader_status = exit_code & 0xFFFFFFFF
        loader_failures = {
            0xC0000135: "missing runtime DLL",
            0xC0000139: "runtime DLL entry point missing",
            0xC000007B: "incompatible executable or runtime DLL",
            0xC0000142: "runtime DLL initialization failed",
        }
        if loader_status in loader_failures:
            failure = {"failure_class": "runtime_environment", "exit_code": exit_code,
                       "windows_status": f"0x{loader_status:08X}",
                       "reason": loader_failures[loader_status], "runner": str(runner)}
            (run_dir / "startup_failure.json").write_text(
                json.dumps(failure, indent=2), encoding="utf-8")
            raise RuntimeError(
                f"Runner startup failed ({failure['windows_status']}: {failure['reason']}); "
                "check the matching toolchain bin directory in PATH. "
                "Voltage bisection cannot repair a runtime environment failure.")
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        def decoded(value: str | bytes | None) -> str:
            return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else value or ""
        (run_dir / "stdout.log").write_text(decoded(exc.stdout), encoding="utf-8")
        (run_dir / "stderr.log").write_text(decoded(exc.stderr), encoding="utf-8")
    return exit_code, timed_out, time.monotonic() - started


def enable_local_update_diagnostics(document: dict, run_dir: Path,
                                    nodes: list[int] | None) -> None:
    if not nodes:
        return
    solver = document["solver"]
    solver["diagnostic_history"] = True
    solver["local_update_diagnostics"] = {
        "enabled": True,
        "csv_file": str((run_dir / "local_updates.csv").resolve()),
        "nodes": nodes,
        "first_iterations": solver.get("max_iter", 10),
        "every_iterations": 1,
    }


def audit_local_candidate(document: dict, state: Path, bias: float,
                          runner: Path, repo: Path, run_dir: Path,
                          timeout_seconds: float, config_directory: Path) -> dict:
    """Require all physical rows to meet the deck's eps_row before continuation."""
    probe_dir = run_dir / "local_audit"
    probe_dir.mkdir(parents=True, exist_ok=True)
    probe = json.loads(json.dumps(document))
    probe["state_format"] = "hdf5"
    probe["simulation_type"] = "newton_carrier_term_probe"
    probe["state_file"] = str(state.resolve())
    probe["output_csv"] = str((probe_dir / "terms.csv").resolve())
    probe["carrier_term_probe"] = {"solved_equation_terms": True}
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        if key in probe:
            probe[key] = str((config_directory / probe[key]).resolve())
    driven = probe["sweep"]["contact"]
    probe.pop("sweep")
    for contact in probe["contacts"]:
        if contact["name"] == driven:
            contact["bias"] = bias
    config_path = probe_dir / "config.json"
    config_path.write_text(json.dumps(probe, indent=2), encoding="utf-8")
    code, timed_out, elapsed = execute_config(runner, repo, config_path, probe_dir, timeout_seconds)
    result = {"passed": False, "exit_code": code, "timed_out": timed_out, "elapsed_s": elapsed}
    if code == 0 and not timed_out:
        spec = importlib.util.spec_from_file_location(
            "slot_local_audit", Path(__file__).with_name("audit_slot_ldmos_local_continuity.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        gate = document.get("solver", {}).get("carrier_row_convergence", {})
        metrics = module.audit(probe_dir / "terms.csv", gate.get("eps_row", 1e-3),
                               gate.get("min_source_scale", 0.0), gate.get("scale_floor", 1e-300))
        result["metrics"] = metrics
        result["passed"] = all(c["rows_above_eps"] == 0 for c in metrics["carriers"].values())
    (probe_dir / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    return result


def restore_poly_candidate(document: dict, source: Path, output: Path,
                           config_directory: Path, bias: float) -> dict:
    """Prepare a guarded equilibrium restart; this never qualifies a state."""
    spec = importlib.util.spec_from_file_location(
        "slot_poly_equilibrium", Path(__file__).with_name("slot_poly_equilibrium.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.write_poly_equilibrium_restart(
        document, source, output, config_directory, bias)


def poisson_stagnation(rows: list, passes: int) -> bool:
    """Only shorten failed recovery; this never grants convergence."""
    if passes < 2 or len(rows) < passes:
        return False
    values = []
    for row in rows[-passes:]:
        if (row.get('failure_reason') != 'line_search_non_decrease' or
                row.get('carrier_row_violations') != '0' or
                row.get('global_continuity_closure_satisfied') != '1'):
            return False
        try:
            value = float(row['final_psi_residual_norm'])
        except (KeyError, ValueError, TypeError):
            return False
        if not math.isfinite(value) or value <= 0:
            return False
        values.append(value)
    return max(values)-min(values) <= 1e-4*max(values)


def correct_local_candidate(document: dict, state: Path, bias: float,
                            runner: Path, repo: Path, run_dir: Path,
                            timeout_seconds: float, config_directory: Path,
                            max_passes: int, mode: str = "carrier_only",
                            restore_poly_equilibrium: bool = False,
                            stagnation_passes: int = 0) -> dict:
    """Try bounded carrier Newton steps, then require the original coupled gates.

    The probe supplies candidates only. It cannot mark an I-V point accepted.
    Both the full coupled validator and the all-active row audit must pass.
    """
    if mode not in ("carrier_only", "carrier_density", "density_then_qf"):
        raise ValueError("Unknown carrier candidate mode")
    if mode == "density_then_qf" and max_passes < 2:
        raise ValueError("density_then_qf requires at least two candidate passes")
    root = run_dir / "local_corrector"
    root.mkdir(parents=True, exist_ok=True)
    base = json.loads(json.dumps(document))
    base["state_format"] = "hdf5"
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        if key in base:
            base[key] = str((config_directory / base[key]).resolve())
    result = {"passed": False, "timed_out": False, "mode": mode, "steps": [], "state_file": "", "validation_row": {},
              "coupled_newton_iterations": 0}
    current = state.resolve()
    validation_rows = []
    for iteration in range(1, max_passes + 1):
        step_mode = ("carrier_density" if iteration <= (max_passes + 1) // 2 else "carrier_only") if mode == "density_then_qf" else mode
        step_dir = root / f"step_{iteration:03d}"
        step_dir.mkdir(exist_ok=True)
        probe = json.loads(json.dumps(base))
        probe["simulation_type"] = "newton_block_step_probe"
        probe["state_file"] = str(current)
        probe["output_csv"] = str((step_dir / "step.csv").resolve())
        probe["output_state_prefix"] = str((step_dir / "trial").resolve())
        probe["block_modes"] = [step_mode]
        driven = probe.pop("sweep")["contact"]
        for contact in probe["contacts"]:
            if contact["name"] == driven:
                contact["bias"] = bias
        # Use the physical carrier Jacobian for the correction. The full
        # validator below retains the original deck's solver configuration.
        probe["solver"].setdefault("carrier_diagonal_floor", {})["enabled"] = False
        redirect_diagnostics(probe, step_dir)
        config = step_dir / "config.json"
        config.write_text(json.dumps(probe, indent=2), encoding="utf-8")
        started = time.monotonic()
        try:
            code, timed_out, elapsed = execute_config(runner, repo, config, step_dir, timeout_seconds)
        except Exception as exc:
            # Adapters may raise before returning an exit tuple. Persist the
            # attempted step without changing exception/retry semantics.
            result["steps"].append({"iteration": iteration, "mode": step_mode,
                "exit_code": getattr(exc, "exit_code", None),
                "elapsed_s": time.monotonic() - started, "local_audit_passed": False,
                "exception_type": type(exc).__name__, "failure": str(exc),
                "failure_status_file": getattr(exc, "status_file", None)})
            result["failure_reason"] = str(exc)
            (root / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
            raise
        trial = step_dir / f"trial_{step_mode}.h5"
        step = {"iteration": iteration, "mode": step_mode, "exit_code": code, "timed_out": timed_out,
                "elapsed_s": elapsed, "local_audit_passed": False}
        result["steps"].append(step)
        if timed_out:
            result.update(timed_out=True, timeout_stage="block_probe")
        if code == 0 and not timed_out and not trial.is_file():
            result["failure_reason"] = "missing_hdf5_candidate"
            (root / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
            raise RuntimeError(f"Successful block probe did not write HDF5 candidate: {trial}")
        if code != 0 or timed_out:
            break
        if restore_poly_equilibrium:
            projected = step_dir / "poly_equilibrium_candidate.h5"
            step["poly_equilibrium"] = restore_poly_candidate(
                base, trial, projected, config_directory, bias)
            trial = projected
        audit = audit_local_candidate(base, trial, bias, runner, repo, step_dir,
                                      timeout_seconds, config_directory)
        step["local_audit_passed"] = audit["passed"]
        if audit.get("timed_out", False):
            result.update(timed_out=True, timeout_stage="candidate_audit")
            break
        current = trial
        if audit["passed"]:
            validation_dir = root / ("validation" if iteration == 1 else f"validation_{iteration:03d}")
            validation_dir.mkdir(exist_ok=True)
            validation = json.loads(json.dumps(base))
            validation["output_csv"] = str((validation_dir / "iv.csv").resolve())
            validation["sweep"]["bias_points"] = [bias]
            validation["sweep"]["initial_state_file"] = str(current.resolve())
            validation["sweep"]["write_state_file"] = str((validation_dir / "final_state.h5").resolve())
            validation["sweep"]["write_state_every_point_prefix"] = str((validation_dir / "states/state").resolve())
            redirect_diagnostics(validation, validation_dir)
            config = validation_dir / "config.json"
            config.write_text(json.dumps(validation, indent=2), encoding="utf-8")
            code, timed_out, elapsed = execute_config(runner, repo, config, validation_dir, timeout_seconds)
            row = load_last_row(validation_dir / "iv.csv")
            candidate = validation_dir / "final_state.h5"
            result["validation_row"] = row
            validation_rows.append(row)
            result["validation_exit_code"] = code
            result["validation_timed_out"] = timed_out
            if timed_out:
                result.update(timed_out=True, timeout_stage="coupled_validation")
            result["validation_elapsed_s"] = elapsed
            result["coupled_newton_iterations"] += int(row.get("newton_iterations", "0"))
            step["validation_directory"] = str(validation_dir.resolve())
            step["validation_converged"] = code == 0 and not timed_out and row.get("converged") == "1"
            step["validation_newton_iterations"] = int(row.get("newton_iterations", "0"))
            if code == 0 and not timed_out and row.get("converged") == "1" and candidate.is_file():
                final_audit = audit_local_candidate(validation, candidate, bias, runner, repo,
                    validation_dir, timeout_seconds, config_directory)
                result["passed"] = final_audit["passed"]
                if final_audit.get("timed_out", False):
                    result.update(timed_out=True, timeout_stage="final_audit")
                    break
                if result["passed"]:
                    result["state_file"] = str(candidate.resolve())
                    break
            # A carrier correction can disturb Poisson; coupled revalidation
            # can in turn disturb carrier rows. Continue the same bounded
            # correction loop from its saved result, never accept either alone.
            if timed_out:
                break
            if poisson_stagnation(validation_rows, stagnation_passes):
                result['failure_reason'] = 'poisson_stagnation_with_qualified_carriers'
                result['stagnation_validation_rows'] = validation_rows[-stagnation_passes:]
                break
            if candidate.is_file() and row.get("converged") == "1":
                current = candidate
            else:
                rejected = load_rejected_final_state(validation_dir / "attempts.csv")
                if rejected is None:
                    break
                current = rejected
        (root / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (root / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def apply_qf_limit_overrides(document: dict, majority_V: float | None,
                             minority_V: float | None) -> None:
    solver = document["solver"]
    if majority_V is not None:
        solver["quasi_fermi_update_limit_V"] = majority_V
    if minority_V is not None:
        solver["quasi_fermi_update_limit_minority_V"] = minority_V


def apply_carrier_diagonal_floor_override(
        document: dict, scale: float | None,
        minority_density_ratio: float | None) -> None:
    if scale is None and minority_density_ratio is None:
        return
    floor = document["solver"].setdefault("carrier_diagonal_floor", {})
    floor["enabled"] = True
    if scale is not None:
        floor["scale"] = scale
    if minority_density_ratio is not None:
        floor["minority_density_ratio"] = minority_density_ratio


def apply_failure_state_policy(document: dict, return_best_feasible: bool) -> None:
    """Make the CLI policy authoritative even when the base deck enables it."""
    document["solver"]["return_best_feasible_on_failure"] = (
        return_best_feasible)


def run(args: argparse.Namespace) -> int:
    base_config = args.base_config.resolve()
    bundle = base_config.parent
    base = json.loads(base_config.read_text(encoding="utf-8"))
    branch = physics_branch(base)
    if getattr(args, "restore_poly_equilibrium", False) and not args.carrier_corrector_passes:
        raise ValueError("poly equilibrium restoration requires carrier corrector passes")
    require_direct_device_bias(base)
    base["state_format"] = "hdf5"
    runner = args.runner.resolve()
    cold_start = getattr(args, "cold_start", False)
    if cold_start:
        if args.batch_predictor or args.secant_previous_state is not None:
            raise ValueError("cold initialization requires individually audited local states")
        if any(t <= 0 for t in args.targets) or args.targets != sorted(set(args.targets)):
            raise ValueError("cold-start targets must be positive and strictly increasing")
        document, config, run_dir = render_cold_equilibrium(base, bundle)
        if run_dir.exists():
            raise ValueError("Refusing to overwrite a cold-equilibrium run")
        config.write_text(json.dumps(document, indent=2), encoding="utf-8")
        code, timed_out, _ = execute_config(runner, args.repo_root, config, run_dir, args.timeout_seconds)
        if code != 0 or timed_out or load_last_row(run_dir / "iv.csv").get("converged") != "1":
            raise RuntimeError("Local equilibrium Poisson initialization failed")
        initial_state = run_dir / "final_state.h5"
        # The Poisson seed is only a candidate. Zero must pass the same full
        # coupled and all-active carrier-row gates as every later point.
        args.initial_bias = 0.0
        args.targets = [0.0, *args.targets]
        args.audit_local_continuity = True
    else:
        initial_state = args.initial_state.resolve()
    if not initial_state.exists():
        raise FileNotFoundError(initial_state)
    states.read(initial_state, base, bundle)
    if args.audit_local_continuity and args.batch_predictor:
        raise ValueError("local qualification requires auditing each state before the next bias")
    if args.carrier_corrector_passes < 0:
        raise ValueError("carrier corrector passes must be nonnegative")
    if args.carrier_corrector_passes and not args.audit_local_continuity:
        raise ValueError("carrier correction requires --audit-local-continuity")

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
        enable_local_update_diagnostics(
            document, run_dir, args.local_update_node)
        apply_qf_limit_overrides(
            document, args.qf_update_limit_V,
            args.qf_update_limit_minority_V)
        apply_carrier_diagonal_floor_override(
            document, args.carrier_diagonal_floor_scale,
            args.carrier_diagonal_floor_minority_density_ratio)
        apply_failure_state_policy(
            document, args.return_best_feasible_on_failure)
        config_path.write_text(json.dumps(document, indent=2), encoding="utf-8")
        exit_code, timed_out, elapsed = execute_config(
            runner, args.repo_root, config_path, run_dir, args.timeout_seconds)
        rows = []
        iv_path = run_dir / "iv.csv"
        if iv_path.exists():
            with iv_path.open(newline="", encoding="utf-8-sig") as stream:
                rows = list(csv.DictReader(stream))
        summary = {
            "schema": f"vela.slot_ldmos.{branch}_batch_predictor.v1",
            "physics_branch": branch,
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
    # Loading a seed (possibly from another model) does not qualify a point in
    # this campaign. Keep its identity separately until a target passes.
    current_qualified = False
    previous_state = args.secant_previous_state.resolve() if args.secant_previous_state is not None else None
    previous_bias = args.secant_previous_bias
    attempts: list[dict[str, object]] = []
    attempt_index = 0
    while requested:
        target = requested[0]
        attempt_index += 1
        if attempt_index > args.max_attempts:
            break
        attempt_initial_state = current_state
        if previous_state is not None and (attempt_index == 1 or args.secant_predictor):
            if not previous_state.exists():
                raise FileNotFoundError(previous_state)
            predicted_state = bundle / "outputs" / "stages" / (
                stage_name(base)) / (
                f"secant_seed_{attempt_index:03d}_v_{tag(target)}.h5")
            write_secant_restart(
                previous_state, current_state, predicted_state,
                previous_bias, current_bias, target,
                args.secant_max_ratio, document=base, config_directory=bundle)
            attempt_initial_state = predicted_state
        prewarm_converged = False
        prewarm_elapsed = 0.0
        prewarm_failure = ""
        coupled_initial_state = attempt_initial_state
        if args.poisson_prewarm:
            prewarm_document, prewarm_config, prewarm_dir = render_poisson_prewarm(
                base, bundle, target, attempt_initial_state, attempt_index,
                args.poisson_prewarm_max_iterations,
                getattr(args, "poisson_prewarm_qf", "preserve_initial"))
            prewarm_config.write_text(
                json.dumps(prewarm_document, indent=2), encoding="utf-8")
            prewarm_exit, prewarm_timed_out, prewarm_elapsed = execute_config(
                runner, args.repo_root, prewarm_config, prewarm_dir,
                args.timeout_seconds)
            prewarm_row = load_last_row(prewarm_dir / "iv.csv")
            prewarm_state = prewarm_dir / "final_state.h5"
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
        if args.save_rejected_states or args.carrier_corrector_passes:
            document["solver"]["diagnostic_history"] = True
            document["sweep"].setdefault("diagnostics", {})["newton_history"] = {
                "enabled": True,
                "attempts_csv_file": str(run_dir.relative_to(bundle) / "attempts.csv"),
                "iterations_csv_file": str(run_dir.relative_to(bundle) / "iterations.csv"),
                "rejected_state_directory": str(
                    run_dir.relative_to(bundle) / "rejected_states"),
            }
        apply_failure_state_policy(
            document, args.return_best_feasible_on_failure)
        if args.closure_feasibility_restoration is not None:
            document["solver"].setdefault("global_continuity_closure", {})[
                "feasibility_restoration"] = args.closure_feasibility_restoration
        if args.closure_feasibility_restoration_factor is not None:
            document["solver"].setdefault("global_continuity_closure", {})[
                "feasibility_restoration_factor"] = (
                    args.closure_feasibility_restoration_factor)
        enable_local_update_diagnostics(
            document, run_dir, args.local_update_node)
        apply_qf_limit_overrides(
            document, args.qf_update_limit_V,
            args.qf_update_limit_minority_V)
        apply_carrier_diagonal_floor_override(
            document, args.carrier_diagonal_floor_scale,
            args.carrier_diagonal_floor_minority_density_ratio)
        config_path.write_text(json.dumps(document, indent=2), encoding="utf-8")
        exit_code, timed_out, elapsed = execute_config(
            runner, args.repo_root, config_path, run_dir, args.timeout_seconds)
        iv_row = load_last_row(run_dir / "iv.csv")
        next_state = run_dir / "final_state.h5"
        converged = exit_code == 0 and iv_row.get("converged") == "1" and next_state.exists()
        local_audit = None
        if converged and args.audit_local_continuity:
            local_audit = audit_local_candidate(
                document, next_state, target, runner, args.repo_root, run_dir,
                args.timeout_seconds, config_path.parent)
            converged = local_audit["passed"] and (
                target == current_bias or int(iv_row.get("newton_iterations", "0")) > 0)
        failure_reason = iv_row.get("failure_reason", "")
        if local_audit is not None and not converged:
            failure_reason = "local_continuity_audit" if not local_audit["passed"] else "zero_newton_new_voltage"
        timeout_stage = "coupled_solve" if timed_out else None
        if local_audit is not None and local_audit.get("timed_out", False):
            timed_out = True
            timeout_stage = "local_audit"
        rejected_final_state = load_rejected_final_state(run_dir / "attempts.csv")
        corrector = None
        coupled_iterations = int(iv_row.get("newton_iterations", "0"))
        if not converged and args.carrier_corrector_passes and not timed_out:
            correction_seed = next_state if next_state.is_file() else rejected_final_state
            if correction_seed is not None:
                corrector = correct_local_candidate(document, correction_seed, target,
                    runner, args.repo_root, run_dir, args.timeout_seconds,
                    config_path.parent, args.carrier_corrector_passes,
                    getattr(args, "carrier_corrector_mode", "carrier_only"),
                    getattr(args, "restore_poly_equilibrium", False))
                coupled_iterations += corrector["coupled_newton_iterations"]
                if corrector["passed"] and (target == current_bias or coupled_iterations > 0):
                    converged = True
                    next_state = Path(corrector["state_file"])
                    iv_row = corrector["validation_row"]
                    failure_reason = ""
                else:
                    failure_reason = "local_corrector_unqualified"
                if corrector.get("timed_out", False):
                    timed_out = True
                    timeout_stage = "corrector_" + corrector["timeout_stage"]
        if timed_out:
            failure_reason = "timeout"
        attempts.append({
            "attempt": attempt_index, "start_bias_V": current_bias,
            "target_bias_V": target, "elapsed_s": elapsed,
            "poisson_prewarm_enabled": args.poisson_prewarm,
            "poisson_prewarm_converged": prewarm_converged,
            "poisson_prewarm_elapsed_s": prewarm_elapsed,
            "poisson_prewarm_failure_reason": prewarm_failure,
            "timed_out": timed_out, "exit_code": exit_code,
            "timeout_stage": timeout_stage,
            "converged": converged, "failure_reason": failure_reason,
            "local_audit_passed": local_audit["passed"] if local_audit is not None else None,
            "local_corrector_passed": corrector["passed"] if corrector is not None else None,
            "local_corrector_steps": len(corrector["steps"]) if corrector is not None else 0,
            "coupled_newton_iterations": coupled_iterations,
            "iterations": iv_row.get("iterations", ""),
            "current_A_per_um": iv_row.get("current_total_A_per_um", ""),
            "state_file": str(next_state.resolve()) if converged else "",
            "rejected_final_state_file": (
                str(rejected_final_state) if rejected_final_state is not None else ""),
            "rejected_state_policy": (
                "best_feasible" if args.return_best_feasible_on_failure
                else "last_accepted_iterate"),
            "config_file": str(config_path.resolve()),
        })
        stop_after_attempt = False
        if converged:
            if target != current_bias:
                previous_state, previous_bias = current_state, current_bias
            current_bias = target
            current_state = next_state
            current_qualified = True
            requested.popleft()
        # A wall-clock budget exhaustion supplies no numerical convergence
        # result. Keep the target pending instead of inferring a smaller step.
        elif timed_out:
            stop_after_attempt = True
        else:
            midpoint = 0.5 * (current_bias + target)
            if target - current_bias <= args.min_step_V or midpoint <= current_bias:
                stop_after_attempt = True
            else:
                requested.appendleft(midpoint)
        write_adaptive_checkpoint(bundle, {
            "schema": "vela.slot_ldmos.adaptive_checkpoint.v1",
            "completed": not requested,
            "base_config_sha256": hashlib.sha256(base_config.read_bytes()).hexdigest(),
            "last_successful_bias_V": current_bias if current_qualified else None,
            "last_successful_state": str(current_state) if current_qualified else None,
            "state_sha256": hashlib.sha256(current_state.read_bytes()).hexdigest() if current_qualified else None,
            "previous_bias_V": previous_bias if previous_state is not None else None,
            "previous_state": str(previous_state) if previous_state is not None else None,
            "previous_state_sha256": hashlib.sha256(previous_state.read_bytes()).hexdigest() if previous_state is not None else None,
            "remaining_targets_V": list(requested),
            "last_finished_attempt": attempts[-1],
        })
        if stop_after_attempt:
            break

    report_dir = bundle / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    csv_path = report_dir / f"adaptive_{branch}_attempts.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(attempts[0]) if attempts else ["attempt"])
        writer.writeheader()
        writer.writerows(attempts)
    completed_all = not requested
    summary = {
        "schema": f"vela.slot_ldmos.{branch}_adaptive.v1",
        "physics_branch": branch,
        "impact_ionization": base.get("solver", {}).get("impact_ionization", "none"),
        "completed": completed_all,
        "termination_reason": "completed" if completed_all else (
            "timeout" if attempts and attempts[-1]["timed_out"] else "unqualified_or_attempt_budget"),
        "initial_bias_V": args.initial_bias,
        "initial_state": str(initial_state),
        "last_successful_bias_V": current_bias if current_qualified else None,
        "last_successful_state": str(current_state) if current_qualified else None,
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
        "poisson_prewarm_qf": getattr(args, "poisson_prewarm_qf", "preserve_initial"),
        "cold_start": cold_start,
        "secant_predictor": args.secant_predictor,
        "audit_local_continuity": args.audit_local_continuity,
        "carrier_corrector_passes": args.carrier_corrector_passes,
        "carrier_corrector_mode": getattr(args, "carrier_corrector_mode", "carrier_only"),
        "restore_poly_equilibrium": getattr(args, "restore_poly_equilibrium", False),
        "poisson_prewarm_max_iterations": args.poisson_prewarm_max_iterations,
        "min_step_V": args.min_step_V,
        "attempts": attempts,
    }
    (report_dir / f"adaptive_{branch}_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({key: summary[key] for key in (
        "completed", "last_successful_bias_V", "remaining_targets_V")}, sort_keys=True))
    return 0 if completed_all else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    seed = parser.add_mutually_exclusive_group(required=True)
    seed.add_argument("--initial-state", type=Path)
    seed.add_argument("--cold-start", action="store_true",
                      help="derive local equilibrium, then qualify 0 V before the requested positive biases")
    parser.add_argument("--initial-bias", type=float, default=0.105)
    parser.add_argument("--targets", type=float, nargs="+", default=DEFAULT_TARGETS)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--max-newton-iterations", type=int, default=8)
    parser.add_argument("--stall-residual-floor", type=float)
    parser.add_argument(
        "--qf-update-limit-mode",
        choices=("componentwise_poisson_recorrect", "uniform_trust_region"))
    parser.add_argument("--qf-update-limit-V", type=float)
    parser.add_argument("--qf-update-limit-minority-V", type=float)
    parser.add_argument("--carrier-diagonal-floor-scale", type=float)
    parser.add_argument(
        "--carrier-diagonal-floor-minority-density-ratio", type=float)
    parser.add_argument(
        "--line-search-mode",
        choices=("merit", "block_filter", "block_closure_filter"))
    parser.add_argument(
        "--closure-feasibility-restoration",
        choices=("off", "linear_qf_projection", "final_linear_qf_projection"))
    parser.add_argument(
        "--closure-feasibility-restoration-factor", type=float)
    parser.add_argument(
        "--solver-method", choices=("newton", "gummel_newton"), default="newton")
    parser.add_argument("--gummel-max-iterations", type=int)
    parser.add_argument("--poisson-prewarm", action="store_true")
    parser.add_argument("--poisson-prewarm-qf", choices=("preserve_initial", "contact_basin"),
                        default="preserve_initial")
    parser.add_argument("--poisson-prewarm-max-iterations", type=int, default=20)
    parser.add_argument("--min-step-V", type=float, default=0.01)
    parser.add_argument("--max-attempts", type=int, default=100)
    parser.add_argument("--batch-predictor", action="store_true")
    parser.add_argument("--audit-local-continuity", action="store_true",
                        help="audit every active carrier row before accepting a continuation checkpoint")
    parser.add_argument("--carrier-corrector-passes", type=int, default=0,
                        help="bounded carrier Newton corrections followed by full coupled validation; requires local audit")
    parser.add_argument("--carrier-corrector-mode", choices=("carrier_only", "carrier_density", "density_then_qf"),
                        default="carrier_only", help="coordinate choice for the candidate correction, preserving full acceptance gates")
    parser.add_argument("--restore-poly-equilibrium", action="store_true",
                        help="restore guarded source-free Poly equilibrium in corrector candidates; retain full validation")
    parser.add_argument("--secant-previous-state", type=Path)
    parser.add_argument("--secant-predictor", action="store_true",
                        help="use the two most recent qualified states for every subsequent bias prediction")
    parser.add_argument("--secant-previous-bias", type=float, default=0.0)
    parser.add_argument("--secant-max-ratio", type=float, default=2.0)
    parser.add_argument("--save-rejected-states", action="store_true")
    parser.add_argument("--return-best-feasible-on-failure", action="store_true")
    parser.add_argument(
        "--local-update-node", type=int, action="append", default=[],
        help="Node id to include in opt-in raw/capped/applied Newton diagnostics; repeatable")
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
