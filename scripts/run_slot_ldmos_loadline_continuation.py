#!/usr/bin/env python3
"""Bounded SLOT independent-device load-line continuation with audited checkpoints.

Only qualified direct device evaluations may seed another solve. External points
also require a closed load line and byte-identical circuit confirmation. Changed
confirmation states remain unqualified; no placeholder circuit row proves the
device equations. This driver does not establish original AvalDensGradQF parity.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time


def _module(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


direct = _module("run_slot_ldmos_noavalanche_continuation")
audit = _module("audit_slot_ldmos_local_continuity")
SCHEMA = "vela.slot_ldmos.loadline_continuation.v1"
RESISTANCE = 1e12  # ohm*um; current is signed A/um throughout
TOLERANCE = 1e-6  # V


class EvidenceError(ValueError):
    """Missing/corrupt/nonfinite proof: stop, never voltage-bisect."""


class NumericalFailure(RuntimeError):
    """A completed candidate did not qualify; bounded shrink is allowed."""


class BudgetExhausted(RuntimeError):
    pass


def finite(value):
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise EvidenceError(f"invalid numeric evidence: {value!r}") from exc
    if not math.isfinite(number):
        raise EvidenceError("nonfinite evidence")
    return number


def sha256(path):
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError as exc:
        raise EvidenceError(f"missing evidence: {path}") from exc


def hashes(paths):
    return {str(Path(path).resolve()): sha256(path) for path in paths}


def verify_hashes(files):
    for path, expected in files.items():
        if sha256(path) != expected:
            raise EvidenceError(f"hash mismatch: {path}")


def clone(value):
    return json.loads(json.dumps(value, allow_nan=False))


def validate_profile(base):
    direct.require_direct_device_bias(base)
    solver = base["solver"]
    if solver.get("method") != "newton":
        raise ValueError("requires explicit full Newton device solver")
    checks = ((solver.get("carrier_row_convergence", {}),
               {"mode": "enforce", "eps_row": .001, "scale_floor": 1e-30}),
              (solver.get("global_continuity_closure", {}),
               {"mode": "enforce", "source_floor": 1e-11,
                "tolerance": .01, "absolute_tolerance": 1e-12}),
              (base["sweep"].get("external_circuit", {}),
               {"enabled": False, "mode": "series_resistor",
                "resistance_ohm_um": RESISTANCE, "residual_tolerance_V": TOLERANCE}))
    for actual, required in checks:
        if any(actual.get(key) != value for key, value in required.items()):
            raise ValueError(f"requires original SLOT qualification profile: {required}")
    if base["sweep"]["external_circuit"].get("current_direction") not in (-1, 1):
        raise ValueError("explicit current_direction must be +1 or -1")


def freeze_inputs(base, directory, runner):
    """Resolve physical resources before any relocation, including nested IALMob."""
    document = clone(base)
    resources = [Path(runner).resolve()]
    def resolve(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in {"mesh_file", "materials_file", "node_doping_file", "geometry_file"}:
                    path = (Path(directory) / item).resolve()
                    value[key] = str(path)
                    resources.append(path)
                else:
                    resolve(item)
        elif isinstance(value, list):
            for item in value:
                resolve(item)
    resolve(document)
    encoded = json.dumps(document, sort_keys=True, allow_nan=False).encode()
    return document, {"config_sha256": hashlib.sha256(encoded).hexdigest(),
                      "files": hashes(resources)}


def residual(inner, current, outer, direction):
    if direction not in (-1, 1):
        raise EvidenceError("invalid current direction")
    return finite(finite(inner) + direction * RESISTANCE * finite(current) - finite(outer))


def halfway(seed, failed, minimum):
    middle = seed + (failed - seed) * .5
    return middle if abs(middle - seed) >= minimum and middle not in (seed, failed) else None


def proposal(samples, seed, maximum, exploration):
    """Safeguard true brackets; same-sign samples are only secant estimates."""
    pairs = [(finite(x), finite(f)) for x, f in samples]
    brackets = [(a, b) for i, a in enumerate(pairs) for b in pairs[i + 1:]
                if a[0] != b[0] and ((a[1] < 0 < b[1]) or (b[1] < 0 < a[1]))]
    bracket = min(brackets, key=lambda ab: abs(ab[0][0] - ab[1][0])) if brackets else None
    pair = bracket or (pairs[-2:] if len(pairs) >= 2 else None)
    candidate = None
    if pair:
        (x0, f0), (x1, f1) = pair
        denominator = f1 - f0
        if math.isfinite(denominator) and abs(denominator) > 1e-12 * max(abs(f0), abs(f1), 1e-300):
            value = x1 - f1 * (x1 - x0) / denominator
            if math.isfinite(value):
                candidate = value
    if bracket:
        lo, hi = sorted([bracket[0][0], bracket[1][0]])
        if candidate is None or not lo + .05 * (hi - lo) < candidate < hi - .05 * (hi - lo):
            candidate = lo + .5 * (hi - lo)
    elif candidate is None:
        # Alternate bounded exploration; do not assume monotonic I(V).
        sign = -1 if pairs[-1][1] > 0 else 1
        candidate = seed + sign * maximum * (1 if exploration % 2 == 0 else -1)
    candidate = min(seed + maximum, max(seed - maximum, candidate))
    inside_bracket = bracket is not None and min(bracket[0][0], bracket[1][0]) < candidate < max(bracket[0][0], bracket[1][0])
    return candidate, inside_bracket


def read_rows(path):
    try:
        with Path(path).open(newline="", encoding="utf-8-sig") as stream:
            rows = list(csv.DictReader(stream))
        if not rows:
            raise EvidenceError(f"empty evidence: {path}")
        return rows
    except OSError as exc:
        raise EvidenceError(f"missing evidence: {path}") from exc


def node_ids(rows, expected):
    try:
        ids = [int(row["node_id"]) for row in rows]
    except (KeyError, ValueError) as exc:
        raise EvidenceError("invalid node IDs") from exc
    if len(set(ids)) != len(ids) or set(ids) != set(expected):
        raise EvidenceError("missing, duplicate or unexpected mesh node IDs")
    return ids


def qualify_files(state, iv, terms, bias, base, mesh_ids, *, expected_count=12634,
                  expected_active=None, updates=0, same_bias=False):
    """Recompute row/closure checks; only the full solver proves final device gates."""
    row = read_rows(iv)[-1]
    try:
        if row["solver_method"] != "newton" or "boundary_seed" in row["newton_convergence_reason"]:
            raise EvidenceError("placeholder or non-full-solver device proof")
        if not math.isclose(finite(row["bias_V"]), bias, rel_tol=0, abs_tol=1e-12):
            raise EvidenceError("actual device bias mismatch")
        sweep = base["sweep"]
        if row["bias_contact"] != sweep["contact"] or row["current_contact"] != sweep.get("current_contact", sweep["contact"]):
            raise EvidenceError("terminal identity mismatch")
        current = finite(row["current_total_A_per_um"])
        for field in ("final_psi_residual_norm", "final_electron_continuity_residual_norm", "final_hole_continuity_residual_norm"):
            if finite(row[field]) < 0:
                raise EvidenceError("negative residual norm")
        if row["converged"] != "1" or row["global_continuity_closure_satisfied"] != "1":
            raise NumericalFailure("device_numerical: full final gates failed")
        if updates <= 0 and not same_bias:
            raise NumericalFailure("zero_coupled_update: new voltage has no coupled Newton update")
        for carrier in ("electron", "hole"):
            prefix = f"global_{carrier}_"
            flux = finite(row[prefix + "contact_flux"])
            source = finite(row[prefix + "integrated_source"])
            mismatch = finite(row[prefix + "continuity_closure_mismatch"])
            scale = max(abs(flux), abs(source), 1e-11)
            if not math.isclose(mismatch, flux - source, rel_tol=1e-8, abs_tol=scale * 1e-12):
                raise EvidenceError("closure mismatch is inconsistent with flux/source")
            ratio = abs(mismatch) / scale
            if not math.isclose(finite(row[prefix + "continuity_closure_ratio"]), ratio, rel_tol=1e-8, abs_tol=1e-12):
                raise EvidenceError("closure ratio inconsistent with physical evidence")
            relative = max(abs(flux), abs(source)) >= 1e-11
            if row[prefix + "continuity_closure_criterion"] != ("relative" if relative else "absolute"):
                raise EvidenceError("closure criterion inconsistent with physical evidence")
            if (relative and ratio > .01) or (not relative and abs(mismatch) > 1e-12):
                raise NumericalFailure("global_closure: actual carrier mismatch exceeds gate")
        states = read_rows(state)
        node_ids(states, mesh_ids)
        for value in states:
            for key in ("psi", "phin", "phip", "electrons_m3", "holes_m3"):
                finite(value[key])
            for key, number in value.items():
                if key.endswith("_qf_reference_V") or key.endswith("_qf_increment_V"):
                    finite(number)
        rows = read_rows(terms)
        node_ids(rows, mesh_ids)
        active = {carrier: sorted(int(r["node_id"]) for r in rows if r[f"{carrier}_continuity_active"] == "1")
                  for carrier in ("electron", "hole")}
        if any(len(ids) != expected_count for ids in active.values()) or (expected_active is not None and active != expected_active):
            raise EvidenceError("all-active carrier row identity/count mismatch")
        metrics = audit.audit(Path(terms), .001, 0., 1e-30)
        if any(c["rows_above_eps"] for c in metrics["carriers"].values()):
            raise NumericalFailure("local_audit: at least one all-active carrier row fails")
    except (KeyError, ValueError) as exc:
        if isinstance(exc, EvidenceError):
            raise
        raise EvidenceError(f"malformed qualification evidence: {exc}") from exc
    return {"inner_V": bias, "current_A_per_um": current, "state_file": str(Path(state).resolve()),
            "state_sha256": sha256(state), "active_ids": active,
            "coupled_newton_updates": updates, "same_bias_requalification": same_bias,
            "iv_file": str(Path(iv).resolve()), "terms_file": str(Path(terms).resolve()),
            "evidence": hashes([state, iv, terms])}


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def solver_environment():
    return {key: os.environ.get(key) for key in (
        "VELA_LINEAR_SOLVER", "VELA_LINEAR_THREADS", "VELA_BLAS_THREADS",
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}


class Campaign:
    """Two durable ledgers plus finite, bidirectional scalar/device recovery."""
    def __init__(self, output, identity, targets, initial_state, initial_inner,
                 direction, backend, *, resume=False, max_attempts=100,
                 max_shots=12, max_inner_retries=8, max_inner_step=.05,
                 min_inner_step=1e-9, min_outer_step=1e-6,
                 max_total_seconds=7200., initial_outer=None):
        self.output = Path(output).resolve()
        self.identity = {**identity, "solver_environment": solver_environment()}
        self.backend = backend
        self.direction = direction
        self.started = time.monotonic()
        self.budgets = dict(max_attempts=max_attempts, max_shots=max_shots,
            max_inner_retries=max_inner_retries, max_inner_step=max_inner_step,
            min_inner_step=min_inner_step, min_outer_step=min_outer_step,
            max_total_seconds=max_total_seconds)
        if any(not math.isfinite(float(v)) or v <= 0 for v in self.budgets.values()):
            raise ValueError("all continuation budgets must be finite and positive")
        if max_inner_step < min_inner_step or not targets:
            raise ValueError("invalid continuation step bounds or empty targets")
        self.request = {"targets": [finite(t) for t in targets], "initial_inner_V": finite(initial_inner),
                        "initial_state": str(Path(initial_state).resolve()), "direction": direction,
                        "initial_outer_V": initial_outer}
        verify_hashes(identity["files"])
        path = self.output / "checkpoint.json"
        if resume:
            try:
                self.snapshot = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                raise EvidenceError("missing/corrupt checkpoint") from exc
            self._validate_resume()
        else:
            self.output.mkdir(parents=True, exist_ok=True)
            if path.exists() or (self.output / "attempts").exists():
                raise ValueError("existing campaign requires --resume or a fresh output directory")
            self.snapshot = {"schema": SCHEMA, "identity": self.identity, "request": self.request,
                "budgets": self.budgets, "attempt_count": 0, "target_shots": 0,
                "elapsed_s": 0., "inner_ledger": [], "outer_ledger": [],
                "pending_targets_V": list(targets), "events": [], "completed": False,
                "failure": "", "active_ids": None}
        self.previous_elapsed = finite(self.snapshot["elapsed_s"])
        self.output.joinpath("attempts").mkdir(exist_ok=True)
        if hasattr(backend, "set_deadline"):
            backend.set_deadline(self.started + max_total_seconds - self.previous_elapsed)
        self.save()

    def _validate_resume(self):
        s = self.snapshot
        if s.get("schema") != SCHEMA or s.get("identity") != self.identity or s.get("request") != self.request or s.get("budgets") != self.budgets:
            raise EvidenceError("checkpoint schema/config/runner/resources/request/budgets identity changed")
        for index, point in enumerate(s["inner_ledger"]):
            if point["id"] != index or point["active_ids"] != s["active_ids"]:
                raise EvidenceError("qualified ledger identity mismatch")
            if point["same_bias_requalification"] != (index == 0) or (index == 0 and point["inner_V"] != self.request["initial_inner_V"]):
                raise EvidenceError("same-bias import exemption does not match initial checkpoint")
            verify_hashes(point["evidence"])
            self.backend.validate(point)
        for outer in s["outer_ledger"]:
            try:
                point = s["inner_ledger"][outer["device_id"]]
            except (KeyError, IndexError) as exc:
                raise EvidenceError("outer ledger references missing device") from exc
            actual = residual(point["inner_V"], point["current_A_per_um"], outer["outer_V"], self.direction)
            if abs(actual) > TOLERANCE or actual != outer["residual_V"]:
                raise EvidenceError("stored load-line residual does not match device proof")
            confirmation = outer["confirmation"]
            verify_hashes(confirmation["evidence"])
            if confirmation["state_sha256"] != point["state_sha256"]:
                raise EvidenceError("confirmation state identity mismatch")
            if hasattr(self.backend, "validate_confirmation"):
                self.backend.validate_confirmation(point, outer["outer_V"], confirmation)
        accepted = [p["outer_V"] for p in s["outer_ledger"]]
        pending = s["pending_targets_V"]
        if any(t not in accepted and t not in pending for t in self.request["targets"]):
            raise EvidenceError("checkpoint discarded requested target")
        if s["completed"] != (not pending):
            raise EvidenceError("checkpoint completion disagrees with pending targets")

    def save(self):
        self.snapshot["elapsed_s"] = self.previous_elapsed + time.monotonic() - self.started
        self.snapshot["completed"] = not self.snapshot["pending_targets_V"]
        atomic_json(self.output / "checkpoint.json", self.snapshot)

    def allocate(self, kind, **details):
        if self.identity["solver_environment"] != solver_environment():
            raise EvidenceError("solver environment changed during campaign")
        verify_hashes(self.identity["files"])
        if self.snapshot["attempt_count"] >= self.budgets["max_attempts"] or self.previous_elapsed + time.monotonic() - self.started >= self.budgets["max_total_seconds"]:
            raise BudgetExhausted("budget_exhausted: total attempts or elapsed time")
        index = self.snapshot["attempt_count"] + 1
        directory = self.output / "attempts" / f"attempt_{index:06d}"
        while directory.exists():  # Orphan created before checkpoint replacement.
            index += 1
            directory = self.output / "attempts" / f"attempt_{index:06d}"
        if index > self.budgets["max_attempts"]:
            raise BudgetExhausted("budget_exhausted: orphan attempts consume remaining budget")
        self.snapshot["attempt_count"] = index
        directory.mkdir()
        atomic_json(directory / "attempt.json", {"kind": kind, **details})
        self.save()  # Interrupted attempts consume budget and cannot be reused.
        return directory

    def evaluate(self, bias, seed, *, same_bias=False):
        directory = self.allocate("device", inner_V=bias, seed=seed.get("id"), same_bias=same_bias)
        try:
            point = self.backend.solve(bias, seed, directory, same_bias=same_bias)
            if self.snapshot["active_ids"] is not None and point["active_ids"] != self.snapshot["active_ids"]:
                raise EvidenceError("all-active row IDs changed from qualified baseline")
        except NumericalFailure as exc:
            atomic_json(directory / "result.json", {"qualified": False, "failure": str(exc)})
            self.snapshot["events"].append({"event": "device_failed", "inner_V": bias, "failure": str(exc), "attempt": str(directory)})
            self.save()
            raise

        point["id"] = len(self.snapshot["inner_ledger"])
        atomic_json(directory / "result.json", {"qualified": True, "device_id": point["id"]})
        self.snapshot["inner_ledger"].append(point)
        self.snapshot["active_ids"] = point["active_ids"]
        self.save()
        return point

    def accept_outer(self, point, target):
        value = residual(point["inner_V"], point["current_A_per_um"], target, self.direction)
        if abs(value) > TOLERANCE:
            raise NumericalFailure("loadline_unclosed")
        directory = self.allocate("confirmation", outer_V=target, device_id=point["id"])
        try:
            confirmation = self.backend.confirm(point, target, directory)
        except NumericalFailure as exc:
            atomic_json(directory / "result.json", {"qualified": False, "failure": str(exc)})
            raise
        if confirmation["state_sha256"] != point["state_sha256"]:
            raise EvidenceError("confirmation inherited proof for a different state")
        checkpoint = {"outer_V": target, "device_id": point["id"], "residual_V": value, "confirmation": confirmation}
        self.snapshot["outer_ledger"].append(checkpoint)
        if self.snapshot["pending_targets_V"] and self.snapshot["pending_targets_V"][0] == target:
            self.snapshot["pending_targets_V"].pop(0)
            self.snapshot["target_shots"] = 0
        self.snapshot["inner_recovery"] = None
        atomic_json(directory / "result.json", {"qualified": True, "outer_V": target})
        self.save()

    def solve_outer(self, target):
        saved = self.snapshot.get("inner_recovery")
        recovery = None
        retries = 0
        if saved and saved["outer_V"] == target:
            if not any(p["inner_V"] == saved["trial_V"] for p in self.snapshot["inner_ledger"]):
                recovery = (saved["trial_V"], self.snapshot["inner_ledger"][saved["seed_id"]])
                retries = saved["retries"]
        while True:
            ledger = self.snapshot["inner_ledger"]
            samples = [(p["inner_V"], residual(p["inner_V"], p["current_A_per_um"], target, self.direction)) for p in ledger]
            best_index = min(range(len(samples)), key=lambda i: abs(samples[i][1]))
            seed = ledger[best_index]
            if abs(samples[best_index][1]) <= TOLERANCE:
                self.accept_outer(seed, target)
                return True
            if self.snapshot["target_shots"] >= self.budgets["max_shots"]:
                return False
            if recovery is None:
                trial, bracket = proposal(samples, seed["inner_V"], self.budgets["max_inner_step"], self.snapshot["target_shots"])
                self.snapshot["last_scalar"] = {"target_outer_V": target, "seed_id": seed["id"], "true_bracket": bracket, "sample_ids": [p["id"] for p in ledger]}
                if any(abs(trial - p["inner_V"]) < self.budgets["min_inner_step"] for p in ledger):
                    # Bounded alternative exploration for a repeated secant estimate.
                    sign = 1 if self.snapshot["target_shots"] % 2 else -1
                    trial = seed["inner_V"] + sign * self.budgets["max_inner_step"] / (self.snapshot["target_shots"] + 1)
            else:
                trial, seed = recovery
            if abs(trial - seed["inner_V"]) < self.budgets["min_inner_step"]:
                return False
            self.snapshot["target_shots"] += 1
            self.snapshot["inner_recovery"] = {"outer_V": target, "trial_V": trial, "seed_id": seed["id"], "retries": retries}
            self.save()
            try:
                self.evaluate(trial, seed)
                recovery = None
                retries = 0
                self.snapshot["inner_recovery"] = None
                self.save()
            except NumericalFailure:
                smaller = halfway(seed["inner_V"], trial, self.budgets["min_inner_step"])
                retries += 1
                if smaller is None or retries > self.budgets["max_inner_retries"]:
                    return False
                recovery = (smaller, seed)
                self.snapshot["inner_recovery"] = {"outer_V": target, "trial_V": smaller, "seed_id": seed["id"], "retries": retries}
                self.save()

    def run(self):
        try:
            if not self.snapshot["inner_ledger"]:
                seed = {"inner_V": self.request["initial_inner_V"], "state_file": self.request["initial_state"]}
                self.evaluate(seed["inner_V"], seed, same_bias=True)
            if not self.snapshot["outer_ledger"]:
                point = self.snapshot["inner_ledger"][0]
                computed_outer = point["inner_V"] + self.direction * RESISTANCE * point["current_A_per_um"]
                outer = self.request["initial_outer_V"]
                self.accept_outer(point, computed_outer if outer is None else outer)
            while self.snapshot["pending_targets_V"]:
                target = self.snapshot["pending_targets_V"][0]
                if self.solve_outer(target):
                    continue
                last_outer = self.snapshot["outer_ledger"][-1]["outer_V"]
                smaller = halfway(last_outer, target, self.budgets["min_outer_step"])
                if smaller is None:
                    raise NumericalFailure("loadline_unclosed: minimum outer step reached")
                self.snapshot["pending_targets_V"].insert(0, smaller)
                self.snapshot["target_shots"] = 0
                self.snapshot["inner_recovery"] = None
                self.snapshot["events"].append({"event": "outer_step_halved", "failed_outer_V": target, "queued_outer_V": smaller})
                self.save()
            self.snapshot["completed"] = True
            self.snapshot["failure"] = ""
            self.save()
            return True
        except (BudgetExhausted, NumericalFailure) as exc:
            self.snapshot["failure"] = str(exc)
            self.save()
            return False
        except (Exception, KeyboardInterrupt) as exc:
            self.snapshot["failure"] = f"{type(exc).__name__}: {exc}"
            self.save()
            raise


class SolverBackend:
    """Adapt existing candidate helpers without changing their solver semantics."""
    def __init__(self, base, repo, runner, *, timeout_seconds=240.,
                 max_newton_iterations=None, poisson_prewarm=False,
                 poisson_prewarm_max_iterations=60, carrier_corrector_passes=0,
                 carrier_corrector_mode="carrier_only", restore_poly_equilibrium=False,
                 expected_node_count=13881, expected_active_count=12634):
        validate_profile(base)
        self.base = clone(base)
        self.repo, self.runner = Path(repo).resolve(), Path(runner).resolve()
        self.timeout = timeout_seconds
        self.max_newton_iterations = max_newton_iterations
        self.prewarm = poisson_prewarm
        self.prewarm_iterations = poisson_prewarm_max_iterations
        self.corrector_passes = carrier_corrector_passes
        self.corrector_mode = carrier_corrector_mode
        self.restore_poly = restore_poly_equilibrium
        self.expected_active_count = expected_active_count
        mesh = json.loads(Path(base["mesh_file"]).read_text(encoding="utf-8"))
        ids = [node["id"] for node in mesh["nodes"]]
        if any(type(i) is not int for i in ids) or len(ids) != expected_node_count or len(set(ids)) != len(ids):
            raise ValueError("unexpected SLOT mesh node count/identity")
        self.mesh_ids = set(ids)
        self.helper = _module("run_slot_ldmos_noavalanche_continuation")
        self.raw_execute = self.helper.execute_config
        # This private module instance also routes nested corrector/audit calls
        # through one deadline and process-failure classifier.
        self.helper.execute_config = self._execute
        self.deadline = math.inf
        self.lineage = []

    def set_deadline(self, deadline):
        self.deadline = deadline

    def _execute(self, runner, repo, config, run_dir, timeout):
        remaining = self.deadline - time.monotonic()
        if remaining <= 0:
            raise BudgetExhausted("budget_exhausted: total wall time")
        result = self.raw_execute(runner, repo, config, run_dir, min(timeout, remaining))
        code, timed_out, _ = result
        if timed_out:
            raise NumericalFailure("timeout: device or probe process exceeded time budget")
        document = json.loads(Path(config).read_text(encoding="utf-8"))
        output = (Path(config).parent / document["output_csv"]).resolve()
        if code != 0 and not output.is_file():
            raise RuntimeError(f"runtime_environment_or_invalid_input: process exit {code} without solver output; see {run_dir}")
        if (document.get("simulation_type", "dc_sweep") in ("dc_sweep", "newton")
                and document["solver"].get("method") == "newton"
                and not document["sweep"]["external_circuit"]["enabled"] and output.is_file()):
            row = read_rows(output)[-1]
            if row.get("solver_method") != "newton":
                raise EvidenceError("direct solve returned placeholder/non-Newton proof")
            self.lineage.append({"iv_file": str(output), "config_file": str(Path(config).resolve()),
                                 "updates": int(row["newton_iterations"])})
        return result

    def _write_run(self, rendered):
        document, config, directory = rendered
        # Existing helpers intentionally use verbose stage paths. Relocate the
        # adapter's run before execution so nested Poly temporary names remain
        # below ordinary Windows path limits for typical workspace roots.
        short = "p" if document["solver"].get("method") == "poisson_only" else "d"
        directory = config.parent / short
        document["output_csv"] = str((directory / "iv.csv").resolve())
        document["sweep"]["write_state_file"] = str((directory / "final_state.csv").resolve())
        document["sweep"]["write_state_every_point_prefix"] = str((directory / "states/s").resolve())
        self.helper.redirect_diagnostics(document, directory)
        config.parent.mkdir(parents=True, exist_ok=True)
        config.write_text(json.dumps(document, indent=2, allow_nan=False), encoding="utf-8")
        code, timed_out, _ = self._execute(self.runner, self.repo, config, directory, self.timeout)
        if code != 0 or timed_out:
            return document, config, directory, False
        row = self.helper.load_last_row(directory / "iv.csv")
        return document, config, directory, row.get("converged") == "1"

    def _qualified(self, document, config, directory, state, bias, same_bias, seed):
        audit_result = self.helper.audit_local_candidate(document, state, bias,
            self.runner, self.repo, directory, self.timeout, config.parent)
        if not audit_result.get("passed"):
            raise NumericalFailure("local_audit: final all-active rows fail")
        terms = directory / "local_audit/terms.csv"
        audit_config = directory / "local_audit/config.json"
        self._check_audit_binding(audit_config, state, bias)
        updates = self._lineage_updates(self.lineage, bias)
        point = qualify_files(state, directory / "iv.csv", terms, bias, self.base, self.mesh_ids,
            expected_count=self.expected_active_count, expected_active=seed.get("active_ids"),
            updates=updates, same_bias=same_bias)
        proof = {"state_sha256": point["state_sha256"], "inner_V": bias,
                 "current_A_per_um": point["current_A_per_um"],
                 "audit_config": str(audit_config.resolve()), "lineage": clone(self.lineage)}
        proof_path = directory / "qualification.json"
        atomic_json(proof_path, proof)
        point["proof_file"] = str(proof_path.resolve())
        files = [config, audit_config, proof_path, Path(seed["state_file"])]
        for step in self.lineage:
            files.extend([step["iv_file"], step["config_file"]])
            step_config = json.loads(Path(step["config_file"]).read_text(encoding="utf-8"))
            files.append(step_config["sweep"]["initial_state_file"])
        point["evidence"].update(hashes(files))
        return point

    def _check_audit_binding(self, path, state, bias):
        probe = json.loads(Path(path).read_text(encoding="utf-8"))
        if Path(probe["state_file"]).resolve() != Path(state).resolve():
            raise EvidenceError("audit state binding mismatch")
        driven = self.base["sweep"]["contact"]
        contacts = [c for c in probe["contacts"] if c["name"] == driven]
        if len(contacts) != 1 or finite(contacts[0]["bias"]) != bias or probe.get("simulation_type") != "newton_carrier_term_probe":
            raise EvidenceError("audit actual bias binding mismatch")

    def _lineage_updates(self, lineage, bias):
        total = 0
        for step in lineage:
            row = read_rows(step["iv_file"])[-1]
            config = json.loads(Path(step["config_file"]).read_text(encoding="utf-8"))
            actual = finite(row["bias_V"])
            updates = int(row["newton_iterations"])
            if (not math.isclose(actual, bias, abs_tol=1e-12, rel_tol=0)
                    or config["sweep"]["bias_points"] != [bias]
                    or config["solver"]["method"] != "newton"
                    or config["sweep"]["external_circuit"]["enabled"]
                    or row["solver_method"] != "newton" or updates < 0 or updates != step["updates"]):
                raise EvidenceError("coupled update lineage mismatch")
            total += updates
        return total

    def solve(self, bias, seed, directory, same_bias=False):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.lineage = []
        initial = Path(seed["state_file"])
        if self.prewarm and bias != seed["inner_V"]:
            _, _, prewarm_dir, passed = self._write_run(self.helper.render_poisson_prewarm(
                self.base, directory, bias, initial, 1, self.prewarm_iterations))
            if passed:
                initial = prewarm_dir / "final_state.csv"
        document, config, run_dir, passed = self._write_run(self.helper.render_attempt(
            self.base, directory, bias, initial, 1, self.max_newton_iterations))
        candidate = run_dir / "final_state.csv"
        failure = NumericalFailure("device_numerical: full device solve did not converge")
        if passed:
            try:
                return self._qualified(document, config, run_dir, candidate, bias, same_bias, seed)
            except NumericalFailure as exc:
                failure = exc
        if not candidate.is_file():
            candidate = self.helper.load_rejected_final_state(run_dir / "attempts.csv")
        if self.corrector_passes and candidate is not None and candidate.is_file():
            result = self.helper.correct_local_candidate(document, candidate, bias, self.runner,
                self.repo, run_dir, self.timeout, config.parent, self.corrector_passes,
                self.corrector_mode, self.restore_poly)
            if result.get("passed"):
                final = Path(result["state_file"])
                # The helper already audited this file. Reuse its bound final
                # probe; never overwrite immutable completed audit evidence.
                validation = final.parent
                config = validation / "config.json"
                audit_config = validation / "local_audit/config.json"
                self._check_audit_binding(audit_config, final, bias)
                point = qualify_files(final, validation / "iv.csv", validation / "local_audit/terms.csv",
                    bias, self.base, self.mesh_ids, expected_count=self.expected_active_count,
                    expected_active=seed.get("active_ids"), updates=self._lineage_updates(self.lineage, bias),
                    same_bias=same_bias)
                proof = {"state_sha256": point["state_sha256"], "inner_V": bias,
                         "current_A_per_um": point["current_A_per_um"],
                         "audit_config": str(audit_config.resolve()), "lineage": clone(self.lineage)}
                proof_path = validation / "qualification.json"
                atomic_json(proof_path, proof)
                point["proof_file"] = str(proof_path.resolve())
                files = [config, audit_config, proof_path, Path(seed["state_file"])]
                for step in self.lineage:
                    files.extend([step["config_file"], step["iv_file"]])
                    files.append(json.loads(Path(step["config_file"]).read_text(encoding="utf-8"))["sweep"]["initial_state_file"])
                point["evidence"].update(hashes(files))
                return point
            failure = NumericalFailure("local_corrector_unqualified: bounded correction exhausted")
        raise failure

    def validate(self, point):
        verify_hashes(point["evidence"])
        proof = json.loads(Path(point["proof_file"]).read_text(encoding="utf-8"))
        self._check_audit_binding(proof["audit_config"], point["state_file"], point["inner_V"])
        updates = self._lineage_updates(proof["lineage"], point["inner_V"])
        checked = qualify_files(Path(point["state_file"]), Path(point["iv_file"]), Path(point["terms_file"]),
            point["inner_V"], self.base, self.mesh_ids, expected_count=self.expected_active_count,
            expected_active=point["active_ids"], updates=updates, same_bias=point["same_bias_requalification"])
        for key in ("inner_V", "current_A_per_um", "state_sha256", "coupled_newton_updates", "active_ids"):
            if checked[key] != point[key]:
                raise EvidenceError(f"qualified point metadata mismatch: {key}")
        for key in ("inner_V", "current_A_per_um", "state_sha256"):
            if proof[key] != point[key]:
                raise EvidenceError(f"qualification proof binding mismatch: {key}")

    def confirm(self, point, outer, directory):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        document, config, run_dir = self.helper.render_attempt(self.base, directory, outer,
            Path(point["state_file"]), 1, self.max_newton_iterations)
        document["sweep"]["external_circuit"].update(enabled=True, initial_inner_voltage_V=point["inner_V"])
        _, _, run_dir, passed = self._write_run((document, config, run_dir))
        if not passed:
            raise NumericalFailure("confirmation_unqualified: circuit solve failed")
        state = run_dir / "final_state.csv"
        if sha256(state) != point["state_sha256"]:
            raise NumericalFailure("confirmation_state_changed: full device requalification is required; target remains pending")
        confirmation = {"state_sha256": sha256(state), "state_file": str(state.resolve()),
                        "iv_file": str((run_dir / "iv.csv").resolve()), "config_file": str(config.resolve()),
                        "evidence": hashes([state, run_dir / "iv.csv", config])}
        self.validate_confirmation(point, outer, confirmation)
        return confirmation

    def validate_confirmation(self, point, outer, confirmation):
        verify_hashes(confirmation["evidence"])
        row = read_rows(confirmation["iv_file"])[-1]
        config = json.loads(Path(confirmation["config_file"]).read_text(encoding="utf-8"))
        circuit = config["sweep"]["external_circuit"]
        if (row.get("converged") != "1" or row.get("boundary_control_mode") != "external_resistor"
                or not circuit["enabled"] or circuit["current_direction"] != self.base["sweep"]["external_circuit"]["current_direction"]
                or finite(row["series_resistance_ohm_um"]) != RESISTANCE
                or finite(row["outer_voltage_V"]) != outer
                or not math.isclose(finite(row["inner_voltage_V"]), point["inner_V"], rel_tol=0, abs_tol=1e-12)
                or finite(row["current_total_A_per_um"]) != point["current_A_per_um"]
                or row["bias_contact"] != self.base["sweep"]["contact"]
                or row["current_contact"] != self.base["sweep"].get("current_contact", self.base["sweep"]["contact"])
                or sha256(confirmation["state_file"]) != point["state_sha256"]):
            raise EvidenceError("confirmation terminal/current/voltage/state identity mismatch")
        value = residual(finite(row["inner_voltage_V"]), finite(row["current_total_A_per_um"]), outer, circuit["current_direction"])
        if abs(value) > TOLERANCE or not math.isclose(finite(row["load_line_residual_V"]), value, rel_tol=0, abs_tol=1e-12):
            raise EvidenceError("confirmation load-line mismatch")


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    for name in ("repo-root", "runner", "base-config", "initial-state", "output-dir"):
        result.add_argument("--" + name, type=Path, required=True)
    result.add_argument("--initial-inner-bias", type=float, required=True)
    result.add_argument("--initial-outer-bias", type=float,
                        help="optional asserted outer coordinate; otherwise derive it from the requalified seed")
    result.add_argument("--outer-targets", type=float, nargs="+", required=True)
    result.add_argument("--resume", action="store_true")
    result.add_argument("--max-attempts", type=int, default=100,
                        help="total device and confirmation attempts, including interrupted attempts")
    result.add_argument("--max-shots", type=int, default=12, help="device attempts per outer target before outer-step halving")
    result.add_argument("--max-inner-retries", type=int, default=8)
    result.add_argument("--max-inner-step", type=float, default=.05)
    result.add_argument("--min-inner-step", type=float, default=1e-9)
    result.add_argument("--min-outer-step", type=float, default=1e-6)
    result.add_argument("--max-total-seconds", type=float, default=7200.)
    result.add_argument("--timeout-seconds", type=float, default=240.)
    result.add_argument("--max-newton-iterations", type=int,
                        help="optional solve budget; final physics and tolerances remain frozen")
    result.add_argument("--poisson-prewarm", action="store_true")
    result.add_argument("--poisson-prewarm-max-iterations", type=int, default=60)
    result.add_argument("--carrier-corrector-passes", type=int, default=0)
    result.add_argument("--carrier-corrector-mode", choices=("carrier_only", "carrier_density"), default="carrier_only")
    result.add_argument("--restore-poly-equilibrium", action="store_true")
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if (not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0
                or args.poisson_prewarm_max_iterations <= 0 or args.carrier_corrector_passes < 0
                or (args.max_newton_iterations is not None and args.max_newton_iterations <= 0)):
            raise ValueError("invalid solver time/iteration budget")
        base = json.loads(args.base_config.read_text(encoding="utf-8"))
        validate_profile(base)
        base, identity = freeze_inputs(base, args.base_config.resolve().parent, args.runner)
        identity["files"].update(hashes([args.base_config, args.initial_state, Path(__file__),
            Path(__file__).with_name("run_slot_ldmos_noavalanche_continuation.py"),
            Path(__file__).with_name("audit_slot_ldmos_local_continuity.py"),
            Path(__file__).with_name("slot_poly_equilibrium.py")]))
        identity["frozen_config"] = base
        identity["thread_environment"] = {key: os.environ.get(key) for key in
                                          ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
        solver_options = {key: getattr(args, key) for key in (
            "timeout_seconds", "max_newton_iterations", "poisson_prewarm",
            "poisson_prewarm_max_iterations", "carrier_corrector_passes",
            "carrier_corrector_mode", "restore_poly_equilibrium")}
        identity["solver_options"] = solver_options
        backend = SolverBackend(base, args.repo_root, args.runner, **solver_options)
        budget_options = {key: getattr(args, key) for key in (
            "max_attempts", "max_shots", "max_inner_retries", "max_inner_step",
            "min_inner_step", "min_outer_step", "max_total_seconds")}
        campaign = Campaign(args.output_dir, identity, args.outer_targets, args.initial_state,
            args.initial_inner_bias, base["sweep"]["external_circuit"]["current_direction"], backend,
            resume=args.resume, initial_outer=args.initial_outer_bias, **budget_options)
        complete = campaign.run()
        print(json.dumps({"completed": complete, "checkpoint": str(campaign.output / "checkpoint.json"),
            "failure": campaign.snapshot["failure"], "pending_targets_V": campaign.snapshot["pending_targets_V"],
            "accepted_outer_V": [p["outer_V"] for p in campaign.snapshot["outer_ledger"]]}, allow_nan=False))
        return 0 if complete else 1
    except KeyboardInterrupt:
        print(json.dumps({"completed": False, "failure_class": "interrupted"}), file=sys.stderr)
        return 130
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        kind = "missing_or_corrupt_evidence" if isinstance(exc, EvidenceError) else "runtime_environment" if isinstance(exc, (OSError, RuntimeError)) else "invalid_input"
        print(json.dumps({"completed": False, "failure_class": kind, "reason": str(exc)}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
