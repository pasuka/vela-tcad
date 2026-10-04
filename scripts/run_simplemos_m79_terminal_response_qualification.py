"""Freeze and qualify full-DD terminal sensitivity against the existing M74 A/B.

Only read-only runner probes and linear solves are dispatched. Source states,
historical evidence, production defaults and remote decks are never written.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
import math
import os
from pathlib import Path
import subprocess

import run_simplemos_m78_electron_volume_localization as m78

m74 = m78.m74
REPO, ROOT = m74.REPO, m74.ROOT
SCRIPT = Path(__file__).resolve()
SLUG = "m79_terminal_response_qualification"
CONTRACT = ROOT / f"simplemos_{SLUG}_contract_v1.json"
FREEZE = ROOT / f"simplemos_{SLUG}_contract_freeze_v1.json"
EVIDENCE = ROOT / f"simplemos_{SLUG}_evidence.json"
OUT = ROOT / "terminal_response_qualification"
LOCAL = REPO / "build-release" / SLUG
DOC = REPO / "docs/validation/simplemos_m79_terminal_response_qualification_2026-09-05.md"
TEST = REPO / "tests/regression/test_simplemos_m79_terminal_response_qualification.py"
read_json, write_json, read_csv, write_csv = m74.read_json, m74.write_json, m74.read_csv, m74.write_csv
sha256, portable = m74.sha256, m74.portable
PILOT = ("n19", "n23", "n20", "n24")


def freeze():
    if CONTRACT.exists() or FREEZE.exists():
        raise ValueError("M79 already frozen; do not overwrite")
    m78.verify()
    inputs, identity = m78.historical_inputs()
    inputs.update((m78.CONTRACT, m78.FREEZE, m78.EVIDENCE, SCRIPT, TEST,
                   m74.RUNNER, REPO / "src/tools/vela_example_runner.cpp",
                   REPO / "src/solver/NewtonSolver.cpp",
                   REPO / "src/equation/CoupledDDAssembler.cpp",
                   REPO / "include/vela/equation/DDAssembler.h",
                   REPO / "docs/validation/simplemos_post_m78_research_debug_plan_2026-09-05.md"))
    contract = {
        "schema": "vela.simplemos.m79_terminal_response_qualification.v1",
        "purpose": "Qualify full-DD terminal response on known M65/M74 intervention before M80 attribution.",
        "supersedes_proposed_action_only": "M78 conditional local-support A/B is not executed; user authorized post-M78 research plan instead.",
        "matrix": read_json(m74.CONTRACT)["matrix"], "pilot_devices": list(PILOT),
        "execution": {"new_nonlinear_solves": 0, "new_sentaurus_runs": 0,
                      "gate_states_only": True, "expand_only_after_pilot_pass": True,
                      "no_production_changes": True},
        "definitions": {
            "reference": "Frozen matched-ni/no-BGN M65 and electron-only M74 gate endpoints",
            "parameter": "alpha from legacy electron Poisson volume (0) to Si barycentric (1)",
            "source": "F_alpha=F_candidate(u)-F_baseline(u) at each endpoint; coordinates and all other settings identical",
            "current": "Same signed drain current A/um as frozen curve; check both residual functional and production extractor",
            "prediction": "-lambda.T F_alpha plus explicit g_alpha; compare dI/(I ln10) at both endpoints against log10|I1/I0|",
            "finite_difference": "Centered +/- directions, max physical component 1e-5 and 5e-6 V; response plus deterministic psi/phin/phip directions with contact exclusion",
            "signal": "Current directional derivative above max(1e-20 A/um, 1e-12*full_gradient_norm); other rows reported but not qualified",
            "failure": "No M80 on failed identity/Jacobian/duality gates. Finite-amplitude failure requires separately frozen small-amplitude calibration. Never relax thresholds retrospectively."
        },
        "thresholds": {"current_replay_relative": 1e-8, "adjoint_relative_residual": 1e-10,
                       "jvp_relative": 1e-3, "current_fd_relative": 1e-3,
                       "duality_relative": 1e-8, "linear_solve_relative": 1e-8,
                       "finite_prediction_relative": 0.1, "finite_prediction_floor_dex": 1e-4,
                       "continuity_source_norm": 0.0, "direct_current_term": 0.0},
    }
    write_json(CONTRACT, contract)
    write_json(FREEZE, {"status": "frozen_before_execution", "contract_sha256": sha256(CONTRACT),
                       "baseline_stage_count": len(identity), "candidate_stage_count": len(identity),
                       "input_hashes": {portable(p): sha256(p) for p in sorted(inputs)}})
    write_csv(OUT / "m79_input_identity_ledger.csv", identity)
    return {"status": "frozen_before_execution", "inputs": len(inputs), "stage_states": 2*len(identity)}


def validate():
    f = read_json(FREEZE)
    if f["status"] != "frozen_before_execution":
        raise ValueError("M79 not frozen")
    m78.check_hash(CONTRACT, f["contract_sha256"])
    for rel, digest in f["input_hashes"].items():
        m78.check_hash(REPO / rel, digest)
    return read_json(CONTRACT)


def run_probe(deck, config):
    write_json(config, deck)
    env = os.environ.copy()
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    result = subprocess.run([str(m74.RUNNER), "--config", str(config), "--log", "off"],
                            cwd=REPO, env=env, capture_output=True, text=True, check=False)
    config.with_suffix(".stdout.txt").write_text(result.stdout, encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"read-only probe failed: {config}: {result.stderr[-1500:]}")
    status = json.loads(result.stdout.strip().splitlines()[-1])
    write_json(config.with_suffix(".status.json"), status)
    return status


def base_deck(source, state, workflow):
    deck = copy.deepcopy(read_json(source))
    for key in ("sweep", "log_file", "output_csv"):
        deck.pop(key, None)
    for c in deck["contacts"]:
        if c["name"] == "gate":
            c["bias"] = workflow["gate_voltage_V"]
        if c["name"] == "drain":
            c["bias"] = workflow["drain_voltage_V"]
    deck["state_file"] = str(state.resolve())
    deck["contact"] = "drain"
    return deck


def qualify_endpoint(status, functional, reference, thresholds):
    response = status["electron_volume_response"]
    currents = [status["current_A_per_um"], functional["current_A_per_um"],
                functional["contact_current_extractor_A_per_um"]]
    replay = max(abs(i/reference-1) for i in currents)
    fd = response["finite_difference_checks"]
    signal_floor = max(1e-20, 1e-12*status["state_derivative_norm"])
    signal_rows = [r for r in fd if abs(r["current_analytic"]) > signal_floor]
    gates = {
        "current_identity": replay <= thresholds["current_replay_relative"],
        "adjoint": status["adjoint_relative_residual"] <= thresholds["adjoint_relative_residual"],
        "jvp": all(r["jvp_relative_error"] <= thresholds["jvp_relative"] for r in fd),
        "current_derivative": bool(signal_rows) and all(r["current_relative_error"] <= thresholds["current_fd_relative"] for r in signal_rows),
        "duality": response["duality_relative_error"] <= thresholds["duality_relative"],
        "linear_solve": response["linear_solve_relative_residual"] <= thresholds["linear_solve_relative"],
        "single_axis": response["continuity_source_norm"] == thresholds["continuity_source_norm"] and
                       response["parameter_direct_current_A_per_um"] == thresholds["direct_current_term"],
    }
    return gates, {"current_replay_relative": replay,
                   "adjoint_relative_residual": status["adjoint_relative_residual"],
                   "jvp_max_relative": max(r["jvp_relative_error"] for r in fd),
                   "current_fd_max_relative": max((r["current_relative_error"] for r in signal_rows), default=0.0),
                   "duality_relative": response["duality_relative_error"],
                   "linear_solve_relative": response["linear_solve_relative_residual"],
                   "qualified_current_fd_rows": len(signal_rows)}


def run_case(workflow):
    phase = workflow["phases"][-1]
    case_dir = LOCAL / workflow["device"] / f"vd_{m74.voltage_tag(workflow['drain_voltage_V'])}"
    rows, fdrows = [], []
    for role in ("baseline", "candidate"):
        source = REPO / phase["source_config" if role == "baseline" else "config"]
        state = source.parent / "state.csv"
        reference = float(read_csv(source.parent / "curve.csv")[-1]["current_total_A_per_um"])
        deck = base_deck(source, state, workflow)
        root = case_dir / role
        adj = dict(deck, simulation_type="terminal_current_adjoint_probe",
                   output_csv=str(root / "adjoint.csv"),
                   electron_volume_response={"output_csv": str(root / "response.csv")})
        status = run_probe(adj, root / "adjoint.json")
        functional = run_probe(dict(deck, simulation_type="terminal_current_functional_probe"), root / "functional.json")
        gates, metrics = qualify_endpoint(status, functional, reference, read_json(CONTRACT)["thresholds"])
        row = {"device": workflow["device"], "drain_voltage_V": workflow["drain_voltage_V"],
               "gate_voltage_V": workflow["gate_voltage_V"], "role": role,
               "reference_current_A_per_um": reference, "replay_current_A_per_um": status["current_A_per_um"],
               "response_A_per_um": status["electron_volume_response"]["adjoint_response_A_per_um"],
               "prediction_dex": status["electron_volume_response"]["adjoint_response_A_per_um"]/(reference*math.log(10)),
               **metrics, **{f"pass_{k}": v for k, v in gates.items()},
               "observer_pass": all(gates.values())}
        rows.append(row)
        for check in status["electron_volume_response"]["finite_difference_checks"]:
            fdrows.append({"device": workflow["device"], "drain_voltage_V": workflow["drain_voltage_V"],
                           "role": role, **{k: v for k, v in check.items() if k != "blocks"}})
    delta = math.log10(abs(rows[1]["reference_current_A_per_um"]/rows[0]["reference_current_A_per_um"]))
    th = read_json(CONTRACT)["thresholds"]
    tolerance = max(th["finite_prediction_relative"]*abs(delta), th["finite_prediction_floor_dex"])
    for row in rows:
        row.update(actual_delta_dex=delta, prediction_error_dex=row["prediction_dex"]-delta,
                   prediction_tolerance_dex=tolerance,
                   finite_prediction_pass=abs(row["prediction_dex"]-delta) <= tolerance)
    print(f"completed {workflow['device']} Vd={workflow['drain_voltage_V']}: "
          f"observer={all(r['observer_pass'] for r in rows)} finite={all(r['finite_prediction_pass'] for r in rows)}", flush=True)
    return rows, fdrows


def analyze_batch(rows, fdrows, matrix_complete):
    write_csv(OUT / "m79_endpoint_ledger.csv", rows)
    write_csv(OUT / "m79_directional_derivative_ledger.csv", fdrows)
    observer_pass = all(r["observer_pass"] for r in rows)
    finite_pass = all(r["finite_prediction_pass"] for r in rows)
    report = {"schema": "vela.simplemos.m79.report.v1", "endpoints": len(rows),
              "cases": len(rows)//2, "matrix_complete": matrix_complete,
              "observer_pass": observer_pass, "finite_prediction_pass": finite_pass,
              "m80_released": matrix_complete and observer_pass and finite_pass,
              "new_nonlinear_solves": 0, "new_sentaurus_runs": 0,
              "failed_gates": sorted({k for r in rows for k,v in r.items() if k.startswith("pass_") and not v}),
              "status": "qualified" if matrix_complete and observer_pass and finite_pass else "stopped_at_qualification_gate",
              "maximum": {k: max(r[k] for r in rows) for k in ("current_replay_relative", "adjoint_relative_residual",
                           "jvp_max_relative", "current_fd_max_relative", "duality_relative", "linear_solve_relative")}}
    write_json(OUT / "m79_report.json", report)
    lines = ["# M79 完整 DD 端口电流响应资格化", "", "日期：2026-09-05。", "",
             "本阶段只读取 M65/M74 冻结状态；没有新增非线性器件求解或 Sentaurus 运行。", "",
             f"状态：`{report['status']}`；完成 {report['cases']} 工况 / {len(rows)} 端点。",
             f"观察器通过：{observer_pass}；有限幅度预测通过：{finite_pass}；M80 放行：{report['m80_released']}。", "",
             "| 器件 | Vd | 端点 | 电流回放相对误差 | Jv 相对误差 | 电流导数相对误差 | 预测/实际 dex |",
             "|---|---:|---|---:|---:|---:|---|"]
    for r in rows:
        lines.append(f"| {r['device']} | {r['drain_voltage_V']} | {r['role']} | {r['current_replay_relative']:.3e} | {r['jvp_max_relative']:.3e} | {r['current_fd_max_relative']:.3e} | {r['prediction_dex']:.6g} / {r['actual_delta_dex']:.6g} |")
    lines += ["", "失败项：" + ", ".join(report["failed_gates"]), "",
              "失败结果保留原判据，不放宽门槛；任何修正后重跑必须另建版本与冻结记录。",
              "不能将未资格化的伴随分项解释为 Vela–Sentaurus 的物理根因。", ""]
    DOC.write_text("\n".join(lines), encoding="utf-8")
    return report


def run(workers):
    validate()
    if (OUT / "m79_report.json").exists():
        raise ValueError("existing M79 results must not be overwritten")
    workflows = read_json(m74.PORTABLE_MANIFEST)["workflows"]
    rows, fdrows = [], []
    for pilot in (True, False):
        batch = [w for w in workflows if (w["device"] in PILOT) == pilot]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for case, fd in pool.map(run_case, batch):
                rows.extend(case)
                fdrows.extend(fd)
        report = analyze_batch(rows, fdrows, not pilot)
        if not report["observer_pass"] or not report["finite_prediction_pass"]:
            break
    validate()
    artifacts = list(OUT.glob("*.csv")) + [OUT / "m79_report.json", DOC]
    artifacts += list(LOCAL.rglob("*.json")) + list(LOCAL.rglob("*.csv"))
    write_json(EVIDENCE, {"status": "frozen", "qualification_status": report["status"],
                         "contract_sha256": sha256(CONTRACT),
                         "artifacts": {portable(p): sha256(p) for p in sorted(artifacts)}})
    return report


def verify():
    validate()
    evidence = read_json(EVIDENCE)
    m78.check_hash(CONTRACT, evidence["contract_sha256"])
    for rel, digest in evidence["artifacts"].items():
        m78.check_hash(REPO / rel, digest)
    return read_json(OUT / "m79_report.json")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--freeze-contract", action="store_true")
    group.add_argument("--run", action="store_true")
    group.add_argument("--verify", action="store_true")
    p.add_argument("--workers", type=int, default=2)
    args = p.parse_args()
    result = freeze() if args.freeze_contract else (run(args.workers) if args.run else verify())
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
