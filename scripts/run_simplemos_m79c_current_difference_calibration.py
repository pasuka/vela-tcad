"""Signal-resolved current FD using Richardson extrapolation; no DC solves."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
import math
from pathlib import Path

import run_simplemos_m79b_numerical_calibration as previous
old = previous.old
REPO, ROOT = old.REPO, old.ROOT
LOCAL = REPO / "build-release/m79c_current_difference_calibration"
OUT = ROOT / "terminal_current_difference_calibration"
CONTRACT = ROOT / "simplemos_m79c_current_difference_calibration_contract_v1.json"
FREEZE = ROOT / "simplemos_m79c_current_difference_calibration_contract_freeze_v1.json"
EVIDENCE = ROOT / "simplemos_m79c_current_difference_calibration_evidence.json"
DOC = REPO / "docs/validation/simplemos_m79c_current_difference_calibration_2026-09-05.md"
SCRIPT = Path(__file__).resolve()
STEPS = (0.008, 0.004, 0.002)


def richardson(coarse, fine):
    return (4*fine-coarse)/3


def perturb(rows, direction, step, potential_scale):
    result = copy.deepcopy(rows)
    for row, tangent in zip(result, direction, strict=True):
        dp, dn = step*tangent[0], step*tangent[1]
        row["psi"] = float(row["psi"])+dp
        row["electron_qf_increment_V"] = float(row["electron_qf_increment_V"])+dn
        row["phin"] = float(row["electron_qf_reference_V"])+row["electron_qf_increment_V"]
        row["electrons_m3"] = float(row["electrons_m3"])*math.exp((dp-dn)/potential_scale)
        row["holes_m3"] = float(row["holes_m3"])*math.exp(-dp/potential_scale)
    return result


def freeze():
    previous.verify()
    if CONTRACT.exists() or FREEZE.exists():
        raise ValueError("M79c already frozen")
    contract = old.read_json(previous.CONTRACT)
    contract.update(schema="vela.simplemos.m79c_current_difference_calibration.v1",
                    purpose="Resolve high-NWell current FD signal; retain previous Jv, adjoint and direct-solve gates.")
    contract["current_difference_calibration"] = {
        "steps_V": STEPS, "direction": "Same normalized psi/phin response tangent as M79b; hole tangent zero only in this FD test.",
        "method": "Centered physical-state current differences then (4*D(h/2)-D(h))/3 for two consecutive pairs.",
        "acceptance": "Both Richardson estimates and their mutual difference relative to analytic derivative <=1e-3. Other measured psi/phin/phip current directions retain M79b thresholds. Full J, adjoint and source unchanged.",
        "state": "Write disposable tangent states preserving reference/increment consistency and Boltzmann n,p; do not alter frozen source states.",
        "stop": "No M80 unless all 16 full-matrix endpoints pairs pass all gates; preserve previous failures."
    }
    paths = [SCRIPT, previous.SCRIPT, previous.CONTRACT, previous.FREEZE, previous.EVIDENCE,
             previous.RUNNER, old.SCRIPT, old.CONTRACT, old.FREEZE, old.EVIDENCE,
             REPO / "tests/regression/test_simplemos_m79c_current_difference_calibration.py"]
    old.write_json(CONTRACT, contract)
    old.write_json(FREEZE, {"status": "frozen_before_execution", "contract_sha256": old.sha256(CONTRACT),
                           "input_hashes": {old.portable(p): old.sha256(p) for p in paths}})
    return {"status": "frozen_before_execution", "steps_V": STEPS}


def validate():
    previous.verify()
    frozen = old.read_json(FREEZE)
    old.m78.check_hash(CONTRACT, frozen["contract_sha256"])
    for rel,digest in frozen["input_hashes"].items():
        old.m78.check_hash(REPO / rel, digest)


def case(workflow):
    pilot = workflow["device"] in old.PILOT
    if not pilot:
        old.run_case(workflow)
    source_root = previous.LOCAL if pilot else LOCAL
    case_rel = Path(workflow["device"]) / f"vd_{old.m74.voltage_tag(workflow['drain_voltage_V'])}"
    phase = workflow["phases"][-1]
    endpoints, differences = [], []
    for role in ("baseline", "candidate"):
        root = LOCAL / case_rel / role
        original = source_root / case_rel / role
        status = old.read_json(original / "adjoint.status.json")
        functional = old.read_json(original / "functional.status.json")
        response = old.read_csv(original / "response.csv")
        adjoint = old.read_csv(original / "adjoint.csv")
        norm = max(abs(float(r["response_"+f+"_scaled"])) for r in response for f in ("psi", "phin"))
        tangent = [[float(r["response_"+f+"_scaled"])/norm for f in ("psi", "phin")] for r in response]
        analytic = math.fsum(float(a["dI_d"+f+"_scaled"])*t[i] for a,t in zip(adjoint,tangent) for i,f in enumerate(("psi", "phin")))
        config = REPO / phase["source_config" if role == "baseline" else "config"]
        state = config.parent / "state.csv"
        state_rows = old.read_csv(state)
        scale = status["potential_scale_V"]
        deck = old.base_deck(config, state, workflow)
        centered = []
        for step in STEPS:
            values = []
            for sign in (1,-1):
                label = f"h_{step:g}_{'plus' if sign>0 else 'minus'}"
                temporary = root / (label + ".csv")
                old.write_csv(temporary, perturb(state_rows, tangent, sign*step, scale))
                value = old.run_probe(dict(deck, simulation_type="terminal_current_functional_probe",
                    state_file=str(temporary)), root / (label+".json"))["current_A_per_um"]
                values.append(value)
            derivative = (values[0]-values[1])/(2*step/scale)
            centered.append(derivative)
            differences.append({"device": workflow["device"], "drain_voltage_V": workflow["drain_voltage_V"],
                                "role": role, "step_V": step, "analytic": analytic,
                                "plus_current_A_per_um": values[0], "minus_current_A_per_um": values[1],
                                "central_derivative": derivative, "relative_error": abs(derivative/analytic-1)})
        extrapolated = [richardson(centered[i], centered[i+1]) for i in (0,1)]
        errors = [abs(x/analytic-1) for x in extrapolated]
        convergence = abs((extrapolated[0]-extrapolated[1])/analytic)
        reference = float(old.read_csv(config.parent / "curve.csv")[-1]["current_total_A_per_um"])
        # Preserve mandatory non-response directions; score the response via
        # the pre-frozen higher-order current difference instead of noisy raw FD.
        scored = copy.deepcopy(status)
        for row in scored["electron_volume_response"]["finite_difference_checks"]:
            if row["direction"] == "response_psi_phin_tangent":
                row["current_relative_error"] = max(*errors, convergence)
        gates, metrics = old.qualify_endpoint(scored, functional, reference, old.read_json(CONTRACT)["thresholds"])
        endpoints.append({"device": workflow["device"], "drain_voltage_V": workflow["drain_voltage_V"],
                          "gate_voltage_V": workflow["gate_voltage_V"], "role": role,
                          "reference_current_A_per_um": reference, "response_analytic": analytic,
                          "richardson_coarse": extrapolated[0], "richardson_fine": extrapolated[1],
                          "richardson_max_relative_error": max(errors), "richardson_convergence_relative": convergence,
                          "prediction_dex": status["electron_volume_response"]["adjoint_response_A_per_um"]/(reference*math.log(10)),
                          **metrics, **{f"pass_{k}": v for k,v in gates.items()}, "observer_pass": all(gates.values())})
    delta = math.log10(abs(endpoints[1]["reference_current_A_per_um"]/endpoints[0]["reference_current_A_per_um"]))
    for row in endpoints:
        row.update(actual_delta_dex=delta, prediction_error_dex=row["prediction_dex"]-delta,
                   finite_prediction_pass=abs(row["prediction_dex"]-delta) <= max(.1*abs(delta),1e-4))
    print(f"M79c {workflow['device']} Vd={workflow['drain_voltage_V']}: "
          f"observer={all(r['observer_pass'] for r in endpoints)}", flush=True)
    return endpoints,differences


def run(workers):
    validate()
    if (OUT / "m79c_report.json").exists():
        raise ValueError("M79c results already exist")
    old.LOCAL, old.CONTRACT = LOCAL, CONTRACT
    old.m74.RUNNER = previous.RUNNER
    rows, fdrows = [], []
    workflows = old.read_json(old.m74.PORTABLE_MANIFEST)["workflows"]
    complete = False
    for pilot in (True,False):
        batch = [w for w in workflows if (w["device"] in old.PILOT) == pilot]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for endpoint,fd in pool.map(case,batch):
                rows.extend(endpoint); fdrows.extend(fd)
        complete = not pilot
        if not all(r["observer_pass"] and r["finite_prediction_pass"] for r in rows):
            break
    old.write_csv(OUT / "m79c_endpoint_ledger.csv", rows)
    old.write_csv(OUT / "m79c_difference_ledger.csv", fdrows)
    report = {"status": "qualified" if complete and all(r["observer_pass"] and r["finite_prediction_pass"] for r in rows) else "stopped_at_qualification_gate",
              "cases": len(rows)//2, "endpoints": len(rows), "matrix_complete": complete,
              "observer_pass": all(r["observer_pass"] for r in rows),
              "finite_prediction_pass": all(r["finite_prediction_pass"] for r in rows),
              "new_nonlinear_solves": 0, "new_sentaurus_runs": 0,
              "failed_gates": sorted({k for r in rows for k,v in r.items() if k.startswith('pass_') and not v}),
              "maximum": {k:max(r[k] for r in rows) for k in ("current_replay_relative", "adjoint_relative_residual",
                  "jvp_max_relative", "current_fd_max_relative", "duality_relative", "richardson_max_relative_error",
                  "richardson_convergence_relative")}}
    report["m80_released"] = report["status"] == "qualified"
    old.write_json(OUT / "m79c_report.json", report)
    lines = ["# M79c 电流差分信号校准", "", f"状态：`{report['status']}`。完成 {report['cases']} 工况、{len(rows)} 端点。",
             f"M80 放行：{report['m80_released']}。无新非线性求解。", "",
             "M79/M79b 的失败记录保留；本阶段采用三档对称扰动及两组 Richardson 外推，导数误差阈值仍为 1e-3。",
             "完整 DD 矩阵、伴随和参数响应保留空穴块；只对 FD 测试方向作 psi/phin 投影，并保留独立 phip 方向检查。", "",
             "| 器件 | Vd | 端点 | 外推最大相对误差 | 两组外推差 | 预测/实际 dex |", "|---|---:|---|---:|---:|---|"]
    for r in rows:
        lines.append(f"| {r['device']} | {r['drain_voltage_V']} | {r['role']} | {r['richardson_max_relative_error']:.3e} | {r['richardson_convergence_relative']:.3e} | {r['prediction_dex']:.6g}/{r['actual_delta_dex']:.6g} |")
    DOC.write_text("\n".join(lines)+"\n",encoding="utf-8")
    artifacts = list(OUT.glob("*.csv")) + [OUT / "m79c_report.json",DOC]
    artifacts += list(LOCAL.rglob("*.json")) + list(LOCAL.rglob("*.csv"))
    old.write_json(EVIDENCE,{"status":"frozen","qualification_status":report["status"],
                            "contract_sha256":old.sha256(CONTRACT),
                            "artifacts":{old.portable(p):old.sha256(p) for p in sorted(artifacts)}})
    return report


def verify():
    validate()
    e = old.read_json(EVIDENCE)
    old.m78.check_hash(CONTRACT,e["contract_sha256"])
    for rel,digest in e["artifacts"].items(): old.m78.check_hash(REPO/rel,digest)
    return old.read_json(OUT / "m79c_report.json")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--freeze-contract",action="store_true");g.add_argument("--run",action="store_true");g.add_argument("--verify",action="store_true")
    p.add_argument("--workers",type=int,default=2);args=p.parse_args()
    print(json.dumps(freeze() if args.freeze_contract else (run(args.workers) if args.run else verify()),indent=2))


if __name__ == "__main__": main()
