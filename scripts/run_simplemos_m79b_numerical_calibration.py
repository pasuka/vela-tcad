"""Separate, frozen numerical calibration of the M79 observer (no DC solves).

Build an isolated diagnostic executable from the hash-qualified M79 runner.
Keep the original executable/source, contract and failed results unchanged.
Only its diagnostic finite differences and linear refinement are adjusted.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess

import run_simplemos_m79_terminal_response_qualification as old

REPO, ROOT = old.REPO, old.ROOT
LOCAL = REPO / "build-release/m79b_numerical_calibration"
OUT = ROOT / "terminal_response_numerical_calibration"
CONTRACT = ROOT / "simplemos_m79b_numerical_calibration_contract_v1.json"
FREEZE = ROOT / "simplemos_m79b_numerical_calibration_contract_freeze_v1.json"
EVIDENCE = ROOT / "simplemos_m79b_numerical_calibration_evidence.json"
DOC = REPO / "docs/validation/simplemos_m79b_numerical_calibration_2026-09-05.md"
SCRIPT = Path(__file__).resolve()
SOURCE = LOCAL / "vela_example_runner_calibrated.cpp"
RUNNER = LOCAL / "vela_example_runner_calibrated.exe"

REFINEMENT = r'''
vela::VectorXd extendedLinearResidual(const vela::SparseMatrixd& matrix,
    const vela::VectorXd& x, const vela::VectorXd& rhs)
{
    std::vector<long double> sums(static_cast<std::size_t>(rhs.size()));
    for (int i = 0; i < rhs.size(); ++i) sums[i] = rhs(i);
    for (int col = 0; col < matrix.outerSize(); ++col)
        for (vela::SparseMatrixd::InnerIterator entry(matrix, col); entry; ++entry)
            sums[entry.row()] -= static_cast<long double>(entry.value()) * x(entry.col());
    vela::VectorXd residual(rhs.size());
    for (int i = 0; i < rhs.size(); ++i) residual(i) = static_cast<double>(sums[i]);
    return residual;
}

void refineLinearSolution(const vela::SparseMatrixd& matrix,
    const vela::VectorXd& rhs, vela::VectorXd& x)
{
    Eigen::SparseLU<vela::SparseMatrixd> factor;
    factor.compute(matrix);
    if (factor.info() != Eigen::Success)
        throw std::runtime_error("M79b refinement factorization failed");
    for (int iteration = 0; iteration < 3; ++iteration) {
        const vela::VectorXd correction = factor.solve(extendedLinearResidual(matrix, x, rhs));
        if (factor.info() != Eigen::Success || !correction.allFinite())
            throw std::runtime_error("M79b refinement correction failed");
        x += correction;
    }
}
'''


def replace_once(text, before, after):
    if text.count(before) != 1:
        raise ValueError(f"unexpected source for bounded diagnostic transform: {before[:60]}")
    return text.replace(before, after, 1)


def build():
    old.verify()
    if SOURCE.exists() or RUNNER.exists():
        raise ValueError("calibration build already exists")
    source = (REPO / "src/tools/vela_example_runner.cpp").read_text(encoding="utf-8")
    source = replace_once(source, 'nlohmann::json qualifyElectronVolumeResponse(',
                          REFINEMENT+'\nnlohmann::json qualifyElectronVolumeResponse(')
    source = replace_once(source, '    const vela::Real directResponse = evaluation.stateDerivative.dot(response);',
                          '    refineLinearSolution(jacobian, -source, response);\n'
                          '    const vela::Real directResponse = evaluation.stateDerivative.dot(response);')
    source = replace_once(source, '            direction = response;',
                          '            direction = response;\n'
                          '            // FD tangent only; full DD solution and adjoint retain holes.\n'
                          '            direction.tail(count).setZero();')
    source = replace_once(source, '        for (const vela::Real stepV : {1.0e-5, 5.0e-6}) {',
                          '        for (const vela::Real stepV : {1.0e-6, 5.0e-7, 2.5e-7}) {')
    source = replace_once(source, 'mode == 0 ? "linear_response" :',
                          'mode == 0 ? "response_psi_phin_tangent" :')
    source = replace_once(source, '    const vela::NewtonTerminalCurrentAdjointEvaluation evaluation =',
                          '    vela::NewtonTerminalCurrentAdjointEvaluation evaluation =')
    marker = '    if (evaluation.nodeCount != static_cast<int>(problem.mesh.numNodes()))'
    source = replace_once(source, marker, '''    {
        const auto packed = solver.packArclengthState(state);
        const auto system = solver.makeArclengthSystem(contact);
        vela::SparseMatrixd transpose = system.jacobian(packed, evaluation.contactBias_V).transpose();
        transpose.makeCompressed();
        refineLinearSolution(transpose, evaluation.stateDerivative, evaluation.adjoint);
        evaluation.adjointNorm = evaluation.adjoint.norm();
        evaluation.adjointResidualNorm = extendedLinearResidual(
            transpose, evaluation.adjoint, evaluation.stateDerivative).norm();
        evaluation.adjointRelativeResidual = evaluation.adjointResidualNorm / evaluation.stateDerivativeNorm;
    }
''' + marker)
    LOCAL.mkdir(parents=True, exist_ok=True)
    SOURCE.write_text('#include <Eigen/SparseLU>\n'+source, encoding="utf-8")
    toolchain = Path("D:/msys64/ucrt64")
    args = [str(toolchain / "bin/g++.exe"), "-O3", "-DNDEBUG", "-std=c++20",
            "-DSPDLOG_COMPILED_LIB", "-DSPDLOG_FMT_EXTERNAL", "-DSPDLOG_SHARED_LIB",
            "-DVELA_HAS_SPQR=1", "-DVELA_HAS_UMFPACK=1", '-DVELA_VERSION="0.1.0"',
            "-I"+str(REPO / "include"), "-I"+str(toolchain / "include/eigen3"),
            "-I"+str(toolchain / "include/suitesparse"), str(SOURCE),
            str(REPO / "build-release/libvela_core.a"),
            *[str(toolchain / "lib" / x) for x in
              ("libspdlog.dll.a", "libfmt.a", "libumfpack.dll.a", "libspqr.dll.a", "libcholmod.dll.a")],
            "-o", str(RUNNER)]
    env = os.environ.copy()
    env["PATH"] = str(toolchain / "bin") + os.pathsep + env.get("PATH", "")
    result = subprocess.run(args, cwd=REPO, env=env, capture_output=True, text=True, check=False)
    (LOCAL / "build.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (LOCAL / "build.stderr.txt").write_text(result.stderr, encoding="utf-8")
    old.write_json(LOCAL / "build_command.json", args)
    if result.returncode:
        raise RuntimeError(result.stderr[-4000:])
    return {"build": "complete", "production_runner_unchanged": True}


def freeze():
    old.verify()
    if CONTRACT.exists() or FREEZE.exists():
        raise ValueError("M79b already frozen")
    contract = old.read_json(old.CONTRACT)
    contract.update(schema="vela.simplemos.m79b_numerical_calibration.v1",
                    purpose="Resolve M79 numerical differentiation/linear solve failures; no threshold relaxation.")
    contract["calibration"] = {
        "finite_difference_steps_V": [1e-6, 5e-7, 2.5e-7],
        "response_direction": "Project only the FD test tangent onto psi/phin to avoid normalization by nearly unconstrained hole increments; full J, adjoint, parameter source and solved response retain all three blocks. Separate hole FD direction is mandatory.",
        "linear_refinement": "Three iterative-refinement corrections with long-double residual accumulation, same double matrix and RHS",
        "acceptance": "All original numerical and finite-response thresholds unchanged; all three FD steps must pass for measured signals.",
        "prior_result": "M79 remains stopped_at_qualification_gate and is not rewritten."
    }
    paths = [SCRIPT, SOURCE, RUNNER, LOCAL / "build_command.json", REPO / "build-release/libvela_core.a",
             old.SCRIPT, old.CONTRACT, old.FREEZE, old.EVIDENCE,
             REPO / "tests/regression/test_simplemos_m79b_numerical_calibration.py"]
    old.write_json(CONTRACT, contract)
    old.write_json(FREEZE, {"status": "frozen_before_execution", "contract_sha256": old.sha256(CONTRACT),
                           "input_hashes": {old.portable(p): old.sha256(p) for p in paths}})
    return {"status": "frozen_before_execution", "thresholds_unchanged": True}


def validate():
    old.verify()
    frozen = old.read_json(FREEZE)
    old.m78.check_hash(CONTRACT, frozen["contract_sha256"])
    for rel, digest in frozen["input_hashes"].items():
        old.m78.check_hash(REPO / rel, digest)


def run(workers):
    validate()
    if (OUT / "m79_report.json").exists():
        raise ValueError("M79b results already exist")
    # Reuse the exact original scoring and orchestration with separate outputs.
    old.LOCAL, old.OUT, old.CONTRACT, old.DOC = LOCAL, OUT, CONTRACT, DOC
    old.m74.RUNNER = RUNNER
    workflows = old.read_json(old.m74.PORTABLE_MANIFEST)["workflows"]
    rows, fdrows = [], []
    for pilot in (True, False):
        batch = [w for w in workflows if (w["device"] in old.PILOT) == pilot]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for case, fd in pool.map(old.run_case, batch):
                rows.extend(case); fdrows.extend(fd)
        report = old.analyze_batch(rows, fdrows, not pilot)
        if not report["observer_pass"] or not report["finite_prediction_pass"]:
            break
    document = DOC.read_text(encoding="utf-8").replace(
        "# M79 完整 DD 端口电流响应资格化", "# M79b 端口观察器数值校准")
    document += ("\n原 M79 失败结果保留。本补充阶段保持验收阈值、完整 DD 矩阵及物理模型不变；"
                 "仅采用独立诊断可执行文件进行线性迭代改进和差分步长校准。"
                 "FD 响应方向只投影测试切向量的 psi/phin 分量，完整伴随与直接响应仍保留空穴块，"
                 "并单独验证 phip 方向。\n")
    DOC.write_text(document, encoding="utf-8")
    artifacts = list(OUT.glob("*.csv")) + [OUT / "m79_report.json", DOC]
    artifacts += list(LOCAL.rglob("*.json")) + list(LOCAL.rglob("*.csv"))
    old.write_json(EVIDENCE, {"status": "frozen", "qualification_status": report["status"],
                            "contract_sha256": old.sha256(CONTRACT),
                            "artifacts": {old.portable(p): old.sha256(p) for p in sorted(artifacts)}})
    return report


def verify():
    validate()
    evidence = old.read_json(EVIDENCE)
    old.m78.check_hash(CONTRACT, evidence["contract_sha256"])
    for rel, digest in evidence["artifacts"].items():
        old.m78.check_hash(REPO / rel, digest)
    return old.read_json(OUT / "m79_report.json")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--build", action="store_true")
    g.add_argument("--freeze-contract", action="store_true")
    g.add_argument("--run", action="store_true")
    g.add_argument("--verify", action="store_true")
    p.add_argument("--workers", type=int, default=2)
    args = p.parse_args()
    result = build() if args.build else (freeze() if args.freeze_contract else
              (run(args.workers) if args.run else verify()))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
