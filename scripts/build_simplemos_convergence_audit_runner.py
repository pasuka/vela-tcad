"""Build an isolated read-only acceptance probe; leave frozen production files intact."""
from pathlib import Path
import json
import os
import subprocess

REPO = Path(__file__).resolve().parents[1]
LOCAL = REPO / "build-release/simplemos_convergence_audit"
RUNNER = LOCAL / "vela_convergence_audit.exe"

METHOD = r'''
std::pair<NewtonCarrierRowConvergenceEvaluation, NewtonGlobalContinuityClosureEvaluation>
NewtonSolver::evaluateAcceptanceReadOnly(const DDSolution& state) const
{
    auto assembler = makeArclengthAssembler();
    restoreElectronQuantumPotential(*assembler, state);
    const auto bcs = buildBoundaryConditions(*assembler);
    const VectorXd x = packReferencedSolution(*assembler, state, bcs);
    const auto local = evaluateCarrierRowConvergence(
        assembler->carrierContinuityEquationTermDiagnostics(x, bcs), cfg_.carrierRowConvergence);
    std::vector<Index> electrons, holes;
    for (const auto& [node, value] : bcs.phin) electrons.push_back(node);
    for (const auto& [node, value] : bcs.phip) holes.push_back(node);
    const auto global = evaluateGlobalContinuityClosure(
        assembler->carrierContinuityEquationTermDiagnostics(x, CoupledDDBoundaryConditions{}),
        electrons, holes, cfg_.globalContinuityClosure);
    return {local, global};
}
'''


def build():
    if RUNNER.exists():
        raise FileExistsError(RUNNER)
    overlay = LOCAL / "include/vela/solver/NewtonSolver.h"
    overlay.parent.mkdir(parents=True, exist_ok=True)
    header = (REPO / "include/vela/solver/NewtonSolver.h").read_text()
    marker = "class NewtonSolver {\npublic:"
    assert header.count(marker) == 1
    header = header.replace(marker, marker + "\n    std::pair<NewtonCarrierRowConvergenceEvaluation, NewtonGlobalContinuityClosureEvaluation>\n    evaluateAcceptanceReadOnly(const DDSolution& state) const;\n")
    overlay.write_text("#include <utility>\n" + header)
    source = (REPO / "src/solver/NewtonSolver.cpp").read_text()
    marker = "NewtonCarrierTermDiagnosticsEvaluation NewtonSolver::evaluateCarrierTermDiagnostics("
    assert source.count(marker) == 1
    (LOCAL / "NewtonSolver.cpp").write_text(source.replace(marker, METHOD + "\n" + marker))
    runner = (REPO / "src/tools/vela_example_runner.cpp").read_text()
    start = runner.index("nlohmann::json runNewtonCarrierTermProbe(")
    end = runner.index("void writeSgEdgeFluxProbeCsv", start)
    part = runner[start:end]
    part = part.replace("    return {", "    const auto acceptance = solver.evaluateAcceptanceReadOnly(state);\n    return {\n        {\"read_only\", true},\n        {\"carrier_row_convergence\", carrierRowConvergenceJson(acceptance.first)},\n        {\"global_continuity_closure\", globalContinuityClosureJson(acceptance.second)},", 1)
    runner = runner[:start] + part + runner[end:]
    (LOCAL / "runner.cpp").write_text(runner)
    tc = Path("D:/msys64/ucrt64")
    args = [str(tc / "bin/g++.exe"), "-O2", "-DNDEBUG", "-std=c++20",
        "-DSPDLOG_COMPILED_LIB", "-DSPDLOG_FMT_EXTERNAL", "-DSPDLOG_SHARED_LIB",
        "-DVELA_HAS_SPQR=1", "-DVELA_HAS_UMFPACK=1", '-DVELA_VERSION="0.1.0"',
        "-I"+str(LOCAL / "include"), "-I"+str(REPO / "include"),
        "-I"+str(tc / "include/eigen3"), "-I"+str(tc / "include/suitesparse"),
        str(LOCAL / "runner.cpp"), str(LOCAL / "NewtonSolver.cpp"), str(REPO / "build-release/libvela_core.a"),
        *[str(tc / "lib" / p) for p in ("libspdlog.dll.a", "libfmt.a", "libumfpack.dll.a", "libspqr.dll.a", "libcholmod.dll.a")],
        "-o", str(RUNNER)]
    (LOCAL / "build_command.json").write_text(json.dumps(args, indent=2))
    env = os.environ.copy(); env["PATH"] = str(tc / "bin") + os.pathsep + env["PATH"]
    result = subprocess.run(args, capture_output=True, text=True, env=env)
    (LOCAL / "build.stdout.txt").write_text(result.stdout)
    (LOCAL / "build.stderr.txt").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(result.stderr[-5000:])
    print("Isolated acceptance probe built; original solver, headers and runner unchanged.")


if __name__ == "__main__":
    build()
