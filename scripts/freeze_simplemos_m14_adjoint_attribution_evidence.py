#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M14 adjoint-attribution evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/adjoint_attribution"
REPORT = ROOT / "m14_adjoint_attribution_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m14_adjoint_attribution_evidence.json")
ARTIFACTS = [
    "m14_adjoint_attribution_report.json",
    "m14_state_summary.csv",
    "m14_substitution_summary.csv",
    "m14_weak_region_80_points.csv",
]
FIGURES = [
    "m14_weak_region_80_points.png",
    "m14_field_cancellation.png",
    "m14_region_response_n21_vg0p8.png",
    "m14_adjoint_quality.png",
]
IMPLEMENTATION = [
    "include/vela/solver/NewtonSolver.h",
    "src/solver/NewtonSolver.cpp",
    "src/tools/vela_example_runner.cpp",
    "tests/test_newton_solver.cpp",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m14_adjoint_attribution_contract_v1.json",
    "scripts/run_simplemos_m14_adjoint_attribution.py",
    "scripts/plot_simplemos_m14_adjoint_attribution.py",
    "scripts/freeze_simplemos_m14_adjoint_attribution_evidence.py",
    "tests/regression/test_simplemos_m14_adjoint_attribution.py",
    "docs/validation/simplemos_m14_adjoint_attribution_2026-08-29.md",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entry(path: Path) -> dict[str, str]:
    return {"path": path.resolve().relative_to(REPO.resolve()).as_posix(),
            "sha256": sha256(path)}


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    execution = report["execution"]
    if execution["weak_region_terminal_point_count"] != 80:
        raise ValueError("M14 weak-region spectrum is incomplete")
    if execution["paired_spatial_state_count"] != 16:
        raise ValueError("M14 paired spatial matrix is incomplete")
    if report["adjoint_quality"]["maximum_relative_residual"] >= 1.0e-12:
        raise ValueError("M14 adjoint solve failed the residual ceiling")
    evidence = {
        "schema": "vela.simplemos.sdevice.m14_adjoint_attribution_evidence.v1",
        "status": "complete_with_declared_field_coverage",
        "scope": "SDevice-only two-layer frozen-state and adjoint terminal-current attribution",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "coverage_boundary": report["coverage_boundary"],
        "adjoint_quality": report["adjoint_quality"],
        "weak_region": report["weak_region"],
        "conclusions": {
            "supported": [
                "All 80 high-NWell low-drain weak-region terminal-error points are covered and have Vela current above Sentaurus.",
                "All 16 production terminal-current adjoint solves close below a 1e-12 relative residual.",
                "Electron quasi-Fermi potential is the dominant frozen direct terminal-current state coordinate in the representative n21 low-drain states.",
                "At n21, Vd=0.05 V, Vg=0.8 V, self-consistent relaxation cancels 97.45 percent of the full-field frozen response and the remaining response is localized to the channel partition.",
                "The M10 mobility-to-current disconnect is mathematically consistent with terminal-low-sensitivity edge changes and self-consistent cancellation."
            ],
            "not_supported": [
                "Treating an independent carrier-density substitution as a production Newton coordinate.",
                "Claiming that phin, HFS, electrostatics, or SRH is the unique remaining cross-code root cause.",
                "Using the sub-fA Vg=0.05 V relative response for headline quantitative attribution.",
                "Extending the 16-state spatial result to all 80 points without additional Sentaurus nodal field exports.",
                "Changing a default physical model from these read-only diagnostics."
            ]
        },
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m14" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name) for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "point_count": 80, "state_count": 16}))


if __name__ == "__main__":
    main()
