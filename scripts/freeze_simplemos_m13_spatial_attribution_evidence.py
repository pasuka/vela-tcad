#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M13 spatial-attribution evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022/spatial_attribution"
          / "m13_spatial_attribution_report.json")
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m13_spatial_attribution_evidence.json")
FIGURE_DIR = REPO / "docs/validation/figures/simplemos_m13"
FIGURES = [
    "m13_drive_discretization_p95.png",
    "m13_barrier_current_attribution.png",
    "m13_surface_profile_n21_vd0p05.png",
    "m13_old_slotboom_replay.png",
    "m13_mobility_stage_chain.png",
]
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m13_spatial_attribution_contract_v1.json",
    "scripts/run_simplemos_m13_spatial_attribution.py",
    "scripts/plot_simplemos_m13_spatial_attribution.py",
    "scripts/freeze_simplemos_m13_spatial_attribution_evidence.py",
    "tests/regression/test_simplemos_m13_spatial_attribution.py",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--output", type=Path, default=EVIDENCE)
    args = parser.parse_args()
    report_path = args.report.resolve()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    findings = report["findings"]
    if report["status"] != "complete" or report["execution"]["state_count"] != 16:
        raise ValueError("M13 report is not a complete 16-state matrix")
    if findings["sentaurus_drive_improves_active_p95_state_count"] != 16:
        raise ValueError("M13 Sentaurus-drive replay did not improve all states")
    if findings["maximum_old_slotboom_p95_difference_meV"] >= 1.0e-10:
        raise ValueError("M13 OldSlotboom replay is not numerically closed")
    artifacts = dict(report["artifacts"])
    artifacts["report"] = {"path": portable(report_path),
                           "sha256": sha256(report_path)}
    figures = [{"path": portable(FIGURE_DIR / name),
                "sha256": sha256(FIGURE_DIR / name)} for name in FIGURES]
    evidence = {
        "schema": "vela.simplemos.sdevice.m13_spatial_attribution_evidence.v1",
        "status": "complete",
        "scope": report["scope"],
        "sentaurus_release": "T-2022.03-SP2",
        "execution": report["execution"],
        "findings": findings,
        "conclusions": {
            "supported": [
                "Sentaurus GradQuasiFermi drive replay improves active-edge mobility P95 in all 16 states.",
                "Reference-aware edge projection improves all eight low-drain states but worsens all eight high-drain states.",
                "The source-barrier proxy has the correct current-error sign in all three high-NWell low-drain representative weak-current states, but does not close their magnitude.",
                "Vela OldSlotboom evaluated on total ionized impurity reproduces the Sentaurus BandgapNarrowing field to floating-point precision.",
                "Surface electrostatic state differences reach a P95 of more than 16 mV and remain a viable self-consistent residual source."
            ],
            "not_supported": [
                "Replacing transport_cell_vector globally with edge_projection.",
                "Attributing the residual platform to an OldSlotboom formula mismatch.",
                "Treating the source-barrier proxy as a complete quantitative terminal-current model.",
                "Inferring a default-model change from these read-only diagnostics."
            ]
        },
        "diagnostic_scope": {
            "edge_projection": "contact-basin-reference-aware drive reconstructed from M10 SG endpoint quasi-Fermi potentials",
            "barrier_proxy": "maximum(phin-psi) on source-side exact Silicon/Oxide interface nodes",
            "bgn": "production total-ionized-impurity OldSlotboom replay"
        },
        "source_integrity": report["inputs"],
        "artifacts": artifacts,
        "figures": figures,
        "implementation_sha256": {
            relative: sha256(REPO / relative) for relative in IMPLEMENTATION
        },
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps({"status": "frozen", "output": portable(output),
                      "state_count": report["execution"]["state_count"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
