#!/usr/bin/env python3
"""Regression coverage for failed Newton reports in diagnose_pn2d_0v_current_balance.py."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"


class DiagnosePn2d0vCurrentBalanceFailureTest(unittest.TestCase):
    def test_line_search_failure_report_is_classified(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_current_balance_failure_") as tmp:
            root = Path(tmp)
            reference = root / "reference"
            vela = reference / "vela"
            vela.mkdir(parents=True)
            (vela / "simulation_0v.json").write_text(json.dumps({
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "solver": {
                    "method": "gummel_newton",
                    "line_search": True,
                    "diagnostics": False,
                },
                "sweep": {
                    "contact": "Anode",
                    "current_contact": "Anode",
                    "start": 0.0,
                    "stop": 0.0,
                    "step": 1.0,
                },
            }, indent=2) + "\n")

            fake_runner = root / "fake_runner.py"
            fake_runner.write_text(textwrap.dedent(r'''
                import csv
                import json
                import sys
                from pathlib import Path

                config_path = Path(sys.argv[sys.argv.index("--config") + 1])
                deck = json.loads(config_path.read_text())
                output_csv = Path(deck["output_csv"])
                output_csv.parent.mkdir(parents=True, exist_ok=True)
                failure_json = output_csv.with_name(output_csv.stem + "_newton_failure_diagnostics.json")
                header = [
                    "mode", "bias_contact", "bias_V", "current_contact",
                    "current_electron", "current_electron_drift", "current_electron_diffusion",
                    "current_hole", "current_hole_drift", "current_hole_diffusion",
                    "current_total", "converged", "iterations", "solver_method",
                    "gummel_iterations", "newton_iterations", "handoff_stage",
                    "step_diagnostics", "validation_diagnostics", "failure_reason",
                    "newton_failure_class", "newton_failure_diagnostics_json",
                ]
                row = [
                    "iv", "Anode", "0", "Anode",
                    "0", "0", "0", "0", "0", "0", "0",
                    "0", "11", "gummel_newton", "8", "11", "newton_failed",
                    "attempted_step=0;accepted_step=0;retry_count=0", "",
                    "line_search_non_decrease", "line_search_non_decrease", str(failure_json),
                ]
                with output_csv.open("w", newline="") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(header)
                    writer.writerow(row)
                failure_json.write_text(json.dumps([{
                    "point_index": 0,
                    "bias_V": 0.0,
                    "bias_contact": "Anode",
                    "solver_method": "gummel_newton",
                    "handoff_stage": "newton_failed",
                    "failure_reason": "line_search_non_decrease",
                    "failed_iteration": 12,
                    "residual_norm": 0.0671994,
                    "step_norm": 1.0,
                    "damping_factor": 0.0,
                    "line_search_attempts": 12,
                    "line_search_failure_reason": "line_search_non_decrease",
                    "block_residuals": {"psi": 17.0661, "phin": 1.0e-9, "phip": 2.0e-9, "combined": 17.0661},
                    "carrier_diagnostics": {"positive_finite": True},
                    "line_search_history": [{
                        "attempt": 11,
                        "damping": 0.00048828125,
                        "residual_norm": 0.0671995,
                        "target_residual_norm": 0.0671994,
                        "finite": True,
                        "carrier_positive_finite": True,
                        "sufficient_decrease": False,
                        "accepted": False,
                        "rejection_reason": "line_search_non_decrease",
                    }],
                    "top_poisson_residual_nodes": [{
                        "node_id": 33,
                        "x_um": 1.0,
                        "y_um": 0.5,
                        "abs_poisson_residual": 17.0661,
                        "donors_cm3": 1.0e17,
                        "acceptors_cm3": 1.0e17,
                        "net_doping_cm3": 0.0,
                        "ni_eff_cm3": 1.0e10,
                    }],
                }], indent=2) + "\n")
                sys.exit(1)
            '''))

            reports = root / "reports"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--reference-root",
                    str(reference),
                    "--runner",
                    f"{sys.executable} {fake_runner}",
                    "--output-dir",
                    str(reports),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads((reports / "pn2d_0v_current_balance.json").read_text())
            self.assertEqual(report["status"], "diagnostic_fail")
            self.assertEqual(report["classification"], "line_search_non_decrease")
            failure = report["newton_failure"]
            self.assertEqual(failure["failure_class"], "line_search_non_decrease")
            latest = failure["diagnostics"][-1]
            self.assertEqual(latest["line_search_failure_reason"], "line_search_non_decrease")
            self.assertEqual(latest["top_poisson_residual_nodes"][0]["x_um"], 1.0)
            markdown = (reports / "pn2d_0v_current_balance.md").read_text()
            self.assertIn("line_search_non_decrease", markdown)
            self.assertIn("Top Poisson Residual Nodes", markdown)


if __name__ == "__main__":
    unittest.main()
