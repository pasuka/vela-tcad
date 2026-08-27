import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "analyze_templates_ldmos_g3_min_transition.py"


def load_module():
    spec = importlib.util.spec_from_file_location("ldmos_g3_alignment", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TemplatesLdmosG3MinTransitionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = load_module()

    def test_distribution_and_signed_ratio_preserve_diagnostic_semantics(self) -> None:
        stats = self.module.distribution([3.0, 1.0, 2.0, math.nan])
        self.assertEqual(stats["count"], 3)
        self.assertEqual(stats["maximum"], 3.0)
        self.assertEqual(stats["median"], 2.0)
        self.assertAlmostEqual(stats["l2"], math.sqrt(14.0))

        ratio = self.module.signed_ratio(-2.0, -4.0)
        self.assertEqual(ratio["signed_ratio"], 0.5)
        self.assertAlmostEqual(ratio["magnitude_error_dex"], math.log10(2.0))

    def test_state_delta_uses_voltage_and_density_metrics(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_ldmos_g3_alignment_") as tmp:
            root = Path(tmp)
            left = root / "left.csv"
            right = root / "right.csv"
            header = "node_id,psi,phin,phip,electrons_m3,holes_m3\n"
            left.write_text(
                header + "0,0,0,0,1e20,1e10\n1,1,2,3,1e22,1e12\n",
                encoding="utf-8",
            )
            right.write_text(
                header + "0,0.1,0.2,0.3,1e21,1e9\n1,1.4,2.5,3.6,1e20,1e13\n",
                encoding="utf-8",
            )
            result = self.module.state_delta(left, right)

        self.assertEqual(result["common_nodes"], 2)
        self.assertAlmostEqual(result["psi_absolute_difference_V"]["maximum"], 0.4)
        self.assertAlmostEqual(result["phip_absolute_difference_V"]["maximum"], 0.6)
        self.assertEqual(
            result["electrons_m3_absolute_difference_dex"]["maximum"], 2.0
        )
        self.assertEqual(
            result["holes_m3_absolute_difference_dex"]["maximum"], 1.0
        )

    def test_state_delta_can_be_restricted_to_silicon_nodes(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_ldmos_g3_alignment_") as tmp:
            root = Path(tmp)
            header = "node_id,psi,phin,phip,electrons_m3,holes_m3\n"
            (root / "left.csv").write_text(
                header + "0,0,0,0,1e20,1e10\n1,1,2,3,1e22,1e12\n",
                encoding="utf-8",
            )
            (root / "right.csv").write_text(
                header + "0,9,9,9,1e2,1e30\n1,1.1,2.2,3.3,1e21,1e13\n",
                encoding="utf-8",
            )
            result = self.module.state_delta(
                root / "left.csv", root / "right.csv", {1}
            )

        self.assertEqual(result["common_nodes"], 1)
        self.assertAlmostEqual(result["psi_absolute_difference_V"]["maximum"], 0.1)
        self.assertEqual(
            result["electrons_m3_absolute_difference_dex"]["maximum"], 1.0
        )

    def test_redirect_csv_files_keeps_control_artifacts_together(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_ldmos_g3_redirect_") as tmp:
            root = Path(tmp)
            config = {
                "write_state_file": "old.csv",
                "diagnostics": {
                    "terminal": {"csv_file": "terminal.csv"},
                    "audit": {"summary_file": "audit.json"},
                },
            }
            self.module.redirect_csv_files(config, root, "sweep")

        self.assertTrue(config["write_state_file"].endswith("sweep_write_state_file.csv"))
        self.assertTrue(
            config["diagnostics"]["terminal"]["csv_file"].endswith(
                "sweep_diagnostics_terminal_csv_file.csv"
            )
        )
        self.assertTrue(
            config["diagnostics"]["audit"]["summary_file"].endswith(
                "sweep_diagnostics_audit_summary_file.json"
            )
        )

    def test_sentaurus_newton_export_reports_native_units_and_locations(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_ldmos_g3_newton_") as tmp:
            root = Path(tmp)
            fields = root / "fields"
            fields.mkdir()
            (root / "nodes.csv").write_text(
                "id,x_um,y_um\n0,1,2\n1,3,4\n", encoding="utf-8"
            )
            data = {
                "PoissonRhs": (1.0, -3.0),
                "eContinuityRhs": (2.0, -4.0),
                "hContinuityRhs": (3.0, -5.0),
                "NewtonStepElectrostaticPotentialUpdate": (0.1, -0.3),
                "NewtonStepEDensityUpdate": (20.0, -40.0),
                "NewtonStepHDensityUpdate": (30.0, -50.0),
            }
            for name, values in data.items():
                (fields / f"{name}_region0.csv").write_text(
                    "node_id,component0\n"
                    f"0,{values[0]}\n1,{values[1]}\n",
                    encoding="utf-8",
                )
            result = self.module.sentaurus_newton_export(root)

        self.assertEqual(result["electron_rhs"]["unit"], "A")
        self.assertEqual(result["electron_rhs"]["top"]["node_id"], 1)
        self.assertEqual(result["electron_rhs"]["top"]["x_um"], 3.0)
        self.assertEqual(result["psi_update"]["absolute"]["maximum"], 0.3)

    def test_constant_high_field_control_uses_tcad_internal_values(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_ldmos_g3_hfs_") as tmp:
            root = Path(tmp)
            baseline = root / "baseline.json"
            output = root / "control.json"
            baseline.write_text(
                json.dumps({"solver": {"mobility": {"model": "masetti_field"}}}),
                encoding="utf-8",
            )
            self.module.constant_high_field_control(baseline, output)
            mobility = json.loads(output.read_text(encoding="utf-8"))["solver"][
                "mobility"
            ]

        self.assertEqual(mobility["model"], "caughey_thomas_field")
        self.assertEqual(mobility["electron_mu_min_m2_V_s"], 1417.0)
        self.assertEqual(mobility["electron_saturation_velocity_m_s"], 1.07e7)
        self.assertEqual(mobility["hole_mu_min_m2_V_s"], 470.5)

    def test_newton_state_delta_separates_predictor_state_fields(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_ldmos_g3_predictor_") as tmp:
            root = Path(tmp)
            for label, values in (
                ("left", ((0.0, 1.0e10, 1.0e8), (1.0, 1.0e12, 1.0e6))),
                ("right", ((0.1, 1.0e11, 1.0e9), (1.4, 1.0e10, 1.0e7))),
            ):
                fields = root / label / "fields"
                fields.mkdir(parents=True)
                for name, index in (
                    ("ElectrostaticPotential", 0),
                    ("eDensity", 1),
                    ("hDensity", 2),
                ):
                    (fields / f"{name}_region0.csv").write_text(
                        "node_id,component0\n"
                        + "".join(
                            f"{node},{row[index]}\n"
                            for node, row in enumerate(values)
                        ),
                        encoding="utf-8",
                    )
            result = self.module.sentaurus_newton_state_delta(
                root / "left", root / "right"
            )

        self.assertEqual(result["common_silicon_nodes"], 2)
        self.assertAlmostEqual(result["psi_absolute_difference_V"]["maximum"], 0.4)
        self.assertEqual(
            result["electron_density_absolute_difference_dex"]["maximum"], 2.0
        )


if __name__ == "__main__":
    unittest.main()
