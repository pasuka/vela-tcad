#!/usr/bin/env python3
"""Regression coverage for the corrected TransportModels material contract."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
FIXTURE = REPO / "reference_tcad" / "transportmodels_sentaurus2022" / "vela"
sys.path.insert(0, str(SCRIPTS))


def load_script(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PARAMETER_SWEEP = load_script(
    "transportmodels_dg_parameter_sweep_test",
    "run_transportmodels_dg_parameter_sweep.py",
)
PHASE7 = load_script(
    "transportmodels_dg_phase7_test",
    "run_transportmodels_dg_phase7_regression.py",
)
DEEP_OFF = load_script(
    "transportmodels_dd_deep_off_self_consistent_test",
    "run_transportmodels_dd_deep_off_self_consistent_srh_fix.py",
)


class TransportModelsSentaurusMaterialContractTest(unittest.TestCase):
    def test_corrected_materials_use_sentaurus_silicon_intrinsic_density(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "materials.json"
            baseline = Path(directory) / "baseline_materials.json"
            source = json.loads((FIXTURE / "materials_sentaurus2022.json").read_text(encoding="utf-8"))
            for material in source["materials"]:
                if material["name"] in {"Si", "PolySilicon"}:
                    material["ni"] = 1.0e10
            baseline.write_text(json.dumps(source), encoding="utf-8")
            with mock.patch.object(PARAMETER_SWEEP, "BASE_MATERIALS", baseline):
                PARAMETER_SWEEP.write_corrected_materials(output)
            payload = json.loads(output.read_text(encoding="utf-8"))
            materials = {row["name"]: row for row in payload["materials"]}

            expected = PARAMETER_SWEEP.SENTAURUS_SILICON_NI_CM3
            self.assertEqual(expected, materials["Si"]["ni"])
            self.assertEqual(expected, materials["PolySilicon"]["ni"])

    def test_phase7_enables_fermi_statistics_bgn_correction(self) -> None:
        curve = dict(PHASE7.CURVES[0])
        curve["points"] = [-1.0]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = json.loads((FIXTURE / "configs" / "09_dg_idvg_curve.json").read_text(encoding="utf-8"))
            source["solver"]["bandgap_narrowing"] = "old_slotboom"
            for carrier in ("electron", "hole"):
                source["solver"]["srh_doping_dependence"][carrier]["reference_doping_m3"] = 1.0e22
            baseline = root / "baseline_config.json"
            baseline.write_text(json.dumps(source), encoding="utf-8")
            curve["config"] = baseline
            curve["initial_state"] = root / "initial_state.csv"
            with mock.patch.object(PHASE7, "OUTPUT_ROOT", root / "output"), mock.patch.object(
                PHASE7, "CORRECTED_MATERIALS", FIXTURE / "materials_sentaurus2022.json"
            ):
                path = PHASE7.make_config(curve)

            config = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(
                {
                    "model": "old_slotboom",
                    "fermi_statistics_correction": True,
                },
                config["solver"]["bandgap_narrowing"],
            )
            self.assertEqual(
                1.0e16,
                config["solver"]["srh_doping_dependence"]["electron"][
                    "reference_doping_m3"
                ],
            )
            self.assertEqual(
                1.0e16,
                config["solver"]["srh_doping_dependence"]["hole"][
                    "reference_doping_m3"
                ],
            )

    def test_deep_off_config_combines_all_three_srh_corrections(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            output = root / "output"
            materials = root / "materials.json"
            source.mkdir()
            materials.write_text('{"materials": []}\n', encoding="utf-8")
            (source / "final_state.csv").write_text("state\n", encoding="utf-8")
            (source / "config.json").write_text(
                json.dumps(
                    {
                        "materials_file": "old.json",
                        "solver": {
                            "bandgap_narrowing": "old_slotboom",
                            "srh_doping_dependence": {
                                "electron": {"reference_doping_m3": 1.0e16},
                                "hole": {"reference_doping_m3": 1.0e16},
                            },
                        },
                        "sweep": {
                            "initial_state_file": "old_state.csv",
                            "diagnostics": {
                                "terminal_balance": {
                                    "csv_file": str(source / "terminal.csv")
                                }
                            },
                        },
                        "output_csv": str(source / "curve.csv"),
                    }
                ),
                encoding="utf-8",
            )
            old_source, old_output, old_materials = (
                DEEP_OFF.SOURCE_RUN,
                DEEP_OFF.OUTPUT,
                DEEP_OFF.MATERIALS,
            )
            try:
                DEEP_OFF.SOURCE_RUN = source
                DEEP_OFF.OUTPUT = output
                DEEP_OFF.MATERIALS = materials
                path = DEEP_OFF.make_config()
            finally:
                DEEP_OFF.SOURCE_RUN = old_source
                DEEP_OFF.OUTPUT = old_output
                DEEP_OFF.MATERIALS = old_materials

            config = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(str(materials.resolve()), config["materials_file"])
            self.assertEqual(
                {"model": "old_slotboom", "fermi_statistics_correction": True},
                config["solver"]["bandgap_narrowing"],
            )
            self.assertEqual(2.0e-11, config["solver"]["stall_residual_floor"])
            self.assertEqual(
                1.0e16,
                config["solver"]["srh_doping_dependence"]["electron"][
                    "reference_doping_m3"
                ],
            )
            self.assertEqual(str(output / "curve.csv"), config["output_csv"])
            self.assertEqual(
                str(source / "final_state.csv"),
                config["sweep"]["initial_state_file"],
            )


if __name__ == "__main__":
    unittest.main()
