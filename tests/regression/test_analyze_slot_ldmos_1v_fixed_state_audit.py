import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "analyze_slot_fixed",
    REPO / "scripts" / "analyze_slot_ldmos_1v_fixed_state_audit.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class AnalyzeSlotLdmosFixedStateAuditTest(unittest.TestCase):
    def test_summarizes_current_recombination_and_state_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fields = root / "import" / "fields"
            fields.mkdir(parents=True)
            manifest = {"fields": []}
            for index, (contact, current) in enumerate(
                    (("source", -1.0), ("drain", 2.0), ("substrate", -1.0))):
                name = f"contact_{index}.csv"
                _write_csv(fields / name, [{"node_id": index, "component0": current}])
                manifest["fields"].append({
                    "name": "ContactCurrentFlux",
                    "region_name": contact,
                    "csv_file": name,
                })
            (root / "import" / "field_manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8")

            terminal = root / "terminal.csv"
            _write_csv(terminal, [{
                "contact": "drain",
                "current_electron_A_per_um": 3.0,
                "current_hole_A_per_um": -0.1,
                "current_total_A_per_um": 4.0,
            }])
            recombination = root / "recombination.csv"
            _write_csv(recombination, [{
                "srh_net_current_A_per_um": -1.0,
                "srh_generation_current_A_per_um": 1.0,
                "srh_recombination_current_A_per_um": 0.0,
                "auger_net_current_A_per_um": -0.01,
                "auger_generation_current_A_per_um": 0.01,
                "auger_recombination_current_A_per_um": 0.0,
                "total_recombination_net_current_A_per_um": -1.01,
            }])
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            rows = [{
                "node_id": 0, "psi": 0.0, "phin": 0.0, "phip": 0.0,
                "electrons_m3": 1.0e20, "holes_m3": 1.0e10,
            }]
            _write_csv(reference, rows)
            changed = [dict(rows[0])]
            changed[0]["psi"] = 0.2
            changed[0]["electrons_m3"] = 1.0e22
            _write_csv(candidate, changed)

            result = MODULE.analyze(
                root / "import", terminal, recombination, reference, candidate)
            self.assertAlmostEqual(
                result["drain_current_comparison"]
                    ["absolute_log10_error_decade"],
                0.3010299956639812)
            self.assertFalse(result["gate"]["passed"])
            self.assertEqual(
                result["vela_fixed_state_recombination_integrals_A_per_um"]
                    ["auger_net_current_A_per_um"],
                -0.01)
            state = result[
                "sentaurus_vs_existing_vela_self_consistent_state"]
            self.assertAlmostEqual(
                state["absolute_voltage_error_V"]["psi"]["maximum"], 0.2)
            self.assertAlmostEqual(
                state["absolute_log10_density_ratio"]["electrons_m3"]
                    ["maximum_decade"],
                2.0)


if __name__ == "__main__":
    unittest.main()
