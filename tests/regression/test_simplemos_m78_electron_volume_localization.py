"""Physical accounting and evidence guards for the read-only M78 observer."""
import csv
import math
from pathlib import Path
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m78_electron_volume_localization as m78


def rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


class M78Accounting(unittest.TestCase):
    def test_current_pair_closure_does_not_claim_absolute_improvement(self):
        upstream = {(r["device"], float(r["drain_voltage_V"])): r for r in rows(m78.m74.CASES)}
        pairs = rows(m78.OUT / "m78_pair_ledger.csv")
        self.assertEqual(len(pairs), 8)
        for pair in pairs:
            low, high = [upstream[(pair[k], float(pair["drain_voltage_V"]))]
                         for k in ("low_device", "high_device")]
            shifts = [math.log10(float(r["candidate_vela_current_A_per_um"])/float(r["baseline_vela_current_A_per_um"])) for r in (low, high)]
            self.assertAlmostEqual(float(pair["pair_reduction_dex"]), shifts[0]-shifts[1], delta=1e-13)
            self.assertGreater(shifts[0], shifts[1])
            for r in (low, high):
                self.assertGreater(float(r["candidate_error_dex"]), float(r["baseline_error_dex"]))
                self.assertGreater(float(r["baseline_error_dex"]), 0)

    def test_stage_increments_telescope_and_charge_terms_close(self):
        stage_rows = rows(m78.OUT / "m78_stage_ledger.csv")
        keys = {(r["device"], r["drain_voltage_V"]) for r in stage_rows}
        self.assertEqual((len(keys), len(stage_rows)), (16, 48))
        for key in keys:
            stages = {r["stage"]: r for r in stage_rows if (r["device"], r["drain_voltage_V"]) == key}
            for field in m78.FIELDS:
                total = sum(float(r[f"control_increment_{field}"]) for r in stages.values())
                self.assertAlmostEqual(total, float(stages["gate"][f"control_delta_{field}"]), delta=1e-13)
            for r in stages.values():
                total = sum(float(r[f"control_{term}_response_V"]) for term in ("forcing", "electron_relaxation", "hole_relaxation"))
                self.assertAlmostEqual(total, float(r["control_delta_psi"]), delta=2e-9)

    def test_spatial_partition_is_conservative_at_all_depths(self):
        spatial = rows(m78.OUT / "m78_partition_sensitivity_ledger.csv")
        groups = {}
        for r in spatial:
            key = (r["device"], r["drain_voltage_V"], r["stage"], r["depth_um"])
            groups.setdefault(key, []).append(r)
        self.assertEqual(len(groups), 16*3*3)
        for group in groups.values():
            self.assertEqual({r["region"] for r in group}, {"source", "channel", "drain", "substrate"})
            self.assertAlmostEqual(sum(float(r["increment_psi_l1_share"]) for r in group), 1.0, delta=1e-12)

    def test_changed_or_missing_input_is_rejected_without_rewrite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state.csv"
            path.write_bytes(b"frozen state\n")
            expected = m78.sha256(path)
            path.write_bytes(b"changed state\n")
            with self.assertRaisesRegex(ValueError, "do not recompute"):
                m78.check_hash(path, expected)
            self.assertEqual(path.read_bytes(), b"changed state\n")
            path.unlink()
            with self.assertRaisesRegex(ValueError, "do not recompute"):
                m78.check_hash(path, expected)
            self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main()
