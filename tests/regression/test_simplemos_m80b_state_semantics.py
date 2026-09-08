"""Validate the source constant inference is kept out of solver parameters."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
import run_simplemos_m80b_state_semantics as m


class StateSemantics(unittest.TestCase):
    def test_input_identity_is_consistent_across_both_carriers_and_devices(self):
        rows=m.audit_rows()
        self.assertEqual(len(rows),16)
        vt=[r["source_inferred_thermal_voltage_V"] for r in rows]
        self.assertLess(max(vt)-min(vt),1e-14)
        self.assertLess(max(r["source_identity_max_error_V"] for r in rows),1e-12)
        self.assertGreater(max(r["raw_density_identity_max_dex"] for r in rows),1e-6)


if __name__=="__main__":unittest.main()
