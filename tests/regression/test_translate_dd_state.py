import csv
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from translate_dd_state import translate_state_csv


class FrameTranslationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.row = dict(node_id="0", psi="0.45", phin="2.24e-17", phip="1.89e-11",
                        electron_qf_reference_V="2.24e-17", electron_qf_increment_V="-1.31e-29",
                        hole_qf_reference_V="1.89e-11", hole_qf_increment_V="4.93e-21",
                        electrons_m3="1e24", holes_m3="600000000", electron_quantum_potential_V="0.001")

    def write(self, row):
        p = self.root / "input.csv"
        with p.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=row.keys()); w.writeheader(); w.writerow(row)
        return p

    def read(self, p):
        with p.open(newline="") as f:
            return next(csv.DictReader(f))

    def test_preserves_split_qf_and_physical_bias(self):
        source = self.write(self.row); result = self.root / "shifted.csv"
        translate_state_csv(source, result, 28., {0})
        row = self.read(result)
        with localcontext() as ctx:
            ctx.prec = 80
            for who in ("electron", "hole"):
                ref, inc = who + "_qf_reference_V", who + "_qf_increment_V"
                exact = Decimal.from_float(float(self.row[ref])) + Decimal.from_float(float(self.row[inc])) - 28
                expected_inc = float(exact - Decimal.from_float(float(row[ref])))
                self.assertEqual(float(row[inc]), expected_inc)
                self.assertNotEqual(float(row[inc]), float(self.row[inc]))
        for key in ("psi", "phin", "phip"):
            self.assertEqual(float(row[key]), float(self.row[key]) - 28.)
        for key in ("electrons_m3", "holes_m3", "electron_quantum_potential_V"):
            self.assertEqual(row[key], self.row[key])
        self.assertEqual(self.read(source), self.row)

    def test_round_trip_preserves_low_bits(self):
        source = self.write(self.row); shifted = self.root / "shifted.csv"; back = self.root / "back.csv"
        translate_state_csv(source, shifted, 28., {0}); translate_state_csv(shifted, back, -28., {0})
        row = self.read(back)
        for who in ("electron", "hole"):
            ref, inc = who + "_qf_reference_V", who + "_qf_increment_V"
            self.assertAlmostEqual(float(row[ref]) + float(row[inc]),
                                   float(self.row[ref]) + float(self.row[inc]), delta=1e-30)

    def test_rounding_tie_uses_exact_binary_sum(self):
        # Actual LDMOS node180: finite decimal precision can round a binary
        # halfway case in the wrong direction after cancellation of 28 V.
        self.row['electron_qf_reference_V'] = '0.027249160822519092'
        self.row['electron_qf_increment_V'] = '-4.5881879959268769e-17'
        out = self.root / 'shifted.csv'
        translate_state_csv(self.write(self.row), out, 28., {0})
        row = self.read(out)
        exact = (Fraction.from_float(float(self.row['electron_qf_reference_V']))
                 + Fraction.from_float(float(self.row['electron_qf_increment_V'])) - 28
                 - Fraction.from_float(float(row['electron_qf_reference_V'])))
        self.assertEqual(float(row['electron_qf_increment_V']), float(exact))

    def test_inactive_qf_stays_zero(self):
        out = self.root / "shifted.csv"
        translate_state_csv(self.write(self.row), out, 28., set())
        row = self.read(out)
        for k in ("phin", "phip", "electron_qf_reference_V", "hole_qf_reference_V",
                  "electron_qf_increment_V", "hole_qf_increment_V"):
            self.assertEqual(float(row[k]), 0.)
        self.assertEqual(float(row["psi"]), float(self.row["psi"]) - 28.)

    def test_refuses_overwrite(self):
        source = self.write(self.row)
        with self.assertRaises(FileExistsError):
            translate_state_csv(source, source, 28., {0})
        self.assertEqual(self.read(source), self.row)

    def test_rejects_incomplete_reference_pair(self):
        row = dict(self.row); del row["electron_qf_increment_V"]
        with self.assertRaises(ValueError):
            translate_state_csv(self.write(row), self.root / "out.csv", 28., {0})

    def test_rejects_nonfinite_shift(self):
        with self.assertRaises(ValueError):
            translate_state_csv(self.write(self.row), self.root / "out.csv", float("nan"), {0})


if __name__ == "__main__":
    unittest.main()
