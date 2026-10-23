"""Physical and topology checks for read-only native box evidence."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from bvmethods_native_box import (read_box, dielectric_entries, node_measures,
                                 box_measures_from_coefficients)

TEXT = """Measure { # um^2
0 8 2 .125 .25 .125
1 9 2 .25 .125 .125
2 -1 1 # boundary, no bulk measure
}
Coefficients { # dimensionless
0 8 2 .5 .5 0
1 9 2 .5 0 .5
2 -1 1
}"""
CELLS = [{"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
         {"id": 1, "region_id": 1, "node_ids": [1, 3, 2]}]


class NativeBoxTests(unittest.TestCase):
    def test_right_triangle_measure_and_coordinate_units(self):
        points = [[0., 0.], [1., 0.], [0., 1.]]
        original = box_measures_from_coefficients(points, [.5, 0., .5])
        self.assertEqual(original, (.25, .125, .125))
        scaled = box_measures_from_coefficients([[2*x, 2*y] for x, y in points],
                                               [.5, 0., .5])
        self.assertEqual(scaled, tuple(4*x for x in original))

    def test_transferred_measure_is_not_renormalized_to_cell_area(self):
        # A cell receiving positive edge support from its neighbor need not
        # retain its isolated geometric area. Source support must preserve it.
        measure = box_measures_from_coefficients([[0, 0], [1, 0], [0, 1]],
                                                [1., 0., 1.])
        self.assertEqual(sum(measure), 1.)
        self.assertNotEqual(sum(measure), .5)

    def test_rejects_invalid_box_reconstruction(self):
        for points, coefficients in [([[0, 0], [1, 0], [2, 0]], [1, 1, 1]),
                                     ([[0, 0], [1, 0], [0, 1]], [1, -1, 1]),
                                     ([[0, 0], [1, 0], [0, float('nan')]], [1, 1, 1]),
                                     ([[0, 0], [1, 0], [0, 1]], [1, 1])]:
            with self.subTest(points=points), self.assertRaises(ValueError):
                box_measures_from_coefficients(points, coefficients)

    def test_constant_potential_and_dielectric_energy(self):
        box = read_box(TEXT)
        matrix = dielectric_entries(CELLS, box["Coefficients"], {0: 11.7, 1: 3.9})
        for node in range(4):
            self.assertAlmostEqual(sum(x for (r, c), x in matrix.items() if r == node), 0)
        potential = [0., 1., 2., 3.]
        energy = sum(potential[r] * x * potential[c] for (r, c), x in matrix.items())
        self.assertGreater(energy, 0)
        for (r, c), value in matrix.items():
            self.assertEqual(value, matrix[c, r])

    def test_semiconductor_support_excludes_dielectric_at_interface(self):
        box = read_box(TEXT)
        total = node_measures(CELLS, box["Measure"])
        semiconductor = node_measures(CELLS, box["Measure"], {0})
        self.assertEqual(total[1], .5)
        self.assertEqual(semiconductor[1], .25)
        self.assertNotIn(3, semiconductor)
        self.assertEqual(sum(total.values()), 1.)

    def test_rejects_incomplete_duplicate_or_nonfinite_evidence(self):
        for text in (TEXT.replace('1 9 2 .5 0 .5', '0 9 2 .5 0 .5'),
                     TEXT.replace('.125 .25 .125', '.125 .25'),
                     TEXT.replace('.125 .25 .125', '.125 nan .125'),
                     TEXT.replace('1 9 2 .5 0 .5', '1 10 2 .5 0 .5')):
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_box(text)


if __name__ == '__main__':
    unittest.main()
