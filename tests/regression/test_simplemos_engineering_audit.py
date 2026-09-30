"""Coverage failures must not silently turn into a passing validation matrix."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from audit_simplemos_engineering_20260929 import index_points


class CoverageTest(unittest.TestCase):
    def setUp(self):
        self.keys={('n23',.05,i) for i in range(2)}
        self.rows=[dict(device='n23',vd=.05,vg=i*.05) for i in range(2)]

    def test_valid_exact_grid(self):
        self.assertEqual(set(index_points(self.rows,self.keys)),self.keys)

    def test_missing_point_is_error(self):
        with self.assertRaises(ValueError):index_points(self.rows[:1],self.keys)

    def test_duplicate_is_error(self):
        with self.assertRaises(ValueError):index_points(self.rows+self.rows[:1],self.keys)

    def test_wrong_device_is_error(self):
        with self.assertRaises(ValueError):index_points([dict(r,device='n24') for r in self.rows],self.keys)

    def test_off_grid_is_error(self):
        with self.assertRaises(ValueError):index_points([dict(self.rows[0],vg=.0501)],self.keys)

    def test_nonfinite_is_error(self):
        for value in (float('nan'),float('inf')):
            with self.assertRaises(ValueError):index_points([dict(self.rows[0],vg=value)],self.keys)


if __name__=='__main__':unittest.main()
