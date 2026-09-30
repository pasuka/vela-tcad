"""Independent acceptance boundaries for full split-state comparisons."""
import copy
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import run_simplemos_engineering_matrix_20260929 as matrix
import supervise_simplemos_engineering_20260929 as supervisor


class DualGateTest(unittest.TestCase):
    def setUp(self):
        fields={k:np.array([0.]) for k in ('packed_psi','packed_psi_low','packed_electron_qf_increment',
            'packed_electron_qf_increment_low','packed_hole_qf_increment','packed_hole_qf_increment_low',
            'electron_qf_reference_V','hole_qf_reference_V')}
        fields.update(electrons_m3=np.array([1e10]),holes_m3=np.array([1e16]))
        self.a=(fields,dict(mesh_sha256='same',packed_potential_scale_V=.025))
        self.b=copy.deepcopy(self.a)
        self.cold=dict(qualified=True,current_A_per_um=1e-16)
        self.warm=copy.deepcopy(self.cold)

    def compare(self):
        with patch.object(matrix,'load',side_effect=[self.a,self.b]):
            return matrix.compare_pair('a','b',self.cold,self.warm,dict(free_si=[0]))

    def test_identical(self):self.assertTrue(self.compare()['qualified'])

    def test_microvolt_boundary(self):
        self.b[0]['packed_psi'][0]=1.1e-6/.025
        self.assertFalse(self.compare()['qualified'])

    def test_density_boundary(self):
        self.b[0]['electrons_m3'][0]*=1.001
        self.assertFalse(self.compare()['qualified'])

    def test_current_boundary(self):
        self.warm['current_A_per_um']*=1.00001
        self.assertFalse(self.compare()['qualified'])

    def test_unqualified_state_is_not_rescued_by_agreement(self):
        self.warm['qualified']=False
        self.assertFalse(self.compare()['qualified'])

    def test_nan_is_not_a_pass(self):
        self.warm['current_A_per_um']=float('nan')
        self.assertFalse(self.compare()['qualified'])

    def test_mesh_mismatch(self):
        self.b[1]['mesh_sha256']='other'
        with self.assertRaises(ValueError):self.compare()


class CompletedCaseTest(unittest.TestCase):
    def check(self,rows,passed=True):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d)
            supervisor.h.write(base/'matrix/n17_vd_0p05/summary.json',dict(points=51,passed=passed,results=rows))
            return supervisor.verified_summary(base,'n17_vd_0p05')

    def setUp(self):
        self.rows=[dict(case='n17_vd_0p05',index=i,passed=True) for i in range(51)]

    def test_complete_case(self):self.assertEqual(len(self.check(self.rows)),51)

    def test_duplicate_cannot_hide_missing_index(self):
        self.rows[-1]['index']=0
        with self.assertRaises(ValueError):self.check(self.rows)

    def test_failed_point_cannot_hide_under_passed_summary(self):
        self.rows[-1]['passed']=False
        with self.assertRaises(ValueError):self.check(self.rows)

    def test_wrong_case(self):
        self.rows[-1]['case']='n23_vd_0p05'
        with self.assertRaises(ValueError):self.check(self.rows)


if __name__=='__main__':unittest.main()
