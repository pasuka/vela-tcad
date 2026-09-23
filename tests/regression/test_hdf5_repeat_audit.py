import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_templates_ldmos_hdf5_repeats import validate_batch


class RepeatAuditTest(unittest.TestCase):
    def fixture(self):
        cases=[dict(repeat=r,mode=m,gate=g,status='completed',trajectory_exact=True,bias_V=40.,
            audit=dict(integrity_pass=True,exact_points=31),verdicts={'curve':{'pass_all':True}})
            for r in range(3) for m in ('binary','hdf5') for g in (4,8)]
        joint=[dict(repeat=r,mode=m,checks={'ratio':True}) for r in range(3) for m in ('binary','hdf5')]
        return dict(status='completed',points=31,control='binary',candidate='hdf5',planned_curves=12,cases=cases,joint=joint)
    def test_complete(self):self.assertEqual(len(validate_batch(self.fixture())),12)
    def test_missing_duplicate_and_partial(self):
        for field in ('cases','joint'):
            b=self.fixture();b[field].pop()
            with self.assertRaises(ValueError):validate_batch(b)
            b=self.fixture();b[field][-1]=copy.deepcopy(b[field][0])
            with self.assertRaises(ValueError):validate_batch(b)
        b=self.fixture();b['points']=8
        with self.assertRaises(ValueError):validate_batch(b)
    def test_failed_gates_and_unfinished_state(self):
        b=self.fixture();b['joint'][0]['checks']['ratio']=False
        with self.assertRaises(ValueError):validate_batch(b)
        b=self.fixture();b['cases'][0]['verdicts']['curve']['pass_all']=False
        with self.assertRaises(ValueError):validate_batch(b)
        b=self.fixture();b['cases'][0]['bias_V']=38.
        with self.assertRaises(ValueError):validate_batch(b)

if __name__=='__main__':unittest.main()
