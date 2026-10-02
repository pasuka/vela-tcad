"""The cloud assignment must never overlap retained T470p work."""
import copy
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from run_simplemos_engineering_codespaces_20260929 import validate_assignment,checked_path


class OwnershipTest(unittest.TestCase):
    def setUp(self):
        self.contract=dict(cases=[dict(case=f'n{n}_vd_{v}') for n in range(17,25) for v in ('0p05','1')])
        keep={'n17_vd_0p05':'done','n17_vd_1':'failed','n18_vd_0p05':'running'}
        self.a=dict(cases=[c['case'] for c in self.contract['cases'] if c['case'] not in keep],retained_T470p=keep,points=663,max_workers=2)
    def test_exact_partition(self):self.assertEqual(len(validate_assignment(self.a,self.contract)),13)
    def test_overlap_rejected(self):
        self.a['cases'][0]='n18_vd_0p05'
        with self.assertRaises(ValueError):validate_assignment(self.a,self.contract)
    def test_duplicate_rejected(self):
        self.a['cases'][0]=self.a['cases'][1]
        with self.assertRaises(ValueError):validate_assignment(self.a,self.contract)
    def test_unknown_case_rejected(self):
        self.a['cases'][0]='n25_vd_1'
        with self.assertRaises(ValueError):validate_assignment(self.a,self.contract)
    def test_path_escape_rejected(self):
        with self.assertRaises(ValueError):checked_path(Path.cwd(),'../escape')


if __name__=='__main__':unittest.main()
