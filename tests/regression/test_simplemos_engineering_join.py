import copy
import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('join_audit',Path(__file__).resolve().parents[2]/'scripts/join_simplemos_engineering_20260930.py')
j=importlib.util.module_from_spec(spec);spec.loader.exec_module(j)

def curve():
    return dict(passed=True,points=51,results=[dict(case='n17_vd_1',index=i,vg=i*.05,passed=True,reconstruction_passed=True,Id_error_percent=dict(cold=0.,native=0.),dual=dict(qualified=True,**{k:0. for k in j.LIMITS})) for i in range(51)])

class JoinTests(unittest.TestCase):
    def test_complete_and_boundaries(self):
        c=curve();c['results'][0]['dual'].update(j.LIMITS);c['results'][0]['Id_error_percent']['cold']=-2.
        self.assertEqual(len(j.validate_curve('n17_vd_1',c)),51)
    def test_duplicate_and_missing_index(self):
        c=curve();c['results'][-1]=copy.deepcopy(c['results'][0])
        with self.assertRaises(ValueError):j.validate_curve('n17_vd_1',c)
    def test_false_summary_and_wrong_case(self):
        for key,value in [('passed',False),('points',50)]:
            c=curve();c[key]=value
            with self.assertRaises(ValueError):j.validate_curve('n17_vd_1',c)
        c=curve();c['results'][0]['case']='n18_vd_1'
        with self.assertRaises(ValueError):j.validate_curve('n17_vd_1',c)
    def test_wrong_bias(self):
        c=curve();c['results'][3]['vg']=.151
        with self.assertRaises(ValueError):j.validate_curve('n17_vd_1',c)
    def test_nonfinite_or_failed_gate(self):
        for val in (float('nan'),float('inf'),1.01e-6):
            c=curve();c['results'][0]['dual']['Id_relative']=val
            with self.assertRaises(ValueError):j.validate_curve('n17_vd_1',c)
        c=curve();c['results'][0]['reconstruction_passed']=False
        with self.assertRaises(ValueError):j.validate_curve('n17_vd_1',c)
    def test_fields_need_all_six_and_unique_keys(self):
        rows=[dict(case='n17_vd_1',index=0,field=f,unit='weighted_relative_L1' if f=='SRH' else 'dex' if f.endswith('_m3') else 'V',max_abs=0.) for f in j.FIELDS]
        j.validate_fields(rows,{('n17_vd_1',0)})
        with self.assertRaises(ValueError):j.validate_fields(rows+rows[:1],{('n17_vd_1',0)})
        with self.assertRaises(ValueError):j.validate_fields(rows[:-1],{('n17_vd_1',0)})

if __name__=='__main__':unittest.main()
