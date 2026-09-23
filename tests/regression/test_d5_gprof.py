import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_d5_gprof import relocate,trajectory_signature,validate_completed_curve


class GprofBaselineTest(unittest.TestCase):
    def test_reuse_rejects_partial_or_changed_trajectory(self):
        row=dict(parent_V=0,target_V=40,cap_V=.2,stage='direct',frame_V=0,Newton_updates=5)
        ledger=dict(status='completed',runs=[row],exact_points=[dict(bias_V=i*40/30) for i in range(31)])
        case=dict(status='completed',trajectory_exact=True,audit=dict(integrity_pass=True),
                  verdicts=dict(final=dict(pass_all=True)))
        validate_completed_curve(case,ledger,ledger)
        with self.assertRaises(ValueError):validate_completed_curve(dict(case,status='failed'),ledger,ledger)
        with self.assertRaises(ValueError):validate_completed_curve(case,dict(ledger,exact_points=ledger['exact_points'][:9]),ledger)
        with self.assertRaises(ValueError):validate_completed_curve(case,ledger,dict(ledger,runs=[dict(row,Newton_updates=6)]))

    def test_relocation_preserves_physics_and_external_inputs(self):
        root=Path('work/case');out=Path('@output')
        cfg={'mesh':'frozen/mesh.json','solver':{'eps':1e-8},'outputs':[str(root/'state.csv')]}
        changed=relocate(cfg,root,out)
        self.assertEqual(changed['mesh'],cfg['mesh'])
        self.assertEqual(changed['solver'],cfg['solver'])
        self.assertEqual(changed['outputs'],[str(out/'state.csv')])
        self.assertNotEqual(cfg['outputs'],changed['outputs'])

    def test_trajectory_check_includes_intermediate_work_and_frame(self):
        row=dict(parent_V=1,target_V=1.2,cap_V=.2,stage='direct',frame_V=0,Newton_updates=5)
        original={'runs':[row]}
        self.assertNotEqual(trajectory_signature(original),trajectory_signature({'runs':[dict(row,Newton_updates=6)]}))
        self.assertNotEqual(trajectory_signature(original),trajectory_signature({'runs':[dict(row,frame_V=28)]}))
        self.assertNotEqual(trajectory_signature(original),trajectory_signature({'runs':[row,row]}))


if __name__=='__main__':unittest.main()
