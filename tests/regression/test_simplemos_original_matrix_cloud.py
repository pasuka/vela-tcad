import unittest,sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import simplemos_original_matrix_cloud_20260926 as m
class ControlFlow(unittest.TestCase):
    def test_failed_equilibrium_never_starts_ramp_or_sweep(self):
        with tempfile.TemporaryDirectory() as td:
            obj=object.__new__(m.Matrix);obj.root=Path(td);calls=[]
            def fail(c,label,*args,**kwargs):
                calls.append(label)
                if label=='poisson':return obj.root/'poisson'
                raise RuntimeError('deliberately unqualified equilibrium')
            obj.sweep=fail
            self.assertEqual(obj.case(dict(case='n23_vd_1',vd=1.)),[])
            self.assertEqual(calls,['poisson','equilibrium'])
            self.assertIn('unqualified',m.h.read(obj.root/'n23_vd_1/failure.json')['error'])
    def test_failed_srh_gate_prevents_any_matrix_execution(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);m.h.write(p/'srh/summary.json',dict(passed=False))
            with self.assertRaises(AssertionError):m.Matrix(p/'matrix',p/'missing_inputs',p/'missing_runner',p/'srh')
    def test_failed_poisson_never_starts_coupled_equilibrium(self):
        with tempfile.TemporaryDirectory() as td:
            obj=object.__new__(m.Matrix);obj.root=Path(td);calls=[]
            def fail(c,label,*args,**kwargs):
                calls.append(label);raise RuntimeError('unqualified Poisson initializer')
            obj.sweep=fail
            self.assertEqual(obj.case(dict(case='n23_vd_1',vd=1.)),[])
            self.assertEqual(calls,['poisson'])
    def test_full_poisson_stage_preserves_coupled_gates(self):
        with tempfile.TemporaryDirectory() as td:
            obj=object.__new__(m.Matrix);obj.root=Path(td)
            obj.config=lambda c:dict(contacts=[dict(name='drain',bias=0.)],solver=dict(method='newton',reltol=1e-7,abstol=1e-12,carrier_row_convergence=dict(eps_row=1e-6)))
            class Executor:
                def execute(self,path):return dict(exit_code=0,converged=True)
            obj.executor=Executor();case=dict(case='n23_vd_1',vd=1.)
            first=obj.sweep(case,'poisson','gate',[0.])
            second=obj.sweep(case,'equilibrium','gate',[0.],first/'accepted.h5')
            p=m.h.read(first/'config.json');e=m.h.read(second/'config.json')
            self.assertEqual(p['solver']['method'],'poisson_only')
            self.assertEqual(p['sweep']['initialization']['mode'],'poisson_block')
            self.assertEqual(p['solver']['reltol'],0.)
            self.assertEqual(e['solver']['method'],'newton')
            self.assertEqual(e['solver']['reltol'],1e-7)
            self.assertEqual(e['solver']['abstol'],1e-12)
            self.assertEqual(e['solver']['carrier_row_convergence']['eps_row'],1e-6)
            self.assertNotIn('initialization',e['sweep'])
if __name__=='__main__':unittest.main()
