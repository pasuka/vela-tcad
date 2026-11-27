"""Conservation/proposal properties of the opt-in local patch recovery."""
import sys
from pathlib import Path
import unittest
import tempfile
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
import slot_local_row_repair as r


class LocalRowRepairTests(unittest.TestCase):
    def test_disabled_has_no_callbacks_or_outputs(self):
        def forbidden(*args):
            self.fail('disabled recovery performed work')
        result = r.repair(None, rows_at=forbidden, global_at=forbidden, trial_at=forbidden,
                          neighbors={}, eps=.001, policy=r.Policy())
        self.assertEqual(result['stop_reason'], 'disabled')
        self.assertFalse(result['device_qualified'])

    def test_patch_is_topological_carrier_specific_and_includes_all_violations(self):
        rows = {(0,'hole'):{'ratio':.0012}, (1,'hole'):{'ratio':.0009},
                (2,'hole'):{'ratio':.0002}, (3,'hole'):{'ratio':.0009},
                (0,'electron'):{'ratio':.0008}, (4,'electron'):{'ratio':.002}}
        self.assertEqual(r.select_patch(rows, {0:{1,2}}, .001, r.Policy()),
                         [(0,'hole'), (1,'hole'), (4,'electron')])

    def experiment(self, global_ok=lambda x: True, active_change=False):
        # Selected row improves at alpha=1, but spills into its neighbor.
        # At alpha=1/2 every physical row qualifies, so backtracking must win.
        def rows(x, name):
            values = [1.2e-3*(1-x[0]), .2e-3+1.0e-3*x[0]]
            out = {(i,'hole'):dict(residual=v, scale=1., ratio=abs(v)) for i,v in enumerate(values)}
            if active_change and 'a' in name and name.startswith('c'):
                out.pop((1,'hole'))
            return out
        def trial(x, patch, delta, name):
            y=x.copy()
            for (i,c),d in zip(patch,delta):y[i]+=d
            return y
        # Dimensionless synthetic field uses a smaller coefficient so voltage
        # step remains inside the real policy cap.
        def scaled_rows(x, name):return rows(x/1e-5, name)
        return r.repair(np.zeros(2), rows_at=scaled_rows,
                        global_at=lambda x,n:dict(passed=global_ok(x)), trial_at=trial,
                        neighbors={0:{1}}, eps=.001, policy=r.Policy(enabled=True))

    def test_damping_checks_neighbor_and_prefers_all_row_solution(self):
        out=self.experiment()
        self.assertTrue(out['rows_passed'])
        self.assertFalse(out['device_qualified'])
        self.assertEqual(out['steps'][0]['selected_alpha'], .5)
        self.assertGreater(out['steps'][0]['trials'][0]['max_ratio'], .001)
        self.assertLessEqual(out['steps'][0]['after_max_ratio'], .001)

    def test_row_improvement_cannot_override_original_global_gate(self):
        out=self.experiment(lambda x: x[0]==0)
        self.assertFalse(out['rows_passed'])
        self.assertEqual(out['stop_reason'], 'no_globally_feasible_row_descent')
        self.assertEqual(out['state_file'], '')

    def test_active_set_changes_are_evidence_errors(self):
        with self.assertRaisesRegex(ValueError,'active rows changed'):
            self.experiment(active_change=True)

    def test_global_bad_seed_cannot_enter_local_repair(self):
        out=self.experiment(lambda x:False)
        self.assertFalse(out['steps'])
        self.assertEqual(out['stop_reason'], 'original_global_screen_failed')

    def test_patch_over_budget_does_not_drop_violations(self):
        out=r.repair(None, rows_at=lambda *a:{(i,'hole'):dict(ratio=.002) for i in range(9)},
                     global_at=lambda *a:dict(passed=True), trial_at=lambda *a:self.fail(),
                     neighbors={}, eps=.001, policy=r.Policy(enabled=True))
        self.assertEqual(out['stop_reason'], 'patch_budget_exceeded')

    def test_bad_policy_fails_before_work(self):
        for p in [r.Policy(alphas=(.5,1)),r.Policy(fd_step_V=float('nan')),
                  r.Policy(max_cycles=7),r.Policy(max_dofs=65),r.Policy(neighbor_fraction=0)]:
            with self.assertRaises(ValueError):p.validate()

    def test_explicit_larger_patch_recovers_all_rows_without_changing_default(self):
        def rows(x, name):
            return {(i,'hole'):dict(residual=v, scale=1., ratio=abs(v))
                    for i,v in enumerate(.002-200*x)}
        def trial(x, patch, delta, name):
            y=x.copy()
            for (i,_),d in zip(patch,delta): y[i]+=d
            return y
        for count, limit in ((23, 32), (47, 64)):
            for global_ok in (True, False):
                out=r.repair(np.zeros(count), rows_at=rows, trial_at=trial,
                             global_at=lambda x,n:dict(passed=bool(global_ok or np.all(x==0))),
                             neighbors={}, eps=.001,
                             policy=r.Policy(enabled=True,max_dofs=limit,max_cycles=6))
                self.assertEqual(out['rows_passed'],global_ok)
                self.assertFalse(out['device_qualified'])
                if global_ok:self.assertEqual(len(out['steps'][0]['patch']),count)
        self.assertEqual(r.Policy().max_dofs,8)
        self.assertEqual(r.Policy().max_cycles,3)

    def test_adapter_is_also_disabled_by_default_without_reading_seed(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'unused'
            result=r.run({},Path(tmp)/'missing.h5',0.,Path('missing.exe'),Path(tmp),out)
            self.assertEqual(result['stop_reason'],'disabled')
            self.assertFalse(out.exists())

    def test_singular_patch_fails_without_accepting_any_state(self):
        rows={(0,'hole'):dict(residual=.002,scale=1.,ratio=.002)}
        result=r.repair(0,rows_at=lambda *a:rows,global_at=lambda *a:dict(passed=True),
                        trial_at=lambda *a:1,neighbors={},eps=.001,policy=r.Policy(enabled=True))
        self.assertEqual(result['stop_reason'],'singular_local_jacobian')
        self.assertFalse(result['rows_passed'])
        self.assertEqual(result['state_file'],'')

    def test_connected_weak_neighbor_is_not_lost_after_it_becomes_qualified(self):
        rows={(i,'hole'):dict(ratio=q) for i,q in enumerate([.0012,.0008,.0009,.0001])}
        selected=r.select_patch(rows,{0:{1},1:{0,2},2:{1,3}},.001,r.Policy())
        self.assertEqual(selected,[(0,'hole'),(1,'hole'),(2,'hole')])

    def test_unscreened_trials_are_not_reported_globally_qualified(self):
        out=self.experiment()
        trials=out['steps'][0]['trials']
        self.assertIsNone(trials[0]['global_screen'])
        self.assertTrue(trials[1]['global_screen']['passed'])

    def test_predictive_patch_detects_clean_neighbor_and_reuses_fd_columns(self):
        ids = [72, 913]  # Deliberately unrelated to any production mesh IDs.
        probes = []
        def rows(x, name):
            v = [.0012-120*x[0]+30*x[1], 100*x[0]-100*x[1]]
            if '_fd' in name: probes.append(name)
            return {(i,'hole'):dict(residual=y,scale=1.,ratio=abs(y)) for i,y in zip(ids,v)}
        def trial(x, patch, delta, name):
            y=x.copy()
            for (i,_),dv in zip(patch,delta):y[ids.index(i)]+=dv
            return y
        out=r.repair(np.zeros(2),rows_at=rows,trial_at=trial,
                     global_at=lambda *a:dict(passed=True),neighbors={72:{913}},eps=.001,
                     policy=r.Policy(enabled=True,predict_boundary=True))
        self.assertTrue(out['rows_passed'])
        self.assertFalse(out['device_qualified'])
        self.assertEqual(out['steps'][0]['boundary_prediction']['added'],[(913,'hole')])
        self.assertEqual(len(probes),4)
        self.assertEqual(out['steps'][0]['selected_alpha'],1.)

    def test_boundary_ranking_is_same_carrier_and_predictive(self):
        rows={(i,c):dict(residual=0.,scale=1.,ratio=0.) for i,c in
              [(0,'hole'),(1,'hole'),(2,'hole'),(1,'electron'),(9,'hole')]}
        col={(0,'hole'):{k:(2. if k==(1,'hole') else 3.) for k in rows}}
        ranked=r.boundary_candidates(rows,[(0,'hole')],col,[.001],{0:{1,2}},.001,r.Policy())
        self.assertEqual([x['key'] for x in ranked],[(2,'hole'),(1,'hole')])

    def test_predictive_mode_keeps_patch_budget_and_original_global_gate(self):
        rows={(0,'hole'):dict(residual=.002,scale=1.,ratio=.002)}
        out=r.repair(None,rows_at=lambda *a:rows,global_at=lambda *a:dict(passed=False),
                     trial_at=lambda *a:self.fail(),neighbors={},eps=.001,
                     policy=r.Policy(enabled=True,predict_boundary=True))
        self.assertEqual(out['stop_reason'],'original_global_screen_failed')
        with self.assertRaises(ValueError):r.Policy(predict_boundary=1).validate()

    def test_boundary_expansion_is_limited_to_configured_initial_cycles(self):
        from unittest.mock import patch
        import math
        def rows(x,name):
            v=.008*math.exp(-x/1e-5)
            return {(4,'hole'):dict(residual=v,scale=1.,ratio=v)}
        for cycles in (1,3):
            with patch.object(r,'boundary_candidates',wraps=r.boundary_candidates) as predictor:
                out=r.repair(0.,rows_at=rows,trial_at=lambda x,p,d,n:x+d[0],
                    global_at=lambda *a:dict(passed=True),neighbors={},eps=.001,
                    policy=r.Policy(enabled=True,predict_boundary=True,predict_boundary_cycles=cycles))
                self.assertTrue(out['rows_passed'])
                self.assertEqual(len(out['steps']),3)
                self.assertEqual(predictor.call_count,cycles)
        with self.assertRaises(ValueError):r.Policy(predict_boundary_cycles=0).validate()


if __name__=='__main__':unittest.main()
