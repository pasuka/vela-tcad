"""Control-flow tests only; mocked solver rows are not numerical evidence."""
import contextlib
import csv
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('idvd', ROOT / 'scripts/run_slot_ldmos_idvd.py')
idvd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(idvd)


class AdaptiveGatePolicy(unittest.TestCase):
    def test_resume_warns_above_point_one_and_halves_actual_failed_step(self):
        parent=.4962890625;calls=[];decisions=[]
        from slot_fermi_poly_predictor import monitor_bias_change
        def solve(v):
            nonlocal parent
            calls.append((parent,v));ok=len(calls)!=1
            monitor_bias_change(parent,v)
            if ok:parent=v
            return dict(qualified=ok,iv={'newton_iterations':'1'})
        with self.assertLogs('slot_fermi_poly_predictor',level='WARNING'):
            idvd.advance_adaptive_gate(solve,decisions.append,start=parent,target=1.2,initial=.25)
        self.assertAlmostEqual(calls[0][1]-calls[0][0],.25)
        self.assertAlmostEqual(calls[1][1]-calls[1][0],.125)
        self.assertEqual(calls[0][0],calls[1][0])
        self.assertEqual(parent,1.2)

    def test_predictor_cap_covers_actual_failure_and_reaches_goal(self):
        from slot_fermi_poly_predictor import GATE_BIAS_WARNING_THRESHOLD_V
        parent=0.;decisions=[]
        def solve(v):
            nonlocal parent
            self.assertLessEqual(v-parent,GATE_BIAS_WARNING_THRESHOLD_V+1e-12)
            parent=v
            return dict(qualified=True,iv={'newton_iterations':'1'})
        idvd.advance_adaptive_gate(solve,decisions.append,predictor_maximum=GATE_BIAS_WARNING_THRESHOLD_V)
        self.assertAlmostEqual(decisions[7]['target_V'],.2462890625)
        self.assertAlmostEqual(decisions[8]['requested_step_V'],.12814453125)
        self.assertAlmostEqual(decisions[8]['target_V'],.3462890625)
        self.assertEqual(decisions[8]['limit_reason'],'predictor_capability')
        self.assertEqual(parent,5.)
        self.assertLessEqual(max(d['next_step_V'] for d in decisions),.1)

    def test_cap_honored_on_first_attempt_and_retry_remains_halved(self):
        decisions=[];calls=[]
        def solve(v):
            calls.append(v)
            return dict(qualified=len(calls)>1,iv={'newton_iterations':'1'})
        idvd.advance_adaptive_gate(solve,decisions.append,target=.4,initial=.25,predictor_maximum=.1)
        self.assertEqual(calls[:2],[.1,.05])
        self.assertEqual(decisions[0]['parent_V'],decisions[1]['parent_V'])
        self.assertTrue(all(d['span_V']<=.1+1e-12 for d in decisions))

    def test_native_first_gate_step_replay(self):
        # Native first gate solve used 3 updates: 5 mV -> 7.233333 mV.
        self.assertAlmostEqual(.005*idvd.adaptive_growth(3), .007233333333333333)
        self.assertAlmostEqual(.006*idvd.adaptive_growth(4,1.2), .007008)
        self.assertEqual(idvd.adaptive_growth(25),1.)
        self.assertEqual(idvd.adaptive_growth(0),1.5)

    def test_failure_keeps_parent_and_discards_old_target(self):
        accepted=0.;calls=[];decisions=[]
        def solve(v):
            nonlocal accepted
            calls.append((accepted,v))
            ok=len(calls)!=2
            if ok:accepted=v
            return dict(qualified=ok,iv={'newton_iterations':'1'})
        idvd.advance_adaptive_gate(solve,decisions.append,target=.025)
        self.assertAlmostEqual(calls[0][1],.005)
        self.assertEqual(calls[1][0],calls[2][0])
        self.assertAlmostEqual(calls[2][1],.00875)
        self.assertAlmostEqual(calls[3][1],.014375) # no forced return to .0125
        self.assertEqual(accepted,.025)

    def test_bounded_failure_and_successful_maximum(self):
        decisions=[]
        with self.assertRaisesRegex(RuntimeError,'failure_limit'):
            idvd.advance_adaptive_gate(lambda v:dict(qualified=False),decisions.append)
        self.assertEqual(len(decisions),8)
        self.assertTrue(all(d['parent_V']==0 for d in decisions))
        decisions=[]
        end=idvd.advance_adaptive_gate(lambda v:dict(qualified=True,iv={'newton_iterations':'2'}),decisions.append)
        self.assertEqual(end,5.)
        self.assertLessEqual(max(d['span_V'] for d in decisions),.25+1e-14)
        self.assertLess(len(decisions),50)


class ProgressiveDrainPolicy(unittest.TestCase):
    def test_failed_trial_does_not_become_parent_and_ratio_stays_bounded(self):
        parent=2.;span=.04;calls=[]
        def solve(v,allow):
            nonlocal parent,span
            calls.append((parent,v,allow))
            ok=v-parent<=.03000001
            if ok: span=v-parent;parent=v
            return dict(qualified=ok,returncode=0 if ok else 1,predictor={} if allow else None)
        self.assertAlmostEqual(idvd.advance_progressive(2.,2.12,.04,solve),2.12)
        self.assertEqual(sum(not a for _,_,a in calls),1)
        self.assertEqual(calls[0][0],calls[1][0])
        self.assertLessEqual(len(calls),32)

    def test_timeout_skips_warm_and_halving_is_bounded(self):
        calls=[]
        def solve(v,allow):
            calls.append((v,allow));return dict(qualified=False,returncode=124,predictor={})
        with self.assertRaisesRegex(RuntimeError,'progressive_halvings'):
            idvd.advance_progressive(2.,3.,.04,solve)
        self.assertEqual(len(calls),4)
        self.assertTrue(all(a for _,a in calls))
        self.assertAlmostEqual(calls[-1][0]-2.,.05/8)

    def test_interval_attempt_cap_survives_interleaved_success(self):
        def solve(v,allow):return dict(qualified=True)
        with self.assertRaisesRegex(RuntimeError,'attempt_limit'):
            idvd.advance_progressive(0.,60.,.001,solve,max_attempts=3)

    def test_soft_budget_failure_is_not_repeated_as_warm_start(self):
        calls=[]
        def solve(v,allow):
            calls.append(allow)
            return dict(qualified=False,returncode=1,predictor={},iv={'failure_reason':'soft_time_limit'})
        with self.assertRaisesRegex(RuntimeError,'progressive_halvings'):
            idvd.advance_progressive(2.,3.,.04,solve)
        self.assertEqual(calls,[True]*4)

    def test_invalid_target_and_attempt_limit_are_rejected(self):
        for target,limit in ((float('nan'),32),(3.,0)):
            with self.assertRaises(ValueError):
                idvd.advance_progressive(2.,target,.04,lambda v,a:None,max_attempts=limit)


class IdVdCampaign(unittest.TestCase):
    def test_adaptive_gate_only_never_launches_drain_and_finishes_summary(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);resource=root/'resource';resource.write_text('{}')
            c=json.loads((ROOT/'configs/regression/slot_ldmos_idvd_fermi_control.json').read_text())
            for key in ('mesh_file','materials_file','node_doping_file'):c[key]=str(resource)
            config=root/'config.json';config.write_text(json.dumps(c))
            native=root/'native.csv';native.write_text('gate_V,drain_V,Id_A_per_um\n5,60,.001\n')
            calls=[];out=root/'run';parents=[]
            def predict(config,parent,output,target):
                from slot_fermi_poly_predictor import GATE_BIAS_WARNING_THRESHOLD_V
                old=json.loads((Path(parent).parent/'config.json').read_text())['sweep']['bias_points'][0]
                self.assertLessEqual(target-old,.25+1e-12)
                parents.append((old,target));output.write_text('unqualified')
                return dict(qualified=False)
            def solve(command,**kwargs):
                cfg=json.loads(Path(command[-1]).read_text());s=cfg['sweep'];v=s['bias_points'][0]
                Path(s['write_state_file']).write_text('accepted')
                row=dict(converged='1',solver_method=cfg['solver']['method'],bias_V=str(v),newton_iterations='3',
                    carrier_row_violations='0',carrier_row_max_ratio='.0001',
                    global_continuity_closure_satisfied='1',qf_bounds_violations='0',current_total_A_per_um='0')
                with Path(cfg['output_csv']).open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=row);w.writeheader();w.writerow(row)
                calls.append(cfg);return SimpleNamespace(returncode=0)
            argv=['idvd','--config',str(config),'--runner',str(resource),'--native',str(native),
                '--output',str(out),'--gate-step-policy','adaptive','--gate-only',
                '--gate-initialization','fermi_poly_poisson']
            with patch.object(sys,'argv',argv),patch.object(idvd.subprocess,'run',side_effect=solve), \
                 patch.object(idvd,'make_gate_seed',side_effect=predict), \
                 patch.dict(idvd.os.environ,{'VELA_LINEAR_SOLVER':''}),contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(idvd.main(),0)
            s=json.loads((out/'summary.json').read_text())
            self.assertTrue(s['completed'] and s['goal_reached'] and s['resources_unchanged'])
            self.assertEqual(s['stop_reason'],'gate_5V_zero_drain_qualified')
            self.assertEqual(calls[1]['sweep']['bias_points'],[.005])
            self.assertEqual(calls[-1]['sweep']['bias_points'],[5.])
            self.assertTrue(all(c['sweep']['contact']=='gate' for c in calls[1:]))
            self.assertFalse(any(r['stage']=='drain' for r in s['records']))
            self.assertEqual(s['policy']['gate_effective_max_step_V'],.25)
            self.assertEqual(s['policy']['gate_predictor_warning_threshold_V'],.1)
            self.assertTrue(any(v-old>.1+1e-12 for old,v in parents))
            self.assertEqual(s['policy']['gate_max_step_V'],.25)
            self.assertEqual(len(parents),len(s['records'])-1)

    def test_poisson_handoff_is_not_a_qualified_parent_and_is_gate_only(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);resource=root/'resource';resource.write_text('{}')
            c=json.loads((ROOT/'configs/regression/slot_ldmos_idvd_fermi_control.json').read_text())
            for key in ('mesh_file','materials_file','node_doping_file'):c[key]=str(resource)
            config=root/'config.json';config.write_text(json.dumps(c))
            native=root/'native.csv';native.write_text('gate_V,drain_V,Id_A_per_um\n5,.1,.001\n')
            out=root/'run';calls=[];parents=[];rejected=False
            def predict(config,parent,output,target):
                parents.append(Path(parent));output.write_text('unqualified predictor')
                return dict(qualified=False)
            def solve(command,**kwargs):
                nonlocal rejected
                c=json.loads(Path(command[-1]).read_text());s=c['sweep'];method=c['solver']['method']
                v=s['bias_points'][0]
                reject=s['contact']=='gate' and method=='newton' and v==.1 and not rejected
                rejected |= reject
                Path(s['write_state_file']).write_text('state')
                row=dict(converged='1',solver_method=method,bias_V=str(v),
                    carrier_row_violations='1' if reject or method=='poisson_only' else '0',
                    carrier_row_max_ratio='.002' if reject else '.0001',
                    global_continuity_closure_satisfied='1',qf_bounds_violations='0',current_total_A_per_um='.001')
                with Path(c['output_csv']).open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=row);w.writeheader();w.writerow(row)
                calls.append(c);return SimpleNamespace(returncode=0)
            argv=['idvd','--config',str(config),'--runner',str(resource),'--native',str(native),
                  '--output',str(out),'--budget','60','--child-cap','3',
                  '--drain-stop','.1',
                  '--gate-initialization','fermi_poly_poisson']
            with patch.object(sys,'argv',argv),patch.object(idvd,'make_gate_seed',side_effect=predict), \
                 patch.object(idvd.subprocess,'run',side_effect=solve), \
                 patch.dict(idvd.os.environ,{'VELA_LINEAR_SOLVER':''}),contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(idvd.main(),0)
            summary=json.loads((out/'summary.json').read_text())
            self.assertEqual(parents[0],parents[1])  # rejected strict state is isolated
            self.assertEqual(parents[0],Path(calls[0]['sweep']['write_state_file']))
            self.assertTrue(all('/poisson/' not in r['state'].replace('\\','/') for r in summary['records']))
            self.assertTrue(all(r['iv']['solver_method']=='newton' for r in summary['records']))
            self.assertEqual(calls[-1]['sweep']['contact'],'drain')
            self.assertEqual(calls[-1]['solver']['method'],'newton')
            self.assertEqual(calls[-1]['sweep']['initial_state_file'],summary['records'][-2]['state'])
            self.assertTrue(summary['resources_unchanged'])

    def test_rejected_state_never_becomes_retry_parent(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            resource = base / 'resource.json'
            resource.write_text('{}')
            config = json.loads((ROOT / 'configs/regression/slot_ldmos_idvd_fermi_control.json').read_text())
            for key in ['mesh_file', 'node_doping_file', 'materials_file']:
                config[key] = str(resource)
            config_path = base / 'config.json'
            config_path.write_text(json.dumps(config))
            runner = base / 'stub.exe'
            runner.write_bytes(b'not-an-executable')
            native = base / 'native.csv'
            native.write_text('gate_V,drain_V,Id_A_per_um\n5,1,1e-4\n')
            out = base / 'run'
            calls = []
            rejected_once = False

            def solve(command, **kwargs):
                nonlocal rejected_once
                c = json.loads(Path(command[-1]).read_text())
                sweep = c['sweep']
                target = sweep['bias_points'][0]
                reject = sweep['contact'] == 'gate' and target == .1 and not rejected_once
                rejected_once |= reject
                state = Path(sweep['write_state_file'])
                # Even a process that exits successfully and leaves a state
                # must not become a parent if one original gate is violated.
                state.write_text('rejected' if reject else f'accepted-{len(calls)}')
                row = dict(converged='1', solver_method='newton', carrier_row_violations='1' if reject else '0',
                           global_continuity_closure_satisfied='1', qf_bounds_violations='0',
                           carrier_row_max_ratio='.002' if reject else '.0001',
                           current_total_A_per_um='1e-4', bias_V=str(target))
                with Path(c['output_csv']).open('w', newline='') as file:
                    writer = csv.DictWriter(file, fieldnames=row)
                    writer.writeheader()
                    writer.writerow(row)
                calls.append(c)
                self.assertLessEqual(kwargs['timeout'], 3.)
                return SimpleNamespace(returncode=0)

            argv = ['idvd', '--config', str(config_path), '--runner', str(runner), '--native', str(native),
                    '--output', str(out), '--budget', '30', '--child-cap', '3']
            argv += ['--drain-stop','1']
            with patch.object(sys, 'argv', argv), patch.object(idvd.subprocess, 'run', side_effect=solve), \
                 patch.dict(idvd.os.environ, {'VELA_LINEAR_SOLVER': ''}), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(idvd.main(), 0)
            summary = json.loads((out / 'summary.json').read_text())
            self.assertTrue(summary['goal_reached'])
            self.assertTrue(summary['resources_unchanged'])
            self.assertFalse(summary['attempts'][1]['qualified'])
            self.assertEqual(calls[2]['sweep']['initial_state_file'], calls[0]['sweep']['write_state_file'])
            self.assertEqual(calls[3]['sweep']['initial_state_file'], calls[2]['sweep']['write_state_file'])
            self.assertEqual(summary['records'][-1]['current_relative_difference'], 0.)
            self.assertTrue(all(r['state'] != calls[1]['sweep']['write_state_file'] for r in summary['records']))

    def test_secant_failure_warm_failure_and_bisection_keep_qualified_history(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);resource=root/'resource';resource.write_text('{}')
            c=json.loads((ROOT/'configs/regression/slot_ldmos_idvd_fermi_control.json').read_text())
            for key in ('mesh_file','materials_file','node_doping_file'):c[key]=str(resource)
            config=root/'config.json';config.write_text(json.dumps(c))
            native=root/'native.csv';native.write_text('gate_V,drain_V,Id_A_per_um\n5,.2,.001\n5,.3,.002\n')
            imports=[]
            for i,v in enumerate((0.,.1)):
                p=root/f'accepted{i}.h5';p.write_text(str(v))
                imports.append(dict(state=str(p),state_sha256=idvd.sha(p),voltage_V=v))
            out=root/'run';seeds=[];calls=[]
            def predictor(config,older,parent,output,target):
                seeds.append((str(older),str(parent),target));output.write_text('unqualified prediction')
                return dict(qualified=False)
            def solve(command,**kwargs):
                cfg=json.loads(Path(command[-1]).read_text());s=cfg['sweep'];v=s['bias_points'][0]
                rejected=len(calls)<2
                Path(s['write_state_file']).write_text('bad' if rejected else str(v))
                row=dict(converged='1',solver_method='newton',bias_V=str(v),
                    carrier_row_violations='1' if rejected else '0',carrier_row_max_ratio='.002' if rejected else '.0001',
                    global_continuity_closure_satisfied='1',qf_bounds_violations='0',current_total_A_per_um='.001')
                with Path(cfg['output_csv']).open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=row);w.writeheader();w.writerow(row)
                calls.append(cfg);return SimpleNamespace(returncode=0)
            argv=['idvd','--config',str(config),'--runner',str(resource),'--native',str(native),
                  '--output',str(out),'--budget','60','--child-cap','3','--drain-stop','.3',
                  '--drain-initialization','secant','--drain-import-summary',str(root/'import.json')]
            with patch.object(sys,'argv',argv),patch.object(idvd,'load_drain_history',return_value=list(imports)), \
                 patch.object(idvd,'make_drain_seed',side_effect=predictor),patch.object(idvd.subprocess,'run',side_effect=solve), \
                 patch.dict(idvd.os.environ,{'VELA_LINEAR_SOLVER':''}),contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(idvd.main(),0)
            s=json.loads((out/'summary.json').read_text())
            self.assertEqual(calls[1]['sweep']['initial_state_file'],imports[-1]['state'])
            self.assertEqual(seeds[0][:2],seeds[1][:2])
            self.assertEqual([a['qualified'] for a in s['attempts']],[False,False,True,True,True])
            rejected={calls[i]['sweep']['write_state_file'] for i in (0,1)}
            self.assertTrue(all(old not in rejected and parent not in rejected for old,parent,_ in seeds))
            self.assertEqual(len(s['imported_records']),2)
            self.assertEqual(len(s['records']),3)  # midpoint is recorded, imports are not new
            self.assertEqual(s['stop_reason'],'requested_drain_segment_qualified')

    def test_import_rejects_tampered_state_before_prediction(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);state=root/'state.h5';state.write_text('modified')
            cfg={'solver':{'carrier_row_convergence':{}},'contacts':[],'sweep':{'contact':'drain'}}
            (root/'config.json').write_text(json.dumps(cfg))
            (root/'iv.csv').write_text('bias_V\n0.1\n')
            (root/'exit.json').write_text(json.dumps(dict(qualified=True,returncode=0,
                config_sha256=idvd.sha(root/'config.json'))))
            record=dict(state=str(state),state_sha256='stale hash',voltage_V=.1,iv={'bias_V':'0.1'})
            manifest=root/'summary.json';manifest.write_text(json.dumps(dict(completed=True,records=[record,record])))
            with self.assertRaisesRegex(ValueError,'evidence/qualification mismatch'):
                idvd.load_drain_history(manifest,cfg,{})


if __name__ == '__main__':
    unittest.main()
