"""Candidate evidence is distinct from final device acceptance."""
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
import slot_candidate_conservation as audit
import run_slot_ldmos_bv_campaign as campaign
from tests.regression.test_run_slot_ldmos_loadline_continuation import deck


class ConservationTests(unittest.TestCase):
    def result(self):
        policy = dict(tolerance=.01, source_floor=1e-11, absolute_tolerance=1e-12)
        row = dict(flux=1., source=.7, mismatch=1.-.7, ratio=1.-.7,
                   relative_criterion=True, satisfied=False)
        return dict(schema='frozen_candidate_conservation_v1', role='initial_guess_only',
                    device_qualified=False, packing='preserve_referenced_contact_low_parts',
                    equations='before_dirichlet_replacement', policy=policy,
                    electron=row, hole=dict(row), satisfied=False), policy

    def test_rejects_nonfinite_wrong_decision_and_wrong_protocol(self):
        result, policy = self.result()
        self.assertFalse(audit.validate_result(result, policy))
        for mutate in [lambda r: r['electron'].update(flux=float('nan')),
                       lambda r: r['electron'].update(ratio=.001),
                       lambda r: r.update(satisfied=True),
                       lambda r: r.update(electron=None),
                       lambda r: r.update(device_qualified=True),
                       lambda r: r.update(packing='ordinary_restart'),
                       lambda r: r.update(policy={**policy, 'tolerance':.5})]:
            r = campaign.ll.clone(result); mutate(r)
            with self.assertRaises(campaign.ll.EvidenceError):
                audit.validate_result(r, policy)

    def test_absolute_branch_does_not_require_relative_qualification(self):
        result, policy = self.result()
        row = dict(flux=0.,source=1e-14,mismatch=-1e-14,ratio=.001,
                   relative_criterion=False,satisfied=True)
        result.update(electron=row,hole=row,satisfied=True)
        self.assertTrue(audit.validate_result(result,policy))

    def test_bad_candidate_is_reclosed_and_raw_current_never_returned(self):
        for satisfied in (False, True):
            with self.subTest(satisfied=satisfied), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);base=deck();directory=root/'d';directory.mkdir()
                config=root/'config.json';config.write_text('{}')
                (directory/'final_state.h5').write_bytes(b'raw')
                (directory/'iv.csv').write_text('candidate')
                calls=[]
                def restore(base,raw,projected,repo,inner):
                    calls.append('poly');projected.write_bytes(b'poly')
                final=dict(inner_V=38.41,current_A_per_um=6e-14,coupled_newton_updates=2)
                backend=SimpleNamespace(base=base,repo=root,
                    helper=SimpleNamespace(render_attempt=lambda *a:(campaign.ll.clone(base),config,directory),
                                           restore_poly_candidate=restore),
                    _write_run=lambda x:(*x,True),
                    solve=lambda *a,**kw:(calls.append('strict') or final),
                    validate=lambda p:calls.append('validate'))
                def inspect(*a):
                    calls.append('audit');folder=root/'conservation';folder.mkdir()
                    (folder/'summary.json').write_text('{}')
                    return dict(conservation_satisfied=satisfied,evidence={})
                with patch.object(audit,'audit',side_effect=inspect), patch.object(campaign.ll,'read_rows',
                        return_value=[dict(inner_voltage_V='38.41',newton_iterations='2',current_A_per_um='9e-14')]):
                    q,proof=campaign.candidate(backend,dict(inner_V=38.4,state_file='seed',active_ids={}),
                        38.46,38.48,root,conservation_runner=root/'runner.exe')
                self.assertEqual(calls,['audit','poly','strict','validate'])
                self.assertIs(q,final)
                self.assertEqual(q['current_A_per_um'],6e-14)
                self.assertEqual(proof['conservation']['conservation_satisfied'],satisfied)

    def test_mutated_state_aborts_as_evidence_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);state=root/'raw.h5';state.write_bytes(b'raw')
            runner=root/'runner.exe';runner.write_bytes(b'runner');base=deck()
            base['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=.01,
                source_floor=1e-11,absolute_tolerance=1e-12)
            def execute(*a):
                state.write_bytes(b'changed');return 0,False,1.
            b=SimpleNamespace(base=base,repo=root,timeout=1.,_execute=execute)
            with patch.object(campaign.ll,'freeze_inputs',return_value=(base,dict(files={}))):
                with self.assertRaises(campaign.ll.EvidenceError):
                    audit.audit(b,state,1.,root/'audit',runner)


if __name__=='__main__': unittest.main()
