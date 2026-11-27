"""Failure protocol and exception recording; no solver invocation."""
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from tests.regression.test_run_slot_ldmos_loadline_continuation import M


class ProbeFailureTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root/'seed.h5').write_bytes(b'fixture-state')
        self.doc = {'simulation_type':'newton_block_step_probe', 'block_modes':['carrier_density'],
                    'state_file':'seed.h5', 'output_csv':'step.csv', 'solver':{}}
        self.config = self.root/'config.json'
        self.config.write_text(json.dumps(self.doc))
        self.backend = object.__new__(M.SolverBackend)
        self.backend.deadline = math.inf
        self.report = {'schema':'vela.block_candidate_failure.v1',
            'failure_class':'numerical_no_candidate','mode':'carrier_density',
            'initial_merit':1., 'config_sha256':hashlib.sha256(self.config.read_bytes()).hexdigest(),
            'input_file_sha256':M.probe_input_hashes(self.doc,self.root),
            'backtracks':[dict(alpha=2.**-i,finite=True,trial_merit=1.,depletion_clips=1,qf_clips=0) for i in range(20)]}

    def run_probe(self, mutate=None, code=3, timeout=False, write=True):
        def execute(*args):
            if mutate: mutate(self.report)
            if write: (self.root/'numerical_failure.json').write_text(json.dumps(self.report))
            return code,timeout,.01
        self.backend.raw_execute = execute
        return self.backend._execute(None,None,self.config,self.root,1.)

    def test_bound_numerical_failure_is_retryable(self):
        with self.assertRaises(M.NumericalFailure) as caught: self.run_probe()
        self.assertEqual(caught.exception.exit_code,3)

    def test_missing_report_is_evidence_failure(self):
        with self.assertRaises(M.EvidenceError): self.run_probe(write=False)

    def test_malformed_reports_stop(self):
        mutations = [lambda d:d.update(config_sha256='wrong'),
                     lambda d:d.update(input_file_sha256={}),
                     lambda d:d.update(initial_merit=float('nan')),
                     lambda d:d['backtracks'].pop(),
                     lambda d:d['backtracks'][0].update(trial_merit=.5),
                     lambda d:d['backtracks'][0].update(finite=False,trial_merit=1.),
                     lambda d:d['backtracks'][0].update(qf_clips=-1)]
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                original=json.loads(json.dumps(self.report))
                (self.root/'numerical_failure.json').unlink(missing_ok=True)
                with self.assertRaises(M.EvidenceError): self.run_probe(mutate)
                self.report=original

    def test_nonfinite_trial_has_explicit_null(self):
        with self.assertRaises(M.NumericalFailure):
            self.run_probe(lambda d:d['backtracks'][0].update(finite=False,trial_merit=None))

    def test_changed_input_stops(self):
        with self.assertRaises(M.EvidenceError):
            self.run_probe(lambda d:(self.root/'seed.h5').write_bytes(b'changed'))

    def test_stale_report_not_reused(self):
        (self.root/'numerical_failure.json').write_text('{}')
        with self.assertRaises(M.EvidenceError): self.run_probe()

    def test_exit_one_is_not_downgraded(self):
        with self.assertRaisesRegex(RuntimeError,'runtime_environment'):
            self.run_probe(code=1)

    def test_success_without_csv_is_evidence_failure(self):
        with self.assertRaises(M.EvidenceError): self.run_probe(code=0,write=False)

    def test_timeout_and_budget_keep_separate_semantics(self):
        with self.assertRaisesRegex(M.NumericalFailure,'timeout'):
            self.run_probe(code=-1,timeout=True,write=False)
        self.backend.deadline=0
        with self.assertRaises(M.BudgetExhausted): self.run_probe(write=False)

    def test_corrector_records_exception_before_propagating(self):
        helper=M.direct
        doc={'solver':{},'contacts':[{'name':'drain','bias':0.}], 'sweep':{'contact':'drain'}}
        failure=M.NumericalFailure('numerical_no_candidate')
        failure.exit_code=3;failure.status_file='bound-status.json'
        with patch.object(helper,'execute_config',side_effect=failure):
            with self.assertRaises(M.NumericalFailure):
                helper.correct_local_candidate(doc,self.root/'seed.h5',.8,None,self.root,
                    self.root,1.,self.root,8,'density_then_qf')
        data=json.loads((self.root/'local_corrector/summary.json').read_text())
        self.assertFalse(data['passed'])
        self.assertEqual(len(data['steps']),1)
        self.assertEqual(data['steps'][0]['exit_code'],3)
        self.assertEqual(data['steps'][0]['failure_status_file'],'bound-status.json')
        self.assertEqual(data['validation_row'],{})


if __name__=='__main__': unittest.main()
