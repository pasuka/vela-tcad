import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_host_pilot import verify_package, schedule, CONFIGS
from analyze_templates_ldmos_host_pilot import verified_point, compare_hosts


class HostPilotTests(unittest.TestCase):
    def test_reject_changed_or_escaping_package(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); payload = root/'runner.exe'; payload.write_bytes(b'original')
            manifest = root/'package_manifest.json'
            sha = hashlib.sha256(b'original').hexdigest()
            manifest.write_text(json.dumps({'files':{'runner.exe':sha}}))
            verify_package(root)
            payload.write_bytes(b'changed')
            with self.assertRaises(ValueError): verify_package(root)
            manifest.write_text(json.dumps({'files':{'../outside.exe':sha}}))
            with self.assertRaises(ValueError): verify_package(root)

    def test_pilot_is_bounded_and_single_threaded(self):
        rows = schedule()
        self.assertEqual(len(rows), 5)
        self.assertEqual(len({r['key'] for r in rows}), 5)
        self.assertEqual({CONFIGS[r['config']]['backend'] for r in rows},
                         {'umfpack','strumpack','sparselu','mumps','superlu_mt'})
        for row in rows:
            self.assertEqual((row['points'],row['gate'],row['round']), (8,8,0))
            self.assertEqual(CONFIGS[row['config']]['threads'], 1)
            self.assertTrue(CONFIGS[row['config']]['reuse'])

    def test_remote_evidence_uses_relative_case_and_verifies_hash(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);case=root/'fixed/child_0';case.mkdir(parents=True)
            state=case/'state.csv';state.write_bytes(b'exact point')
            point=dict(case='fixed/child_0',state='Z:/remote-only/state.csv',
                       sha256=hashlib.sha256(b'exact point').hexdigest())
            self.assertEqual(verified_point(root,point),state.resolve())
            state.write_bytes(b'changed')
            with self.assertRaises(ValueError): verified_point(root,point)

    def test_reject_different_binary_before_comparing_results(self):
        with tempfile.TemporaryDirectory() as temp:
            a,b=Path(temp)/'a',Path(temp)/'b';a.mkdir();b.mkdir()
            for root,sha in ((a,'one'),(b,'two')):
                (root/'summary.json').write_text(json.dumps(dict(status='pass',
                    package_sha256='same',runner_sha256=sha)))
            with self.assertRaisesRegex(ValueError,'runner_sha256'):compare_hosts(a,b)


if __name__ == '__main__': unittest.main()
