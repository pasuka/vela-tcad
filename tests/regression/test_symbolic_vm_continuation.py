"""Continuation preserves interrupted evidence and rejects changed identities."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_symbolic_vm import carry_completed, sha


class Continuation(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'vela_example_runner').write_bytes(b'frozen test identity')
        self.digest = sha(self.root/'vela_example_runner')
        self.native = {'cases_sha256': {}}
        self.directory = self.root/'native'
        self.directory.mkdir()
        for name in ('IdVd.cmd', 'n1_fps.tdr', 'sdevice.par', 'Siliconc100.par'):
            (self.directory/name).write_text(name)
            self.native['cases_sha256']['native_vg4/'+name] = sha(self.directory/name)
        for k in range(31):
            (self.directory/('field_vg4_%04d_des.tdr' % k)).touch()
        self.row = dict(repeat=0, gate_V=4, variant='native', exit_code=0,
                        directory=str(self.directory), external_wall_seconds=123.)
        self.prior = dict(status='paused_at_deadline', runner_sha256=self.digest,
                          profile_sha256='profile', repeats=2, runs=[self.row])

    def invoke(self, digest=None):
        path = self.root/'summary.json'
        path.write_text(json.dumps(self.prior))
        before = path.read_bytes()
        result = carry_completed(self.root, digest or self.digest, 'profile', 2,
                                 self.root, {}, self.native)
        self.assertEqual(before, path.read_bytes())
        return result

    def test_carries_complete_keeps_interrupted_cost(self):
        interrupted = dict(self.row, gate_V=8, exit_code=-15,
                           interrupted_at_deadline=True, external_wall_seconds=215.)
        self.prior['runs'].append(interrupted)
        complete, excluded = self.invoke()
        self.assertEqual(complete, [self.row])
        self.assertEqual(excluded, [interrupted])

    def test_changed_binary_rejected(self):
        with self.assertRaises(AssertionError):
            self.invoke('changed')

    def test_changed_input_rejected(self):
        (self.directory/'IdVd.cmd').write_text('changed')
        with self.assertRaises(AssertionError):
            self.invoke()

    def test_missing_native_fields_rejected(self):
        (self.directory/'field_vg4_0000_des.tdr').unlink()
        with self.assertRaises(AssertionError):
            self.invoke()

    def test_unfinished_audit_rejected(self):
        self.row['variant'] = 'high'
        with self.assertRaises(AssertionError):
            self.invoke()


if __name__ == '__main__':
    unittest.main()
