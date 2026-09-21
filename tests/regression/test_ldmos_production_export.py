"""Synthetic provenance/path tests, not device calibration evidence."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts.export_templates_ldmos_production import child, export_bundle


class ProductionExportTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.evidence = self.root / 'evidence'
        self.evidence.mkdir()
        self.output = self.root / 'new bundle with spaces'
        self.cfg = dict(mesh_file='/frozen/mesh.json',
                        mobility_SI=dict(geometry_file='/frozen/geometry.json'),
                        electrothermal_linear_solver='umfpack',
                        boundaries=[dict(node=1, kind='temperature', value=300.)],
                        electrical_gate_solver=dict(mode='enforce', eps_row=.001),
                        diagnostic_csv='@output/rows.csv')
        self.deck = dict(input_file='/frozen/input.json', output_directory='/frozen/results',
                         simulation_type='electrothermal_dc_sweep',
                         initialization=dict(mode='neutral_300K', gate_voltage_V=4),
                         sweep=dict(bias_points_V=[0., 1.3333333333333333, 40.]))
        for name, value in [('input.json', self.cfg), ('deck.json', self.deck),
                            ('mesh.json', {'nodes': [1]}), ('geometry.json', {'weights': [.7]})]:
            (self.evidence / name).write_text(json.dumps(value))
        profile = dict(files_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in self.evidence.iterdir()},
                       dependencies=[dict(recorded_path='/frozen/' + n, evidence_path=n,
                                          bundle_path='data/' + n)
                                     for n in ('mesh.json', 'geometry.json')],
                       cases=[dict(gate_V=4, input='input.json', deck='deck.json')],
                       qualified_linux_runner_sha256='synthetic-test-only')
        self.profile = self.root / 'profile.json'
        self.profile.write_text(json.dumps(profile))

    def export(self):
        return export_bundle(self.evidence, self.output, self.profile)

    def test_only_declared_paths_change_and_manifest_covers_all_generated_inputs(self):
        manifest = self.export()
        actual = json.loads((self.output / 'input_vg4.json').read_text())
        expected = copy.deepcopy(self.cfg)
        expected['mesh_file'] = str(self.output / 'data/mesh.json')
        expected['mobility_SI']['geometry_file'] = str(self.output / 'data/geometry.json')
        self.assertEqual(actual, expected)
        deck = json.loads((self.output / 'vg4.json').read_text())
        self.assertEqual(deck['sweep'], self.deck['sweep'])
        self.assertEqual(deck['initialization'], self.deck['initialization'])
        self.assertEqual(deck['input_file'], 'input_vg4.json')
        for name, digest in manifest['files_sha256'].items():
            self.assertEqual(hashlib.sha256((self.output / name).read_bytes()).hexdigest(), digest)
        self.assertEqual(len(manifest['files_sha256']), 4)

    def test_cli_explicit_profile_exports_its_hashed_configuration(self):
        script = Path(__file__).resolve().parents[2] / 'scripts/export_templates_ldmos_production.py'
        result = subprocess.run([sys.executable, str(script), '--evidence-root', str(self.evidence),
                                 '--output', str(self.output), '--profile', str(self.profile)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = json.loads((self.output / 'manifest.json').read_text())
        self.assertEqual(manifest['profile_sha256'], hashlib.sha256(self.profile.read_bytes()).hexdigest())
        self.assertEqual(manifest['qualified_linux_runner_sha256'], 'synthetic-test-only')
        actual = json.loads((self.output / 'input_vg4.json').read_text())
        self.assertEqual(actual['electrical_gate_solver'], self.cfg['electrical_gate_solver'])

    def test_changed_geometry_or_gate_is_rejected_before_writing(self):
        for name in ('geometry.json', 'input.json'):
            with self.subTest(name=name):
                p = self.evidence / name
                original = p.read_bytes()
                p.write_text('{}')
                with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
                    self.export()
                self.assertFalse(self.output.exists())
                p.write_bytes(original)

    def test_existing_output_is_never_overwritten(self):
        self.output.mkdir()
        marker = self.output / 'keep.txt'
        marker.write_text('existing run')
        with self.assertRaisesRegex(ValueError, 'already exists'):
            self.export()
        self.assertEqual(marker.read_text(), 'existing run')

    def test_missing_dependency_is_rejected_before_writing(self):
        (self.evidence / 'geometry.json').unlink()
        with self.assertRaises(FileNotFoundError):
            self.export()
        self.assertFalse(self.output.exists())

    def test_traversal_and_foreign_absolute_paths_are_rejected(self):
        for value in ('../outside', '/absolute', 'data/../../outside', 'C:/absolute', 'data\\file'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                child(self.evidence, value)


if __name__ == '__main__':
    unittest.main()
