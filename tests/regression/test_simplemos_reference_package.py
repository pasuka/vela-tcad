"""Portable delivery audit, with no external states or numerical solver."""
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'reference_tcad/simplemos_sentaurus2022/engineering'
SPEC = importlib.util.spec_from_file_location('simplemos_reference_package', PACKAGE / 'verify.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class ReferencePackageTests(unittest.TestCase):
    def test_shipped_evidence(self):
        report = audit.verify(PACKAGE)
        self.assertTrue(report['passed'])
        self.assertFalse(report['external_inputs_checked'])

    def test_tampered_payload(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name) / 'engineering'
            shutil.copytree(PACKAGE, root, ignore=shutil.ignore_patterns('__pycache__'))
            with (root / 'summary.json').open('a') as stream:
                stream.write(' ')
            with self.assertRaisesRegex(ValueError, 'Missing or changed file'):
                audit.verify(root)

    def test_duplicate_point_even_after_reseal(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name) / 'engineering'
            shutil.copytree(PACKAGE, root, ignore=shutil.ignore_patterns('__pycache__'))
            csvfile = root / 'evidence/joined/comparison.csv'
            lines = csvfile.read_text().splitlines()
            lines[-1] = lines[-2]
            csvfile.write_text('\n'.join(lines) + '\n')
            sealfile = root / 'evidence/joined/seal.json'
            seal = audit.read(sealfile); seal['comparison.csv'] = audit.sha(csvfile)
            sealfile.write_text(json.dumps(seal))
            manifest = audit.read(root / 'sha256.json')
            for rel in ('evidence/joined/comparison.csv', 'evidence/joined/seal.json'):
                manifest[rel] = audit.sha(root / rel)
            (root / 'sha256.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'Missing/duplicate comparison points'):
                audit.verify(root)

    def test_missing_external_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            with self.assertRaisesRegex(ValueError, 'Missing or changed file'):
                audit.verify(PACKAGE, Path(name))

    def test_manifest_cannot_escape_root(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name); (root / 'inside').mkdir()
            path = root / 'outside'; path.write_text('untrusted')
            with self.assertRaisesRegex(ValueError, 'Missing or changed file'):
                audit.verify_hashes(root / 'inside', {'../outside': audit.sha(path)})


if __name__ == '__main__':
    unittest.main()
