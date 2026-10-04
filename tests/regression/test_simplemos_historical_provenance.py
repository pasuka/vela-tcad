"""Negative controls for exact historical provenance, not live-source hashes."""
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
import zipfile

from tests.regression.simplemos_evidence_chain import (
    ARCHIVE, MANIFEST, REPO, validate_historical_provenance)


class HistoricalProvenanceTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.manifest = self.root / 'manifest.json'
        self.archive = self.root / 'sources.zip'
        shutil.copyfile(MANIFEST, self.manifest)
        shutil.copyfile(ARCHIVE, self.archive)
        self.data = json.loads(self.manifest.read_text())
        name = next(n for n in self.data['evidence'] if n.startswith('simplemos_m46_'))
        self.evidence = self.root / name
        shutil.copyfile(REPO / 'reference_tcad/simplemos_sentaurus2022' / name,
                        self.evidence)

    def check(self):
        return validate_historical_provenance(self.evidence, self.manifest, self.archive)

    def test_all_evidence_and_explicit_gap(self):
        gaps = {}
        for name in self.data['evidence']:
            gaps[name] = validate_historical_provenance(
                REPO / 'reference_tcad/simplemos_sentaurus2022' / name)
        self.assertEqual(len(gaps), 16)
        self.assertEqual([n.split('_')[1] for n, g in gaps.items() if g],
                         ['m33', 'm35', 'm37'])
        self.assertEqual(sum(v['status'] == 'recovered'
                             for v in self.data['sources'].values()), 132)

    def test_changed_evidence_bytes_rejected(self):
        with self.evidence.open('ab') as stream:
            stream.write(b' ')
        with self.assertRaisesRegex(ValueError, 'evidence bytes'):
            self.check()

    def test_unknown_evidence_rejected(self):
        del self.data['evidence'][self.evidence.name]
        self.manifest.write_text(json.dumps(self.data))
        with self.assertRaises(KeyError):
            self.check()

    def test_source_map_drift_rejected_even_with_refreshed_evidence_digest(self):
        data = json.loads(self.evidence.read_text())
        data['source_hashes']['invented.cpp'] = '0' * 64
        self.evidence.write_text(json.dumps(data))
        self.data['evidence'][self.evidence.name]['sha256'] = hashlib.sha256(
            self.evidence.read_bytes()).hexdigest()
        self.manifest.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(ValueError, 'source map'):
            self.check()

    def test_missing_source_bytes_rejected(self):
        with zipfile.ZipFile(self.archive, 'w'):
            pass
        with self.assertRaises(KeyError):
            self.check()

    def test_corrupt_source_bytes_rejected(self):
        with zipfile.ZipFile(self.archive) as z:
            blobs = {n: z.read(n) for n in z.namelist()}
        digest = next(iter(self.data['evidence'][self.evidence.name]['sources'].values()))
        blobs[digest] += b'changed'
        with zipfile.ZipFile(self.archive, 'w') as z:
            for name, raw in blobs.items():
                z.writestr(name, raw)
        with self.assertRaisesRegex(ValueError, 'source bytes'):
            self.check()

    def test_new_gap_cannot_be_authorized_by_manifest(self):
        path, digest = next(iter(self.data['evidence'][self.evidence.name]['sources'].items()))
        key = path + ':' + digest
        self.data['sources'][key] = self.data['known_gap']
        self.data['known_gap_key'] = key
        self.manifest.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(ValueError, 'Unrecognized archival gap'):
            self.check()


if __name__ == '__main__':
    unittest.main()
