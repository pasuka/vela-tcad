"""Historical source integrity; these checks do not qualify the live solver."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from unittest import TestCase
import zipfile

REPO = Path(__file__).resolve().parents[2]
ARCHIVE = REPO / 'tests/fixtures/simplemos_historical_sources/sources.zip'
MANIFEST = ARCHIVE.with_name('manifest.json')
KNOWN_GAP = ('src/equation/CoupledDDAssembler.cpp:'
             'e525d828d0199dafdea947377dc6faed7cea419b5cb76fc1f51ae05d3dfe64ed')


def validate_historical_provenance(evidence_path: Path, manifest_path=MANIFEST,
                                 archive_path=ARCHIVE) -> list[dict]:
    """Check exact evidence/source bytes; return explicitly unrecovered sources.

    A recorded gap is metadata only, never a recovered blob or a live-source
    qualification. Unknown evidence, altered maps and absent members fail.
    """
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    if manifest['schema'] != 'vela.simplemos.historical-source-provenance/1':
        raise ValueError('Unsupported historical source manifest')
    raw = evidence_path.read_bytes()
    record = manifest['evidence'][evidence_path.name]
    if hashlib.sha256(raw).hexdigest() != record['sha256']:
        raise ValueError('Historical evidence bytes changed')
    evidence = json.loads(raw)
    sources = evidence.get('source_hashes', evidence.get('implementation_sha256'))
    if not sources or sources != record['sources']:
        raise ValueError('Historical source map changed')
    gaps = []
    with zipfile.ZipFile(archive_path) as archive:
        for path, digest in sources.items():
            key = path + ':' + digest
            entry = manifest['sources'][key]
            if entry['status'] == 'unavailable':
                if (entry != manifest['known_gap'] or key != KNOWN_GAP
                        or key != manifest['known_gap_key']):
                    raise ValueError('Unrecognized archival gap')
                gaps.append(dict(path=path, sha256=digest, **entry))
                continue
            if entry['status'] != 'recovered':
                raise ValueError('Unknown source status')
            if hashlib.sha256(archive.read(digest)).hexdigest() != digest:
                raise ValueError('Historical source bytes changed: ' + path)
    return gaps


def assert_historical_source_provenance(testcase: TestCase, evidence_path: Path):
    """Validate the archive, with an explicit record of unavailable provenance."""
    try:
        return validate_historical_provenance(evidence_path)
    except (KeyError, ValueError, OSError, zipfile.BadZipFile) as exc:
        testcase.fail(str(exc))
