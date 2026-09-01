"""Helpers for validating superseded SimpleMOS frozen source hashes."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from unittest import TestCase


REPO = Path(__file__).resolve().parents[2]
M44_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m44_qf_coordinate_consistency_evidence.json")
M45_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m45_post_qf_rebaseline_evidence.json")
M46_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m46_full_matrix_requalification_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assert_source_hashes_current_or_m44(
        testcase: TestCase, evidence: dict) -> None:
    """Accept an old hash when frozen M44, M45, or M46 pins its successor."""
    assert_hash_mapping_current_or_superseded(
        testcase, evidence["source_hashes"])


def assert_hash_mapping_current_or_superseded(
        testcase: TestCase, expected_hashes: dict[str, str]) -> None:
    """Validate a source-hash map against current or superseding evidence."""
    m44 = json.loads(M44_EVIDENCE.read_text(encoding="utf-8"))
    testcase.assertEqual(m44["status"], "frozen")
    superseding = [m44]
    if M46_EVIDENCE.is_file():
        m46 = json.loads(M46_EVIDENCE.read_text(encoding="utf-8"))
        testcase.assertEqual(m46["status"], "frozen")
        superseding.insert(0, m46)
    if M45_EVIDENCE.is_file():
        m45 = json.loads(M45_EVIDENCE.read_text(encoding="utf-8"))
        testcase.assertEqual(m45["status"], "frozen")
        superseding.insert(0, m45)
    for relative, expected in expected_hashes.items():
        current = sha256(REPO / relative)
        if expected == current:
            continue
        testcase.assertTrue(
            any(item["source_hashes"].get(relative) == current
                for item in superseding),
            f"changed source is not frozen by superseding M44/M45/M46: {relative}")
