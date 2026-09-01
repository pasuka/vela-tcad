"""Helpers for validating superseded SimpleMOS frozen source hashes."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from unittest import TestCase


REPO = Path(__file__).resolve().parents[2]
M44_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m44_qf_coordinate_consistency_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assert_source_hashes_current_or_m44(
        testcase: TestCase, evidence: dict) -> None:
    """Accept an old hash only when frozen M44 pins its current successor."""
    m44 = json.loads(M44_EVIDENCE.read_text(encoding="utf-8"))
    testcase.assertEqual(m44["status"], "frozen")
    for relative, expected in evidence["source_hashes"].items():
        current = sha256(REPO / relative)
        if expected == current:
            continue
        testcase.assertEqual(
            m44["source_hashes"].get(relative), current,
            f"changed source is not frozen by superseding M44: {relative}")
