"""Verify the accepted state that carries extended80 controls between stages."""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from hashlib import sha256
import json
import math
import re
import subprocess

from ..cmd_parser.syntax import SourceNumber
from .support import RunError


def control_decimal(value):
    if isinstance(value, SourceNumber):
        return Decimal(value.literal)
    return Decimal.from_float(value) if isinstance(value, float) else Decimal(value)


def _source_control(text):
    if not isinstance(text, str) or not text or len(text) > 128:
        raise ValueError("expected a decimal control string")
    value = Decimal(text)
    if not value.is_finite():
        raise ValueError("nonfinite control")
    view = SourceNumber(text)
    if not math.isfinite(view):
        raise ValueError("control is outside the output coordinate range")
    return view


def _file_sha256(path):
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass
class AcceptedExtendedStage:
    record: dict
    voltages: dict
    outer: SourceNumber
    currents: dict


def verify_extended_stage_parent(work, archive, parent):
    """Check the selected inter-stage archive immediately before starting the runner."""
    try:
        if parent is None or parent.record["archive"] != archive:
            raise ValueError("missing verified accepted parent")
        if parent.record["archive_sha256"] != _file_sha256(work / archive):
            raise ValueError("parent archive changed after acceptance")
    except (KeyError, OSError, ValueError) as error:
        raise RunError("E202", f"Cannot use extended80 stage parent '{archive}': {error}") from error


def read_accepted_extended_stage(work, stage, summary, records, sweep_contact, runner):
    """Use the exact round-trip decimal fields bound to an accepted archive."""
    try:
        record = summary["accepted_extended80_state"]
        if not isinstance(record, dict) or record.get("schema") != "vela.extended80-stage-state.v1":
            raise ValueError("missing extended80 accepted state")
        if (summary["arithmetic_precision"] != "extended80" or record["role"] != "accepted" or
                record["archive"] != stage["output_state"] or record["accepted_points"] != len(records)):
            raise ValueError("accepted state identity mismatch")
        for name, key in ((stage["deck"], "deck_sha256"), (stage["output_state"], "archive_sha256"),
                          (stage["csv"], "csv_sha256")):
            if record[key] != _file_sha256(work / name):
                raise ValueError(f"{key} mismatch")
        voltages = {name: _source_control(value) for name, value in record["voltages_V"].items()}
        currents = {name: _source_control(value) for name, value in record["currents_A_per_um"].items()}
        if (set(voltages) != set(stage["electrode_voltages_start_V"]) or
                sweep_contact not in currents or (set(currents) - set(voltages))):
            raise ValueError("accepted control electrode map mismatch")
        for criterion in stage.get("break", {}).get("current", []):
            if criterion["contact"] not in currents:
                raise ValueError("accepted state lacks a configured break current")
        outer = _source_control(record["outer_voltage_V"])
        identities = {key: record[key] for key in ("mesh_sha256", "input_sha256")}
        if not all(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value)
                   for value in identities.values()):
            raise ValueError("invalid archive identity")
        node_count = record["node_count"]
        if type(node_count) is not int or node_count <= 0:
            raise ValueError("invalid archive mesh size")
        request = {"schema": "vela.verify-extended80-stage-state.v1", **identities,
                   "archive": stage["output_state"], "node_count": node_count,
                   "device_voltage_V": voltages[sweep_contact].literal,
                   "supply_voltage_V": outer.literal,
                   "current_A_per_um": currents[sweep_contact].literal}
        result = subprocess.run(runner + ["--verify-extended80-stage-state"],
                                cwd=work, input=json.dumps(request), capture_output=True,
                                encoding="utf-8", timeout=45)
        if result.returncode:
            raise ValueError("archive verification failed: " + result.stderr.strip()[:400])
        if json.loads(result.stdout) != {"verified": True, "arithmetic_precision": "extended80"}:
            raise ValueError("runner did not verify the accepted archive")
        return AcceptedExtendedStage(record, voltages, outer, currents)
    except (AttributeError, InvalidOperation, KeyError, OSError, subprocess.SubprocessError,
            TypeError, ValueError) as error:
        raise RunError("E202", f"Cannot verify accepted extended80 stage state: {error}") from error
