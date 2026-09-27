"""Build ``.plt`` rows from the C++ runner's sweep CSV files (contract C5 -> C6).

Column names come from ``src/simulation/DCSweep.cpp``.  Main sweep CSV
(header assembled near ``std::vector<std::string> header = {"mode", ...``):

``bias_V``
    Device (inner) voltage of the swept contact for the point.
``bias_contact`` / ``current_contact``
    Swept contact and the contact whose currents the row reports.  The main
    CSV currents are attributed to ``current_contact`` (falls back to
    ``sweep_contact`` when the column is absent).
``converged``
    ``1``/``0``; only ``1`` rows are accepted steps and become ``.plt`` rows.
``current_total_A_per_um``, ``current_electron_A_per_um``,
``current_hole_A_per_um``
    Present when the deck uses ``unit_scaling`` (``writeUnitScaledColumns``);
    the total is required, the carrier splits are optional.
``outer_voltage_V``, ``inner_voltage_V``
    Present only with ``sweep.external_resistor`` / ``voltage_to_current``
    (``boundaryControlEnabled``).  When present the swept electrode's
    OuterVoltage/InnerVoltage are taken from them and ``t`` is computed from
    the outer voltage; otherwise both equal ``bias_V``.

Terminal-balance CSV (``sweep.diagnostics.terminal_balance.enabled``,
default ``<csv stem>_terminal_balance.csv``): long format with one row per
(point, contact) -- ``point_index``, ``bias_V``, ``bias_contact``,
``contact``, ``current_electron_A_per_um``, ``current_hole_A_per_um``,
``current_total_A_per_um``, ``converged``.  Rows are joined to the main CSV
by ``bias_V`` and fill the currents of the electrodes the main CSV does not
report (missing electrodes keep ``0.0``).

Every physical column name can be overridden through ``column_map`` keyed by
the logical names in :data:`DEFAULT_COLUMN_MAP`.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Union

from .document import TIME_DATASET

PathLike = Union[str, "os.PathLike[str]"]

DEFAULT_COLUMN_MAP: dict[str, str] = {
    "bias_V": "bias_V",
    "outer_voltage_V": "outer_voltage_V",
    "inner_voltage_V": "inner_voltage_V",
    "current_contact": "current_contact",
    "converged": "converged",
    "total_current_A_per_um": "current_total_A_per_um",
    "e_current_A_per_um": "current_electron_A_per_um",
    "h_current_A_per_um": "current_hole_A_per_um",
    "terminal_balance_contact": "contact",
}

ALL_ACCEPTED_STEPS = "all_accepted_steps"
T_MATCH_TOLERANCE = 1.0e-9


def filter_rows(rows: Iterable[Mapping[str, Any]],
                plt_write: Union[str, Sequence[float]]) -> list[Mapping[str, Any]]:
    """Apply a stage's ``plt_write`` selection (``"all_accepted_steps"`` or t list)."""
    if isinstance(plt_write, str):
        if plt_write != ALL_ACCEPTED_STEPS:
            raise ValueError(f"unsupported plt_write value {plt_write!r}")
        return list(rows)
    requested = [float(t) for t in plt_write]
    return [row for row in rows
            if any(abs(float(row[TIME_DATASET]) - t) <= T_MATCH_TOLERANCE for t in requested)]


def rows_from_runner_csv(main_csv_path: PathLike,
                         *,
                         sweep_contact: str,
                         electrodes: Sequence[str],
                         t_to_bias: tuple[float, float],
                         terminal_balance_csv_path: Optional[PathLike] = None,
                         column_map: Optional[Mapping[str, str]] = None,
                         electrode_voltages: Optional[Mapping[str, float]] = None) -> list[dict[str, Any]]:
    """Convert accepted runner points into ``.plt`` rows for one stage.

    ``time = (V_outer - v0) / (v1 - v0)`` with ``(v0, v1) = t_to_bias``;
    ``v1 == v0`` (equilibrium stage) gives ``time = 0.0``.  ``electrode_voltages``
    supplies the fixed voltages of the non-swept electrodes (stage plan
    ``electrode_voltages_start_V``); missing electrodes get ``0.0``.  Carrier
    currents that no CSV provides are ``None`` (written as ``0.0`` by the
    document builder); other electrodes without terminal-balance data are ``0.0``.
    """
    names = list(electrodes)
    if sweep_contact not in names:
        raise ValueError(f"sweep contact {sweep_contact!r} is not one of the electrodes {names}")
    columns = dict(DEFAULT_COLUMN_MAP)
    columns.update(column_map or {})
    v0, v1 = (float(t_to_bias[0]), float(t_to_bias[1]))
    fixed = {name: float(value) for name, value in (electrode_voltages or {}).items()}

    main_rows, main_fields = _read_csv(main_csv_path)
    for logical in ("bias_V", "total_current_A_per_um"):
        if columns[logical] not in main_fields:
            raise ValueError(
                f"main sweep CSV {os.fspath(main_csv_path)} lacks column {columns[logical]!r}; "
                f"available: {main_fields}")
    balance = _read_terminal_balance(terminal_balance_csv_path, columns) if terminal_balance_csv_path else {}

    rows: list[dict[str, Any]] = []
    for record in main_rows:
        if not _accepted(record, columns["converged"]):
            continue
        bias = float(record[columns["bias_V"]])
        outer = _optional_float(record, columns["outer_voltage_V"])
        inner = _optional_float(record, columns["inner_voltage_V"])
        outer = bias if outer is None else outer
        inner = bias if inner is None else inner
        time = 0.0 if v1 == v0 else (outer - v0) / (v1 - v0)
        if time == 0.0:
            time = 0.0  # (0 - v0) / negative span yields -0.0

        current_contact = (record.get(columns["current_contact"]) or "").strip() or sweep_contact
        if current_contact not in names:
            raise ValueError(f"current contact {current_contact!r} is not one of the electrodes {names}")

        row: dict[str, Any] = {TIME_DATASET: time}
        for name in names:
            voltage = (outer, inner) if name == sweep_contact else (fixed.get(name, 0.0),) * 2
            row[name] = {
                "outer_voltage_V": voltage[0],
                "inner_voltage_V": voltage[1],
                "e_current_A_per_um": 0.0,
                "h_current_A_per_um": 0.0,
                "total_current_A_per_um": 0.0,
            }
        row[current_contact].update({
            "e_current_A_per_um": _optional_float(record, columns["e_current_A_per_um"]),
            "h_current_A_per_um": _optional_float(record, columns["h_current_A_per_um"]),
            "total_current_A_per_um": float(record[columns["total_current_A_per_um"]]),
        })
        for name, currents in _balance_at(balance, bias).items():
            if name not in names:
                continue
            if name == current_contact:
                for key in ("e_current_A_per_um", "h_current_A_per_um"):
                    if row[name][key] is None:
                        row[name][key] = currents[key]
            else:
                row[name].update(currents)
        rows.append(row)
    return rows


def _read_csv(path: PathLike) -> tuple[list[dict[str, str]], list[str]]:
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        records = list(reader)
        return records, list(reader.fieldnames or [])


def _read_terminal_balance(path: PathLike,
                           columns: Mapping[str, str]) -> dict[float, dict[str, dict[str, Optional[float]]]]:
    records, fields = _read_csv(path)
    for logical in ("bias_V", "terminal_balance_contact", "total_current_A_per_um"):
        if columns[logical] not in fields:
            raise ValueError(
                f"terminal balance CSV {os.fspath(path)} lacks column {columns[logical]!r}; available: {fields}")
    balance: dict[float, dict[str, dict[str, Optional[float]]]] = {}
    for record in records:
        if not _accepted(record, columns["converged"]):
            continue
        bias = float(record[columns["bias_V"]])
        contact = record[columns["terminal_balance_contact"]].strip()
        balance.setdefault(bias, {})[contact] = {
            "e_current_A_per_um": _optional_float(record, columns["e_current_A_per_um"]),
            "h_current_A_per_um": _optional_float(record, columns["h_current_A_per_um"]),
            "total_current_A_per_um": float(record[columns["total_current_A_per_um"]]),
        }
    return balance


def _balance_at(balance: Mapping[float, dict[str, dict[str, Optional[float]]]],
                bias: float) -> dict[str, dict[str, Optional[float]]]:
    if bias in balance:
        return balance[bias]
    for key, value in balance.items():
        if math.isclose(key, bias, rel_tol=1.0e-12, abs_tol=1.0e-12):
            return value
    return {}


def _accepted(record: Mapping[str, str], converged_column: str) -> bool:
    value = record.get(converged_column)
    if value is None:
        return True
    return value.strip().lower() in {"1", "true", "yes"}


def _optional_float(record: Mapping[str, str], column: str) -> Optional[float]:
    value = record.get(column)
    if value is None or not value.strip():
        return None
    return float(value)


__all__ = ["ALL_ACCEPTED_STEPS", "DEFAULT_COLUMN_MAP", "T_MATCH_TOLERANCE",
           "filter_rows", "rows_from_runner_csv"]
