"""DF-ISE ``.plt`` and ParaView ``.pvd`` writer library (WP8, contract C6).

Pure functions, no CLI and no solver access.  The orchestrator (WP5) feeds
them from the stage plan (contract C4) and the runner CSV files (C5):

* :func:`rows_from_runner_csv` -- accepted runner points -> ``.plt`` rows with
  ``time = (V - v0) / (v1 - v0)`` (``0.0`` for equilibrium stages);
* :func:`filter_rows` -- stage ``plt_write`` selection;
* :func:`build_plt_document` / :func:`read_plt_document` -- DF-ISE xyplot text
  (datasets ``time`` then ``<electrode> OuterVoltage InnerVoltage eCurrent
  hCurrent TotalCurrent`` per electrode, currents times
  :data:`CURRENT_SIGN_FACTOR`);
* :func:`write_plt_atomic`, :func:`append_plt`, :func:`write_pvd_atomic` --
  ``<name>.tmp`` + ``os.replace`` publication.

Runner CSV columns consumed (``src/simulation/DCSweep.cpp``): main sweep CSV
``bias_V``, ``current_contact``, ``converged``, ``current_total_A_per_um``
(required) plus optional ``current_electron_A_per_um``,
``current_hole_A_per_um``, ``outer_voltage_V``, ``inner_voltage_V``; terminal
balance CSV (long format, one row per point and contact) ``bias_V``,
``contact``, ``current_electron_A_per_um``, ``current_hole_A_per_um``,
``current_total_A_per_um``, ``converged``.  See :mod:`vela_sdevice.plt.runner_csv`
for the exact semantics and the ``column_map`` override.

Caveats: carrier currents the runner does not report are written as ``0.0``;
electrodes without terminal-balance rows get ``0.0`` currents; the ``time``
origin of multi-stage sweeps (reset per Quasistationary or cumulative) is
decided by WP10, this library only computes the per-stage ``t``.
"""

from __future__ import annotations

from .document import (CURRENT_SIGN_FACTOR, ELECTRODE_FUNCTIONS, FUNCTION_KEYS, ROW_KEYS, TIME_DATASET,
                       build_plt_document, read_plt_document)
from .files import append_plt, build_pvd_document, write_plt_atomic, write_pvd_atomic
from .runner_csv import ALL_ACCEPTED_STEPS, DEFAULT_COLUMN_MAP, filter_rows, rows_from_runner_csv

__all__ = [
    "ALL_ACCEPTED_STEPS",
    "CURRENT_SIGN_FACTOR",
    "DEFAULT_COLUMN_MAP",
    "ELECTRODE_FUNCTIONS",
    "FUNCTION_KEYS",
    "ROW_KEYS",
    "TIME_DATASET",
    "append_plt",
    "build_plt_document",
    "build_pvd_document",
    "filter_rows",
    "read_plt_document",
    "rows_from_runner_csv",
    "write_plt_atomic",
    "write_pvd_atomic",
]
