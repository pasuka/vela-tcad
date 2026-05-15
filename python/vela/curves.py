"""Lightweight Python helpers for Vela device curve sweeps.

These helpers intentionally keep simulation work in the C++ core. They only
validate the requested sweep mode and forward the configuration file to the
pybind11-backed :func:`vela.run_dc_sweep` binding.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

from ._core import run_dc_sweep

CurvePoint = Mapping[str, Any]


def _config_path(config: str | Path | Mapping[str, Any]) -> tuple[str, tempfile.TemporaryDirectory[str] | None]:
    if isinstance(config, Mapping):
        tmpdir = tempfile.TemporaryDirectory(prefix="vela_curve_config_")
        path = Path(tmpdir.name) / "curve_config.json"
        path.write_text(json.dumps(config), encoding="utf-8")
        return str(path), tmpdir
    return str(config), None


def _sweep_mode(config_file: str) -> str:
    with open(config_file, encoding="utf-8") as handle:
        cfg = json.load(handle)
    mode = cfg.get("sweep", {}).get("mode", "iv")
    return str(mode).lower().replace("-", "_")


def _run_curve(config: str | Path | Mapping[str, Any], accepted_modes: Sequence[str]) -> list[CurvePoint]:
    config_file, tmpdir = _config_path(config)
    try:
        mode = _sweep_mode(config_file)
        normalized_modes = {item.lower().replace("-", "_") for item in accepted_modes}
        if mode not in normalized_modes:
            expected = ", ".join(sorted(normalized_modes))
            raise ValueError(f"Expected sweep mode {expected}, got {mode!r}.")
        return run_dc_sweep(config_file)
    finally:
        if tmpdir is not None:
            tmpdir.cleanup()


def run_iv_curve(config: str | Path | Mapping[str, Any]) -> list[CurvePoint]:
    """Run an IV curve through the C++ DC sweep implementation."""

    return _run_curve(config, ("iv", ""))


def run_cv_curve(config: str | Path | Mapping[str, Any]) -> list[CurvePoint]:
    """Run a quasi-static CV curve through the C++ DC sweep implementation."""

    return _run_curve(config, ("cv", "cv_quasistatic"))


def run_bv_curve(config: str | Path | Mapping[str, Any]) -> list[CurvePoint]:
    """Run a reverse-bias BV diagnostic curve through the C++ DC sweep implementation."""

    return _run_curve(config, ("bv", "bv_reverse", "reverse_breakdown"))
