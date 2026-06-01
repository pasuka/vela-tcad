"""Python front-end for Vela TCAD.

The Python package intentionally exposes a small API surface and delegates mesh
I/O, simulation, and export work to the C++ Vela core through pybind11.
"""

import os


def _configure_windows_dll_search_path() -> None:
    # Python 3.8+ on Windows tightens DLL search rules for extension modules.
    # Add common MSYS2 UCRT64 runtime paths so _core can resolve libstdc++.
    if os.name != "nt" or not hasattr(os, "add_dll_directory"):
        return

    seen = set()
    candidates = [
        os.environ.get("MSYS2_UCRT64_BIN", ""),
        r"D:\msys64\ucrt64\bin",
        r"D:\msys64\usr\bin",
    ]
    for path in candidates:
        if not path:
            continue
        norm = os.path.normcase(os.path.abspath(path))
        if norm in seen or not os.path.isdir(path):
            continue
        seen.add(norm)
        os.add_dll_directory(path)


_configure_windows_dll_search_path()

from ._core import (  # noqa: F401
    DCSweep,
    DeviceMesh,
    MaterialDatabase,
    PoissonSimulation,
    load_mesh,
    run_dc_sweep,
    run_poisson,
    write_vtk,
)
from .curves import run_bv_curve, run_cv_curve, run_iv_curve  # noqa: F401


__all__ = [
    "DCSweep",
    "DeviceMesh",
    "MaterialDatabase",
    "PoissonSimulation",
    "load_mesh",
    "run_bv_curve",
    "run_cv_curve",
    "run_dc_sweep",
    "run_iv_curve",
    "run_poisson",
    "write_vtk",
]
