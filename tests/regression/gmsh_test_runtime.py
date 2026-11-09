"""Shared cross-platform Gmsh runtime discovery for integration tests."""

from __future__ import annotations

import unittest

from gmsh_runtime import GmshRuntimeError, RuntimeSources, discover_runtime_sources


def require_gmsh_runtime() -> RuntimeSources:
    """Return the configured runtime or skip a developer machine without it."""
    try:
        return discover_runtime_sources()
    except GmshRuntimeError as error:
        raise unittest.SkipTest(
            f"Gmsh 4.15.2 runtime is not installed or configured: {error}"
        ) from error
