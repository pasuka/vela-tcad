#!/usr/bin/env python3
"""Helpers for PN2D fixed-state diagnostic manifests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def resolve_manifest_path(raw: str | Path, manifest_path: Path | None = None) -> Path:
    path = Path(raw)
    if path.is_absolute() or manifest_path is None:
        return path
    return manifest_path.parent / path


def read_config_for_manifest_item(item: dict[str, Any], manifest_path: Path | None = None) -> tuple[Path, dict[str, Any]] | None:
    raw = item.get("config")
    if not raw:
        return None
    config_path = resolve_manifest_path(str(raw), manifest_path)
    return config_path, json.loads(config_path.read_text(encoding="utf-8-sig"))


def require_materials_file(item: dict[str, Any], manifest_path: Path | None = None) -> Path | None:
    loaded = read_config_for_manifest_item(item, manifest_path)
    if loaded is None:
        return None
    config_path, config = loaded
    raw = config.get("materials_file")
    if not raw:
        raise ValueError(
            f"fixed-state config '{config_path}' is missing materials_file; "
            "PN2D Sentaurus diagnostics must not fall back to Vela's built-in Si ni"
        )
    materials_path = Path(str(raw))
    if not materials_path.is_absolute():
        materials_path = config_path.parent / materials_path
    if not materials_path.exists():
        raise FileNotFoundError(
            f"fixed-state config '{config_path}' materials_file does not exist: {materials_path}"
        )
    return materials_path
