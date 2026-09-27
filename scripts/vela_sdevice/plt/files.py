"""Atomic file writers for ``.plt`` and ``.pvd`` products (contract C6, 3.1/3.2).

Every write goes to ``<name>.tmp`` in the target directory and is published
with :func:`os.replace`, so readers never observe a partially written file and
a failed write leaves the previous product untouched.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence, Union

from .document import CURRENT_SIGN_FACTOR, build_plt_document, read_plt_document

PathLike = Union[str, "os.PathLike[str]"]


def write_plt_atomic(path: PathLike, text: str) -> None:
    _write_text_atomic(Path(path), text)


def append_plt(path: PathLike,
               electrodes: Sequence[str],
               new_rows: Iterable[Mapping[str, Any]],
               *,
               sign_factor: float = CURRENT_SIGN_FACTOR) -> list[Mapping[str, Any]]:
    """Append ``new_rows`` to the ``.plt`` at ``path`` (created when absent).

    The existing document is parsed with the same ``sign_factor`` so its rows
    come back in Vela convention, extended, and rewritten atomically.  Returns
    the complete row list that was written.
    """
    target = Path(path)
    names = list(electrodes)
    rows: list[Mapping[str, Any]] = []
    if target.exists():
        existing_electrodes, existing_rows = read_plt_document(
            target.read_text(encoding="utf-8"), sign_factor=sign_factor)
        if existing_electrodes != names:
            raise ValueError(
                f"{target} lists electrodes {existing_electrodes}, cannot append rows for {names}")
        rows.extend(existing_rows)
    rows.extend(new_rows)
    _write_text_atomic(target, build_plt_document(names, rows, sign_factor=sign_factor))
    return rows


def build_pvd_document(entries: Iterable[Mapping[str, Any]]) -> str:
    """ParaView collection: ``[{"timestep": int, "file": "name.vtu"}, ...]``."""
    root = ET.Element("VTKFile", {"type": "Collection", "version": "0.1", "byte_order": "LittleEndian"})
    collection = ET.SubElement(root, "Collection")
    seen: set[int] = set()
    for index, entry in enumerate(entries):
        timestep = entry["timestep"]
        file_name = entry["file"]
        if isinstance(timestep, bool) or not isinstance(timestep, int):
            raise ValueError(f"pvd entry {index}: timestep must be an int, got {timestep!r}")
        if timestep in seen:
            raise ValueError(f"pvd entry {index}: duplicate timestep {timestep}")
        if not isinstance(file_name, str) or not file_name:
            raise ValueError(f"pvd entry {index}: file must be a non-empty string, got {file_name!r}")
        seen.add(timestep)
        ET.SubElement(collection, "DataSet",
                      {"timestep": str(timestep), "group": "", "part": "0", "file": file_name})
    ET.indent(root, space="  ")
    return '<?xml version="1.0"?>\n' + ET.tostring(root, encoding="unicode") + "\n"


def write_pvd_atomic(path: PathLike, entries: Iterable[Mapping[str, Any]]) -> None:
    _write_text_atomic(Path(path), build_pvd_document(entries))


def _write_text_atomic(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


__all__ = ["append_plt", "build_pvd_document", "write_plt_atomic", "write_pvd_atomic"]
