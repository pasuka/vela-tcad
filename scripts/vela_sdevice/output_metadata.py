"""Read complete, published solver records and VTU metadata."""
from __future__ import annotations

import csv
import io
from pathlib import Path
import re


def accepted_records(path):
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except FileNotFoundError:
        return []
    # A runner may be halfway through a row. Never publish that row.
    if not text.endswith("\n"):
        text = text[:text.rfind("\n") + 1]
    return [r for r in csv.DictReader(io.StringIO(text)) if None not in r and
            all(v is not None for v in r.values()) and r.get("converged", "1") in ("1", "true")]


def vtu_metadata(path):
    data = Path(path).read_bytes()
    marker = data.find(b"<AppendedData")
    prefix = data if marker < 0 else data[:marker]
    if not data.rstrip().endswith(b"</VTKFile>"):
        return None
    text = prefix.decode("utf-8")
    result = {}
    for dtype, name, value in re.findall(r'<DataArray type="([^"]+)" Name="([^"]+)"[^>]*>\s*([^<]*)</DataArray>', text):
        if dtype == "String":
            result[name] = bytes(int(v) for v in value.split() if int(v)).decode("utf-8")
        elif value.strip():
            try:
                result[name] = float(value.strip())
            except ValueError:
                pass  # ordinary mesh arrays contain multiple values
    return result
