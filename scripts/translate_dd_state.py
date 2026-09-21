"""Translate a DD checkpoint while retaining split quasi-Fermi low bits.

Potential columns are volts. Carrier densities and quantum potentials are
unchanged. The caller supplies transport node IDs; inactive QF placeholders
remain zero. This utility does not qualify the translated state for a sweep.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable


def translate_state_csv(
    source: Path, destination: Path, subtract_V: float,
    transport_nodes: Iterable[int],
) -> None:
    if not math.isfinite(subtract_V):
        raise ValueError("Frame translation must be finite")
    active = set(transport_nodes)
    with Path(source).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        if columns is None or not {"node_id", "psi", "phin", "phip"}.issubset(columns):
            raise ValueError("Checkpoint requires node_id, psi, phin and phip")
        data = list(reader)
    for carrier in ("electron", "hole"):
        pair = {f"{carrier}_qf_reference_V", f"{carrier}_qf_increment_V"}
        if len(pair.intersection(columns)) == 1:
            raise ValueError("Referenced QF checkpoints require both reference and increment")
    for row in data:
        live = int(row["node_id"]) in active
        for carrier in ("electron", "hole"):
            ref_key = f"{carrier}_qf_reference_V"
            inc_key = f"{carrier}_qf_increment_V"
            if ref_key not in row:
                continue
            old_ref, old_inc = float(row[ref_key]), float(row[inc_key])
            new_ref = old_ref - subtract_V
            if not all(map(math.isfinite, (old_ref, old_inc, new_ref))):
                raise ValueError("QF reference and increment must be finite")
            # fsum retains the rounding remainder of old_ref - subtract_V
            # in the explicit increment, including when new_ref == old_ref.
            new_inc = math.fsum((old_ref, -subtract_V, -new_ref, old_inc))
            row[ref_key] = format(new_ref, ".17g") if live else "0"
            row[inc_key] = format(new_inc, ".17g") if live else "0"
        for key in ("psi", "phin", "phip"):
            value = float(row[key]) - subtract_V
            if not math.isfinite(value):
                raise ValueError("Checkpoint potentials must be finite")
            row[key] = format(value, ".17g") if key == "psi" or live else "0"
    # Refuse to overwrite an input or an existing experiment checkpoint.
    with Path(destination).open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(data)
