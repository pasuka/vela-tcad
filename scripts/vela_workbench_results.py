"""Read the published curve and final two-dimensional Tri3 field."""
from __future__ import annotations

from dataclasses import dataclass
import xml.etree.ElementTree as ET

import numpy as np

from vela_sdevice.output_metadata import accepted_records, vtu_metadata
from vela_sdevice.plt import read_plt_document


FIELDS = {
    "ElectrostaticPotential": "Potential",
    "ElectricField": "Electric field magnitude",
    "eDensity": "Electron density",
    "hDensity": "Hole density",
    "DopingConcentration": "Net doping",
}


@dataclass
class Field:
    points: np.ndarray
    triangles: np.ndarray
    quantities: dict[str, np.ndarray]
    metadata: dict


def read_curve(path, voltage="inner"):
    _, rows = read_plt_document(path.read_text(encoding="utf-8"))
    column = f"{voltage}_voltage_V"
    return [(row["Anode"].get(column), row["Anode"]["total_current_A_per_um"])
            for row in rows]


def read_stage_curves(work, stages, voltage="inner"):
    """Keep separate swept terminals and bias axes for multistage decks."""
    curves = []
    column = f"{voltage}_voltage_V"
    for index, stage in enumerate(stages, 1):
        if "csv" not in stage:
            continue
        rows = accepted_records(work / stage["csv"])
        if not rows:
            continue
        sweep = stage.get("sweep") or {}
        terminal = sweep.get("contact") or rows[0]["current_contact"]
        # Equilibrium rows can omit voltage columns because the contacts are fixed.
        fixed_voltage = None
        if stage.get("kind") == "equilibrium" and not sweep:
            fixed_voltage = (stage.get("electrode_voltages_goal_V") or {}).get(terminal)
        points = []
        for row in rows:
            current = row.get(f"contact_{terminal}_total_A_per_um")
            if not current and row["current_contact"] == terminal:
                current = row.get("current_total_A_per_um")
            if current:
                value = row.get(column)
                if value in (None, ""):
                    value = row.get("bias_V") if voltage == "outer" else fixed_voltage
                points.append((float(value) if value not in (None, "") else None, float(current)))
        if points:
            curves.append((index, terminal, points))
    return curves


def display_coordinates(points, quarter_turns):
    x, y = points[:, :2].T
    if quarter_turns == 0:
        return x, y
    if quarter_turns == 1:
        return -y, x
    if quarter_turns == 2:
        return -x, -y
    return y, -x


def _array(parent, name=None, dtype=float):
    expression = "DataArray" if name is None else f"DataArray[@Name='{name}']"
    item = parent.find(expression) if parent is not None else None
    if item is None or item.get("format", "ascii") != "ascii":
        raise ValueError(f"Expected an ASCII VTU array: {name or 'Points'}")
    return np.fromstring(item.text or "", sep=" ", dtype=dtype)


def read_field(path):
    try:
        piece = ET.parse(path).getroot().find("./UnstructuredGrid/Piece")
    except ET.ParseError as error:
        raise ValueError("Incomplete or invalid VTU file") from error
    if piece is None:
        raise ValueError("VTU is missing UnstructuredGrid/Piece")
    points = _array(piece.find("Points")).reshape(-1, 3)
    cells = piece.find("Cells")
    types = _array(cells, "types", int)
    if not len(types) or not np.all(types == 5):
        raise ValueError("Expected a 2D Tri3 mesh")
    triangles = _array(cells, "connectivity", int).reshape(-1, 3)
    if len(triangles) != len(types) or triangles.min() < 0 or triangles.max() >= len(points):
        raise ValueError("Invalid VTU triangle connectivity")

    data = piece.find("PointData")
    quantities = {}
    if data is not None:
        for name in FIELDS:
            if data.find(f"DataArray[@Name='{name}']") is None:
                continue
            values = _array(data, name)
            if name == "ElectricField":
                values = np.linalg.norm(values.reshape(-1, 3), axis=1)
            if len(values) != len(points):
                raise ValueError(f"VTU point data length mismatch: {name}")
            quantities[name] = values
    return Field(points, triangles, quantities, vtu_metadata(path) or {})
