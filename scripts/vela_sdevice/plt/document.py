"""DF-ISE ``xyplot`` (``.plt``) text document builder and reader.

The layout mirrors the Sentaurus current files archived in the repository
fixtures (``pn2d_0v.plt`` in ``tests/regression/test_sentaurus_import_tools.py``):

* banner ``DF-ISE text``;
* ``Info { version = 1.0  type = xyplot  datasets = [...]  functions = [...] }``
  where ``datasets`` are the quoted ``"<electrode> <function>"`` names and
  ``functions`` the bare function names, ``"time"`` first and then five per
  line;
* ``Data { ... }`` with, per row, the ``time`` value on its own line followed
  by five values per line.

Numbers are written with :func:`repr`, the shortest text that reproduces the
binary64 value exactly, so a file read back (by :func:`read_plt_document` or by
``scripts/sentaurus_import.py``) yields bit-identical values.
"""

from __future__ import annotations

import math
import re
from typing import Any, Iterable, Mapping, Sequence

# Global factor applied to eCurrent/hCurrent/TotalCurrent on output (Vela
# terminal-current sign -> Sentaurus .plt sign). Frozen by WP10 from the
# verified T-2022.03-SP2 v2 pn2d_iv.plt (first three positive Anode points).
# SHA-256 a5f712e3ba0bfeea221e0de681ce380b3ad4ab14eb42c63b1324eac431709b25.
# Evidence and per-stage Time semantics: docs/sentaurus_sdevice_workflow.md 3.3.
CURRENT_SIGN_FACTOR = 1.0

TIME_DATASET = "time"
ELECTRODE_FUNCTIONS = ("OuterVoltage", "InnerVoltage", "eCurrent", "hCurrent", "TotalCurrent")
FUNCTION_KEYS = {
    "OuterVoltage": "outer_voltage_V",
    "InnerVoltage": "inner_voltage_V",
    "eCurrent": "e_current_A_per_um",
    "hCurrent": "h_current_A_per_um",
    "TotalCurrent": "total_current_A_per_um",
}
CURRENT_FUNCTIONS = frozenset({"eCurrent", "hCurrent", "TotalCurrent"})
ROW_KEYS = tuple(FUNCTION_KEYS[function] for function in ELECTRODE_FUNCTIONS)

_VALUES_PER_LINE = 5
_DATASETS_RE = re.compile(r"datasets\s*=\s*\[(.*?)\]", re.DOTALL)
_DATA_RE = re.compile(r"\bData\s*\{(.*?)\}", re.DOTALL)
_NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def build_plt_document(electrodes: Sequence[str],
                       rows: Iterable[Mapping[str, Any]],
                       *,
                       sign_factor: float = CURRENT_SIGN_FACTOR) -> str:
    """Render ``rows`` as a DF-ISE xyplot document.

    ``rows`` are ``{"time": t, "<electrode>": {<ROW_KEYS>}}`` in Vela units
    (V, A/um).  A ``None`` value declares a quantity as unavailable and is
    written as ``0.0``; missing keys are errors.  Currents are multiplied by
    ``sign_factor``.
    """
    names = validate_electrodes(electrodes)
    factor = _validated_factor(sign_factor)
    datasets = [TIME_DATASET] + [f"{name} {function}" for name in names for function in ELECTRODE_FUNCTIONS]
    functions = [TIME_DATASET] + [function for _ in names for function in ELECTRODE_FUNCTIONS]

    lines = ["DF-ISE text", "", "Info {", "  version   = 1.0", "  type      = xyplot", "  datasets  = ["]
    lines.extend(_bracket_list([f'"{dataset}"' for dataset in datasets]))
    lines.append("  functions = [")
    lines.extend(_bracket_list(functions))
    lines.extend(["}", "", "Data {"])
    for index, row in enumerate(rows):
        lines.extend(_row_lines(row, names, factor, index))
    lines.append("}")
    return "\n".join(lines) + "\n"


def read_plt_document(text: str,
                      *,
                      sign_factor: float = CURRENT_SIGN_FACTOR) -> tuple[list[str], list[dict[str, Any]]]:
    """Inverse of :func:`build_plt_document`: ``(electrodes, rows)``.

    Currents are divided by ``sign_factor`` so the rows come back in Vela
    convention.  Datasets beyond the five written functions (for example
    ``QuasiFermiPotential`` or ``Charge`` in a Sentaurus file) are kept under
    their Sentaurus function name, unscaled; the writer ignores them.
    """
    factor = _validated_factor(sign_factor)
    datasets_match = _DATASETS_RE.search(text)
    if datasets_match is None:
        raise ValueError("DF-ISE xyplot text does not declare datasets")
    datasets = re.findall(r'"([^"]*)"', datasets_match.group(1))
    if not datasets or datasets[0] != TIME_DATASET:
        raise ValueError(f"DF-ISE xyplot datasets must start with \"{TIME_DATASET}\", got {datasets[:1]}")
    data_match = _DATA_RE.search(text)
    if data_match is None:
        raise ValueError("DF-ISE xyplot text does not contain a Data block")
    body = re.sub(r"Values\s*=\s*\[", " ", data_match.group(1)).replace("]", " ")
    numbers = [float(token) for token in _NUMBER_RE.findall(body)]
    if len(numbers) % len(datasets) != 0:
        raise ValueError(
            f"DF-ISE Data block has {len(numbers)} values, not a multiple of {len(datasets)} datasets")

    electrodes: list[str] = []
    columns: list[tuple[str, str]] = []
    for dataset in datasets[1:]:
        electrode, separator, function = dataset.rpartition(" ")
        if not separator or not electrode or not function:
            raise ValueError(f"unrecognised DF-ISE dataset name {dataset!r}")
        if electrode not in electrodes:
            electrodes.append(electrode)
        columns.append((electrode, function))

    rows: list[dict[str, Any]] = []
    width = len(datasets)
    for start in range(0, len(numbers), width):
        chunk = numbers[start:start + width]
        row: dict[str, Any] = {TIME_DATASET: chunk[0]}
        for electrode in electrodes:
            row[electrode] = {}
        for (electrode, function), value in zip(columns, chunk[1:]):
            if function in CURRENT_FUNCTIONS:
                value = value / factor
                if value == 0.0:
                    value = 0.0
            row[electrode][FUNCTION_KEYS.get(function, function)] = value
        rows.append(row)
    return electrodes, rows


def validate_electrodes(electrodes: Sequence[str]) -> list[str]:
    names = list(electrodes)
    if not names:
        raise ValueError("at least one electrode is required")
    seen: set[str] = set()
    for name in names:
        if not isinstance(name, str) or not name or name == TIME_DATASET or '"' in name \
                or any(character.isspace() for character in name):
            raise ValueError(f"invalid electrode name {name!r}")
        if name in seen:
            raise ValueError(f"duplicate electrode name {name!r}")
        seen.add(name)
    return names


def format_number(value: Any) -> str:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite value {value!r} cannot be written to a DF-ISE file")
    return repr(number)


def _validated_factor(sign_factor: float) -> float:
    factor = float(sign_factor)
    if not math.isfinite(factor) or factor == 0.0:
        raise ValueError(f"sign_factor must be finite and non-zero, got {sign_factor!r}")
    return factor


def _bracket_list(tokens: Sequence[str]) -> list[str]:
    lines = ["    " + tokens[0]]
    for start in range(1, len(tokens), _VALUES_PER_LINE):
        lines.append("    " + " ".join(tokens[start:start + _VALUES_PER_LINE]))
    lines[-1] += " ]"
    return lines


def _row_lines(row: Mapping[str, Any], names: Sequence[str], factor: float, index: int) -> list[str]:
    try:
        values = [float(row[TIME_DATASET])]
        for name in names:
            entry = row[name]
            for function in ELECTRODE_FUNCTIONS:
                value = entry[FUNCTION_KEYS[function]]
                number = 0.0 if value is None else float(value)
                values.append(number * factor if function in CURRENT_FUNCTIONS else number)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"plt row {index} is incomplete or not numeric: {exc!r}") from exc
    values = [0.0 if value == 0.0 else value for value in values]  # never emit -0.0
    lines = ["      " + format_number(values[0])]
    for start in range(1, len(values), _VALUES_PER_LINE):
        lines.append("      " + "   ".join(format_number(value) for value in values[start:start + _VALUES_PER_LINE]))
    return lines
