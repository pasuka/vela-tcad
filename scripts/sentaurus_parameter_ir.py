#!/usr/bin/env python3
"""Sentaurus ``.par`` -> Vela parameter IR (``vela.sentaurus_parameter_ir.v1``).

This is the third IR in the Sentaurus frontend, alongside
``sentaurus_device_ir`` (SDE geometry) and ``sentaurus_execution_ir`` (the
models an SDevice deck actually activates).  It answers the remaining question:
*what numeric values do those activated models use?*

Scope discipline
================

The module is a **syntax** frontend only.  It performs no unit conversion and
no physics interpretation: every value is kept as its original lexeme, its
original unit annotation, and a parsed numeric form.  Deciding which sections
matter, which formula variant is active, and how a value maps onto a Vela field
is the emitter's job, and the emitter must be driven by the execution IR rather
than by the ``.par`` file.  A ``.par`` file is a *candidate library*, not a list
of enabled models.

Losslessness
============

The IR keeps an ordered list of blocks rather than a dictionary keyed by
section name.  A dictionary would silently drop repeated sections, the relative
order of definitions, and the override order introduced by ``#include``.  Each
block records its scope, its source span, and the include stack that reached
it; each parameter records its raw name, raw lexeme, raw unit, and value kind.
Later definitions do not overwrite earlier ones -- both are retained and the
earlier one is annotated with ``shadowed_by``.

Fail-closed
===========

Valid ``.par`` syntax that this module cannot classify is an error, not a
warning: :func:`parse_parameter_ir` raises :class:`ParParseError` and no IR is
produced.  Sections that parse cleanly but that Vela has no model for are kept
as ordinary blocks and additionally listed in ``unsupported_sections`` so an
emitter can reject them by name.  Only comments and blank lines are ignored.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence


PARAMETER_IR_SCHEMA = "vela.sentaurus_parameter_ir.v1"

#: Maximum ``#include`` nesting depth before the reader gives up.  A cycle is
#: detected exactly by the include stack, so this only guards pathological but
#: acyclic fan-out.
MAX_INCLUDE_DEPTH = 32

_NUMBER_RE = re.compile(
    r"^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$")

# ``Material = "Silicon" {`` / ``Region = "drift" {``.  The trailing brace may
# sit on the following line.
_SCOPE_HEADER_RE = re.compile(
    r"""^(?P<kind>Material|Region|MaterialInterface|RegionInterface)
         \s*=\s*"(?P<name>[^"]*)"\s*(?P<brace>\{)?\s*$""",
    re.VERBOSE)

# ``Bandgap``/``ConstantMobility:``/``Auger * coefficients:``/
# ``QuantumPotentialParameters "theta_zero" {``.
_SECTION_HEADER_RE = re.compile(
    r"""^(?P<name>[A-Za-z_][A-Za-z0-9_]*)
         (?:\s*"(?P<variant>[^"]*)")?
         (?P<comment>\s*\*[^{]*?)?
         \s*(?P<colon>:)?
         \s*(?P<brace>\{)?\s*$""",
    re.VERBOSE)

_INCLUDE_RE = re.compile(r'^#include\s+"(?P<path>[^"]+)"\s*$')

# ``dEg0(OldSlotboom)``, ``Dop[1][1]``, ``efit( 0)``, ``(tau_w)_ele``,
# ``1/kappa``.  The name is captured whole and decomposed separately so that no
# information is lost when a shape is not recognised.
_PAREN_VARIANT_RE = re.compile(
    r"^(?P<base>[A-Za-z_][A-Za-z0-9_/]*)\(\s*(?P<variant>[^()]*?)\s*\)$")
_INDEXED_RE = re.compile(
    r"^(?P<base>[A-Za-z_][A-Za-z0-9_/]*)(?P<indices>(?:\[\s*-?\d+\s*\])+)$")
_WRAPPED_RE = re.compile(
    r"^\(\s*(?P<base>[A-Za-z_][A-Za-z0-9_]*)\s*\)(?P<suffix>_[A-Za-z0-9_]+)$")
_PLAIN_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(?:/[A-Za-z_][A-Za-z0-9_]*)?$")
_RECIPROCAL_RE = re.compile(r"^1/(?P<base>[A-Za-z_][A-Za-z0-9_]*)$")

# ``eValley"Delta1"(1,0,0)(ml=0.914 ...)`` -- a bare structured record with no
# top-level ``=``.
_RECORD_RE = re.compile(
    r"""^(?P<base>[A-Za-z_][A-Za-z0-9_]*)
         \s*"(?P<label>[^"]*)"
         \s*(?P<groups>\(.*\))\s*$""",
    re.VERBOSE)

# A ``NumericalTable`` written without a left-hand side.
_BARE_TABLE_RE = re.compile(r"^NumericalTable\s*\(.*\)$",
                            re.DOTALL | re.IGNORECASE)


class ParParseError(ValueError):
    """Raised when a ``.par`` file cannot be translated fail-closed.

    The message always carries ``file:line`` plus the offending source text so
    that a caller never has to reconstruct the location itself.
    """


@dataclass(frozen=True)
class SourceRef:
    """A resolved location inside a concrete ``.par`` file."""

    file: str
    line: int

    def to_json(self) -> dict[str, Any]:
        return {"file": self.file, "line": self.line}

    def __str__(self) -> str:  # pragma: no cover - diagnostics only
        return f"{self.file}:{self.line}"


@dataclass
class PhysicalLine:
    """One logical statement together with everything needed to report it.

    ``text`` is the statement with comments removed but with the original
    spelling otherwise intact.  ``raw`` keeps the untouched source so error
    messages can quote exactly what the author wrote.
    """

    text: str
    raw: str
    file: str
    line: int
    end_line: int
    include_stack: tuple[SourceRef, ...]


def _strip_comments(line: str) -> tuple[str, str | None]:
    """Split a physical line into code and its trailing unit annotation.

    ``.par`` uses two comment markers with different meanings:

    * ``*`` introduces prose documentation, and
    * ``#`` introduces the unit annotation such as ``# [cm^2/(Vs)]``.

    The unit is data, not decoration, so it is returned instead of discarded.
    ``#include`` is not a comment and is handled before this function runs.
    """
    unit: str | None = None
    hash_index = line.find("#")
    if hash_index >= 0:
        comment = line[hash_index + 1:].strip()
        line = line[:hash_index]
        match = re.match(r"^\[(?P<unit>.*)\]$", comment)
        unit = match.group("unit") if match else (comment or None)

    star_index = _find_prose_comment(line)
    if star_index >= 0:
        line = line[:star_index]
    return line.rstrip(), unit


def _find_prose_comment(line: str) -> int:
    """Return the index of a prose ``*`` comment, or ``-1``.

    A ``*`` only starts a comment at the beginning of a statement or after
    whitespace; ``a*b`` inside an expression must not be truncated.  Anything
    inside double quotes is literal.
    """
    in_string = False
    for index, char in enumerate(line):
        if char == '"':
            in_string = not in_string
            continue
        if in_string or char != "*":
            continue
        if index == 0 or line[index - 1] in " \t{":
            return index
    return -1


def _read_physical_lines(path: Path,
                         include_stack: tuple[SourceRef, ...],
                         seen: tuple[Path, ...],
                         root: Path) -> list[PhysicalLine]:
    """Read ``path`` and splice in every ``#include`` at its point of use.

    Include semantics that matter for losslessness:

    * the included text is spliced in place, so it inherits whatever scope is
      open at the ``#include`` site (a fragment included inside
      ``Material = "Silicon" {`` belongs to that material);
    * every spliced line keeps the *including* file's location in its include
      stack, so a value can always be traced back; and
    * a cycle is detected against the active stack rather than a global visited
      set, because including the same fragment twice from different places is
      legal while including it from itself is not.
    """
    resolved = path.resolve()
    if resolved in seen:
        chain = " -> ".join(item.name for item in seen) + f" -> {resolved.name}"
        raise ParParseError(
            f"{_stack_prefix(include_stack)}#include cycle detected: {chain}")
    if len(seen) > MAX_INCLUDE_DEPTH:
        raise ParParseError(
            f"{_stack_prefix(include_stack)}#include nesting deeper than "
            f"{MAX_INCLUDE_DEPTH} levels: {resolved.name}")
    if not resolved.is_file():
        raise ParParseError(
            f"{_stack_prefix(include_stack)}#include target does not exist: "
            f"{path}")

    display = _display_path(resolved, root)
    text = resolved.read_text(encoding="utf-8", errors="strict")
    lines: list[PhysicalLine] = []
    for number, raw in enumerate(text.splitlines(), start=1):
        stripped = raw.strip()
        include = _INCLUDE_RE.match(stripped)
        if include:
            target = (resolved.parent / include.group("path"))
            _reject_escaping_include(target, root, display, number,
                                     include_stack)
            here = SourceRef(display, number)
            lines.extend(_read_physical_lines(
                target, include_stack + (here,), seen + (resolved,), root))
            continue
        lines.append(PhysicalLine(
            text=raw, raw=raw, file=display, line=number, end_line=number,
            include_stack=include_stack))
    return lines


def _reject_escaping_include(target: Path,
                             root: Path,
                             display: str,
                             number: int,
                             include_stack: tuple[SourceRef, ...]) -> None:
    """Refuse an ``#include`` that resolves outside the import root.

    ``.par`` files are third-party input.  Allowing ``../../..`` traversal
    would let an input file pull arbitrary host content into a generated deck,
    so the reader confines resolution to the declared root.
    """
    try:
        resolved_root = root.resolve()
        candidate = target.resolve() if target.exists() else target
        candidate.relative_to(resolved_root)
    except ValueError:
        raise ParParseError(
            f"{display}:{number}: #include escapes the import root "
            f"{resolved_root}: {target}") from None


def _display_path(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.name


def _stack_prefix(include_stack: Sequence[SourceRef]) -> str:
    if not include_stack:
        return ""
    return " <- ".join(str(item) for item in include_stack) + ": "


def _join_continuations(lines: Sequence[PhysicalLine]) -> list[PhysicalLine]:
    """Merge statements that span several physical lines.

    Only bracketed constructs continue across a newline in this grammar:
    ``NumericalTable ( ... )`` tables and the multi-group ``eValley"..."(..)(..)``
    records.  Joining is driven by bracket depth so that a table row containing
    a stray token cannot silently terminate the statement early.
    """
    joined: list[PhysicalLine] = []
    index = 0
    while index < len(lines):
        current = lines[index]
        code, _unit = _strip_comments(current.text)
        depth = _bracket_depth(code)
        if depth <= 0:
            joined.append(current)
            index += 1
            continue

        parts = [current.text]
        end_line = current.line
        cursor = index + 1
        while cursor < len(lines) and depth > 0:
            nxt = lines[cursor]
            parts.append(nxt.text)
            end_line = nxt.line
            nxt_code, _ = _strip_comments(nxt.text)
            depth += _bracket_depth(nxt_code)
            cursor += 1
        if depth > 0:
            raise ParParseError(
                f"{current.file}:{current.line}: unterminated '(' in "
                f"{current.raw.strip()!r}")
        joined.append(PhysicalLine(
            text="\n".join(parts), raw="\n".join(parts), file=current.file,
            line=current.line, end_line=end_line,
            include_stack=current.include_stack))
        index = cursor
    return joined


def _bracket_depth(code: str) -> int:
    depth = 0
    in_string = False
    for char in code:
        if char == '"':
            in_string = not in_string
        elif in_string:
            continue
        elif char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
    return depth


def _parse_number(token: str) -> float | None:
    return float(token) if _NUMBER_RE.match(token) else None


def _decompose_name(raw_name: str) -> dict[str, Any]:
    """Split a parameter name into base, variant, indices, and suffix.

    Recognised shapes, all present in the shipped corpus::

        epsilon              plain
        dEg0(OldSlotboom)    named model variant
        efit( 0)             numeric index in parentheses
        Dop[1][1]            bracketed indices
        (tau_w)_ele          wrapped base with a carrier suffix
        1/kappa              reciprocal-valued name

    An unrecognised shape is reported so the caller can fail closed rather than
    guessing at a decomposition.
    """
    name = raw_name.strip()
    match = _PAREN_VARIANT_RE.match(name)
    if match:
        variant = match.group("variant")
        index = _parse_number(variant)
        if index is not None and float(index).is_integer():
            return {"base_name": match.group("base"), "variant": None,
                    "indices": [int(index)], "recognised": True}
        return {"base_name": match.group("base"), "variant": variant,
                "indices": [], "recognised": True}

    match = _INDEXED_RE.match(name)
    if match:
        indices = [int(value) for value in
                   re.findall(r"-?\d+", match.group("indices"))]
        return {"base_name": match.group("base"), "variant": None,
                "indices": indices, "recognised": True}

    match = _WRAPPED_RE.match(name)
    if match:
        return {"base_name": match.group("base") + match.group("suffix"),
                "variant": None, "indices": [], "recognised": True}

    match = _RECIPROCAL_RE.match(name)
    if match:
        return {"base_name": name, "variant": None, "indices": [],
                "recognised": True}

    if _PLAIN_NAME_RE.match(name):
        return {"base_name": name, "variant": None, "indices": [],
                "recognised": True}

    return {"base_name": name, "variant": None, "indices": [],
            "recognised": False}


def _classify_value(lexeme: str) -> dict[str, Any]:
    """Classify a right-hand side without converting units.

    ``value_kind`` distinguishes the shapes that a downstream emitter must
    treat differently:

    ``scalar``           one number
    ``carrier_pair``     ``electron , hole`` -- the dominant ``.par`` shape
    ``number_list``      three or more numbers, comma or whitespace separated
    ``string``           a bare or quoted identifier such as ``PositiveSpline``
    ``string_list``      comma-separated identifiers
    ``numerical_table``  ``NumericalTable ( r ; r ; ... )``
    ``mixed_list``       anything else that still tokenises cleanly

    A pair of numbers is reported as ``carrier_pair`` because in this grammar
    two comma-separated values always mean electron/hole.  Whitespace-separated
    numbers (``efit( 0) = 1.20698  2.63089``) are a list, never a pair, since
    the file uses commas for the carrier split.
    """
    text = lexeme.strip()
    if not text:
        return {"value_kind": "empty", "values": []}

    # ``eoffset = ()`` is how the corpus spells an explicitly unset value.  It
    # is meaningful (the parameter exists but carries no number) and so is kept
    # rather than treated as a parse failure.
    if re.match(r"^\(\s*\)$", text):
        return {"value_kind": "empty_list", "values": []}

    table = re.match(r"^NumericalTable\s*\((?P<body>.*)\)$", text,
                     re.DOTALL | re.IGNORECASE)
    if table:
        rows = []
        for row in table.group("body").split(";"):
            tokens = row.split()
            if not tokens:
                continue
            numbers = [_parse_number(token) for token in tokens]
            if any(number is None for number in numbers):
                return {"value_kind": "unrecognised", "values": []}
            rows.append(numbers)
        return {"value_kind": "numerical_table", "values": rows}

    if "," in text:
        parts = [part.strip() for part in text.split(",")]
        parts = [part for part in parts if part != ""]
        numbers = [_parse_number(part) for part in parts]
        if all(number is not None for number in numbers):
            kind = "carrier_pair" if len(numbers) == 2 else (
                "scalar" if len(numbers) == 1 else "number_list")
            return {"value_kind": kind, "values": numbers}
        if all(_is_bare_string(part) for part in parts):
            return {"value_kind": "string_list",
                    "values": [_unquote(part) for part in parts]}
        return {"value_kind": "mixed_list", "values": [_unquote(p) for p in parts]}

    tokens = text.split()
    numbers = [_parse_number(token) for token in tokens]
    if all(number is not None for number in numbers):
        if len(numbers) == 1:
            return {"value_kind": "scalar", "values": numbers}
        return {"value_kind": "number_list", "values": numbers}

    if len(tokens) == 1 and _is_bare_string(tokens[0]):
        return {"value_kind": "string", "values": [_unquote(tokens[0])]}

    if all(_is_bare_string(token) or _parse_number(token) is not None
           for token in tokens):
        return {"value_kind": "mixed_list",
                "values": [_unquote(token) for token in tokens]}

    return {"value_kind": "unrecognised", "values": []}


def _is_bare_string(token: str) -> bool:
    if token.startswith('"') and token.endswith('"') and len(token) >= 2:
        return True
    return bool(re.match(r"^[A-Za-z_][A-Za-z0-9_.+-]*$", token))


def _unquote(token: str) -> str:
    if token.startswith('"') and token.endswith('"') and len(token) >= 2:
        return token[1:-1]
    return token


@dataclass
class _Block:
    index: int
    scope_kind: str
    scope_name: str | None
    section: str
    variant: str | None
    file: str
    line: int
    include_stack: tuple[SourceRef, ...]
    end_line: int = 0
    parameters: list[dict[str, Any]] = field(default_factory=list)

    def key(self) -> tuple[str, str | None, str, str | None]:
        return (self.scope_kind, self.scope_name, self.section, self.variant)

    def to_json(self) -> dict[str, Any]:
        scope: dict[str, Any] = {"kind": self.scope_kind}
        if self.scope_name is not None:
            scope["name"] = self.scope_name
        payload: dict[str, Any] = {
            "index": self.index,
            "scope": scope,
            "section": self.section,
            "order": self.index,
            "source_span": {
                "file": self.file,
                "line": self.line,
                "end_line": self.end_line or self.line,
            },
            "include_stack": [item.to_json() for item in self.include_stack],
            "parameters": self.parameters,
        }
        if self.variant is not None:
            payload["variant"] = self.variant
        return payload


class _Reader:
    """Walks joined physical lines and builds ordered blocks."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.blocks: list[_Block] = []
        self.errors: list[dict[str, Any]] = []
        self._scope_stack: list[tuple[str, str | None]] = []
        self._block: _Block | None = None
        self._pending: tuple[str, str | None] | tuple[str, str, str | None] | None = None
        self._pending_kind: str | None = None
        self._pending_line: PhysicalLine | None = None

    # -- diagnostics -----------------------------------------------------
    def _error(self, line: PhysicalLine, reason: str) -> None:
        self.errors.append({
            "file": line.file,
            "line": line.line,
            "text": line.raw.strip()[:200],
            "reason": reason,
            "include_stack": [item.to_json() for item in line.include_stack],
        })

    # -- scope helpers ---------------------------------------------------
    def _current_scope(self) -> tuple[str, str | None]:
        return self._scope_stack[-1] if self._scope_stack else ("global", None)

    def _open_block(self, line: PhysicalLine, section: str,
                    variant: str | None) -> None:
        scope_kind, scope_name = self._current_scope()
        self._block = _Block(
            index=len(self.blocks),
            scope_kind=scope_kind,
            scope_name=scope_name,
            section=section,
            variant=variant,
            file=line.file,
            line=line.line,
            include_stack=line.include_stack,
        )

    def _close_block(self, line: PhysicalLine) -> None:
        assert self._block is not None
        self._block.end_line = line.end_line
        self.blocks.append(self._block)
        self._block = None

    # -- main loop -------------------------------------------------------
    def read(self, lines: Sequence[PhysicalLine]) -> None:
        for line in lines:
            code, unit = _strip_comments(line.text)
            stripped = code.strip()
            if not stripped:
                continue
            self._read_statement(line, stripped, unit)

        if self._block is not None:
            self._error(
                PhysicalLine(self._block.section, self._block.section,
                             self._block.file, self._block.line,
                             self._block.line, self._block.include_stack),
                "section block is never closed by '}'")
        if self._scope_stack:
            kind, name = self._scope_stack[-1]
            self.errors.append({
                "file": self.blocks[-1].file if self.blocks else "<input>",
                "line": 0,
                "text": f'{kind} = "{name}"',
                "reason": "scope block is never closed by '}'",
                "include_stack": [],
            })

    def _read_statement(self, line: PhysicalLine, stripped: str,
                        unit: str | None) -> None:
        # A brace that completes a header seen on the previous line.
        if self._pending_kind is not None:
            if stripped.startswith("{"):
                self._commit_pending(line)
                remainder = stripped[1:].strip()
                if remainder:
                    self._read_statement(line, remainder, unit)
                return
            self._error(line,
                        "section or scope header is not followed by '{'")
            self._pending_kind = None
            self._pending_line = None
            return

        if stripped == "}":
            self._close_scope_or_block(line)
            return

        if stripped.startswith("{"):
            remainder = stripped[1:].strip()
            if remainder:
                self._read_statement(line, remainder, unit)
            return

        if stripped.endswith("}") and stripped != "}":
            self._read_statement(line, stripped[:-1].strip(), unit)
            self._read_statement(line, "}", None)
            return

        scope = _SCOPE_HEADER_RE.match(stripped)
        if scope and self._block is None:
            kind = _SCOPE_KINDS[scope.group("kind")]
            self._pending_kind = "scope"
            self._pending = (kind, scope.group("name"))
            self._pending_line = line
            if scope.group("brace"):
                self._commit_pending(line)
            return

        # Records must be tested before assignments: their group bodies carry
        # inner ``=`` signs (``hValley"LH"(m=0.16 ...)``) that would otherwise
        # be mistaken for a top-level assignment.
        record = _RECORD_RE.match(" ".join(stripped.split()))
        if record and self._block is not None:
            self._append_record(line, record, stripped)
            return

        # ``NumericalTable ( ... )`` also appears as a bare statement, without
        # a name and without ``=``, when the enclosing section has exactly one
        # table.  It is data, so it is kept as a parameter rather than skipped.
        if _BARE_TABLE_RE.match(" ".join(stripped.split())) and self._block is not None:
            self._read_assignment(line, f"NumericalTable = {stripped}", unit)
            return

        if "=" in stripped:
            self._read_assignment(line, stripped, unit)
            return

        header = _SECTION_HEADER_RE.match(stripped)
        if header and self._block is None:
            self._pending_kind = "section"
            self._pending = ("section", header.group("name"),
                             header.group("variant"))
            self._pending_line = line
            if header.group("brace"):
                self._commit_pending(line)
            return

        self._error(line, "line is valid .par syntax but was not classified")

    def _commit_pending(self, line: PhysicalLine) -> None:
        assert self._pending is not None
        if self._pending_kind == "scope":
            kind, name = self._pending  # type: ignore[misc]
            self._scope_stack.append((kind, name))
        else:
            _tag, section, variant = self._pending  # type: ignore[misc]
            self._open_block(self._pending_line or line, section, variant)
        self._pending_kind = None
        self._pending = None
        self._pending_line = None

    def _close_scope_or_block(self, line: PhysicalLine) -> None:
        if self._block is not None:
            self._close_block(line)
        elif self._scope_stack:
            self._scope_stack.pop()
        else:
            self._error(line, "'}' does not close any open block or scope")

    def _read_assignment(self, line: PhysicalLine, stripped: str,
                         unit: str | None) -> None:
        raw_name, _, lexeme = stripped.partition("=")
        raw_name = raw_name.strip()
        lexeme = lexeme.strip()
        if self._block is None:
            self._error(line, "assignment appears outside any section block")
            return

        decomposed = _decompose_name(raw_name)
        if not decomposed["recognised"]:
            self._error(line, f"unrecognised parameter name {raw_name!r}")
            return

        classified = _classify_value(lexeme)
        if classified["value_kind"] in ("unrecognised", "empty"):
            self._error(line, f"unrecognised value for {raw_name!r}")
            return

        entry: dict[str, Any] = {
            "raw_name": raw_name,
            "base_name": decomposed["base_name"],
            "raw_lexeme": " ".join(lexeme.split()),
            "value_kind": classified["value_kind"],
            "values": classified["values"],
            "source": {"file": line.file, "line": line.line},
        }
        if decomposed["variant"] is not None:
            entry["variant"] = decomposed["variant"]
        if decomposed["indices"]:
            entry["indices"] = decomposed["indices"]
        if unit is not None:
            entry["raw_unit"] = unit
        self._block.parameters.append(entry)

    def _append_record(self, line: PhysicalLine, record: "re.Match[str]",
                       stripped: str) -> None:
        assert self._block is not None
        self._block.parameters.append({
            "raw_name": record.group("base"),
            "base_name": record.group("base"),
            "raw_lexeme": " ".join(stripped.split()),
            "value_kind": "record",
            "values": [],
            "variant": record.group("label"),
            "source": {"file": line.file, "line": line.line},
        })


_SCOPE_KINDS = {
    "Material": "material",
    "Region": "region",
    "MaterialInterface": "material_interface",
    "RegionInterface": "region_interface",
}

#: Sections whose parameters describe a model Vela does not implement at all.
#: These still parse and are still kept in the IR -- an emitter needs to see
#: them to report coverage honestly -- but they are flagged so that activating
#: one is a hard error rather than a silent omission.
UNSUPPORTED_SECTIONS: dict[str, str] = {
    "Aniso": "anisotropic transport is not implemented",
    "EnergyRelaxationTime": "carrier energy balance is not solved",
    "HydroHighFieldDependence": "hydrodynamic transport is not implemented",
    "LatticeHeatCapacity": "the lattice heat equation is not solved",
    "LatticeParameters": "lattice/band-structure detail is not modelled",
    "OkutoCrowell": "Okuto-Crowell ionization coefficients are not implemented",
    "Lackner": "Lackner ionization coefficients are not implemented",
    "UniBo": "University of Bologna ionization is not implemented",
    "UniBo2": "University of Bologna ionization is not implemented",
    "UniBoDopingDependence": "UniBo mobility is not implemented",
    "UniBoEnormalDependence": "UniBo mobility is not implemented",
    "NBTI": "bias-temperature instability kinetics are not solved",
    "eNMP": "non-radiative multi-phonon trap kinetics are not solved",
    "HCSDegradation": "hot-carrier stress degradation is not solved",
    "SHEDistribution": "spherical-harmonics BTE is not solved",
    "ThermalConductivity": "the lattice heat equation is not solved",
    "RefractiveIndex": "optics are not modelled",
    "ComplexRefractiveIndex": "optics are not modelled",
}


def parse_parameter_ir(path: Path,
                       root: Path | None = None) -> dict[str, Any]:
    """Parse ``path`` into ``vela.sentaurus_parameter_ir.v1``.

    ``root`` bounds ``#include`` resolution and defaults to the directory
    holding ``path``.  Raises :class:`ParParseError` if any statement cannot be
    classified; in that case no IR is returned, because a partially understood
    parameter file is more dangerous than none at all.
    """
    path = Path(path)
    root = Path(root) if root is not None else path.parent

    lines = _read_physical_lines(path, (), (), root)
    joined = _join_continuations(lines)
    reader = _Reader(root)
    reader.read(joined)

    if reader.errors:
        first = reader.errors[0]
        raise ParParseError(
            f"{first['file']}:{first['line']}: {first['reason']}: "
            f"{first['text']!r} "
            f"({len(reader.errors)} unclassified statement(s) total)")

    blocks = [block.to_json() for block in reader.blocks]
    _annotate_shadowing(reader.blocks, blocks)

    files: list[str] = [_display_path(path, root)]
    for line in joined:
        if line.file not in files:
            files.append(line.file)

    unsupported = [
        {
            "index": block["index"],
            "section": block["section"],
            "scope": block["scope"],
            "reason": UNSUPPORTED_SECTIONS[block["section"]],
        }
        for block in blocks
        if block["section"] in UNSUPPORTED_SECTIONS
    ]

    return {
        "schema": PARAMETER_IR_SCHEMA,
        "source": _display_path(path, root),
        "files": files,
        "blocks": blocks,
        "unsupported_sections": unsupported,
        "parse_errors": [],
    }


def _annotate_shadowing(blocks: Sequence[_Block],
                        payload: list[dict[str, Any]]) -> None:
    """Mark earlier definitions that a later one overrides.

    Sentaurus applies last-definition-wins for a repeated
    ``(scope, section, variant)``.  Both definitions stay in the IR so the
    override is auditable; the earlier one simply records which block replaced
    it, and each shadowed parameter records the same at parameter granularity.
    """
    last_by_key: dict[tuple[Any, ...], int] = {}
    for block in blocks:
        key = block.key()
        previous = last_by_key.get(key)
        if previous is not None:
            payload[previous]["shadowed_by"] = block.index
            _shadow_parameters(payload[previous], payload[block.index])
        last_by_key[key] = block.index


def _shadow_parameters(earlier: dict[str, Any],
                       later: dict[str, Any]) -> None:
    later_names = {
        (item["raw_name"], tuple(item.get("indices", [])))
        for item in later["parameters"]
    }
    for item in earlier["parameters"]:
        if (item["raw_name"], tuple(item.get("indices", []))) in later_names:
            item["shadowed_by"] = later["index"]


def active_blocks(ir: dict[str, Any]) -> list[dict[str, Any]]:
    """Return only the blocks that survive last-definition-wins resolution."""
    return [block for block in ir["blocks"] if "shadowed_by" not in block]


def section_index(ir: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Group active blocks by section name, preserving source order."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for block in active_blocks(ir):
        grouped.setdefault(block["section"], []).append(block)
    return grouped


def parameter_count(ir: dict[str, Any]) -> int:
    return sum(len(block["parameters"]) for block in ir["blocks"])


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Parse a Sentaurus .par file into the Vela parameter IR.")
    parser.add_argument("par", type=Path, help="input .par file")
    parser.add_argument("--root", type=Path, default=None,
                        help="directory that bounds #include resolution "
                             "(defaults to the input file's directory)")
    parser.add_argument("--out", type=Path, default=None,
                        help="write the IR JSON here instead of stdout")
    parser.add_argument("--summary", action="store_true",
                        help="print a block/parameter summary instead of IR")
    args = parser.parse_args(argv)

    try:
        ir = parse_parameter_ir(args.par, args.root)
    except ParParseError as error:
        print(f"error: {error}")
        return 1

    if args.summary:
        print(f"source           : {ir['source']}")
        print(f"files            : {len(ir['files'])}")
        print(f"blocks           : {len(ir['blocks'])}")
        print(f"active blocks    : {len(active_blocks(ir))}")
        print(f"parameters       : {parameter_count(ir)}")
        print(f"unsupported      : {len(ir['unsupported_sections'])}")
        return 0

    if args.out:
        write_json(args.out, ir)
    else:
        print(json.dumps(ir, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
