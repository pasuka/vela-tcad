#!/usr/bin/env python3
"""S-expression tokenizer, reader, and evaluator for Sentaurus SDE ``.cmd`` files.

SDE command files are Scheme sources. Parsing them with regular expressions
loses data silently whenever formatting deviates from the assumed shape, so
this module provides a real reader instead.

Every parsed form keeps its source file name, line, and column so that
diagnostics can point at the exact command that failed, and so that a
fail-closed frontend can report unsupported commands precisely.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Sequence


NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


class SexpError(ValueError):
    """Raised for any tokenizer/reader/evaluator failure.

    The message always identifies the source location so that callers can
    surface a fail-closed diagnostic without inventing their own context.
    """


@dataclass(frozen=True)
class Atom:
    """A non-list token: symbol, string literal, or number."""

    value: Any
    kind: str  # "symbol" | "string" | "number"
    line: int
    column: int
    file: str

    @property
    def is_symbol(self) -> bool:
        return self.kind == "symbol"

    @property
    def is_string(self) -> bool:
        return self.kind == "string"

    @property
    def is_number(self) -> bool:
        return self.kind == "number"

    def location(self) -> str:
        return f"{self.file}:{self.line}:{self.column}"


@dataclass(frozen=True)
class Sexp:
    """A parenthesised list form together with its source location."""

    items: tuple[Any, ...]
    line: int
    column: int
    end_line: int
    file: str
    text: str = ""

    def __len__(self) -> int:
        return len(self.items)

    def __iter__(self) -> Iterator[Any]:
        return iter(self.items)

    def __getitem__(self, index: int) -> Any:
        return self.items[index]

    def location(self) -> str:
        return f"{self.file}:{self.line}:{self.column}"

    def head(self) -> str | None:
        """Return the leading symbol name, or ``None`` when there is none."""
        if not self.items:
            return None
        first = self.items[0]
        if isinstance(first, Atom) and first.is_symbol:
            return str(first.value)
        return None

    def command_text(self, max_chars: int = 200) -> str:
        """Return a single-line excerpt of the original source text."""
        collapsed = " ".join(self.text.split())
        if len(collapsed) <= max_chars:
            return collapsed
        return collapsed[: max_chars - 3] + "..."


@dataclass
class _Token:
    kind: str  # "(" | ")" | "atom"
    text: str
    line: int
    column: int
    offset: int


def _tokenize(text: str, filename: str) -> list[_Token]:
    tokens: list[_Token] = []
    line = 1
    column = 1
    index = 0
    length = len(text)

    while index < length:
        char = text[index]

        if char == "\n":
            line += 1
            column = 1
            index += 1
            continue
        if char.isspace():
            index += 1
            column += 1
            continue
        if char == ";":
            # Scheme line comment: skip to end of line, but leave the newline
            # so that the line counter above stays correct.
            while index < length and text[index] != "\n":
                index += 1
            continue
        if char == "#" and text.startswith("#|", index):
            end = text.find("|#", index + 2)
            if end < 0:
                raise SexpError(
                    f"{filename}:{line}:{column}: unterminated block comment '#|'")
            consumed = text[index:end + 2]
            line += consumed.count("\n")
            trailing = consumed.rsplit("\n", 1)
            column = len(trailing[-1]) + 1 if len(trailing) > 1 else column + len(consumed)
            index = end + 2
            continue
        if char in "()":
            tokens.append(_Token(char, char, line, column, index))
            index += 1
            column += 1
            continue
        if char == '"':
            start_line, start_column, start_index = line, column, index
            index += 1
            column += 1
            buffer: list[str] = []
            while True:
                if index >= length:
                    raise SexpError(
                        f"{filename}:{start_line}:{start_column}: unterminated string literal")
                current = text[index]
                if current == "\\" and index + 1 < length:
                    buffer.append(text[index + 1])
                    index += 2
                    column += 2
                    continue
                if current == '"':
                    index += 1
                    column += 1
                    break
                if current == "\n":
                    line += 1
                    column = 1
                else:
                    column += 1
                buffer.append(current)
                index += 1
            tokens.append(
                _Token("string", "".join(buffer), start_line, start_column, start_index))
            continue

        start_line, start_column, start_index = line, column, index
        while index < length and not text[index].isspace() and text[index] not in "();\"":
            index += 1
            column += 1
        tokens.append(
            _Token("atom", text[start_index:index], start_line, start_column, start_index))

    return tokens


def _atom_from_token(token: _Token, filename: str) -> Atom:
    if token.kind == "string":
        return Atom(token.text, "string", token.line, token.column, filename)
    if NUMBER_RE.fullmatch(token.text):
        return Atom(float(token.text), "number", token.line, token.column, filename)
    return Atom(token.text, "symbol", token.line, token.column, filename)


def parse_forms(text: str, filename: str) -> list[Sexp]:
    """Parse ``text`` into a list of top-level forms.

    Bare top-level atoms are rejected: every SDE statement is a command list.
    """
    tokens = _tokenize(text, filename)
    forms: list[Sexp] = []
    stack: list[tuple[list[Any], int, int, int]] = []

    for token in tokens:
        if token.kind == "(":
            stack.append(([], token.line, token.column, token.offset))
            continue
        if token.kind == ")":
            if not stack:
                raise SexpError(
                    f"{filename}:{token.line}:{token.column}: unbalanced ')'")
            items, line, column, offset = stack.pop()
            form = Sexp(
                items=tuple(items),
                line=line,
                column=column,
                end_line=token.line,
                file=filename,
                text=text[offset:token.offset + 1],
            )
            if stack:
                stack[-1][0].append(form)
            else:
                forms.append(form)
            continue

        atom = _atom_from_token(token, filename)
        if not stack:
            raise SexpError(
                f"{filename}:{atom.line}:{atom.column}: unexpected top-level atom "
                f"'{token.text}'; expected a '(' command form")
        stack[-1][0].append(atom)

    if stack:
        _, line, column, _ = stack[-1]
        raise SexpError(f"{filename}:{line}:{column}: unterminated '(' form")
    return forms


def parse_file(path: Path) -> list[Sexp]:
    """Read and parse an SDE command file."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    return parse_forms(text, str(path))


# ---------------------------------------------------------------------------
# Expression evaluation
# ---------------------------------------------------------------------------

_BINARY_OPS = {
    "+": lambda a, b: a + b,
    "-": lambda a, b: a - b,
    "*": lambda a, b: a * b,
    "/": lambda a, b: a / b,
}

_UNARY_FUNCS = {
    "abs": abs,
    "sqrt": math.sqrt,
    "exp": math.exp,
    "log": math.log,
    "min": min,
    "max": max,
}


@dataclass
class Environment:
    """Symbol table produced by ``(define name value)`` forms."""

    values: dict[str, float | str] = field(default_factory=dict)

    def define(self, name: str, value: float | str) -> None:
        self.values[name] = value

    def __contains__(self, name: object) -> bool:
        return name in self.values

    def get(self, name: str) -> float | str:
        return self.values[name]


def _describe(node: Any) -> str:
    if isinstance(node, Atom):
        return str(node.value)
    if isinstance(node, Sexp):
        return node.command_text()
    return str(node)


def eval_number(node: Any, env: Environment) -> float:
    """Evaluate ``node`` to a float.

    Unknown symbols and unsupported operators raise ``SexpError`` instead of
    falling back to a default value, so that a partially understood SDE file
    can never be mistaken for a fully understood one.
    """
    if isinstance(node, Atom):
        if node.is_number:
            return float(node.value)
        if node.is_symbol:
            name = str(node.value)
            if name not in env:
                raise SexpError(
                    f"{node.location()}: undefined SDE symbol '{name}'")
            value = env.get(name)
            if isinstance(value, (int, float)):
                return float(value)
            raise SexpError(
                f"{node.location()}: SDE symbol '{name}' is not numeric "
                f"(value {value!r})")
        raise SexpError(
            f"{node.location()}: expected a number, got string \"{node.value}\"")

    if isinstance(node, Sexp):
        head = node.head()
        if head is None:
            raise SexpError(
                f"{node.location()}: expected a numeric expression, got "
                f"'{node.command_text()}'")
        args = node.items[1:]
        if head in _BINARY_OPS:
            if not args:
                raise SexpError(
                    f"{node.location()}: operator '{head}' requires at least one operand")
            values = [eval_number(arg, env) for arg in args]
            if len(values) == 1:
                if head == "-":
                    return -values[0]
                if head == "/":
                    return 1.0 / values[0]
                return values[0]
            result = values[0]
            for value in values[1:]:
                result = _BINARY_OPS[head](result, value)
            return result
        if head in _UNARY_FUNCS:
            values = [eval_number(arg, env) for arg in args]
            if head in ("min", "max"):
                if not values:
                    raise SexpError(
                        f"{node.location()}: '{head}' requires at least one operand")
                return float(_UNARY_FUNCS[head](*values))
            if len(values) != 1:
                raise SexpError(
                    f"{node.location()}: '{head}' expects exactly one operand")
            return float(_UNARY_FUNCS[head](values[0]))
        raise SexpError(
            f"{node.location()}: unsupported SDE expression operator '{head}' in "
            f"'{node.command_text()}'")

    raise SexpError(f"unsupported SDE expression node: {_describe(node)}")


def eval_scalar(node: Any, env: Environment) -> float | str:
    """Evaluate ``node`` to a number when possible, otherwise a string."""
    if isinstance(node, Atom) and node.is_string:
        return str(node.value)
    return eval_number(node, env)


def expect_string(node: Any, env: Environment, what: str) -> str:
    """Return a string literal, or a symbol bound to a string by ``define``."""
    if isinstance(node, Atom):
        if node.is_string:
            return str(node.value)
        if node.is_symbol:
            name = str(node.value)
            if name in env:
                value = env.get(name)
                if isinstance(value, str):
                    return value
            raise SexpError(
                f"{node.location()}: expected a quoted {what}, got symbol '{name}'")
    location = node.location() if isinstance(node, (Atom, Sexp)) else "<unknown>"
    raise SexpError(f"{location}: expected a quoted {what}, got '{_describe(node)}'")


def parse_position(node: Any, env: Environment) -> tuple[float, float, float]:
    """Evaluate a ``(position x y z)`` form to a 3-tuple."""
    if not isinstance(node, Sexp) or node.head() != "position":
        location = node.location() if isinstance(node, (Atom, Sexp)) else "<unknown>"
        raise SexpError(
            f"{location}: expected a (position x y z) form, got '{_describe(node)}'")
    args = node.items[1:]
    if len(args) != 3:
        raise SexpError(
            f"{node.location()}: (position ...) requires exactly 3 coordinates, "
            f"got {len(args)}")
    return tuple(eval_number(arg, env) for arg in args)  # type: ignore[return-value]


def collect_positions(nodes: Sequence[Any], env: Environment) -> list[tuple[float, float, float]]:
    return [parse_position(node, env) for node in nodes]
