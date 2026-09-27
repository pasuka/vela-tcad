"""Token-stream utilities shared by the section parsers."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple, Union

from .catalogue import finding
from .lexer import Source, Token

_NUMBER_RE = re.compile(r"^[+-]?(?:\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)$")


class SyntaxIssue(Exception):
    """A structural problem that becomes an E104 finding."""

    def __init__(self, detail: str, line: Optional[int], text: str) -> None:
        super().__init__(detail)
        self.detail = detail
        self.line = line
        self.text = text


class ParseContext:
    """Collects findings while the section parsers run."""

    def __init__(self, source: Source) -> None:
        self.source = source
        self.errors: List[dict] = []
        self.warnings: List[dict] = []

    def error(self, message_id: str, line: Optional[int], section: Optional[str], source_text: str,
              **fields: object) -> dict:
        item = finding(message_id, line, section, source_text, **fields)
        self.errors.append(item)
        return item

    def warn(self, message_id: str, line: Optional[int], section: Optional[str], source_text: str,
             **fields: object) -> dict:
        item = finding(message_id, line, section, source_text, **fields)
        self.warnings.append(item)
        return item

    def syntax_error(self, issue: SyntaxIssue, section: Optional[str]) -> None:
        self.error("E104", issue.line, section, issue.text, detail=issue.detail)


class TokenStream:
    def __init__(self, tokens: Sequence[Token], source: Source) -> None:
        self.tokens = list(tokens)
        self.source = source
        self.i = 0

    def at_end(self) -> bool:
        return self.i >= len(self.tokens)

    def peek(self, offset: int = 0) -> Optional[Token]:
        index = self.i + offset
        if 0 <= index < len(self.tokens):
            return self.tokens[index]
        return None

    def next(self) -> Token:
        if self.at_end():
            last = self.tokens[-1] if self.tokens else None
            raise SyntaxIssue("unexpected end of block", last.line if last else None,
                              last.value if last else "")
        token = self.tokens[self.i]
        self.i += 1
        return token

    def at_punct(self, char: str, offset: int = 0) -> bool:
        token = self.peek(offset)
        return token is not None and token.kind == "punct" and token.value == char

    def expect_punct(self, char: str) -> Token:
        token = self.peek()
        if token is None or token.kind != "punct" or token.value != char:
            found = token.value if token is not None else "end of block"
            line = token.line if token is not None else (self.tokens[-1].line if self.tokens else None)
            raise SyntaxIssue(f"expected '{char}' but found '{found}'", line, found)
        return self.next()

    def take_group(self) -> Tuple[Token, List[Token], Token]:
        """Consume a balanced ``(...)`` or ``{...}`` group starting at the current token."""
        opener = self.next()
        if opener.kind != "punct" or opener.value not in "({":
            raise SyntaxIssue(f"expected '(' or '{{' but found '{opener.value}'", opener.line, opener.value)
        stack = [opener.value]
        inner: List[Token] = []
        while True:
            token = self.peek()
            if token is None:
                raise SyntaxIssue(f"unbalanced '{opener.value}'", opener.line, opener.value)
            self.next()
            if token.kind == "punct" and token.value in "({":
                stack.append(token.value)
            elif token.kind == "punct" and token.value in ")}":
                expected = ")" if stack[-1] == "(" else "}"
                if token.value != expected:
                    raise SyntaxIssue(f"mismatched '{token.value}'", token.line, token.value)
                stack.pop()
                if not stack:
                    return opener, inner, token
            inner.append(token)

    def take_value(self) -> Token:
        """Consume a word or string token used as a value."""
        token = self.next()
        if token.kind not in ("word", "string"):
            raise SyntaxIssue(f"expected a value but found '{token.value}'", token.line, token.value)
        return token

    def text_between(self, first: Token, last: Token) -> str:
        return self.source.slice(first.start, last.end)


@dataclass
class Block:
    name: Token
    scope: Optional[List[Token]]
    body: List[Token]
    open_brace: Token
    close_brace: Token
    scope_tokens: List[Token] = field(default_factory=list)  # including the parentheses


def group_blocks(tokens: Sequence[Token], source: Source) -> Tuple[List[Block], Optional[SyntaxIssue]]:
    """Split the top level into ``Name [ (scope) ] { body }`` blocks.

    Returns the blocks parsed so far and the first structural issue (if any); the
    remainder of the deck after the issue is not parsed.
    """
    stream = TokenStream(tokens, source)
    blocks: List[Block] = []
    while not stream.at_end():
        token = stream.next()
        if token.kind != "word":
            return blocks, SyntaxIssue(f"unexpected '{token.value}' at top level", token.line, token.value)
        scope = None
        scope_tokens: List[Token] = []
        try:
            if stream.at_punct("("):
                opener, scope, closer = stream.take_group()
                scope_tokens = [opener, *scope, closer]
            if not stream.at_punct("{"):
                found = stream.peek()
                found_text = found.value if found is not None else "end of file"
                return blocks, SyntaxIssue(
                    f"expected '{{' after section '{token.value}' but found '{found_text}'",
                    token.line, token.value)
            opener, body, closer = stream.take_group()
        except SyntaxIssue as issue:
            return blocks, issue
        blocks.append(Block(token, scope, body, opener, closer, scope_tokens))
    return blocks, None


def is_number(text: str) -> bool:
    return bool(_NUMBER_RE.match(text))


def parse_number(token: Token, what: str = "number") -> float:
    if token.kind != "word" or not is_number(token.value):
        raise SyntaxIssue(f"invalid {what} '{token.value}'", token.line, token.value)
    return float(token.value)


def parse_int(token: Token, what: str = "integer") -> int:
    value = parse_number(token, what)
    if value != int(value):
        raise SyntaxIssue(f"invalid {what} '{token.value}'", token.line, token.value)
    return int(value)


def coerce_value(token: Token) -> Union[float, int, str]:
    """Number when the token looks numeric, otherwise the raw text."""
    if token.kind == "word" and is_number(token.value):
        value = float(token.value)
        if value == int(value) and "." not in token.value and "e" not in token.value.lower():
            return int(value)
        return value
    return token.value


def assignment_text(key: Token, value: Optional[Token]) -> str:
    """``Key=value`` with the value token as written (quotes kept for strings)."""
    if value is None:
        return key.value
    literal = f'"{value.value}"' if value.kind == "string" else value.value
    return f"{key.value}={literal}"
