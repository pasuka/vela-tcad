"""Source normalisation, Workbench preprocessor scan (gate 1), comment stripping and tokenizer.

Pipeline (docs/sentaurus_sdevice_workflow.md section 5, items 1-2):

1. bytes -> text: strip a UTF-8 BOM, decode UTF-8, CRLF -> LF.  ``source.sha256``
   is taken over this normalised text so BOM/CRLF variants yield identical IR.
2. gate 1 scans every line *before* comment stripping: a ``#`` that starts the
   line and is immediately followed by a Workbench directive word is E101; Tcl
   ``!( ... )!`` blocks are E103; ``@...@`` placeholders are E102.  Placeholders
   are only looked for in code (and on directive lines), so ``#Plot(... @node@ ...)``
   comments are not reported while ``#setdep @previous@`` reports both.
3. comment stripping blanks ``*``/``#`` comments (token-start rule, quoted strings
   protected) and Tcl blocks with spaces, keeping every offset and line number.
4. the tokenizer produces words, strings and the punctuation ``{ } ( ) = ;``.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from typing import List, Optional, Tuple

from .catalogue import finding

DIRECTIVES = ("if", "elif", "else", "endif", "setdep", "define", "include", "noprint",
              "header", "endheader", "split", "endsplit", "exit")
_DIRECTIVE_RE = re.compile(r"#(?:" + "|".join(DIRECTIVES) + r")\b")
PLACEHOLDER_RE = re.compile(r"@[^@\s]+@")
# A comment marker must start a token: line start or after one of these characters.
_COMMENT_BOUNDARY = frozenset(" \t{}();,\"")
_WORD_STOP = frozenset(" \t\r\n\x0b\x0c{}()=;\",")
_PUNCT = frozenset("{}()=;")


class LexError(Exception):
    def __init__(self, detail: str, line: int, text: str) -> None:
        super().__init__(detail)
        self.detail = detail
        self.line = line
        self.text = text


@dataclass(frozen=True)
class Token:
    kind: str      # "word" | "string" | "punct"
    value: str     # word text, string content without quotes, or the punctuation char
    line: int      # 1-based source line
    start: int     # offset into the normalised text
    end: int


def normalise_text(text: str) -> str:
    if text.startswith("\ufeff"):
        text = text[1:]
    return text.replace("\r\n", "\n")


def decode_bytes(data: bytes) -> str:
    """Decode command-file bytes; raises ``UnicodeDecodeError`` for non UTF-8 input."""
    if data.startswith(b"\xef\xbb\xbf"):
        data = data[3:]
    return normalise_text(data.decode("utf-8"))


def sha256_of_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def find_comment_start(line: str, start: int = 0) -> Optional[int]:
    """Index of the first ``*``/``#`` that starts a token outside double quotes."""
    in_string = False
    i = start
    while i < len(line):
        char = line[i]
        if char == '"':
            in_string = not in_string
        elif not in_string and char in "*#" and (i == 0 or line[i - 1] in _COMMENT_BOUNDARY):
            return i
        i += 1
    return None


def _blank(chars: List[str], start: int, end: int) -> None:
    for i in range(start, min(end, len(chars))):
        if chars[i] != "\n":
            chars[i] = " "


def collapse(text: str) -> str:
    return " ".join(text.split())


def scan_and_strip(text: str) -> Tuple[List[dict], str]:
    """Gate 1 scan plus comment stripping.

    Returns ``(preprocessor_findings, code)`` where ``code`` has the same length
    as ``text`` with comments, directive lines and Tcl blocks blanked.
    """
    findings: List[dict] = []
    code_lines: List[str] = []
    in_tcl = False
    for index, raw in enumerate(text.split("\n")):
        lineno = index + 1
        chars = list(raw)
        pos = 0
        if in_tcl:
            close = raw.find(")!")
            if close < 0:
                code_lines.append(" " * len(raw))
                continue
            _blank(chars, 0, close + 2)
            pos = close + 2
            in_tcl = False
        cut = find_comment_start(raw, pos)
        directive = None
        if cut is not None and raw[cut] == "#" and raw[pos:cut].strip() == "":
            match = _DIRECTIVE_RE.match(raw, cut)
            if match:
                directive = match.group(0)
                findings.append(finding("E101", lineno, None, directive, text=directive))
        if cut is not None:
            _blank(chars, cut, len(raw))
        scan_end = len(raw) if (cut is None or directive) else cut
        segments: List[Tuple[int, int]] = []
        seek = pos
        while True:
            open_at = raw.find("!(", seek, scan_end)
            if open_at < 0:
                segments.append((seek, scan_end))
                break
            segments.append((seek, open_at))
            close_at = raw.find(")!", open_at + 2, scan_end)
            if close_at < 0:
                block_text = collapse(raw[open_at:scan_end])
                findings.append(finding("E103", lineno, None, block_text))
                _blank(chars, open_at, len(raw))
                in_tcl = True
                break
            block_text = collapse(raw[open_at:close_at + 2])
            findings.append(finding("E103", lineno, None, block_text))
            _blank(chars, open_at, close_at + 2)
            seek = close_at + 2
        for seg_start, seg_end in segments:
            for match in PLACEHOLDER_RE.finditer(raw, seg_start, seg_end):
                findings.append(finding("E102", lineno, None, match.group(0), text=match.group(0)))
        code_lines.append("".join(chars))
    return findings, "\n".join(code_lines)


def strip_comments(text: str) -> str:
    """Comment-stripped view of ``text`` (same length; quoted strings untouched)."""
    return scan_and_strip(normalise_text(text))[1]


def tokenize(code: str) -> List[Token]:
    tokens: List[Token] = []
    i = 0
    n = len(code)
    line = 1
    while i < n:
        char = code[i]
        if char == "\n":
            line += 1
            i += 1
            continue
        if char.isspace() or char == ",":
            i += 1
            continue
        if char in _PUNCT:
            tokens.append(Token("punct", char, line, i, i + 1))
            i += 1
            continue
        if char == '"':
            close = code.find('"', i + 1)
            newline = code.find("\n", i + 1)
            if close < 0 or (0 <= newline < close):
                end = newline if newline >= 0 else n
                raise LexError("unterminated string literal", line, code[i:end].strip())
            tokens.append(Token("string", code[i + 1:close], line, i, close + 1))
            i = close + 1
            continue
        j = i
        while j < n and code[j] not in _WORD_STOP:
            j += 1
        tokens.append(Token("word", code[i:j], line, i, j))
        i = j
    return tokens


class Source:
    """Normalised command text plus its comment-stripped ``code`` view."""

    def __init__(self, text: str, path_name: str) -> None:
        self.text = text
        self.path_name = path_name
        self.sha256 = sha256_of_text(text)
        self.line_count = len(text.splitlines())
        self.preprocessor_findings, self.code = scan_and_strip(text)
        self._code_lines = self.code.split("\n")

    def line_text(self, line: int) -> str:
        """Comment-stripped source line, trimmed."""
        if 1 <= line <= len(self._code_lines):
            return collapse(self._code_lines[line - 1])
        return ""

    def slice(self, start: int, end: int) -> str:
        """Comment-stripped source between two offsets with whitespace collapsed."""
        return collapse(self.code[start:end])

    def tokens(self) -> List[Token]:
        return tokenize(self.code)
