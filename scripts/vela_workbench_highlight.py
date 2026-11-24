"""Display-only highlighting for original Sentaurus SDevice command decks."""

import re

from vela_sdevice.cmd_parser.lexer import DIRECTIVES, find_comment_start


SECTIONS = frozenset(("File", "Electrode", "Physics", "Plot", "Math", "Solve", "CurrentPlot",
                      "Device", "System", "Quasistationary", "Coupled", "Poisson", "NewCurrentPrefix"))
DIRECTIVE = re.compile(r"\s*#(?:" + "|".join(DIRECTIVES) + r")\b")
TOKEN = re.compile(
    r'(?P<string>"(?:[^"\\]|\\.)*")|(?P<placeholder>@[^@\s]+@)|'
    r'(?P<number>(?<![\w.])[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?![\w.]))|'
    r'(?P<keyword>[A-Za-z_][\w.-]*)'
)
PLACEHOLDER = re.compile(r"@[^@\s]+@")


def script_spans(source):
    """Yield (tag, line, start column, end column) without parsing or changing the deck."""
    for row, line in enumerate(source.split("\n"), 1):
        start = 0
        directive = DIRECTIVE.match(line)
        if directive:
            start = directive.end()
            yield "directive", row, line.index("#"), start
        comment = find_comment_start(line, start)
        if comment is not None:
            yield "comment", row, comment, len(line)
        for match in TOKEN.finditer(line, start, comment if comment is not None else len(line)):
            tag = match.lastgroup
            if tag == "keyword" and match.group() not in SECTIONS:
                continue
            yield tag, row, match.start(), match.end()
            if tag == "string":
                for placeholder in PLACEHOLDER.finditer(line, match.start(), match.end()):
                    yield "placeholder", row, placeholder.start(), placeholder.end()
