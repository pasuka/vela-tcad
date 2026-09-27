"""Message catalogue for the SDevice cmd parser (contract C8).

Every id emitted by the parser must exist in the table of
``docs/sentaurus_sdevice_workflow.md`` section 4 with exactly the template
recorded here; ``tests/regression/test_sdevice_cmd_parser.py`` cross-checks the
two.  Templates are copied verbatim; ``{...}`` placeholders are substituted with
:func:`render`, which only touches the placeholder names it is given so that
literal braces such as ``File{Parameter}`` survive.
"""

from __future__ import annotations

import re
from typing import Dict, Optional, Tuple

WARNING = "warning"
FATAL = "fatal"

TEMPLATES: Dict[str, Tuple[str, str]] = {
    "W101": (WARNING, "Math keyword '{kw}' is unknown and ignored"),
    "W102": (WARNING, "Math keyword '{kw}' is recorded and ignored (no Vela equivalent)"),
    "W103": (WARNING, "Coupled option '{kw}' is recorded and ignored"),
    "W104": (WARNING, "Plot field '{name}' is not a known data name or alias; skipped"),
    "W105": (WARNING, "CurrentPlot coordinate sample '{text}' is ignored"),
    "W106": (WARNING, "CurrentPlot integral or extremum '{text}' is ignored"),
    "W107": (WARNING, "File keyword '{kw}' is recorded and ignored (extra output only)"),
    "W108": (WARNING, "Plot option '{kw}' is ignored"),
    "W109": (WARNING, "Keyword '{kw}' appears more than once; the last occurrence is used"),
    "W110": (WARNING, "Quasistationary step parameter '{kw}' is absent; sdevice default assumed"),
    "W111": (WARNING, "Extrapolate option '{kw}' is ignored"),
    "W112": (WARNING, "Electrode attribute 'Resistor' accepted as an alias of 'Resist' (Table 206)"),
    "W113": (WARNING, "Electrode attribute '{kw}={value}' equals the default and is dropped"),
    "W114": (WARNING, "Physics model '{key}' is recorded as metadata only"),
    "E101": (FATAL, "Workbench preprocessor directive '{text}' is not supported; preprocess the deck first"),
    "E102": (FATAL, "Workbench placeholder '{text}' is not expanded"),
    "E103": (FATAL, "Tcl block '!( ... )!' is not supported"),
    "E104": (FATAL, "Syntax error: {detail}"),
    "E105": (FATAL, "Section '{name}' is not supported"),
    "E106": (FATAL, "Physics model '{key}' is not supported"),
    "E107": (FATAL, "Physics token '{text}' is not recognised"),
    "E108": (FATAL, "Electrode attribute '{kw}' is not supported"),
    "E109": (FATAL, "Solve statement '{text}' is not supported"),
    "E110": (FATAL, "Goal form '{text}' is not supported"),
    "E111": (FATAL, "Electrode attribute '{kw}' is not defined in Table 206"),
    "E112": (FATAL, "File keyword '{kw}' changes inputs or physics and is not supported"),
    "E113": (FATAL, "Effective Physics model sets differ between semiconductor regions: {detail}"),
    "E114": (FATAL, "Grid file '{name}' has an unsupported extension (expected .tdr or .grd)"),
    "E115": (FATAL, "Command file cannot be decoded as UTF-8/ASCII text or is empty"),
    "E116": (FATAL, "Coupled equation '{name}' is not supported"),
    "E117": (FATAL, "Electrode Name=Regexp(...) is not supported"),
    "E301": (FATAL, "File{Parameter} is not declared: vela device requires a parameter file "
                    "(add Parameter=\"<name>.par\" or export the defaults with 'vela device -P')"),
}

_PLACEHOLDER = re.compile(r"\{(kw|text|name|key|detail|value)\}")


def severity_of(message_id: str) -> str:
    return TEMPLATES[message_id][0]


def exit_code_of(message_id: str) -> int:
    """Exit-code segment encoded in the id (E1xx -> 1, E3xx -> 3); warnings map to 0."""
    if message_id.startswith("W"):
        return 0
    return int(message_id[1])


def render(message_id: str, **fields: object) -> str:
    template = TEMPLATES[message_id][1]

    def substitute(match: "re.Match[str]") -> str:
        name = match.group(1)
        if name in fields:
            return str(fields[name])
        return match.group(0)

    return _PLACEHOLDER.sub(substitute, template)


def finding(message_id: str, line: Optional[int], section: Optional[str], source_text: str,
            **fields: object) -> dict:
    """Build a finding object in the C1 ``finding`` shape.

    ``source_text`` is the offending source token; it also fills the ``{text}``
    placeholder unless the caller passes an explicit ``text=`` field.
    """
    fields.setdefault("text", source_text)
    return {
        "id": message_id,
        "severity": severity_of(message_id),
        "line": line,
        "section": section,
        "text": source_text,
        "message": render(message_id, **fields),
    }
