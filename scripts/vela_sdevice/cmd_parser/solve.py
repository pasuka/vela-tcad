"""Solve section: statements, execution stages and the analysis classification."""

from __future__ import annotations

from typing import Dict, List, Optional

from .lexer import Token
from .syntax import (Block, ParseContext, SyntaxIssue, TokenStream, assignment_text, is_number, parse_int,
                     parse_number)

SECTION = "Solve"
_SUPPORTED_EQUATIONS = {"poisson": "Poisson", "electron": "Electron", "hole": "Hole",
                        "equantumpotential": "eQuantumPotential"}
_KNOWN_EQUATIONS = {**_SUPPORTED_EQUATIONS, "contact": "Contact", "circuit": "Circuit",
                    "temperature": "Temperature", "etemperature": "eTemperature", "htemperature": "hTemperature",
                    "hquantumpotential": "hQuantumPotential"}
_UNSUPPORTED_STATEMENTS = {"transient", "accoupled", "plugin", "continuation", "unset", "system",
                           "optics", "cyclic", "harmonicbalance", "acplot", "quasistationary_dg"}
_T_CONTROL_KEYS = {"initialstep": "initial_step", "minstep": "min_step", "maxstep": "max_step",
                   "increment": "increment", "decrement": "decrement"}
_REASONS = {
    "R1": "R1: no Quasistationary statement",
    "R2": "R2: one Quasistationary with a voltage goal and an avalanche model, ionization integral or current break criterion",
    "R3": "R3: one Quasistationary with a voltage goal and no avalanche model, ionization integral or current break criterion",
    "R3c": "R3: one Quasistationary with a current goal",
    "R4": "R4: more than one Quasistationary statement",
    "R5": "R5: deck rejected before its sections were evaluated (preprocessor constructs or undecodable file)",
}


def _parse_time_spec(tokens: List[Token], ctx: ParseContext) -> Optional[dict]:
    """``Range=(a b) Intervals=n`` or ``a; b; c`` inside ``Time=( ... )``."""
    stream = TokenStream(tokens, ctx.source)
    values: List[float] = []
    spec: dict = {}
    while not stream.at_end():
        token = stream.next()
        if token.kind == "punct" and token.value == ";":
            continue
        if token.kind == "word" and token.value.lower() == "range":
            stream.expect_punct("=")
            _, inner, _ = stream.take_group()
            numbers = [parse_number(t, "Range") for t in inner if t.kind == "word"]
            if len(numbers) != 2:
                raise SyntaxIssue("Range expects two values", token.line, token.value)
            spec["range"] = numbers
        elif token.kind == "word" and token.value.lower() == "intervals":
            stream.expect_punct("=")
            spec["intervals"] = parse_int(stream.take_value(), "Intervals")
        elif token.kind == "word" and is_number(token.value):
            values.append(float(token.value))
        else:
            raise SyntaxIssue(f"unexpected '{token.value}' in Time specification", token.line, token.value)
    if "range" in spec or "intervals" in spec:
        if "range" not in spec or "intervals" not in spec:
            raise SyntaxIssue("Time=(Range=(a b) Intervals=n) needs both Range and Intervals",
                              tokens[0].line if tokens else None, "Time")
        return {"range": spec["range"], "intervals": spec["intervals"]}
    if values:
        return {"list": values}
    return None


def _parse_plot_statement(head: Token, args: List[Token], ctx: ParseContext) -> dict:
    stream = TokenStream(args, ctx.source)
    statement: dict = {"kind": "plot", "line": head.line, "file_prefix": "", "time": None, "no_overwrite": False}
    while not stream.at_end():
        token = stream.next()
        if token.kind != "word":
            raise SyntaxIssue(f"unexpected '{token.value}' in Plot", token.line, token.value)
        lower = token.value.lower()
        if lower == "fileprefix":
            stream.expect_punct("=")
            statement["file_prefix"] = stream.take_value().value
        elif lower == "time":
            stream.expect_punct("=")
            _, inner, _ = stream.take_group()
            statement["time"] = _parse_time_spec(inner, ctx)
        elif lower == "nooverwrite":
            statement["no_overwrite"] = True
        elif lower in ("-loadable", "loadable", "-plotloadable", "plotloadable", "datasetsfromgrid", "-explicit",
                       "explicit", "plotexplicit"):
            ctx.warn("W108", token.line, SECTION, token.value, kw=token.value)
        else:
            raise SyntaxIssue(f"unknown Plot option '{token.value}'", token.line, token.value)
    return statement


def _parse_current_plot_statement(head: Token, args: List[Token], ctx: ParseContext) -> dict:
    stream = TokenStream(args, ctx.source)
    time_spec = None
    while not stream.at_end():
        token = stream.next()
        if token.kind == "word" and token.value.lower() == "time":
            stream.expect_punct("=")
            _, inner, _ = stream.take_group()
            time_spec = _parse_time_spec(inner, ctx)
        else:
            raise SyntaxIssue(f"unknown CurrentPlot option '{token.value}'", token.line, token.value)
    return {"kind": "current_plot", "line": head.line, "time": time_spec}


def _parse_coupled(head: Token, stream: TokenStream, ctx: ParseContext) -> dict:
    statement: dict = {"kind": "coupled", "line": head.line, "equations": [], "iterations": None}
    options: Dict[str, str] = {}
    if stream.at_punct("("):
        _, inner, _ = stream.take_group()
        option_stream = TokenStream(inner, ctx.source)
        while not option_stream.at_end():
            key = option_stream.next()
            value: Optional[Token] = None
            if option_stream.at_punct("="):
                option_stream.next()
                value = option_stream.take_value()
            if key.value.lower() == "iterations" and value is not None:
                statement["iterations"] = parse_int(value, "Iterations")
            else:
                options[key.value] = value.value if value is not None else ""
                ctx.warn("W103", key.line, SECTION, assignment_text(key, value), kw=key.value)
    _, body, _ = stream.take_group()
    for token in body:
        if token.kind != "word":
            raise SyntaxIssue(f"unexpected '{token.value}' in Coupled", token.line, token.value)
        canonical = _KNOWN_EQUATIONS.get(token.value.lower(), token.value)
        if token.value.lower() not in _SUPPORTED_EQUATIONS:
            ctx.error("E116", token.line, SECTION, token.value, name=token.value)
        statement["equations"].append(canonical)
    if options:
        statement["options"] = options
    return statement


def _parse_goal(inner: List[Token], ctx: ParseContext, line: int, text: str) -> Optional[dict]:
    goal: dict = {}
    stream = TokenStream(inner, ctx.source)
    unsupported = False
    while not stream.at_end():
        key = stream.next()
        if key.kind != "word":
            raise SyntaxIssue(f"unexpected '{key.value}' in Goal", key.line, key.value)
        stream.expect_punct("=")
        lower = key.value.lower()
        if stream.peek() is not None and stream.peek().kind == "word" and stream.peek().value.lower() == "regexp":
            stream.next()
            stream.take_group()
            unsupported = True
            continue
        value = stream.take_value()
        if lower == "name":
            goal["contact"] = value.value
        elif lower == "voltage":
            goal["voltage_V"] = parse_number(value, "Voltage")
        elif lower == "current":
            goal["current_A_per_um"] = parse_number(value, "Current")
        else:
            unsupported = True
    if unsupported or "contact" not in goal or ("voltage_V" not in goal and "current_A_per_um" not in goal):
        ctx.error("E110", line, SECTION, text, text=text)
        return None
    return goal


def _parse_quasistationary(head: Token, stream: TokenStream, ctx: ParseContext) -> dict:
    statement: dict = {"kind": "quasistationary", "line": head.line, "goal": None,
                       "t_control": {"initial_step": None, "min_step": None, "max_step": None,
                                     "increment": None, "decrement": None, "do_zero": False},
                       "body": []}
    _, params, _ = stream.take_group()
    param_stream = TokenStream(params, ctx.source)
    goals = 0
    while not param_stream.at_end():
        key = param_stream.next()
        if key.kind != "word":
            raise SyntaxIssue(f"unexpected '{key.value}' in Quasistationary", key.line, key.value)
        lower = key.value.lower()
        if lower == "goal":
            opener, inner, closer = param_stream.take_group()
            goals += 1
            text = ctx.source.slice(key.start, closer.end)
            if goals > 1:
                ctx.error("E110", key.line, SECTION, text, text=text)
                continue
            statement["goal"] = _parse_goal(inner, ctx, key.line, text)
        elif lower == "dozero":
            statement["t_control"]["do_zero"] = True
        elif lower in _T_CONTROL_KEYS:
            param_stream.expect_punct("=")
            statement["t_control"][_T_CONTROL_KEYS[lower]] = parse_number(param_stream.take_value(), key.value)
        else:
            value: Optional[Token] = None
            if param_stream.at_punct("="):
                param_stream.next()
                if param_stream.at_punct("("):
                    param_stream.take_group()
                else:
                    value = param_stream.take_value()
            elif param_stream.at_punct("("):
                param_stream.take_group()
            ctx.warn("W103", key.line, SECTION, assignment_text(key, value), kw=key.value)
    if goals == 0:
        ctx.error("E110", head.line, SECTION, ctx.source.line_text(head.line), text="Quasistationary without Goal")
    _, body, _ = stream.take_group()
    statement["body"] = parse_statements(body, ctx, nested=True)
    if statement["goal"] is None:
        del statement["goal"]
    return statement


def parse_statements(tokens: List[Token], ctx: ParseContext, nested: bool = False) -> List[dict]:
    statements: List[dict] = []
    stream = TokenStream(tokens, ctx.source)
    try:
        while not stream.at_end():
            head = stream.next()
            if head.kind != "word":
                raise SyntaxIssue(f"unexpected '{head.value}' in Solve", head.line, head.value)
            lower = head.value.lower()
            if lower == "poisson" and not stream.at_punct("{") and not stream.at_punct("("):
                statements.append({"kind": "poisson", "line": head.line, "equations": ["Poisson"], "iterations": None})
            elif lower == "coupled":
                statements.append(_parse_coupled(head, stream, ctx))
            elif lower == "quasistationary" and not nested:
                statements.append(_parse_quasistationary(head, stream, ctx))
            elif lower == "plot":
                args: List[Token] = []
                last = head
                if stream.at_punct("("):
                    _, args, last = stream.take_group()
                statements.append(_parse_plot_statement(head, args, ctx))
            elif lower == "currentplot":
                args = []
                if stream.at_punct("("):
                    _, args, _ = stream.take_group()
                statements.append(_parse_current_plot_statement(head, args, ctx))
            elif lower == "newcurrentprefix":
                stream.expect_punct("=")
                value = stream.take_value()
                statements.append({"kind": "new_current_prefix", "line": head.line, "prefix": value.value})
            elif lower in ("save", "load"):
                _, args, _ = stream.take_group()
                prefix = ""
                arg_stream = TokenStream(args, ctx.source)
                while not arg_stream.at_end():
                    key = arg_stream.next()
                    arg_stream.expect_punct("=")
                    value = arg_stream.take_value()
                    if key.value.lower() == "fileprefix":
                        prefix = value.value
                statements.append({"kind": lower, "line": head.line, "file_prefix": prefix})
            elif lower == "set":
                opener, args, closer = stream.take_group()
                text = ctx.source.slice(head.start, closer.end)
                words = [t for t in args if t.kind in ("word", "string")]
                if len(words) == 3 and words[0].kind == "string" and words[1].value.lower() == "mode" \
                        and words[2].value.lower() in ("current", "voltage"):
                    statements.append({"kind": "set", "line": head.line, "contact": words[0].value,
                                       "mode": words[2].value.lower()})
                else:
                    ctx.error("E109", head.line, SECTION, text, text=text)
                    statements.append({"kind": "unsupported", "line": head.line, "text": text, "finding_id": "E109"})
            else:
                last = head
                if stream.at_punct("="):
                    stream.next()
                    last = stream.take_value()
                while stream.at_punct("(") or stream.at_punct("{"):
                    _, _, last = stream.take_group()
                text = ctx.source.slice(head.start, last.end)
                ctx.error("E109", head.line, SECTION, text, text=text)
                statements.append({"kind": "unsupported", "line": head.line, "text": text, "finding_id": "E109"})
    except SyntaxIssue as issue:
        ctx.syntax_error(issue, SECTION)
    return statements


def parse_solve_block(block: Block, ctx: ParseContext) -> List[dict]:
    return parse_statements(block.body, ctx)


def _plot_file(statement: dict, numbered: bool) -> dict:
    suffix = "_NNNN_des.tdr" if numbered else "_des.tdr"
    pattern = statement["file_prefix"] + suffix
    return {"file_prefix": statement["file_prefix"], "time": statement["time"],
            "no_overwrite": statement["no_overwrite"], "sentaurus_pattern": pattern,
            "vela_pattern": pattern[: -len(".tdr")] + ".vtu"}


def build_stages(statements: List[dict], current_resolved: str) -> List[dict]:
    stages: List[dict] = []
    prefix = ""
    current_mode: List[str] = []

    def new_stage(kind: str, index: int, statement: dict, equations: List[str], iterations: Optional[int]) -> dict:
        stage: dict = {"index": len(stages), "kind": kind, "statement_index": index, "equations": equations,
                       "iterations": iterations}
        return stage

    for index, statement in enumerate(statements):
        kind = statement["kind"]
        stage: Optional[dict] = None
        if kind in ("poisson", "coupled"):
            stage = new_stage("initial", index, statement, list(statement["equations"]), statement["iterations"])
        elif kind == "quasistationary":
            inner = [s for s in statement["body"] if s["kind"] == "coupled"]
            equations = list(inner[0]["equations"]) if inner else []
            iterations = inner[0]["iterations"] if inner else None
            stage = new_stage("sweep", index, statement, equations, iterations)
            if statement.get("goal") is not None:
                stage["goal"] = statement["goal"]
            stage["t_control"] = statement["t_control"]
        elif kind in ("save", "load"):
            stage = new_stage(kind, index, statement, [], None)
        elif kind == "new_current_prefix":
            prefix = statement["prefix"]
            continue
        elif kind == "set":
            if statement["mode"] == "current" and statement["contact"] not in current_mode:
                current_mode.append(statement["contact"])
            elif statement["mode"] == "voltage" and statement["contact"] in current_mode:
                current_mode.remove(statement["contact"])
            continue
        elif kind == "plot":
            if stages:
                stages[-1]["outputs"]["plot_files"].append(_plot_file(statement, numbered=False))
            continue
        else:
            continue
        if current_mode:
            stage["current_mode_contacts"] = list(current_mode)
        stage["current_prefix"] = prefix
        outputs = {"current_file": prefix + current_resolved, "plot_files": [], "current_plot_times": None}
        if kind == "quasistationary":
            for inner_statement in statement["body"]:
                if inner_statement["kind"] == "plot":
                    outputs["plot_files"].append(_plot_file(inner_statement, numbered=True))
                elif inner_statement["kind"] == "current_plot":
                    outputs["current_plot_times"] = inner_statement["time"]
        stage["outputs"] = outputs
        stage["depends_on"] = [stages[-1]["index"]] if stages else []
        stage["line"] = statement["line"]
        stages.append(stage)
    return stages


def build_analysis(stages: List[dict], physics: dict, math: dict, rejected_at_gate1: bool = False) -> dict:
    sweeps = [s for s in stages if s["kind"] == "sweep" and "goal" in s]
    analysis = {"kind": "equilibrium", "reason": _REASONS["R1"],
                "sweep_contacts": [s["goal"]["contact"] for s in sweeps],
                "final_contact": sweeps[-1]["goal"]["contact"] if sweeps else None,
                "final_voltage": sweeps[-1]["goal"].get("voltage_V") if sweeps else None}
    if rejected_at_gate1:
        analysis["kind"], analysis["reason"] = "rejected", _REASONS["R5"]
        return analysis
    if not sweeps:
        return analysis
    if len(sweeps) > 1:
        analysis["kind"], analysis["reason"] = "multi_stage", _REASONS["R4"]
        return analysis
    goal = sweeps[0]["goal"]
    if "voltage_V" not in goal:
        analysis["kind"], analysis["reason"] = "iv", _REASONS["R3c"]
        return analysis
    has_bv_trigger = any(entry["key"] == "Recombination.Avalanche" for entry in physics["effective"]) \
        or math["compute_ionization_integrals"] or math["break_at_ion_integral"] is not None \
        or bool(math["break_criteria"])
    if has_bv_trigger:
        analysis["kind"], analysis["reason"] = "bv", _REASONS["R2"]
    else:
        analysis["kind"], analysis["reason"] = "iv", _REASONS["R3"]
    return analysis
