"""SDevice command-file parser producing the execution IR v2 (contract C1).

Entry points::

    from vela_sdevice.cmd_parser import parse_cmd_file, parse_cmd_text
    ir = parse_cmd_file(Path("pn2d_iv_sdevice.cmd"))

The parser never reads the file system beyond the command file itself: mesh,
parameter-file and existence checks belong to the orchestrator gates.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional, Union

from . import catalogue
from .lexer import LexError, Source, decode_bytes, normalise_text
from .sections import (PhysicsScope, combine_physics, default_math, parse_current_plot_block,
                       parse_electrode_block, parse_file_block, parse_math_block, parse_physics_block,
                       parse_plot_block, skeleton_files)
from .solve import build_analysis, build_stages, parse_solve_block
from .syntax import ParseContext, SyntaxIssue, group_blocks

SCHEMA_ID = "vela.sentaurus_execution_ir.v2"
SYNTAX_VERSION = "T-2022.03"
_UNSUPPORTED_SECTIONS = {"thermode", "system", "device", "optics", "transient", "acsource", "accoupled",
                         "montecarlo", "ldp", "plugin"}

__all__ = ["parse_cmd_file", "parse_cmd_text", "SCHEMA_ID", "SYNTAX_VERSION"]


def parse_cmd_file(path: Union[str, "os.PathLike[str]"]) -> dict:
    path = Path(path)
    data = path.read_bytes()
    try:
        text = decode_bytes(data)
    except UnicodeDecodeError:
        return _undecodable_ir(path.name)
    return parse_cmd_text(text, path_name=path.name)


def parse_cmd_text(text: str, *, path_name: str) -> dict:
    source = Source(normalise_text(text), path_name)
    if not source.text.strip():
        return _undecodable_ir(path_name, source)
    ctx = ParseContext(source)
    files = None
    electrodes: List[dict] = []
    physics_scopes: List[PhysicsScope] = []
    plot = {"explicit": False, "fields": []}
    current_plot = {"entries": []}
    math = default_math()
    statements: List[dict] = []
    seen_sections = set()

    if not source.preprocessor_findings:
        try:
            tokens = source.tokens()
        except LexError as error:
            ctx.error("E104", error.line, None, error.text, detail=error.detail)
            tokens = []
        blocks, issue = group_blocks(tokens, source)
        for block in blocks:
            name = block.name.value.lower()
            header = source.line_text(block.name.line)
            if name == "file":
                files = parse_file_block(block, ctx, path_name)
            elif name == "electrode":
                electrodes.extend(parse_electrode_block(block, ctx))
            elif name == "physics":
                physics_scopes.append(parse_physics_block(block, ctx))
            elif name == "plot":
                plot = parse_plot_block(block, ctx)
            elif name == "currentplot":
                current_plot = parse_current_plot_block(block, ctx)
            elif name == "math":
                parse_math_block(block, ctx, math)
            elif name == "solve":
                statements.extend(parse_solve_block(block, ctx))
            elif name in _UNSUPPORTED_SECTIONS:
                ctx.error("E105", block.name.line, None, header, name=block.name.value)
            else:
                ctx.error("E104", block.name.line, None, header, detail=f"unknown section '{block.name.value}'")
            seen_sections.add(name)
        if issue is not None:
            ctx.syntax_error(issue, None)
        for required in ("file", "electrode", "solve"):
            if required not in seen_sections:
                ctx.error("E104", None, None, "", detail=f"missing {required.capitalize()} section")

    if files is None:
        files = skeleton_files(path_name)
    physics = combine_physics(physics_scopes, ctx)
    stages = build_stages(statements, files["current"]["resolved"])
    rejected_at_gate1 = bool(source.preprocessor_findings)
    analysis = build_analysis(stages, physics, math, rejected_at_gate1=rejected_at_gate1)
    verdict = _verdict(source.preprocessor_findings, ctx.errors)
    return {
        "schema": SCHEMA_ID,
        "source": {"path": path_name, "sha256": source.sha256, "line_count": source.line_count,
                   "syntax_version": SYNTAX_VERSION},
        "verdict": verdict,
        "preprocessor_findings": source.preprocessor_findings,
        "errors": ctx.errors,
        "warnings": ctx.warnings,
        "files": files,
        "electrodes": electrodes,
        "physics": physics,
        "plot": plot,
        "current_plot": current_plot,
        "math": math,
        "solve": {"statements": statements},
        "stages": stages,
        "analysis": analysis,
    }


def _verdict(preprocessor_findings: List[dict], errors: List[dict]) -> dict:
    fatal_ids: List[str] = []
    for finding in list(preprocessor_findings) + list(errors):
        if finding["id"] not in fatal_ids:
            fatal_ids.append(finding["id"])
    exit_code = 0 if not fatal_ids else min(catalogue.exit_code_of(i) for i in fatal_ids)
    return {"accepted": not fatal_ids, "exit_code": exit_code, "fatal_ids": fatal_ids}


def _undecodable_ir(path_name: str, source: Optional[Source] = None) -> dict:
    error = catalogue.finding("E115", None, None, "")
    text = source.text if source is not None else ""
    sha = Source(text, path_name).sha256
    return {
        "schema": SCHEMA_ID,
        "source": {"path": path_name, "sha256": sha, "line_count": max(1, len(text.splitlines())),
                   "syntax_version": SYNTAX_VERSION},
        "verdict": {"accepted": False, "exit_code": 1, "fatal_ids": ["E115"]},
        "preprocessor_findings": [],
        "errors": [error],
        "warnings": [],
        "files": skeleton_files(path_name),
        "electrodes": [],
        "physics": combine_physics([], ParseContext(Source("", path_name))),
        "plot": {"explicit": False, "fields": []},
        "current_plot": {"entries": []},
        "math": default_math(),
        "solve": {"statements": []},
        "stages": [],
        "analysis": build_analysis([], {"effective": []}, default_math(), rejected_at_gate1=True),
    }
