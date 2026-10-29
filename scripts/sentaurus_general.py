"""Bounded General-expression AST. No Python eval/exec or attribute access."""
from __future__ import annotations
import ast
from functools import lru_cache
import math
import operator
import re

OPS = {ast.Add: ("add", operator.add), ast.Sub: ("sub", operator.sub),
       ast.Mult: ("mul", operator.mul), ast.Div: ("div", operator.truediv),
       ast.Pow: ("pow", operator.pow)}
FUNCS = {"abs": abs, "exp": math.exp, "sqrt": math.sqrt, "log": math.log}


def expression_tree(text: str, names: set[str]) -> dict:
    if len(text) > 4096:
        raise ValueError("General expression exceeds 4096 characters")
    try:
        root = ast.parse(text.replace("^", "**"), mode="eval")
    except (SyntaxError, RecursionError) as error:
        raise ValueError("invalid General expression") from error
    if sum(1 for _ in ast.walk(root)) > 256:
        raise ValueError("General expression exceeds AST budget")

    def visit(node, depth=0):
        if depth > 32:
            raise ValueError("General expression exceeds nesting budget")
        child = lambda n: visit(n, depth+1)
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            value = float(node.value)
            if not math.isfinite(value):
                raise ValueError("nonfinite General constant")
            return {"op": "number", "value": value}
        if isinstance(node, ast.Name) and node.id in names:
            return {"op": "symbol", "name": node.id}
        if isinstance(node, ast.BinOp) and type(node.op) in OPS:
            return {"op": OPS[type(node.op)][0], "args": [child(node.left), child(node.right)]}
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            return {"op": "neg" if isinstance(node.op, ast.USub) else "pos", "args": [child(node.operand)]}
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and
                node.func.id in FUNCS and len(node.args) == 1 and not node.keywords):
            return {"op": node.func.id, "args": [child(node.args[0])]}
        raise ValueError(f"unsupported General AST node: {type(node).__name__}")
    return visit(root.body)


def _compile(tree, parameters):
    op = tree["op"]
    if op == "number":
        return lambda xyz: tree["value"]
    if op == "symbol":
        name = tree["name"]
        if name in "xyz" and len(name) == 1:
            return lambda xyz: xyz["xyz".index(name)]
        return lambda xyz: parameters[name]
    args = [_compile(a, parameters) for a in tree["args"]]
    if op in FUNCS:
        return lambda xyz: FUNCS[op](args[0](xyz))
    if op in ("neg", "pos"):
        return lambda xyz: (-1 if op == "neg" else 1) * args[0](xyz)
    binary = {value[0]: value[1] for value in OPS.values()}
    def result(xyz):
        a, b = args[0](xyz), args[1](xyz)
        if op == "pow" and abs(b) > 64:
            raise ValueError("General exponent exceeds 64")
        return binary[op](a, b)
    return result


def initialise(text: str, globals_: dict[str, float]) -> dict[str, float]:
    """First assignment wins in the General global namespace (source order)."""
    if len(text) > 4096:
        raise ValueError("General init exceeds budget")
    for assignment in text.split(";"):
        if not assignment.strip():
            continue
        match = re.fullmatch(r"\s*([A-Za-z_][A-Za-z_0-9]*)\s*=\s*(.+)\s*", assignment)
        if not match or match[1] in {"x", "y", "z"}:
            raise ValueError("General init requires scalar named assignments")
        name, expr = match[1], match[2]
        tree = expression_tree(expr, set(globals_))
        if name not in globals_:
            try:
                value = float(_compile(tree, globals_)((0.0, 0.0, 0.0)))
            except (ValueError, ArithmeticError, TypeError) as error:
                raise ValueError(f"invalid General initializer: {name}") from error
            if not math.isfinite(value):
                raise ValueError("nonfinite General initializer")
            globals_[name] = value
    return dict(globals_)


@lru_cache(maxsize=256)
def _program(expression: str, parameters: tuple):
    values = dict(parameters)
    if not all(type(v) in (float, int) and math.isfinite(v) for v in values.values()):
        raise ValueError("General parameters must be finite scalars")
    tree = expression_tree(expression, set(values) | {"x", "y", "z"})
    return tree, _compile(tree, values)


def validate_profile(profile):
    tree, _ = _program(profile["expression"], tuple(sorted(profile["parameters"].items())))
    if tree != profile["expression_tree"]:
        raise ValueError("General AST does not match the source expression")
    value = profile["default_value_cm3"]
    if type(value) not in (float, int) or not math.isfinite(value) or value < 0:
        raise ValueError("General default concentration must be finite and nonnegative")


def evaluate(profile: dict, x: float, y: float) -> tuple[float, bool]:
    _, program = _program(profile["expression"], tuple(sorted(profile["parameters"].items())))
    try:
        value = float(program((x, y, 0.0)))
        if not math.isfinite(value):
            raise ValueError("nonfinite General evaluation")
    except (ArithmeticError, ValueError, TypeError):
        return float(profile["default_value_cm3"]), True
    if value < 0:
        raise ValueError("General doping concentration is negative")
    return value, False
