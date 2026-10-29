"""Shared synthetic inputs and subprocess helpers for WP11 acceptance."""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

from tests.regression.test_sdevice_orchestrator import inputs

def synthetic_cmd(root, stem="pn2d_iv"):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    cmd = inputs(root, fixture=True)
    text = cmd.read_text(encoding="utf-8").replace("pn2d_iv", stem)
    target = root / (stem + "_des.cmd")
    if target != cmd:
        cmd.unlink()
    target.write_text(text, encoding="utf-8")
    return target

def isolated_environment():
    environment = {k: v for k, v in os.environ.items()
                   if not k.startswith(("VELA_", "PYTHON", "FAKE_RUNNER_"))}
    # Absolute executable paths are used throughout. No system Python on PATH.
    environment["PATH"] = str(Path(os.environ.get("SystemRoot", "C:/Windows")) / "System32") if os.name == "nt" else "/nonexistent"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "EIGEN_NUM_THREADS"):
        environment[key] = "1"
    return environment

def command(launcher, cmd, *args):
    return [str(launcher), "device", str(cmd), "--allow-lossy", *args]

def run(launcher, cmd, *args, env=None):
    return subprocess.run(command(launcher, cmd, *args), cwd=cmd.parent,
        env=env or isolated_environment(), capture_output=True, text=True, encoding="utf-8", timeout=90)

def numeric_outputs(root, prefix):
    result = {}
    for path in sorted(Path(root).glob(prefix + "*.vtu")):
        data = ET.parse(path)
        arrays = []
        for node in data.findall(".//DataArray"):
            if node.get("format") != "ascii":
                raise AssertionError("Default single-thread acceptance expects ASCII VTU")
            tokens = (node.text or "").split()
            arrays.append((tuple(sorted(node.attrib.items())), tuple(float(t) for t in tokens)))
        result[path.name] = arrays
    if not result:
        raise AssertionError("No VTU output published")
    return result

def manifest(cmd):
    return json.loads((cmd.parent / (cmd.stem + "_vela") / "run_manifest.json").read_text(encoding="utf-8"))

