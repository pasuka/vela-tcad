#!/usr/bin/env python3
"""Reference artifact tools: ``register`` / ``fetch`` / ``verify`` (WP9, plan section 12).

Sentaurus reference artifacts stay outside the public repository (storage
scheme S3: ``/home/tcad/vela_reference_artifacts/<case>/<version>/`` on the
Sentaurus VM).  The repository only carries ``reference_tcad/<case>/
artifact_manifest.json`` (contract C7, ``schemas/vela.reference_artifact_manifest.v1.schema.json``).

* ``register`` lists the archive through a transport backend, hashes the files
  on the store side, assigns categories A-E from the directory, derives the
  producer heuristically from the archived command files and writes a
  C7-valid manifest.
* ``fetch`` copies every listed artifact to ``reference_staging/<case>/<version>/``
  and re-hashes the copies.
* ``verify`` runs the six checks of plan section 12.4 (existence/hash, source
  consistency, structure, completeness, version, report) and writes
  ``verify_report.json`` plus a Markdown summary.  Exit code 0 = all checks
  pass, 1 = a check failed, 2 = the store is unreachable and no local copy
  exists ("unable to verify").  ``--offline`` evaluates the manifest-only
  checks (2 and 5) without touching any file.

Transport backends: :class:`LocalDirBackend` (``--backend local:<root>``, used
by the tests) and :class:`SshBackend` (``ssh``/``scp`` subprocesses, never a
shell string).  ``.tdr`` files are HDF5 and cannot be parsed with the standard
library; their structural inspection is optional and uses
``build/sentaurus_import`` when present (``--tdr-importer``).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import posixpath
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Optional, Sequence

SCRIPTS = Path(__file__).resolve().parent
REPO = SCRIPTS.parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import compare_sentaurus_tdr_tdx as dfise  # noqa: E402
from vela_sdevice import schema_validator as sv  # noqa: E402
from vela_sdevice.cmd_parser.lexer import strip_comments  # noqa: E402
from vela_sdevice.plt import read_plt_document  # noqa: E402

SCHEMA_ID = "vela.reference_artifact_manifest.v1"
SCHEMA_PATH = REPO / "schemas" / "vela.reference_artifact_manifest.v1.schema.json"
REPORT_SCHEMA_ID = "vela.reference_artifact_verify_report.v1"

CATEGORY_DIRS = ("source", "mesh", "plt", "tdr", "log")
CATEGORY_BY_DIR = {"source": "A", "plt": "B", "mesh": "C", "tdr": "D", "log": "E"}
DEFAULT_HOST = "sentaurus"
DEFAULT_STORE_ROOT = "/home/tcad/vela_reference_artifacts"
STAGING_ENV = "VELA_REFERENCE_ARTIFACTS_DIR"
DEFAULT_STAGING = REPO / "reference_staging"
DEFAULT_REPORT_DIR = REPO / "docs" / "validation"
EXPECTED_VERSION_PREFIX = "T-2022.03"

EXIT_OK = 0
EXIT_FAIL = 1
EXIT_UNAVAILABLE = 2

BINARY_SUFFIXES = frozenset({".tdr", ".sav"})
MAX_TEXT_FETCH_BYTES = 64 * 1024 * 1024
MAX_TDR_INSPECT_BYTES = 256 * 1024 * 1024
GOAL_TOLERANCE_ABS = 1e-9
GOAL_TOLERANCE_REL = 1e-6

_SDEVICE_BANNER_RE = re.compile(r"Sentaurus Device")
_VERSION_RE = re.compile(r"Version\s+([A-Z]-\d{4}\.\d{2}(?:-SP\d+)?)")
_BUILD_RE = re.compile(r"\(\s*(\d+(?:\.\d+)?)\s*,")
_GOOD_BYE_RE = re.compile(r"Good Bye")
_HEX64_RE = re.compile(r"^[0-9a-f]{64}$")
_NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


# --------------------------------------------------------------------------- helpers

def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def utc_from_epoch(epoch: float) -> str:
    return dt.datetime.fromtimestamp(int(epoch), dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode_text(data: bytes) -> str:
    """BOM-free, LF-normalised text (UTF-8, latin-1 fallback for stray bytes)."""
    if data.startswith(b"\xef\xbb\xbf"):
        data = data[3:]
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        text = data.decode("latin-1")
    return text.replace("\r\n", "\n")


def lf_sha256(data: bytes) -> str:
    """Hash rule shared with the execution IR ``source.sha256`` and ``producer.cmd_sha256``."""
    return sha256_bytes(decode_text(data).encode("utf-8"))


def deck_tokens(data: bytes) -> list[str]:
    """Comment-stripped (``*``/``#``), whitespace-collapsed token view of a command or parameter file."""
    return strip_comments(decode_text(data)).split()


def decks_equivalent(left: bytes, right: bytes) -> bool:
    return deck_tokens(left) == deck_tokens(right)


def split_artifact_path(rel: str) -> tuple[str, Optional[str], str]:
    """``plt/m0/x.plt`` -> ``("plt", "m0", "x.plt")``; ``plt/x.plt`` -> ``("plt", None, "x.plt")``."""
    parts = rel.split("/")
    variant = "/".join(parts[1:-1]) or None
    return parts[0], variant, parts[-1]


def scratch_directory(prefix: str) -> tempfile.TemporaryDirectory:
    build = REPO / "build"
    if build.is_dir():
        return tempfile.TemporaryDirectory(prefix=prefix, dir=str(build))
    return tempfile.TemporaryDirectory(prefix=prefix)


def load_schema() -> dict:
    return sv.load_schema(SCHEMA_PATH)


def manifest_path_for(case: str) -> Path:
    return REPO / "reference_tcad" / case / "artifact_manifest.json"


def load_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError(f"{path}: manifest is not a JSON object")
    return manifest


def write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


def default_staging_root() -> Path:
    value = os.environ.get(STAGING_ENV)
    return Path(value) if value else DEFAULT_STAGING


def parse_version_banner(text: str) -> Optional[str]:
    """``T-2022.03-SP2 (0.7745337)`` from a Sentaurus log header, or ``None``."""
    match = _VERSION_RE.search(text)
    if match is None:
        return None
    version = match.group(1)
    build = _BUILD_RE.search(text, match.end(), match.end() + 400)
    return f"{version} ({build.group(1)})" if build else version


def log_exit_status(text: str) -> str:
    if re.search(r"(?mi)^\s*instance\b[^\r\n]*\bnot loaded:", text):
        return "failed: Load instance not loaded"
    if failure := re.search(r"Step-size (?:less than MinStep|is too small)", text, re.IGNORECASE):
        return "failed: " + failure.group(0)
    return "Good Bye" if _GOOD_BYE_RE.search(text) else "no Good Bye marker"


# --------------------------------------------------------------------------- transport backends

class TransportError(RuntimeError):
    """The store could not be listed, hashed or copied."""


class TransportBackend(ABC):
    kind = "abstract"
    root = ""

    def describe(self) -> str:
        return f"{self.kind}:{self.root}"

    @abstractmethod
    def store_descriptor(self) -> dict:
        """``store`` field of the manifest (C7)."""

    @abstractmethod
    def is_reachable(self) -> bool:
        ...

    @abstractmethod
    def list_files(self, rel_dir: str) -> list[str]:
        """Sorted POSIX paths (relative to the root) of all files below ``rel_dir``; ``[]`` if absent."""

    @abstractmethod
    def sha256(self, rel_path: str) -> str:
        ...

    @abstractmethod
    def size(self, rel_path: str) -> int:
        ...

    @abstractmethod
    def fetch(self, rel_path: str, local_path: Path) -> None:
        ...

    def mtime(self, rel_path: str) -> Optional[float]:
        return None

    def stat_and_hash(self, rel_path: str) -> tuple[int, str]:
        return self.size(rel_path), self.sha256(rel_path)

    def read_bytes(self, rel_path: str) -> bytes:
        with scratch_directory("vela_ref_read_") as scratch:
            target = Path(scratch) / PurePosixPath(rel_path).name
            self.fetch(rel_path, target)
            return target.read_bytes()


class LocalDirBackend(TransportBackend):
    """A local directory that mirrors the VM layout (tests, offline mirrors, external disks)."""

    kind = "local"

    def __init__(self, root: str | os.PathLike) -> None:
        self.root_path = Path(root)
        self.root = str(root)

    def store_descriptor(self) -> dict:
        return {"kind": "local", "root": self.root}

    def _path(self, rel_path: str) -> Path:
        return self.root_path.joinpath(*PurePosixPath(rel_path).parts) if rel_path else self.root_path

    def _require(self, rel_path: str) -> Path:
        path = self._path(rel_path)
        if not path.is_file():
            raise TransportError(f"{self.describe()}: no such file {rel_path!r}")
        return path

    def is_reachable(self) -> bool:
        return self.root_path.is_dir()

    def list_files(self, rel_dir: str) -> list[str]:
        base = self._path(rel_dir)
        if not base.is_dir():
            return []
        return sorted(path.relative_to(self.root_path).as_posix() for path in base.rglob("*") if path.is_file())

    def sha256(self, rel_path: str) -> str:
        return sha256_file(self._require(rel_path))

    def size(self, rel_path: str) -> int:
        return self._require(rel_path).stat().st_size

    def mtime(self, rel_path: str) -> Optional[float]:
        return self._require(rel_path).stat().st_mtime

    def fetch(self, rel_path: str, local_path: Path) -> None:
        source = self._require(rel_path)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, local_path)

    def read_bytes(self, rel_path: str) -> bytes:
        return self._require(rel_path).read_bytes()


class SshBackend(TransportBackend):
    """``ssh <host> find/sha256sum/stat`` and ``scp`` through ``subprocess.run`` argument lists."""

    kind = "ssh"

    def __init__(self, host: str = DEFAULT_HOST, root: str = DEFAULT_STORE_ROOT, *,
                 ssh: Optional[str] = None, scp: Optional[str] = None,
                 timeout: float = 900.0, connect_timeout: int = 10) -> None:
        self.host = host
        self.root = root
        self.ssh = ssh or self.default_executable("ssh")
        self.scp = scp or self.default_executable("scp")
        self.timeout = timeout
        self.connect_timeout = connect_timeout

    @staticmethod
    def default_executable(name: str) -> str:
        # The MSYS ssh on PATH ignores the user's ssh config; Windows OpenSSH reads it.
        if os.name == "nt":
            system_root = Path(os.environ.get("SystemRoot", r"C:\Windows"))
            candidate = system_root / "System32" / "OpenSSH" / f"{name}.exe"
            if candidate.is_file():
                return str(candidate)
        return shutil.which(name) or name

    def describe(self) -> str:
        return f"ssh:{self.host}:{self.root}"

    def store_descriptor(self) -> dict:
        return {"kind": "ssh", "host": self.host, "root": self.root}

    def _batch_options(self) -> list[str]:
        return ["-o", "BatchMode=yes", "-o", f"ConnectTimeout={self.connect_timeout}"]

    def _run(self, args: Sequence[str], *, timeout: Optional[float] = None,
             cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
        try:
            return subprocess.run(list(args), capture_output=True, text=True, check=False,
                                  timeout=timeout or self.timeout, cwd=str(cwd) if cwd else None)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise TransportError(f"{args[0]} failed: {exc}") from exc

    def _remote(self, command: str) -> str:
        proc = self._run([self.ssh, *self._batch_options(), self.host, command])
        if proc.returncode != 0:
            raise TransportError(f"ssh {self.host} {command!r} exited {proc.returncode}: {proc.stderr.strip()}")
        return proc.stdout

    def _remote_path(self, rel_path: str) -> str:
        return posixpath.join(self.root, rel_path) if rel_path else self.root

    def is_reachable(self) -> bool:
        try:
            proc = self._run([self.ssh, *self._batch_options(), self.host, "true"],
                             timeout=self.connect_timeout + 30)
        except TransportError:
            return False
        return proc.returncode == 0

    def list_files(self, rel_dir: str) -> list[str]:
        quoted = shlex.quote(self._remote_path(rel_dir))
        output = self._remote(f"if [ -d {quoted} ]; then find {quoted} -type f -printf '%P\\n'; fi")
        return sorted(posixpath.join(rel_dir, line.strip()) for line in output.splitlines() if line.strip())

    def sha256(self, rel_path: str) -> str:
        output = self._remote(f"sha256sum -- {shlex.quote(self._remote_path(rel_path))}")
        token = output.split()[0] if output.split() else ""
        if not _HEX64_RE.match(token):
            raise TransportError(f"unexpected sha256sum output for {rel_path!r}: {output.strip()!r}")
        return token

    def size(self, rel_path: str) -> int:
        output = self._remote(f"stat -c %s -- {shlex.quote(self._remote_path(rel_path))}")
        return int(output.split()[0])

    def mtime(self, rel_path: str) -> Optional[float]:
        output = self._remote(f"stat -c %Y -- {shlex.quote(self._remote_path(rel_path))}")
        return float(output.split()[0])

    def stat_and_hash(self, rel_path: str) -> tuple[int, str]:
        quoted = shlex.quote(self._remote_path(rel_path))
        output = self._remote(f"stat -c %s -- {quoted} && sha256sum -- {quoted}")
        lines = [line for line in output.splitlines() if line.strip()]
        if len(lines) < 2 or not _HEX64_RE.match(lines[1].split()[0]):
            raise TransportError(f"unexpected stat/sha256sum output for {rel_path!r}: {output.strip()!r}")
        return int(lines[0].split()[0]), lines[1].split()[0]

    def fetch(self, rel_path: str, local_path: Path) -> None:
        local_path.parent.mkdir(parents=True, exist_ok=True)
        remote = f"{self.host}:{shlex.quote(self._remote_path(rel_path))}"
        # Relative destination + cwd keeps Windows drive letters away from scp's host:path parsing.
        proc = self._run([self.scp, *self._batch_options(), "-q", remote, local_path.name], cwd=local_path.parent)
        if proc.returncode != 0 or not local_path.is_file():
            raise TransportError(f"scp {remote} exited {proc.returncode}: {proc.stderr.strip()}")


def make_backend(spec: Optional[str], store: Optional[dict] = None) -> Optional[TransportBackend]:
    """``--backend local:<root>`` / ``ssh:<host>:<root>`` or the manifest ``store`` field."""
    if spec:
        kind, _, rest = spec.partition(":")
        if kind == "local":
            if not rest:
                raise ValueError("--backend local:<root> needs a directory")
            return LocalDirBackend(rest)
        if kind == "ssh":
            host, _, root = rest.partition(":")
            return SshBackend(host or DEFAULT_HOST, root or DEFAULT_STORE_ROOT)
        raise ValueError(f"unsupported backend spec {spec!r} (use local:<root> or ssh:<host>:<root>)")
    if store:
        kind = store.get("kind")
        if kind == "ssh":
            return SshBackend(store.get("host", DEFAULT_HOST), store.get("root", DEFAULT_STORE_ROOT))
        if kind == "local":
            return LocalDirBackend(store.get("root", ""))
    return None


# --------------------------------------------------------------------------- command file analysis

@dataclass(frozen=True)
class Goal:
    prefix: str
    name: str
    voltage: float
    # Retain the positional target field for existing voltage callers.
    quantity: str = "voltage"


def _find_blocks(code: str, keyword: str) -> list[str]:
    """Bodies of ``keyword { ... }`` blocks (brace matched) in comment-stripped code."""
    bodies: list[str] = []
    for match in re.finditer(rf"\b{keyword}\s*\{{", code):
        depth, index = 1, match.end()
        while index < len(code) and depth:
            if code[index] == "{":
                depth += 1
            elif code[index] == "}":
                depth -= 1
            index += 1
        bodies.append(code[match.end():index - 1])
    return bodies


def _sub_blocks(body: str) -> list[str]:
    blocks: list[str] = []
    depth, start = 0, -1
    for index, char in enumerate(body):
        if char == "{":
            if depth == 0:
                start = index + 1
            depth += 1
        elif char == "}" and depth:
            depth -= 1
            if depth == 0:
                blocks.append(body[start:index])
    return blocks


def sdevice_output_name(value: str, extension: str) -> str:
    """Appendix C Table 194: a basename containing ``.`` is literal, else ``_des.<ext>`` is appended."""
    base = value.replace("\\", "/").rsplit("/", 1)[-1]
    return base if "." in base else f"{base}_des.{extension}"


def sdevice_log_name(value: str) -> str:
    """Observed T-2022.03 rule: ``_des.log`` is appended unless the name already ends with it."""
    base = value.replace("\\", "/").rsplit("/", 1)[-1]
    return base if base.endswith("_des.log") else f"{base}_des.log"


def classify_command_file(text: str, code: str) -> str:
    if re.search(r"\(sde(?:geo|dr|io|pe)?:", text) or re.search(r"\(define\s", text):
        return "sde"
    if re.search(r"\bElectrode\s*\{", code) or (re.search(r"\bFile\s*\{", code) and re.search(r"\bSolve\s*\{", code)):
        return "sdevice"
    if re.search(r"\b(?:Definitions|Placements)\s*\{", code) or re.search(r"\bTitle\s*\"", code):
        return "snmesh"
    if re.search(r"(?m)^\s*(?:line\s+[xyz]\b|init\b|struct\b|deposit\b|etch\b|implant\b|diffuse\b)", text):
        return "sprocess"
    return "other"


class Deck:
    """Regex-level view of an archived command file (enough for producer mapping and check 3)."""

    def __init__(self, name: str, data: bytes) -> None:
        self.name = name
        self.data = data
        self.text = decode_text(data)
        self.code = strip_comments(self.text)
        self.lf_sha256 = lf_sha256(data)
        self.kind = classify_command_file(self.text, self.code)
        self.file_names: dict[str, str] = {}
        self.electrodes: list[str] = []
        self.goals: list[Goal] = []
        self.current_prefixes: list[str] = []
        self.plot_prefixes: list[str] = []
        if self.kind == "sdevice":
            self._analyse_sdevice()

    def _analyse_sdevice(self) -> None:
        for body in _find_blocks(self.code, "File"):
            for key, value in re.findall(r'\b(\w+)\s*=\s*"([^"]*)"', body):
                self.file_names.setdefault(key, value)
        for body in _find_blocks(self.code, "Electrode"):
            for block in _sub_blocks(body):
                match = re.search(r'\bName\s*=\s*"([^"]+)"', block)
                if match and match.group(1) not in self.electrodes:
                    self.electrodes.append(match.group(1))
        pattern = re.compile(r'\bNewCurrentPrefix\s*=\s*"([^"]*)"|\bGoal\s*\{([^}]*)\}|\bFilePrefix\s*=\s*"([^"]*)"')
        for body in _find_blocks(self.code, "Solve"):
            prefix = ""
            for match in pattern.finditer(body):
                if match.group(1) is not None:
                    prefix = match.group(1)
                    if prefix not in self.current_prefixes:
                        self.current_prefixes.append(prefix)
                elif match.group(2) is not None:
                    name = re.search(r'\bName\s*=\s*"([^"]+)"', match.group(2))
                    voltage = re.search(rf"\bVoltage\s*=\s*({_NUMBER})", match.group(2))
                    current = re.search(rf"\bCurrent\s*=\s*({_NUMBER})", match.group(2))
                    if name and voltage:
                        self.goals.append(Goal(prefix, name.group(1), float(voltage.group(1))))
                    elif name and current:
                        self.goals.append(Goal(prefix, name.group(1), float(current.group(1)), "current"))
                elif match.group(3) is not None and match.group(3) not in self.plot_prefixes:
                    self.plot_prefixes.append(match.group(3))

    @property
    def stem(self) -> str:
        stem = self.name.rsplit(".", 1)[0]
        return stem[:-4] if stem.endswith("_des") else stem

    @property
    def parameter_file(self) -> Optional[str]:
        value = self.file_names.get("Parameter") or self.file_names.get("Parameters")
        return value.replace("\\", "/").rsplit("/", 1)[-1] if value else None

    def expected_outputs(self) -> dict[str, str]:
        names: dict[str, str] = {}
        current = self.file_names.get("Current")
        if current:
            for prefix in [""] + self.current_prefixes:
                names[prefix + sdevice_output_name(current, "plt")] = "plt"
        plot = self.file_names.get("Plot")
        if plot:
            names[sdevice_output_name(plot, "tdr")] = "tdr"
        for prefix in self.plot_prefixes:
            names[sdevice_output_name(prefix, "tdr")] = "tdr"
        for key in ("Save", "Load"):
            if self.file_names.get(key):
                names[sdevice_output_name(self.file_names[key], "sav")] = "sav"
        output = self.file_names.get("Output")
        if output:
            names[sdevice_log_name(output)] = "log"
        return names

    def output_stems(self) -> list[str]:
        stems = []
        for key in ("Plot", "Current", "Output", "Save"):
            value = self.file_names.get(key)
            if value:
                base = value.replace("\\", "/").rsplit("/", 1)[-1]
                stems.append(base.split(".")[0])
        return stems

    def match_score(self, filename: str, content: str = "") -> int:
        """How plausibly this sdevice deck produced ``filename`` (0 = no relation)."""
        if filename in self.expected_outputs():
            return 4
        for stem in self.output_stems():
            if stem and (filename.startswith(stem + "_") or filename.startswith(stem + ".")):
                return 3
        for prefix in self.current_prefixes:
            if prefix and filename.startswith(prefix):
                return 3
        workbench = re.match(r"^pp(\d+)_", self.name)
        if workbench and re.search(rf"(?:^|_)n{workbench.group(1)}_", filename):
            return 3
        if self.stem and (filename.startswith(self.stem + "_") or filename.startswith(self.stem + ".")):
            return 2
        if content:
            mentions = sum(content.count(name) for name in self.expected_outputs())
            mentions += sum(content.count(value) for value in self.file_names.values())
            if mentions:
                return 1
        return 0

    def goal_for_output(self, filename: str) -> Optional[Goal]:
        prefix = ""
        for candidate in self.current_prefixes:
            if candidate and filename.startswith(candidate) and len(candidate) > len(prefix):
                prefix = candidate
        matching = [goal for goal in self.goals if goal.prefix == prefix]
        return matching[-1] if matching else None


# --------------------------------------------------------------------------- structural summaries

def plt_dataset_names(text: str) -> list[str]:
    match = re.search(r"datasets\s*=\s*\[(.*?)\]", text, re.DOTALL)
    return re.findall(r'"([^"]*)"', match.group(1)) if match else []


def summarise_plt(text: str, goal: Optional[Goal] = None) -> dict:
    electrodes, rows = read_plt_document(text)
    summary: dict[str, Any] = {"datasets": plt_dataset_names(text), "rows": len(rows)}
    if not rows or not electrodes:
        return summary
    if goal is not None and goal.name in electrodes:
        swept = goal.name
    else:
        def span(name: str) -> float:
            values = [row[name].get("outer_voltage_V", 0.0) for row in rows]
            return max(values) - min(values)
        swept = max(electrodes, key=span)
    first = rows[0][swept].get("outer_voltage_V", 0.0)
    last = rows[-1][swept].get("outer_voltage_V", 0.0)
    summary["first_bias_V"] = first
    summary["last_bias_V"] = last
    if goal is not None and goal.name in electrodes:
        if goal.quantity == "current":
            summary["first_current_A_per_um"] = rows[0][goal.name]["total_current_A_per_um"]
            summary["last_current_A_per_um"] = rows[-1][goal.name]["total_current_A_per_um"]
            summary["goal_reached"] = current_goal_reached(summary["last_current_A_per_um"], goal.voltage)
        else:
            summary["goal_reached"] = goal_reached(last, goal.voltage)
    return summary


def goal_reached(last_bias: float, goal_voltage: float) -> bool:
    return abs(last_bias - goal_voltage) <= GOAL_TOLERANCE_ABS + GOAL_TOLERANCE_REL * abs(goal_voltage)


def current_goal_reached(current: float, target: float) -> bool:
    return abs(current - target) <= 1e-20 + GOAL_TOLERANCE_REL * abs(target)


def declared_break_evidence(deck: Deck, rows: list, logs: list[str]) -> Optional[str]:
    """Accept a declared early stop only with its final native evidence."""
    successful = [text for text in logs if log_exit_status(text) == "Good Bye"]
    if not successful or not rows:
        return None
    for body in _find_blocks(deck.code, "BreakCriteria"):
        for criterion in re.findall(r"\bCurrent\s*\(([^)]*)\)", body):
            contact = re.search(r'\bContact\s*=\s*"([^"]+)"', criterion)
            limit = re.search(rf"\bAbsVal\s*=\s*({_NUMBER})", criterion)
            if contact and limit and contact.group(1) in rows[-1]:
                current = abs(rows[-1][contact.group(1)]["total_current_A_per_um"])
                threshold = float(limit.group(1))
                if threshold > 0 and current >= threshold:
                    return f"declared Current AbsVal {threshold:g} A/um reached: {current:g} A/um"
    integral = re.search(rf"\bBreakAtIonIntegral\s*\(\s*([1-9]\d*)\s+({_NUMBER})\s*\)", deck.code)
    if integral:
        rank = int(integral.group(1))
        threshold = float(integral.group(2))
        if threshold <= 0:
            return None
        final_integrals = [rows[-1].get(name) for name in ("PhiElectron", "PhiHole")]
        if any(value is None or not math.isfinite(value) or value < threshold for value in final_integrals):
            return None
        for text in successful:
            values = re.findall(rf"Electron:\s*({_NUMBER})\s+Hole:\s*({_NUMBER})", text)
            if ("Finished quasi-stat., because ionization integrals meet specifications." in text
                    and values and min(map(float, values[-1])) >= threshold):
                return (f"native confirms declared ionization-integral stop (rank {rank}, threshold {threshold:g}); "
                        f"final PLT maximum integrals {final_integrals}, native Electron/Hole {values[-1]}; "
                        "PLT does not independently enumerate ranked paths")
    return None


def summarise_grd(path: Path) -> dict:
    grid = dfise.parse_grd(path)
    summary: dict[str, Any] = {}
    if grid.get("vertex_count") is not None:
        summary["nodes"] = int(grid["vertex_count"])
    if grid.get("element_count") is not None:
        summary["elements"] = int(grid["element_count"])
    materials = sorted({str(m) for m in grid.get("materials", []) if str(m).lower() != "contact"})
    if materials:
        summary["materials"] = materials
    return summary


def summarise_tdr_inventory(inventory: dict) -> dict:
    """Structural summary from ``sentaurus_import --inventory-json`` output."""
    summary: dict[str, Any] = {}
    if isinstance(inventory.get("vertex_count"), int):
        summary["nodes"] = inventory["vertex_count"]
    regions = inventory.get("regions") or []
    triangles = sum(int(region.get("triangles", 0)) for region in regions if isinstance(region, dict))
    if regions:
        summary["elements"] = triangles
    contacts = []
    materials = []
    geometry_regions = {region.get("name"): region for region in (inventory.get("geometry") or {}).get("regions", [])
                        if isinstance(region, dict)}
    for region in regions:
        if not isinstance(region, dict):
            continue
        if region.get("type") == 1:
            geometry = geometry_regions.get(region.get("name"), {})
            nodes = set()
            for edge in geometry.get("edges", []):
                nodes.update(int(v) for v in edge)
            nodes.update(int(v) for v in geometry.get("points", []))
            node_count = len(nodes) if nodes else int(region.get("points", 0))
            contacts.append({"name": str(region.get("name")), "node_count": node_count})
        elif region.get("type") == 0 and region.get("material") and region["material"] not in materials:
            materials.append(str(region["material"]))
    if contacts:
        summary["contacts"] = contacts
    if materials:
        summary["materials"] = materials
    datasets = sorted({str(f.get("name")) for f in inventory.get("fields", []) if isinstance(f, dict) and f.get("name")})
    if datasets:
        summary["datasets"] = datasets
    return summary


def resolve_tdr_importer(spec: Optional[str]) -> Optional[list[str]]:
    """``None`` disables TDR inspection; the default is ``build/sentaurus_import[.exe]`` when it exists."""
    if spec is not None:
        if spec.strip().lower() in {"", "none", "off"}:
            return None
        return shlex.split(spec, posix=os.name != "nt")
    if os.environ.get("VELA_TDR_IMPORTER"):
        return resolve_tdr_importer(os.environ["VELA_TDR_IMPORTER"])
    candidate = REPO / "build" / ("sentaurus_import.exe" if os.name == "nt" else "sentaurus_import")
    return [str(candidate)] if candidate.is_file() else None


def inspect_tdr(importer: Sequence[str], tdr_path: Path) -> tuple[bool, str, dict]:
    """Run the C++ importer; ``(ok, detail, summary)``."""
    with scratch_directory("vela_ref_tdr_") as scratch:
        inventory_path = Path(scratch) / "inventory.json"
        try:
            proc = subprocess.run([*importer, "--tdr", str(tdr_path), "--inventory-json", str(inventory_path)],
                                  capture_output=True, text=True, check=False, timeout=1800)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return False, f"importer failed to run: {exc}", {}
        if proc.returncode != 0:
            return False, f"importer exited {proc.returncode}: {proc.stderr.strip()[:300]}", {}
        try:
            inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            return True, f"importer succeeded but inventory unreadable: {exc}", {}
    summary = summarise_tdr_inventory(inventory)
    return True, f"readable by sentaurus_import: {summary.get('nodes', '?')} nodes, " \
                 f"{len(summary.get('datasets', []))} datasets", summary


# --------------------------------------------------------------------------- register

@dataclass
class StoreFile:
    path: str            # artifact path (relative to <case>/<version>)
    rel: str             # path relative to the store root
    size: int
    sha256: str
    mtime: Optional[float]
    local: Optional[Path] = None   # fetched copy for text inspection

    @property
    def name(self) -> str:
        return split_artifact_path(self.path)[2]

    @property
    def variant(self) -> Optional[str]:
        return split_artifact_path(self.path)[1]

    @property
    def category_dir(self) -> str:
        return split_artifact_path(self.path)[0]

    def text(self) -> Optional[str]:
        if self.local is None or not self.local.is_file():
            return None
        return decode_text(self.local.read_bytes())


@dataclass
class RegisterOptions:
    repo_source_dir: Optional[Path] = None
    origin: Optional[str] = None
    note: Optional[str] = None
    deck_matches: str = "auto"          # auto | true | false
    tdr_importer: Optional[Sequence[str]] = None
    warnings: list[str] = field(default_factory=list)


def _pick_deck(decks: Sequence[Deck], filename: str, content: str = "") -> Optional[Deck]:
    # Native console logs name the actual input explicitly. Shared mesh/par
    # echoes otherwise tie across every materialized branch.
    sources = re.findall(r"(?m)^\s*Input source:\s*(\S+)\s*$", content)
    if sources:
        named = [deck for deck in decks if deck.name in {s.replace('\\', '/').rsplit('/', 1)[-1] for s in sources}]
        return named[0] if len(named) == 1 else None
    best: list[Deck] = []
    best_score = (0, 0)
    for deck in decks:
        score = deck.match_score(filename, content)
        specificity = max([0] + [len(stem) for stem in deck.output_stems() + deck.current_prefixes
            if stem and (filename.startswith(stem + '_') or filename.startswith(stem + '.'))]) if score == 3 else 0
        ranked_score = (score, specificity)
        if ranked_score > best_score:
            best, best_score = [deck], ranked_score
        elif ranked_score == best_score and score > 0:
            best.append(deck)
    if len(best) == 1:
        return best[0]
    if not best and len(decks) == 1:
        return decks[0]
    return None


def _mesh_producer(name: str, content: str, decks_by_kind: dict[str, list[Deck]]) -> tuple[str, Optional[Deck]]:
    lower = name.lower()
    sde_decks, msh_decks, sprocess_decks = decks_by_kind["sde"], decks_by_kind["snmesh"], decks_by_kind["sprocess"]
    reflected = re.match(r"^n(\d+)_msh\.tdr$", name)
    if reflected:
        for deck in msh_decks:
            if deck.name == f"n{reflected.group(1)}_half_msh.cmd":
                return "other", deck
    if lower.endswith("_fps.tdr") or lower.startswith("sprocess") or lower == "run_sprocess.out" \
            or "Sentaurus Process" in content:
        return "sprocess", sprocess_decks[0] if len(sprocess_decks) == 1 else None
    if re.search(r"_bnd\.tdr$|sde.*\.(?:log|out)$|^sde\.std(?:out|err)$|_dvs\.log$|^run_sde\.out$", lower) \
            or "Sentaurus Structure Editor" in content:
        stem = re.sub(r"_(?:bnd|dvs|sde).*$", "", name)
        for deck in sde_decks:
            if stem and stem in deck.text:
                return "sde", deck
        return "sde", sde_decks[0] if len(sde_decks) == 1 else None
    if re.search(r"_msh\.(?:tdr|log)$", lower) or "Sentaurus Mesh" in content or lower == "run_tdx.out":
        stem = re.sub(r"\.(?:tdr|log)$", "", name)
        for deck in msh_decks:
            if deck.name == f"{stem}.cmd":
                return "snmesh", deck
        return ("other" if lower == "run_tdx.out" else "snmesh"), msh_decks[0] if len(msh_decks) == 1 else None
    return "other", None


def _producer_for(store_file: StoreFile, decks: Sequence[Deck], parameters: dict[str, str],
                  options: RegisterOptions) -> tuple[dict, Optional[Deck]]:
    producer: dict[str, Any] = {"tool": "other"}
    if options.origin:
        producer["run_directory"] = options.origin
    content = store_file.text() or ""
    by_kind: dict[str, list[Deck]] = {kind: [d for d in decks if d.kind == kind]
                                      for kind in ("sdevice", "sde", "snmesh", "sprocess")}
    tool, deck = "other", None
    category = store_file.category_dir
    is_sdevice_log = bool(_SDEVICE_BANNER_RE.search(content))
    if category in ("plt", "tdr") or (category == "log" and is_sdevice_log):
        deck = _pick_deck(by_kind["sdevice"], store_file.name, content)
        tool = "sdevice" if (deck is not None or category != "log") else "other"
        if deck is None and by_kind["sdevice"]:
            options.warnings.append(f"{store_file.path}: could not attribute to one sdevice command file")
    elif category in ("mesh", "log"):
        tool, deck = _mesh_producer(store_file.name, content, by_kind)
        if deck is None and category == "log" and tool == "other":
            deck = _pick_deck(by_kind["sdevice"], store_file.name, content) if by_kind["sdevice"] else None
            if deck is not None and deck.match_score(store_file.name, content) >= 3:
                tool = "sdevice"
            else:
                deck = None
    if deck is None and tool in ("sde", "snmesh", "sprocess", "sdevice"):
        options.warnings.append(f"{store_file.path}: {tool} product without an archived command file")
        tool = "other"
    producer["tool"] = tool
    if deck is not None:
        producer["cmd"] = deck.name
        producer["cmd_sha256"] = deck.lf_sha256
        if deck.kind == "sdevice" and deck.parameter_file:
            producer["parameter_file"] = deck.parameter_file
            sha = parameters.get(deck.parameter_file)
            if sha:
                producer["parameter_sha256"] = sha
            else:
                options.warnings.append(f"{store_file.path}: parameter file {deck.parameter_file} not archived")
    if store_file.mtime is not None:
        producer["finished_utc"] = utc_from_epoch(store_file.mtime)
    return producer, deck


def _source_tool(store_file: StoreFile, repo_source_dir: Optional[Path]) -> str:
    if repo_source_dir is None or store_file.local is None:
        return "other"
    repo_file = repo_source_dir / store_file.name
    if repo_file.is_file() and lf_sha256(repo_file.read_bytes()) == lf_sha256(store_file.local.read_bytes()):
        return "repository"
    return "other"


def register_case(case: str, version: str, backend: TransportBackend,
                  options: Optional[RegisterOptions] = None) -> dict:
    """Build a C7 manifest for ``<root>/<case>/<version>`` on ``backend`` (raises TransportError when unreachable)."""
    options = options or RegisterOptions()
    if not re.match(r"^[a-z0-9_]+$", case):
        raise ValueError(f"case name {case!r} must match ^[a-z0-9_]+$")
    if not re.match(r"^v[0-9]+$", version):
        raise ValueError(f"version {version!r} must match ^v[0-9]+$")
    if not backend.is_reachable():
        raise TransportError(f"store {backend.describe()} is unreachable")
    base = f"{case}/{version}"
    files: list[StoreFile] = []
    for directory in CATEGORY_DIRS:
        for rel in backend.list_files(f"{base}/{directory}"):
            size, digest = backend.stat_and_hash(rel)
            files.append(StoreFile(path=rel[len(base) + 1:], rel=rel, size=size, sha256=digest,
                                   mtime=backend.mtime(rel)))
    if not files:
        raise FileNotFoundError(f"no artifacts under {backend.describe()}/{base}/{{{','.join(CATEGORY_DIRS)}}}")
    repo_source_dir = options.repo_source_dir
    if repo_source_dir is None:
        default_dir = REPO / "reference_tcad" / case / "source"
        repo_source_dir = default_dir if default_dir.is_dir() else None

    with scratch_directory("vela_ref_register_") as scratch_name:
        scratch = Path(scratch_name)
        for store_file in files:
            suffix = Path(store_file.name).suffix.lower()
            wanted = suffix not in BINARY_SUFFIXES and store_file.size <= MAX_TEXT_FETCH_BYTES
            if suffix == ".tdr" and options.tdr_importer and store_file.size <= MAX_TDR_INSPECT_BYTES:
                wanted = True
            if not wanted:
                continue
            target = scratch.joinpath(*PurePosixPath(store_file.path).parts)
            backend.fetch(store_file.rel, target)
            if sha256_file(target) != store_file.sha256:
                raise TransportError(f"{store_file.rel}: fetched copy does not match the store hash")
            store_file.local = target

        decks_by_variant: dict[Optional[str], list[Deck]] = {}
        parameters_by_variant: dict[Optional[str], dict[str, str]] = {}
        for store_file in files:
            if store_file.category_dir != "source" or store_file.local is None:
                continue
            data = store_file.local.read_bytes()
            if store_file.name.lower().endswith(".cmd"):
                decks_by_variant.setdefault(store_file.variant, []).append(Deck(store_file.name, data))
            elif store_file.name.lower().endswith(".par"):
                parameters_by_variant.setdefault(store_file.variant, {})[store_file.name] = lf_sha256(data)

        versions: list[str] = []
        artifacts: list[dict] = []
        deck_matches_repository: Optional[bool] = None
        now = utc_now()
        for store_file in files:
            variant = store_file.variant
            decks = decks_by_variant.get(variant) or decks_by_variant.get(None, [])
            parameters = parameters_by_variant.get(variant) or parameters_by_variant.get(None, {})
            summary: dict[str, Any] = {}
            if store_file.category_dir == "source":
                producer: dict[str, Any] = {"tool": _source_tool(store_file, repo_source_dir)}
                if options.origin:
                    producer["run_directory"] = options.origin
            else:
                producer, deck = _producer_for(store_file, decks, parameters, options)
                text = store_file.text()
                suffix = Path(store_file.name).suffix.lower()
                if text is not None and producer["tool"] == "sdevice":
                    if suffix == ".plt":
                        try:
                            summary = summarise_plt(text, deck.goal_for_output(store_file.name) if deck else None)
                        except ValueError as exc:
                            options.warnings.append(f"{store_file.path}: unreadable .plt ({exc})")
                    elif _SDEVICE_BANNER_RE.search(text):
                        summary = {"exit_status": log_exit_status(text)}
                        banner = parse_version_banner(text)
                        if banner and banner not in versions:
                            versions.append(banner)
                if suffix == ".grd" and store_file.local is not None:
                    try:
                        summary = summarise_grd(store_file.local)
                    except (OSError, ValueError) as exc:
                        options.warnings.append(f"{store_file.path}: unreadable .grd ({exc})")
                if suffix == ".tdr" and store_file.local is not None and options.tdr_importer:
                    ok, detail, tdr_summary = inspect_tdr(options.tdr_importer, store_file.local)
                    if ok:
                        summary = tdr_summary
                    else:
                        options.warnings.append(f"{store_file.path}: {detail}")
            artifact = {
                "path": store_file.path,
                "category": CATEGORY_BY_DIR[store_file.category_dir],
                "sha256": store_file.sha256,
                "bytes": store_file.size,
                "producer": producer,
                "verify": {"status": "registered", "checked_utc": now,
                           "note": f"sha256 computed on the store copy at registration ({backend.describe()})"},
            }
            if summary:
                artifact["summary"] = summary
            artifacts.append(artifact)

        if options.deck_matches in ("true", "false"):
            deck_matches_repository = options.deck_matches == "true"
        else:
            sdevice_decks = [deck for decks in decks_by_variant.values() for deck in decks if deck.kind == "sdevice"]
            deck_matches_repository = bool(sdevice_decks) and repo_source_dir is not None and all(
                (repo_source_dir / deck.name).is_file()
                and decks_equivalent(deck.data, (repo_source_dir / deck.name).read_bytes())
                for deck in sdevice_decks)

    if len(versions) > 1:
        options.warnings.append(f"several Sentaurus versions in the logs: {versions}; recording {versions[0]}")
    manifest: dict[str, Any] = {
        "schema": SCHEMA_ID,
        "case": case,
        "version": version,
        "created_utc": now,
        "store": backend.store_descriptor(),
        "sentaurus_version": versions[0] if versions else "unknown",
        "deck_matches_repository": bool(deck_matches_repository),
        "artifacts": artifacts,
    }
    origin: dict[str, str] = {}
    if options.origin:
        origin["vm_run_directory"] = options.origin
    if options.note:
        origin["note"] = options.note
    if origin:
        manifest["origin"] = origin
    errors = sv.collect_errors(load_schema(), manifest)
    if errors:
        raise ValueError("registered manifest violates C7:\n" + "\n".join(errors))
    return manifest


# --------------------------------------------------------------------------- fetch

def fetch_case(manifest: dict, backend: TransportBackend, staging_root: Path, *, force: bool = False) -> list[dict]:
    """Copy every artifact to ``<staging_root>/<case>/<version>/<path>`` and re-hash it."""
    if not backend.is_reachable():
        raise TransportError(f"store {backend.describe()} is unreachable")
    base = f"{manifest['case']}/{manifest['version']}"
    local_base = staging_root / manifest["case"] / manifest["version"]
    results: list[dict] = []
    for artifact in manifest["artifacts"]:
        path = artifact["path"]
        local = local_base.joinpath(*PurePosixPath(path).parts)
        entry = {"path": path, "local": str(local), "expected": artifact["sha256"]}
        if local.is_file() and not force:
            actual = sha256_file(local)
            if actual == artifact["sha256"]:
                entry.update(status="present", actual=actual)
                results.append(entry)
                continue
        try:
            backend.fetch(f"{base}/{path}", local)
        except TransportError as exc:
            entry.update(status="error", detail=str(exc))
            results.append(entry)
            continue
        actual = sha256_file(local)
        entry.update(status="fetched" if actual == artifact["sha256"] else "mismatch", actual=actual)
        results.append(entry)
    return results


# --------------------------------------------------------------------------- verify

class Check:
    def __init__(self, number: int, name: str) -> None:
        self.number = number
        self.name = name
        self.items: list[dict] = []
        self.skipped: Optional[str] = None

    def add(self, subject: str, status: str, detail: str) -> None:
        self.items.append({"subject": subject, "status": status, "detail": detail})

    def skip(self, reason: str) -> None:
        self.skipped = reason

    @property
    def status(self) -> str:
        if self.skipped is not None:
            return "skip"
        if any(item["status"] == "fail" for item in self.items):
            return "fail"
        return "pass"

    def counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for item in self.items:
            counts[item["status"]] = counts.get(item["status"], 0) + 1
        return counts

    def to_json(self) -> dict:
        data = {"number": self.number, "name": self.name, "status": self.status, "counts": self.counts(),
                "items": self.items}
        if self.skipped is not None:
            data["skipped"] = self.skipped
        return data


class ArtifactAccess:
    """Bytes/paths of archived artifacts for checks 2-4: verified staged copy > verified store copy > unverified copy."""

    def __init__(self, manifest: dict, backend: Optional[TransportBackend], reachable: bool,
                 staged_base: Path, scratch: Path, staged_ok: dict[str, bool]) -> None:
        self.manifest = manifest
        self.backend = backend if reachable else None
        self.staged_base = staged_base
        self.scratch = scratch
        self.staged_ok = staged_ok
        self.expected = {a["path"]: a["sha256"] for a in manifest["artifacts"]}
        self.notes: dict[str, str] = {}
        self._paths: dict[str, Optional[Path]] = {}

    def local_path(self, path: str) -> Optional[Path]:
        if path in self._paths:
            return self._paths[path]
        result: Optional[Path] = None
        staged = self.staged_base.joinpath(*PurePosixPath(path).parts)
        if staged.is_file() and self.staged_ok.get(path):
            result = staged
        elif self.backend is not None and path in self.expected:
            target = self.scratch.joinpath(*PurePosixPath(path).parts)
            try:
                self.backend.fetch(f"{self.manifest['case']}/{self.manifest['version']}/{path}", target)
                if sha256_file(target) == self.expected[path]:
                    result = target
                else:
                    self.notes[path] = "store copy does not match the manifest hash"
            except TransportError as exc:
                self.notes[path] = f"store copy unavailable: {exc}"
        if result is None and staged.is_file():
            self.notes[path] = "using an unverified local copy"
            result = staged
        self._paths[path] = result
        return result

    def read(self, path: str) -> Optional[bytes]:
        local = self.local_path(path)
        return local.read_bytes() if local is not None else None


def find_source_artifact(manifest: dict, variant: Optional[str], name: str) -> Optional[dict]:
    by_path = {a["path"]: a for a in manifest["artifacts"] if a["path"].startswith("source/")}
    if variant and f"source/{variant}/{name}" in by_path:
        return by_path[f"source/{variant}/{name}"]
    if f"source/{name}" in by_path:
        return by_path[f"source/{name}"]
    matches = [a for p, a in by_path.items() if p.rsplit("/", 1)[-1] == name]
    return matches[0] if len(matches) == 1 else None


def _mode_label(mode: str) -> str:
    return {"online": "online (compared with the store)",
            "local_copy": "local copy (store unreachable; not compared with source)",
            "offline": "offline (manifest only; checks 1, 3 and 4 skipped)"}[mode]


def _check_hashes(check: Check, manifest: dict, backend: Optional[TransportBackend], reachable: bool,
                  staged_base: Path, staged_ok: dict[str, bool]) -> None:
    base = f"{manifest['case']}/{manifest['version']}"
    for artifact in manifest["artifacts"]:
        path = artifact["path"]
        staged = staged_base.joinpath(*PurePosixPath(path).parts)
        staged_digest = sha256_file(staged) if staged.is_file() else None
        staged_ok[path] = staged_digest == artifact["sha256"]
        if reachable and backend is not None:
            try:
                size, digest = backend.stat_and_hash(f"{base}/{path}")
            except TransportError as exc:
                check.add(path, "fail", f"missing on store: {exc}")
                continue
            if digest != artifact["sha256"]:
                check.add(path, "fail", f"store sha256 {digest[:16]}... differs from manifest {artifact['sha256'][:16]}...")
            elif size != artifact["bytes"]:
                check.add(path, "fail", f"store size {size} differs from manifest {artifact['bytes']}")
            elif staged_digest is not None and not staged_ok[path]:
                check.add(path, "fail", f"staged copy {staged} differs from the manifest hash (store copy intact)")
            else:
                check.add(path, "pass", "store hash matches" + ("; staged copy matches" if staged_digest else ""))
        else:
            if staged_digest is None:
                check.add(path, "fail", "no local copy and the store is unreachable")
            elif staged_ok[path]:
                check.add(path, "pass", "local copy hash matches (not compared with source)")
            else:
                check.add(path, "fail", f"local copy {staged} sha256 {staged_digest[:16]}... differs from manifest")


def _check_sources(check: Check, manifest: dict, access: Optional[ArtifactAccess], repo_source_dir: Path) -> None:
    deck_matches = bool(manifest.get("deck_matches_repository"))
    seen: set[tuple] = set()
    for artifact in manifest["artifacts"]:
        producer = artifact["producer"]
        cmd = producer.get("cmd")
        tool = producer["tool"]
        if not cmd and tool == "sdevice":
            check.add(artifact["path"], "fail", "sdevice artifact without producer.cmd")
            continue
        if not cmd:
            continue
        _, variant, _ = split_artifact_path(artifact["path"])
        key = (variant, cmd, producer.get("cmd_sha256"), producer.get("parameter_file"), producer.get("parameter_sha256"), tool)
        if key in seen:
            continue
        seen.add(key)
        subject = f"source/{variant + '/' if variant else ''}{cmd}"
        source = find_source_artifact(manifest, variant, cmd)
        if source is None:
            check.add(subject, "fail", f"producer cmd of {artifact['path']} is not archived under source/")
            continue
        expected = producer.get("cmd_sha256")
        if tool == "sdevice" and not expected:
            check.add(subject, "fail", "sdevice producer without cmd_sha256")
        data = access.read(source["path"]) if access is not None else None
        if expected:
            if data is not None:
                actual = lf_sha256(data)
                if actual == expected:
                    check.add(subject, "pass", "LF-normalised sha256 of the archived cmd equals producer.cmd_sha256")
                else:
                    check.add(subject, "fail", f"LF-normalised sha256 {actual[:16]}... of the archived cmd differs "
                                               f"from producer.cmd_sha256 {expected[:16]}...")
            elif source["sha256"] == expected:
                check.add(subject, "pass", "archived raw sha256 equals producer.cmd_sha256 (LF file)")
            elif access is None:
                check.add(subject, "info", "cmd_sha256 is LF-normalised and the archived copy has a different raw "
                                           "hash (CRLF/BOM); comparison needs the file (run without --offline)")
            else:
                check.add(subject, "skip", f"archived cmd unavailable: {access.notes.get(source['path'], 'no copy')}")
        if tool == "sdevice" and deck_matches:
            repo_file = repo_source_dir / cmd
            if not repo_file.is_file():
                check.add(subject, "fail", f"deck_matches_repository is true but {repo_file} does not exist")
            elif data is not None:
                if decks_equivalent(data, repo_file.read_bytes()):
                    check.add(subject, "pass", f"archived cmd equals repository {repo_file.name} "
                                               "after comment stripping and whitespace collapsing")
                else:
                    check.add(subject, "fail", f"archived cmd differs from repository {repo_file.name} "
                                               "after comment stripping and whitespace collapsing")
            else:
                check.add(subject, "info", f"repository deck {repo_file.name} present; content comparison deferred")
        parameter = producer.get("parameter_file")
        if parameter:
            par_subject = f"source/{variant + '/' if variant else ''}{parameter}"
            par_source = find_source_artifact(manifest, variant, parameter)
            if par_source is None:
                check.add(par_subject, "fail", f"parameter file of {cmd} is not archived under source/")
                continue
            par_expected = producer.get("parameter_sha256")
            par_data = access.read(par_source["path"]) if access is not None else None
            if par_expected:
                if par_data is not None:
                    actual = lf_sha256(par_data)
                    if actual == par_expected:
                        check.add(par_subject, "pass", "LF-normalised sha256 equals producer.parameter_sha256")
                    else:
                        check.add(par_subject, "fail", f"LF-normalised sha256 {actual[:16]}... differs from "
                                                       f"producer.parameter_sha256 {par_expected[:16]}...")
                elif par_source["sha256"] == par_expected:
                    check.add(par_subject, "pass", "archived raw sha256 equals producer.parameter_sha256 (LF file)")
                elif access is None:
                    check.add(par_subject, "info", "parameter_sha256 comparison needs the file (run without --offline)")
                else:
                    check.add(par_subject, "skip", "archived parameter file unavailable")
            if tool == "sdevice" and deck_matches:
                repo_par = repo_source_dir / parameter
                if not repo_par.is_file():
                    check.add(par_subject, "warn", f"repository has no {parameter}; parameter comparison skipped")
                elif par_data is not None:
                    if decks_equivalent(par_data, repo_par.read_bytes()):
                        check.add(par_subject, "pass", f"archived parameter file equals repository {parameter} "
                                                       "after normalisation")
                    else:
                        check.add(par_subject, "fail", f"archived parameter file differs from repository {parameter} "
                                                       "after normalisation")
    if not check.items:
        check.add("(manifest)", "info", "no artifacts declare a producer cmd")


def _check_structure(check: Check, manifest: dict, access: ArtifactAccess) -> None:
    deck_cache: dict[str, Optional[Deck]] = {}

    def deck_for(artifact: dict) -> Optional[Deck]:
        producer = artifact["producer"]
        if producer.get("tool") != "sdevice" or not producer.get("cmd"):
            return None
        _, variant, _ = split_artifact_path(artifact["path"])
        source = find_source_artifact(manifest, variant, producer["cmd"])
        if source is None:
            return None
        if source["path"] not in deck_cache:
            data = access.read(source["path"])
            deck_cache[source["path"]] = Deck(producer["cmd"], data) if data is not None else None
        return deck_cache[source["path"]]

    grd_contacts: dict[Optional[str], Optional[set[str]]] = {}
    for artifact in manifest["artifacts"]:
        _, variant, name = split_artifact_path(artifact["path"])
        if artifact["category"] == "C" and name.lower().endswith(".grd"):
            local = access.local_path(artifact["path"])
            if local is not None:
                try:
                    grd_contacts[variant] = set(dfise.parse_grd(local)["contact_element_counts"])
                except (OSError, ValueError):
                    grd_contacts[variant] = None

    plt_seen = 0
    for artifact in manifest["artifacts"]:
        path = artifact["path"]
        _, variant, name = split_artifact_path(path)
        if artifact["category"] != "B" or not name.lower().endswith(".plt"):
            continue
        plt_seen += 1
        deck = deck_for(artifact)
        if deck is None:
            check.add(path, "skip", "producer cmd unavailable; electrode/goal structure not checked")
            continue
        data = access.read(path)
        if data is None:
            check.add(path, "skip", f"current file unavailable: {access.notes.get(path, 'no copy')}")
            continue
        text = decode_text(data)
        try:
            electrodes, rows = read_plt_document(text)
        except ValueError as exc:
            check.add(path, "fail", f"unreadable .plt: {exc}")
            continue
        datasets = set(plt_dataset_names(text))
        missing = [e for e in deck.electrodes if f"{e} OuterVoltage" not in datasets]
        if not deck.electrodes:
            check.add(path, "fail", f"no Electrode{{Name=}} found in {deck.name}")
        elif missing:
            check.add(path, "fail", f"electrodes {missing} of {deck.name} have no '<name> OuterVoltage' dataset")
        else:
            check.add(path, "pass", f"all {len(deck.electrodes)} electrodes of {deck.name} have OuterVoltage datasets")
        goal = deck.goal_for_output(name)
        if goal is None:
            check.add(path, "info", f"{deck.name} has no Quasistationary Goal for this current file (equilibrium only)")
        elif goal.name not in electrodes or not rows:
            check.add(path, "fail", f"Goal electrode {goal.name!r} has no rows in the current file")
        else:
            first = rows[0][goal.name].get("outer_voltage_V", 0.0)
            last = rows[-1][goal.name].get("outer_voltage_V", 0.0)
            if goal.quantity == "current":
                current = rows[-1][goal.name]["total_current_A_per_um"]
                if current_goal_reached(current, goal.voltage):
                    check.add(path, "pass", f"{goal.name} current {current:g} A/um reaches Goal {goal.voltage:g} A/um")
                else:
                    check.add(path, "partial", f"partial reference: {goal.name} current {current:g} A/um, Goal {goal.voltage:g} A/um not reached")
                for earlier in deck.goals:
                    if earlier.prefix == goal.prefix and earlier.quantity == "voltage":
                        reached = any(earlier.name in row and goal_reached(row[earlier.name]["outer_voltage_V"], earlier.voltage) for row in rows)
                        check.add(path, "pass" if reached else "partial",
                                  f"preceding {earlier.name} Voltage Goal {earlier.voltage:g} V {'reached' if reached else 'not reached'} in the shared current file")
            elif goal_reached(last, goal.voltage):
                check.add(path, "pass", f"{goal.name} bias {first:g} V -> {last:g} V reaches Goal {goal.voltage:g} V")
            else:
                logs = []
                for log_artifact in manifest["artifacts"]:
                    producer = log_artifact["producer"]
                    if log_artifact["category"] == "E" and producer.get("tool") == "sdevice" and producer.get("cmd") == deck.name:
                        data = access.read(log_artifact["path"])
                        if data is not None:
                            logs.append(decode_text(data))
                reason = declared_break_evidence(deck, rows, logs)
                if reason:
                    check.add(path, "pass", f"{goal.name} bias {first:g} V -> {last:g} V; {reason}")
                else:
                    check.add(path, "partial", f"partial reference: {goal.name} bias {first:g} V -> {last:g} V, "
                                               f"Goal {goal.voltage:g} V not reached and no declared stop is proven")
        contacts = grd_contacts.get(variant, "absent")
        if contacts == "absent":
            check.add(path, "info", "mesh contact names not inspected (TDR mesh is HDF5)")
        elif contacts is None:
            check.add(path, "fail", "DF-ISE .grd mesh of this variant is unreadable")
        else:
            absent = [e for e in deck.electrodes if e not in contacts]
            if absent:
                check.add(path, "fail", f"electrodes {absent} are not contact regions of the .grd mesh")
            else:
                check.add(path, "pass", "mesh contact names cover the cmd electrodes")
    if not plt_seen:
        check.add("(manifest)", "info", "no .plt artifacts registered")


def _check_completeness(check: Check, manifest: dict, access: ArtifactAccess,
                        importer: Optional[Sequence[str]]) -> None:
    for artifact in manifest["artifacts"]:
        path = artifact["path"]
        _, _, name = split_artifact_path(path)
        suffix = Path(name).suffix.lower()
        tool = artifact["producer"].get("tool")
        if suffix == ".plt":
            data = access.read(path)
            if data is None:
                check.add(path, "skip", f"unavailable: {access.notes.get(path, 'no copy')}")
                continue
            try:
                electrodes, rows = read_plt_document(decode_text(data))
                check.add(path, "pass", f"DF-ISE xyplot parsed: {len(rows)} rows, {len(electrodes)} electrodes")
            except ValueError as exc:
                check.add(path, "fail", f"unreadable .plt: {exc}")
        elif suffix in (".grd", ".dat"):
            local = access.local_path(path)
            if local is None:
                check.add(path, "skip", f"unavailable: {access.notes.get(path, 'no copy')}")
                continue
            try:
                if suffix == ".grd":
                    grid = dfise.parse_grd(local)
                    if not grid.get("vertex_count"):
                        raise ValueError("no nb_vertices header")
                    check.add(path, "pass", f"DF-ISE grid parsed: {grid['vertex_count']} vertices, "
                                            f"{grid.get('element_count')} elements")
                else:
                    datasets = dfise.parse_dat(local)
                    if not datasets:
                        raise ValueError("no Dataset blocks")
                    check.add(path, "pass", f"DF-ISE data parsed: {len(datasets)} datasets")
            except (OSError, ValueError) as exc:
                check.add(path, "fail", f"unreadable DF-ISE file: {exc}")
        elif suffix == ".tdr":
            if tool == "sde" and artifact["category"] == "C" and name.lower().endswith("_bnd.tdr"):
                check.add(path, "info", "not inspected (SDE boundary input): volume-mesh importer requires Tri3 cells; hash/source checks still apply")
            elif importer is None:
                check.add(path, "info", "not inspected (HDF5): sentaurus_import importer not available")
            elif artifact["bytes"] > MAX_TDR_INSPECT_BYTES:
                check.add(path, "info", f"not inspected (HDF5): {artifact['bytes']} bytes exceeds the inspection limit")
            else:
                local = access.local_path(path)
                if local is None:
                    check.add(path, "skip", f"unavailable: {access.notes.get(path, 'no copy')}")
                    continue
                ok, detail, _ = inspect_tdr(importer, local)
                check.add(path, "pass" if ok else "fail", detail)
        elif suffix == ".sav":
            check.add(path, "info", "not inspected (sdevice save state)")
        elif tool == "sdevice" and artifact["category"] == "E":
            data = access.read(path)
            if data is None:
                check.add(path, "skip", f"unavailable: {access.notes.get(path, 'no copy')}")
                continue
            status = log_exit_status(decode_text(data))
            recorded = (artifact.get("summary") or {}).get("exit_status")
            if status.startswith("failed:"):
                check.add(path, "fail", status)
            elif status != "Good Bye":
                check.add(path, "warn", "log has no Good Bye marker")
            elif recorded and recorded != status:
                check.add(path, "warn", f"log ends with Good Bye but the manifest records {recorded!r}")
            else:
                check.add(path, "pass", "log ends with the Good Bye marker")
        else:
            check.add(path, "info", "no structural check for this file type")


def _check_version(check: Check, manifest: dict, access: Optional[ArtifactAccess]) -> None:
    version = str(manifest.get("sentaurus_version", ""))
    if version.startswith(EXPECTED_VERSION_PREFIX):
        check.add("sentaurus_version", "pass", f"{version} is a {EXPECTED_VERSION_PREFIX} release")
    else:
        declared = re.search(r"sentaurus(\d{4})", manifest.get("case", ""))
        banner_year = re.search(r"[A-Z]-(\d{4})\.", version)
        if declared and banner_year and int(banner_year.group(1)) <= int(declared.group(1)) < 2022:
            check.add("sentaurus_version", "pass",
                      f"{version}: syntax-compatible case (declares Sentaurus {declared.group(1)})")
        else:
            check.add("sentaurus_version", "fail",
                      f"{version!r} is neither {EXPECTED_VERSION_PREFIX} nor the version declared by the case name")
    if access is None:
        return
    for artifact in manifest["artifacts"]:
        if artifact["producer"].get("tool") != "sdevice" or artifact["category"] != "E":
            continue
        data = access.read(artifact["path"])
        if data is None:
            continue
        text = decode_text(data)
        if not _SDEVICE_BANNER_RE.search(text):
            continue
        banner = parse_version_banner(text)
        if banner is None:
            check.add(artifact["path"], "warn", "no Sentaurus version banner found in the log header")
        elif version.startswith(banner) or banner.startswith(version):
            check.add(artifact["path"], "pass", f"log header version {banner} matches the manifest")
        else:
            check.add(artifact["path"], "fail", f"log header version {banner} differs from manifest {version}")


def render_markdown(report: dict) -> str:
    lines = [f"# Reference artifact verification: `{report['case']}` {report['version']}", "",
             f"- Checked (UTC): {report['checked_utc']}",
             f"- Manifest: `{report['manifest_path']}`",
             f"- Store: `{report['store_label']}` (reachable: {report['store_reachable']})",
             f"- Mode: {_mode_label(report['mode'])}",
             f"- Compared with source: {report['compared_with_source']}",
             f"- Status: **{report['status'].upper().replace('_', ' ')}** (exit code {report['exit_code']})", "",
             "| # | Check | Status | pass | fail | partial | warn | info | skip |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for check in report["checks"]:
        counts = check["counts"]
        lines.append(f"| {check['number']} | {check['name']} | {check['status']} | "
                     + " | ".join(str(counts.get(key, 0)) for key in ("pass", "fail", "partial", "warn", "info", "skip"))
                     + " |")
    for title, key in (("Failures", "failures"), ("Partial references", "partial_references"),
                       ("Warnings", "warnings"), ("Notes", "notes")):
        entries = report.get(key) or []
        if entries:
            lines.extend(["", f"## {title}", ""])
            lines.extend(f"- {entry}" for entry in entries)
    return "\n".join(lines) + "\n"


def verify_case(manifest: dict, *, manifest_path: Optional[Path], backend: Optional[TransportBackend],
                staging_root: Path, report_dir: Path, offline: bool = False,
                repo_source_dir: Optional[Path] = None, tdr_importer: Optional[Sequence[str]] = None) -> dict:
    """Run the six checks of plan section 12.4; returns the report (``exit_code`` inside)."""
    case, version = manifest["case"], manifest["version"]
    repo_source_dir = repo_source_dir or (REPO / "reference_tcad" / case / "source")
    staged_base = staging_root / case / version
    checks = [Check(1, "existence and sha256"), Check(2, "source consistency"), Check(3, "structure"),
              Check(4, "completeness"), Check(5, "Sentaurus version"), Check(6, "report")]
    notes: list[str] = []
    schema_errors = sv.collect_errors(load_schema(), manifest)
    if schema_errors:
        notes.append("manifest violates C7: " + "; ".join(schema_errors[:5]))
    store = manifest.get("store", {})
    store_label = backend.describe() if backend else f"{store.get('kind')}:{store.get('host', '')}:{store.get('root', '')}"

    reachable: Optional[bool] = None
    compared_with_source = False
    status = "pass"
    if offline:
        mode = "offline"
        for check in checks[:1] + checks[2:4]:
            check.skip("offline mode")
    else:
        reachable = bool(backend and backend.is_reachable())
        if backend is None:
            notes.append(f"store kind {store.get('kind')!r} is not supported by this tool; using local copies only")
        local_copies = [a for a in manifest["artifacts"]
                        if staged_base.joinpath(*PurePosixPath(a["path"]).parts).is_file()]
        if not reachable and not local_copies:
            mode = "local_copy"
            status = "unable_to_verify"
            notes.append(f"store {store_label} is unreachable and {staged_base} holds no copy of any artifact")
            for check in checks[:5]:
                check.skip("unable to verify: store unreachable and no local copy")
        else:
            mode = "online" if reachable else "local_copy"
            compared_with_source = reachable
            if not reachable:
                notes.append("store unreachable: hashes checked against the local staging copies only; "
                             "not compared with source")

    staged_ok: dict[str, bool] = {}
    with scratch_directory("vela_ref_verify_") as scratch_name:
        scratch = Path(scratch_name)
        access: Optional[ArtifactAccess] = None
        if status != "unable_to_verify":
            if not offline:
                _check_hashes(checks[0], manifest, backend, bool(reachable), staged_base, staged_ok)
                access = ArtifactAccess(manifest, backend, bool(reachable), staged_base, scratch, staged_ok)
            _check_sources(checks[1], manifest, access, repo_source_dir)
            if access is not None:
                _check_structure(checks[2], manifest, access)
                _check_completeness(checks[3], manifest, access, tdr_importer)
            _check_version(checks[4], manifest, access)
            if access is not None:
                for path, note in sorted(access.notes.items()):
                    notes.append(f"{path}: {note}")
        if status != "unable_to_verify":
            status = "fail" if any(check.status == "fail" for check in checks[:5]) or schema_errors else "pass"
    exit_code = {"pass": EXIT_OK, "fail": EXIT_FAIL, "unable_to_verify": EXIT_UNAVAILABLE}[status]
    checked = utc_now()
    report_json = report_dir / "verify_report.json"
    report_md = report_dir / f"verify_{case}_{checked[:10]}.md"
    checks[5].add(str(report_json), "pass", "verify_report.json written")
    checks[5].add(str(report_md), "pass", "Markdown summary written")
    report: dict[str, Any] = {
        "schema": REPORT_SCHEMA_ID,
        "case": case,
        "version": version,
        "manifest_path": str(manifest_path) if manifest_path else "",
        "checked_utc": checked,
        "mode": mode,
        "store": store,
        "store_label": store_label,
        "store_reachable": reachable,
        "staging_dir": str(staged_base),
        "compared_with_source": compared_with_source,
        "status": status,
        "exit_code": exit_code,
        "checks": [check.to_json() for check in checks],
        "failures": [f"check {c.number} ({c.name}): {i['subject']}: {i['detail']}"
                     for c in checks for i in c.items if i["status"] == "fail"],
        "partial_references": [f"check {c.number}: {i['subject']}: {i['detail']}"
                               for c in checks for i in c.items if i["status"] == "partial"],
        "warnings": [f"check {c.number}: {i['subject']}: {i['detail']}"
                     for c in checks for i in c.items if i["status"] == "warn"],
        "notes": notes,
        "report_json": str(report_json),
        "report_markdown": str(report_md),
    }
    try:
        write_json(report_json, report)
        report_md.write_text(render_markdown(report), encoding="utf-8", newline="\n")
    except OSError as exc:
        checks[5].items.clear()
        checks[5].add(str(report_dir), "fail", f"could not write the report: {exc}")
        report["checks"] = [check.to_json() for check in checks]
        report["failures"].append(f"check 6 (report): {report_dir}: could not write the report: {exc}")
        if report["status"] == "pass":
            report["status"], report["exit_code"] = "fail", EXIT_FAIL
    return report


# --------------------------------------------------------------------------- CLI

def _resolve_manifest(args: argparse.Namespace) -> tuple[dict, Path]:
    path = Path(args.manifest) if args.manifest else manifest_path_for(args.case)
    if not path.is_file():
        raise FileNotFoundError(f"manifest not found: {path}")
    manifest = load_manifest(path)
    if manifest.get("case") != args.case:
        raise ValueError(f"{path}: manifest case {manifest.get('case')!r} does not match --case {args.case!r}")
    return manifest, path


def _staging_root(args: argparse.Namespace) -> Path:
    return Path(args.dest) if args.dest else default_staging_root()


def cmd_register(args: argparse.Namespace) -> int:
    if args.backend:
        backend = make_backend(args.backend)
    else:
        backend = SshBackend(args.host, args.store_root or DEFAULT_STORE_ROOT)
    assert backend is not None
    options = RegisterOptions(
        repo_source_dir=Path(args.repo_source_dir) if args.repo_source_dir else None,
        origin=args.origin, note=args.note, deck_matches=args.deck_matches_repository,
        tdr_importer=resolve_tdr_importer(args.tdr_importer))
    try:
        manifest = register_case(args.case, args.version, backend, options)
    except TransportError as exc:
        print(f"unable to register: {exc}", file=sys.stderr)
        return EXIT_UNAVAILABLE
    except (FileNotFoundError, ValueError) as exc:
        print(f"register failed: {exc}", file=sys.stderr)
        return EXIT_FAIL
    out = Path(args.manifest_out) if args.manifest_out else manifest_path_for(args.case)
    write_json(out, manifest)
    for warning in options.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    tools: dict[str, int] = {}
    for artifact in manifest["artifacts"]:
        tools[artifact["producer"]["tool"]] = tools.get(artifact["producer"]["tool"], 0) + 1
    print(f"registered {len(manifest['artifacts'])} artifacts of {args.case}/{args.version} from {backend.describe()} "
          f"-> {out} (sentaurus_version={manifest['sentaurus_version']}, "
          f"deck_matches_repository={manifest['deck_matches_repository']}, tools={tools})")
    return EXIT_OK


def cmd_fetch(args: argparse.Namespace) -> int:
    try:
        manifest, _ = _resolve_manifest(args)
        backend = make_backend(args.backend, manifest.get("store"))
        if backend is None:
            print(f"fetch failed: store kind {manifest.get('store', {}).get('kind')!r} is not supported", file=sys.stderr)
            return EXIT_FAIL
        results = fetch_case(manifest, backend, _staging_root(args), force=args.force)
    except TransportError as exc:
        print(f"unable to fetch: {exc}", file=sys.stderr)
        return EXIT_UNAVAILABLE
    except (FileNotFoundError, ValueError) as exc:
        print(f"fetch failed: {exc}", file=sys.stderr)
        return EXIT_FAIL
    failed = [r for r in results if r["status"] in ("mismatch", "error")]
    counts: dict[str, int] = {}
    for result in results:
        counts[result["status"]] = counts.get(result["status"], 0) + 1
    for result in failed:
        detail = result.get("detail") or f"sha256 {result.get('actual', '')[:16]}... != manifest {result['expected'][:16]}..."
        print(f"{result['status']}: {result['path']} -> {result['local']}: {detail}", file=sys.stderr)
    print(f"fetch {manifest['case']}/{manifest['version']} from {backend.describe()}: {counts}")
    return EXIT_FAIL if failed else EXIT_OK


def cmd_verify(args: argparse.Namespace) -> int:
    try:
        manifest, manifest_path = _resolve_manifest(args)
        backend = None if args.offline else make_backend(args.backend, manifest.get("store"))
    except (FileNotFoundError, ValueError) as exc:
        print(f"verify failed: {exc}", file=sys.stderr)
        return EXIT_FAIL
    report = verify_case(
        manifest, manifest_path=manifest_path, backend=backend, staging_root=_staging_root(args),
        report_dir=Path(args.report_dir) if args.report_dir else DEFAULT_REPORT_DIR, offline=args.offline,
        repo_source_dir=Path(args.repo_source_dir) if args.repo_source_dir else None,
        tdr_importer=resolve_tdr_importer(args.tdr_importer))
    for failure in report["failures"]:
        print(f"FAIL {failure}", file=sys.stderr)
    for partial in report["partial_references"]:
        print(f"PARTIAL {partial}")
    for note in report["notes"]:
        print(f"note: {note}")
    summary = ", ".join(f"{c['number']}:{c['status']}" for c in report["checks"])
    print(f"verify {manifest['case']}/{manifest['version']} [{report['mode']}]: {report['status']} ({summary}); "
          f"report {report['report_markdown']}")
    return int(report["exit_code"])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--case", required=True, help="reference_tcad case directory name")
        p.add_argument("--backend", help="transport override: local:<root> or ssh:<host>:<root>")
        p.add_argument("--repo-source-dir", help="directory of the repository decks (default reference_tcad/<case>/source)")
        p.add_argument("--tdr-importer", help="sentaurus_import executable for .tdr inspection; 'none' disables "
                                              "(default: build/sentaurus_import when present)")

    register = sub.add_parser("register", help="list, hash and describe a store directory into a C7 manifest")
    common(register)
    register.add_argument("--version", default="v1", help="archive version sub-directory (default v1)")
    register.add_argument("--store-root", help=f"ssh store root (default {DEFAULT_STORE_ROOT})")
    register.add_argument("--host", default=DEFAULT_HOST, help=f"ssh alias (default {DEFAULT_HOST})")
    register.add_argument("--manifest-out", help="output path (default reference_tcad/<case>/artifact_manifest.json)")
    register.add_argument("--origin", help="provenance run directory recorded as producer.run_directory")
    register.add_argument("--note", help="free-text origin note")
    register.add_argument("--deck-matches-repository", choices=("auto", "true", "false"), default="auto")
    register.set_defaults(func=cmd_register)

    fetch = sub.add_parser("fetch", help="copy the manifest artifacts to the staging directory and re-hash them")
    common(fetch)
    fetch.add_argument("--manifest", help="manifest path override")
    fetch.add_argument("--dest", help=f"staging root (default ${STAGING_ENV} or reference_staging/)")
    fetch.add_argument("--force", action="store_true", help="re-download files that already match")
    fetch.set_defaults(func=cmd_fetch)

    verify = sub.add_parser("verify", help="run the six checks of plan section 12.4")
    common(verify)
    verify.add_argument("--manifest", help="manifest path override")
    verify.add_argument("--dest", help=f"staging root (default ${STAGING_ENV} or reference_staging/)")
    verify.add_argument("--report-dir", help="where verify_report.json and the Markdown summary go (default docs/validation)")
    verify.add_argument("--offline", action="store_true", help="manifest-only checks 2 and 5; no store or file access")
    verify.set_defaults(func=cmd_verify)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args))
    except ValueError as exc:
        print(f"{args.command} failed: {exc}", file=sys.stderr)
        return EXIT_FAIL


if __name__ == "__main__":
    raise SystemExit(main())
