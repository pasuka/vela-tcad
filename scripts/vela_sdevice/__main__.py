"""vela device command-line entry point (stdlib only)."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import sys
from .materials import default_models_par_path
from .orchestrator import DeviceRun, version_info
from .orchestrator.support import RunError


class Parser(argparse.ArgumentParser):
    def error(self, message):
        self.exit(2, f"E202: Option is not supported by vela device: {message}\n")


def resolve_cmd(name):
    path = Path(name)
    for candidate in (path, Path(str(path) + ".cmd"), Path(str(path) + "_des.cmd")):
        if candidate.is_file():
            return candidate.resolve()
    raise RunError("E305", f"Command file '{name}' not found")


def main(argv=None):
    # Windows' default code page also encodes localized OS error details.
    # Keep captured CLI output in the same UTF-8 encoding as the run log.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="backslashreplace")
    parser = Parser(prog="vela", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    device = commands.add_parser("device", help="Run a preprocessed SDevice command file")
    device.add_argument("cmd", nargs="?")
    device.add_argument("-v", action="store_true", dest="version")
    device.add_argument("-P", action="store_true", dest="parameters")
    device.add_argument("-q", action="store_true", dest="quiet")
    device.add_argument("-n", "--check", action="store_true", dest="check")
    device.add_argument("-i", action="store_true", dest="initial_only")
    device.add_argument("--exit-on-failure", action="store_true")
    device.add_argument("--threads", "--max_threads", type=int)
    device.add_argument("--verbose", action="store_true")
    device.add_argument("--allow-lossy", action="store_true")
    keep = device.add_mutually_exclusive_group()
    keep.add_argument("--keep-intermediate", action="store_true", dest="keep_intermediate", default=True)
    keep.add_argument("--no-keep-intermediate", action="store_false", dest="keep_intermediate")
    device.add_argument("--vtu-binary", action="store_true")
    device.add_argument("--plot-final-only", action="store_true")
    options = parser.parse_args(argv)
    if options.version:
        version = version_info()
        print(f"Vela {version['version']} ({version['commit']}); SDevice syntax T-2022.03")
        return 0
    if options.parameters:
        try:
            shutil.copyfile(default_models_par_path(), Path.cwd() / "models.par")
        except OSError as error:
            print(f"E601: Cannot write output file 'models.par': {error}", file=sys.stderr)
            return 6
        return 0
    if not options.cmd:
        parser.error("No command file given (reading from stdin is not supported)")
    if options.threads is not None and options.threads < 1:
        parser.error("--threads must be positive")
    try:
        return DeviceRun(resolve_cmd(options.cmd), options).run()
    except RunError as error:
        print(str(error), file=sys.stderr)
        return error.exit_code
    except OSError as error:
        print(f"E601: Cannot read/write frontend files: {error}", file=sys.stderr)
        return 6


if __name__ == "__main__":
    raise SystemExit(main())
