#!/usr/bin/env python3
"""Prepare an isolated PN2D HFS diagnostic; never patch production sources.

This is an investigation candidate, not a qualified solver profile. The
diagnostic legacy constants, material oracle and tangent choices must be
reported separately from production validation.
"""
import argparse
from pathlib import Path
import shutil
import subprocess


def prepare(output):
    root = Path(__file__).resolve().parents[1]
    output = (root / output).resolve()
    if not output.is_relative_to(root / "build-release") or output.exists():
        raise ValueError("a fresh directory under this worktree's build-release is required")
    cpp = root / "src/equation/CoupledDDAssembler.cpp"
    constants = (root / "include/vela/core/PhysicalConstants.h").read_text(encoding="utf-8")
    for old, new in [("1.602176634e-19", "1.602192e-19"), ("1.380649e-23", "1.380662e-23"),
                     ("8.8541878128e-12", "8.8542e-12"), ("6.62607015e-34", "6.626176e-34"),
                     ("9.1093837015e-31", "9.109534e-31")]:
        if old not in constants:
            raise ValueError("production constants changed; review the diagnostic overlay")
        constants = constants.replace(old, new)
    source = output / "src/equation/CoupledDDAssembler.cpp"
    source.parent.mkdir(parents=True)
    shutil.copy2(cpp, source)
    patch = root / "scripts/diagnostics/pn2d_element_hfs_candidate.patch"
    command = ["git", "-c", "core.fsmonitor=false", "apply",
               "--directory=" + output.relative_to(root).as_posix()]
    subprocess.run(command + ["--check", str(patch)], cwd=root, check=True)
    subprocess.run(command + [str(patch)], cwd=root, check=True)
    include = output / "include/vela/core/PhysicalConstants.h"
    include.parent.mkdir(parents=True)
    include.write_text("// Diagnostic legacy constants only.\n" + constants, encoding="utf-8")
    probe = root / "scripts/diagnostics/pn2d_zero_hfs_probe.cpp"
    cmake = f'''cmake_minimum_required(VERSION 3.21)
project(pn2d_hfs_diagnostic LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
add_subdirectory("{root.as_posix()}" vela)
target_include_directories(vela_core BEFORE PUBLIC "${{CMAKE_CURRENT_SOURCE_DIR}}/include")
get_target_property(pn2d_sources vela_core SOURCES)
list(REMOVE_ITEM pn2d_sources src/equation/CoupledDDAssembler.cpp)
list(APPEND pn2d_sources "{source.as_posix()}")
set_property(TARGET vela_core PROPERTY SOURCES "${{pn2d_sources}}")
add_executable(pn2d_zero_hfs_probe "{probe.as_posix()}")
target_link_libraries(pn2d_zero_hfs_probe PRIVATE vela_core)
'''
    (output / "CMakeLists.txt").write_text(cmake, encoding="utf-8")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    print(prepare(parser.parse_args().output_dir))
