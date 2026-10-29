# Linux relative RUNPATH validation, 2026-10-06

## Result and scope

The Ubuntu 26.04 x86_64 validation package contains 66 dynamic ELF objects.
Every object has package-relative DT_RUNPATH; none retains DT_RPATH or a
path-bearing DT_NEEDED entry. All non-glibc dependencies resolve inside the
relocated package without LD_LIBRARY_PATH. Package fingerprints remain unchanged
after the bundled-Python concurrency matrix.

This is a Release build of the working tree based on HEAD
`2e06db6d1141f99dd402f369c07a2a3cac0ebd1e`, including uncommitted packaging
changes. It is a validation artifact, not a newly published release. The
previous 8dc5498c archives remain separate.

## Changes

- CMake enables new ELF dtags and installs executable RUNPATH as
  `$ORIGIN/../lib`, with automatic link-directory additions disabled.
- `package_vela.py` completes the Python runtime's non-glibc dependency closure
  and normalizes every dynamic ELF, including copied libraries. It verifies
  the resulting tags before computing manifest hashes.
- A verified package-local ORIGIN-qualified dependency in `libpython3.so` is
  converted to `libpython3.12.so.1.0`. Absolute or escaping dependencies fail
  packaging. The manifest records this rewrite.
- Python's `_crypt` extension now includes `libcrypt.so.1` and the `libcrypt1`
  copyright text. `libutil.so.1` is provided by Ubuntu's `libc6` package and
  remains a platform dependency; its ownership was checked with dpkg-query.
- The concurrency controller identifies solver executables through
  `/proc/<pid>/exe`, which supports spaces and Unicode in package paths, and
  disables cache writes when explicitly selecting an interpreter.

| Object | RUNPATH |
| --- | --- |
| Root `vela` | `$ORIGIN/lib` |
| `bin/vela`, `bin/vela_example_runner`, `bin/sentaurus_import` | `$ORIGIN/../lib` |
| Libraries in `lib/` | `$ORIGIN` |
| Python runtime and extensions | Object directory, private Python lib, then Vela lib; every entry uses `$ORIGIN` |

## Verification

Ubuntu environment: 26.04.1 LTS, x86_64, existing GCC/Ninja Release build.
HDF5/TDR import and SuiteSparse UMFPACK/SPQR support were detected during
configuration. The packaging editor is patchelf 0.18.0, extracted from the
Ubuntu package into the ignored build directory without a system installation.

| Check | Result |
| --- | --- |
| Release build of three shipped executables | Passed |
| `sdevice_package_smoke` | 15 tests passed |
| `sdevice_concurrency` via CTest/build launcher | 4 tests passed |
| Fresh relocated package with bundled CPython 3.12.12 | 4 concurrency tests passed |
| Independent executable -> library -> library relocation probe | Passed; detects failure with executable-only RUNPATH |
| All 66 dynamic ELF tags and ldd closures | Passed; no untagged objects or missing dependencies |
| Relocated root/bin launchers, default material export and synthetic IV | Passed; IV reaches 10 V |
| Manifest fingerprints before and after relocated concurrency | Passed |
| Windows offline packaging contracts | 11 passed; Linux ELF probe skipped |
| Existing Windows archive smoke | 13 passed; 2 Linux checks skipped |
| Windows lock conflict after controller change | Passed; 3 Linux matrix checks skipped |
| `git diff --check` | Passed |

The relocated package directory contains spaces and Chinese characters. Its
controller uses the bundled Python with `-B`, PATH is `/nonexistent`, and
LD_LIBRARY_PATH is unset. Four commands in one directory and four directories
match serial PLT bytes and all ASCII VTU numerical arrays exactly. Observed
solver thread limits are 1 for the default runs and at most 2 for the explicit
`[1, 2, 1, 2]` matrix. Duplicate-command locking returns 7 and preserves the
first owner's lock.

The synthetic mesh uses four nodes/two triangles, coordinates in um, donors
at 1e15 cm^-3 and terminal currents in A/um. Full native PN2D 0V/IV/BV spatial
acceptance was not rerun for this runtime-path change.

## Retained artifact and evidence

- Archive: `build-release/runpath_20261006/packages/vela-0.1.0-ubuntu26.04-x86_64.tar.gz`
- SHA-256: `306e533a8d802cec282bdf25b670b7a985dba89e92bccae5fcbb360cf6e6d3c1`
- Adjacent `.manifest.json` and `.sha256` files record the package inventory
  and fingerprints; `elf_runpaths` contains each object's verified tags.
- Ubuntu logs and controllers: `build-linux26/runpath_20261006/`, including
  `configure.log`, `build.log`, `package.log`, `package-tests.log`,
  `elf-inventory.json`, `relocated-concurrency.log` and `relocated-concurrency.json`.
- Windows logs: `build-release/runpath_20261006/windows-*.log`.

These generated artifacts remain in ignored build directories. See the
[packaging guide](../sdevice_packaging.md) for the supported assembly workflow.
