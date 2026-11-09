#!/usr/bin/env bash
set -euo pipefail

gmsh_version="${GMSH_VERSION:-4.15.2}"
wheel_sha256="${GMSH_LINUX_WHEEL_SHA256:-4076a948ce22625330d1413d4982e22b5c69fc2f0f7951f5df64c778cf54108c}"
temporary_base="${RUNNER_TEMP:-/tmp}"
mkdir -p "${temporary_base}"
temporary_root="$(mktemp -d "${temporary_base%/}/vela-gmsh-${gmsh_version}.XXXXXX")"
download_dir="${temporary_root}/download"
runtime_dir="${temporary_root}/runtime"

mkdir -p "${download_dir}" "${runtime_dir}"
python3 -m venv --system-site-packages "${runtime_dir}"
runtime_python="${runtime_dir}/bin/python"
"${runtime_python}" -m pip download \
  --disable-pip-version-check \
  --timeout 120 \
  --no-deps \
  --only-binary=:all: \
  --dest "${download_dir}" \
  "gmsh==${gmsh_version}"

wheel="$(find "${download_dir}" -maxdepth 1 -type f -name 'gmsh-*.whl' -print -quit)"
test -n "${wheel}"
printf '%s  %s\n' "${wheel_sha256}" "${wheel}" | sha256sum --check --strict
"${runtime_python}" -m pip install \
  --disable-pip-version-check \
  --no-deps \
  "${wheel}"

gmsh_python="$(find "${runtime_dir}" -type f -name gmsh.py -print -quit)"
gmsh_library="$(find "${runtime_dir}" -type f -name 'libgmsh.so*' -print -quit)"
gmsh_executable="$(find "${runtime_dir}" -type f -path '*/bin/gmsh' -print -quit)"
test -f "${gmsh_python}"
test -f "${gmsh_library}"
test -f "${gmsh_executable}"
chmod +x "${gmsh_executable}"

if [[ -n "${GITHUB_ENV:-}" ]]; then
  printf 'VELA_GMSH_PYTHON=%s\n' "${gmsh_python}" >> "${GITHUB_ENV}"
  printf 'VELA_GMSH_LIBRARY=%s\n' "${gmsh_library}" >> "${GITHUB_ENV}"
  printf 'VELA_GMSH_EXECUTABLE=%s\n' "${gmsh_executable}" >> "${GITHUB_ENV}"
fi
if [[ -n "${GITHUB_PATH:-}" ]]; then
  dirname "${gmsh_executable}" >> "${GITHUB_PATH}"
fi

printf 'VELA_GMSH_PYTHON=%s\n' "${gmsh_python}"
printf 'VELA_GMSH_LIBRARY=%s\n' "${gmsh_library}"
printf 'VELA_GMSH_EXECUTABLE=%s\n' "${gmsh_executable}"
