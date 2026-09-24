#!/usr/bin/env bash
set -euo pipefail
script_dir=$(dirname "${BASH_SOURCE[0]}")
cd "$script_dir/../.."

action=${1:-build}
jobs=${2:-4}
if [[ ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "Jobs must be a positive integer" >&2
    exit 2
fi
# Keep the previous GCC 13 cache and binaries separate during migration.
build_dir=build-codespaces-gcc16-release
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export VELA_LINEAR_THREADS=${VELA_LINEAR_THREADS:-1}

configure() {
    for compiler in /usr/bin/gcc-16 /usr/bin/g++-16; do
        if [[ ! -x "$compiler" ]]; then
            echo "Missing $compiler; rebuild the Codespaces container after syncing its configuration." >&2
            exit 1
        fi
        compiler_version=$("$compiler" -dumpversion)
        if [[ "${compiler_version%%.*}" != 16 ]]; then
            echo "Expected GCC 16, got $compiler_version from $compiler" >&2
            exit 1
        fi
        "$compiler" --version | head -n 1
    done
    /usr/bin/python3 -c 'import h5py, numpy, PIL; print(h5py.version.info)'
    mkdir -p "$build_dir"
    cmake -S . -B "$build_dir" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DCMAKE_C_COMPILER=/usr/bin/gcc-16 -DCMAKE_CXX_COMPILER=/usr/bin/g++-16 \
        -DPython3_EXECUTABLE=/usr/bin/python3 -DVELA_ENABLE_PYTHON=OFF \
        -DVELA_ENABLE_HDF5=ON -DVELA_ENABLE_HDF5_STATE=ON \
        -DVELA_ENABLE_UMFPACK=ON -DVELA_ENABLE_METIS=ON \
        -DVELA_ENABLE_OPENBLAS_THREAD_CONTROL=ON \
        -DVELA_ENABLE_MUMPS=OFF -DVELA_ENABLE_SUPERLU_MT=OFF \
        -DVELA_ENABLE_STRUMPACK=OFF 2>&1 | tee "$build_dir/configure.log"
    for feature in \
        'SuiteSparse UMFPACK bordered solver enabled' \
        'SuiteSparse SPQR bordered solver enabled' \
        'METIS ordering enabled' \
        'Independent HDF5/HighFive state archive enabled'; do
        if ! grep -Fq "$feature" "$build_dir/configure.log"; then
            echo "Required feature missing: $feature" >&2
            exit 1
        fi
    done
}

case "$action" in
    configure) configure ;;
    build)
        configure
        cmake --build "$build_dir" --parallel "$jobs" ;;
    test)
        # An unmatched filter must fail, rather than report a false success.
        ctest --test-dir "$build_dir" --output-on-failure --no-tests=error \
            --parallel "$jobs" -R "${3:-.}" ;;
    *) echo "Usage: bash scripts/codespaces/build.sh {configure|build|test} [jobs] [test-regex]" >&2; exit 2 ;;
esac
