#!/usr/bin/env bash
set -euo pipefail
action=${1:-configure}
jobs=${2:-8}
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || exit 2
: "${CONDA_PREFIX:?Activate the Vela micromamba environment first}"
build_dir=build-mamba-gcc16-eigen5-release
export CMAKE_PREFIX_PATH="$CONDA_PREFIX"
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
mkdir -p "$build_dir"
case "$action" in
configure)
    [[ $("$CXX" -dumpversion) == 16* ]] || { echo 'GCC 16 required' >&2; exit 1; }
    cmake -S . -B "$build_dir" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER="$(command -v "$CC")" \
        -DCMAKE_CXX_COMPILER="$(command -v "$CXX")" \
        -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
        -DPython3_EXECUTABLE="$CONDA_PREFIX/bin/python" \
        -DVELA_ENABLE_PYTHON=OFF \
        -DVELA_ENABLE_HDF5=ON -DVELA_ENABLE_HDF5_STATE=ON \
        -DVELA_ENABLE_UMFPACK=ON -DVELA_ENABLE_METIS=ON \
        -DVELA_ENABLE_OPENBLAS_THREAD_CONTROL=OFF \
        -DVELA_ENABLE_MUMPS=OFF -DVELA_ENABLE_SUPERLU_MT=OFF \
        -DVELA_ENABLE_STRUMPACK=OFF 2>&1 | tee "$build_dir/configure.log"
    for feature in \
        'SuiteSparse UMFPACK bordered solver enabled' \
        'SuiteSparse SPQR bordered solver enabled' \
        'METIS ordering enabled' \
        'Independent HDF5/HighFive state archive enabled'; do
        grep -Fq "$feature" "$build_dir/configure.log" || {
            echo "Required feature missing: $feature" >&2; exit 1;
        }
    done
    ;;
build)
    cmake --build "$build_dir" --parallel "$jobs" 2>&1 | tee "$build_dir/build.log"
    ;;
test)
    ctest --test-dir "$build_dir" --output-on-failure --no-tests=error \
        --parallel "$jobs" -R "${3:-.}" 2>&1 | tee "$build_dir/ctest.log"
    ;;
*) echo 'Expected configure, build or test' >&2; exit 2 ;;
esac
