FROM mcr.microsoft.com/devcontainers/base:3.0.8-ubuntu24.04

# Ubuntu 24.04 supplies GCC 13 and Catch2 3. Keep dependencies aligned with CI;
# CI's separately installed GCC 16 is not required by the project's C++20 API.
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build pkg-config git ca-certificates \
    libsuitesparse-dev libopenblas-pthread-dev libmetis-dev libhdf5-dev \
    libboost-dev libspdlog-dev libeigen3-dev nlohmann-json3-dev catch2 \
    python3 python3-pil python3-numpy python3-h5py \
    && rm -rf /var/lib/apt/lists/*

ENV CC=/usr/bin/gcc CXX=/usr/bin/g++
