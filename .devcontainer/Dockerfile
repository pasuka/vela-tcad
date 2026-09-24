FROM mcr.microsoft.com/devcontainers/base:3.0.8-ubuntu24.04

# Match egor-tensin/setup-gcc@v2 (version: 16) in .github/workflows/ci.yml:
# it enables ubuntu-toolchain-r/test and installs the versioned compiler.
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    software-properties-common ca-certificates \
    && add-apt-repository --yes --update ppa:ubuntu-toolchain-r/test \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential gcc-16 g++-16 cmake ninja-build pkg-config git \
    libsuitesparse-dev libopenblas-pthread-dev libmetis-dev libhdf5-dev \
    libboost-dev libspdlog-dev libeigen3-dev nlohmann-json3-dev catch2 \
    python3 python3-pil python3-numpy python3-h5py \
    && ln -sfn /usr/bin/gcc-16 /usr/local/bin/gcc \
    && ln -sfn /usr/bin/gcc-16 /usr/local/bin/cc \
    && ln -sfn /usr/bin/g++-16 /usr/local/bin/g++ \
    && ln -sfn /usr/bin/g++-16 /usr/local/bin/c++ \
    && gcc --version && g++ --version \
    && rm -rf /var/lib/apt/lists/*

ENV CC=/usr/bin/gcc-16 CXX=/usr/bin/g++-16
