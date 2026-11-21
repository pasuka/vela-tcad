#!/usr/bin/env bash
set -euo pipefail

packages=(
    build-essential ccache cmake ninja-build pkg-config patchelf
    libsuitesparse-dev libopenblas-pthread-dev libmetis-dev libhdf5-dev
    libboost-dev libspdlog-dev libeigen3-dev nlohmann-json3-dev catch2
    libglu1-mesa python3-venv python3 python3-pil python3-numpy python3-h5py
    python3-jsonschema
)
cache="$RUNNER_TEMP/vela-apt-debs"
manifest="$RUNNER_TEMP/vela-apt-debs.sha256"
archives="$RUNNER_TEMP/vela-apt-install"
apt_archives=(-o "Dir::Cache::archives=$archives")

case "${1:-}" in
    plan)
        mkdir -p "$cache" "$archives"
        sudo apt-get update
        # Resolve dependencies against the current runner and signed APT indexes.
        LC_ALL=C apt-get "${apt_archives[@]}" -o Acquire::ForceHash=sha256 \
            --print-uris --download-only --yes install "${packages[@]}" |
            awk -v quote="'" '
                substr($0, 1, 1) == quote {
                    if ($2 !~ /^[A-Za-z0-9][A-Za-z0-9_.+%:~-]*\.deb$/ ||
                        $4 !~ /^SHA256:[0-9a-f]+$/ || length($4) != 71) exit 1
                    print substr($4, 8) "  " $2
                }
            ' | sort -u > "$manifest"
        echo "key=$(sha256sum "$manifest" | cut -d' ' -f1)" >> "$GITHUB_OUTPUT"
        ;;
    install)
        while read -r digest name; do
            from="$cache/$name"
            if [[ -f $from && ! -L $from ]] &&
               printf '%s  %s\n' "$digest" "$from" | sha256sum --check --status; then
                cp -- "$from" "$archives/$name"
            fi
        done < "$manifest"
        sudo DEBIAN_FRONTEND=noninteractive apt-get "${apt_archives[@]}" \
            -o APT::Keep-Downloaded-Packages=true --yes install "${packages[@]}"
        # A partial cache hit can include older package versions. Save only the
        # archives verified against the manifest for this resolver transaction.
        find "$cache" -maxdepth 1 -name '*.deb' -type f -delete
        while read -r digest name; do
            from="$archives/$name"
            if [[ -f $from ]] &&
               printf '%s  %s\n' "$digest" "$from" | sha256sum --check --status; then
                cp -fT --remove-destination -- "$from" "$cache/$name"
            fi
        done < "$manifest"
        ;;
    *)
        echo "Usage: $0 plan|install" >&2
        exit 2
        ;;
esac
