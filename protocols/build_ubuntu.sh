#!/usr/bin/env bash
set -euo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
pkg-config --exists libsbml clp || {
    printf '%s\n' 'Missing development libraries. Follow protocols/REPRODUCE_LOBULE.md.' >&2
    exit 1
}
# A dedicated temporary build avoids stale objects without changing local sources.
build="$(mktemp -d -t hepatomodel-build-XXXXXXXX)"
printf 'Temporary build: %s\n' "$build"
for dir in BioFVM core modules custom_modules addons/dFBA/src; do
    mkdir -p "$build/$dir"
    cp -a "$root/$dir/." "$build/$dir/"
done
cp "$root/Makefile" "$root/main.cpp" "$root/VERSION.txt" "$build/"
make -C "$build" -j "${BUILD_JOBS:-2}" \
    DFBA_INC="-I./addons/dFBA/src $(pkg-config --cflags libsbml clp)" \
    DFBA_LD= DFBA_RPATH= \
    DFBA_LIBS="$(pkg-config --libs libsbml clp) -lxml2"
if [[ -e "$root/project" ]]; then
    cp -p "$root/project" "$root/project.before-build.$(date +%Y%m%d_%H%M%S)"
fi
cp "$build/project" "$root/project"
printf 'Built %s/project\nBuild files retained at %s\n' "$root" "$build"
