#!/bin/bash
# runs up to 3 arms in parallel (about 1.2 GB RAM each); usage: scripts/run_family.sh arm1 arm2 ...
cd "$(dirname "$0")/.."
N=${JOBS:-3}; i=0
for a in "$@"; do ./scripts/run_arm.sh $a & i=$((i+1)); if [ $((i % N)) -eq 0 ]; then wait; fi; done; wait
