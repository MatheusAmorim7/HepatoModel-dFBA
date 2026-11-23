#!/bin/bash
# usage: scripts/run_arm.sh <arm>   (run from the project directory; uses config/arms/<arm>/control.xml)
set -e
cd "$(dirname "$0")/.."
ARM=$1; [ -n "$ARM" ] || { echo "usage: $0 <arm>"; exit 1; }
ATPM=$(python3 -c "import json,sys;print(json.load(open('config/arms.json'))['arms']['$ARM']['atpm'])")
rm -rf runs/$ARM; mkdir -p runs/$ARM
echo "[run] $ARM (ATPM=$ATPM)"
ATPM_EXPECTED=$ATPM ./hepatic_lobule_dfba config/arms/$ARM/control.xml $ARM runs/$ARM replay warm gzip > runs/${ARM}_stdout.txt 2> runs/${ARM}_stderr.txt || { rc=$?; echo "[run] $ARM exited rc=$rc (rc=42 = atomic rejection, see runs/$ARM/rejection.json)"; exit 0; }
echo "[run] $ARM done"
