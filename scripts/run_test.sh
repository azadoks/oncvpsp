#!/bin/bash
# Runs one test input once for the tests that read its results (the CTest
# fixtures of tests/CMakeLists.txt): run_REL_INPUT, the output of oncvpsp,
# which the reference, TOML, pipe and split tests of the input read from
# $ONCVPSP_RUN; run_ae_REL_INPUT, the all-electron file of oncvpsp_ae.x, which
# the split and all-electron file tests read from $ONCVPSP_AE_FILE.
#
# usage: run_test.sh EXE INPUT.dat OUT [AE_FILE]
#   AE_FILE  EXE is oncvpsp_ae.x, which writes the all-electron file AE_FILE
set -euo pipefail
exe=$(realpath "$1"); dat=$(realpath "$2"); out=$3; ae_file=${4:-}
mkdir -p "$(dirname "$out")"
cd "$(dirname "$out")"
args=()
if [ -n "$ae_file" ]; then
    args=(-a "$ae_file.tmp")
fi
if ! "$exe" "${args[@]}" < "$dat" > "$out.tmp"; then
    echo "$(basename "$exe") failed on $(basename "$dat")" >&2
    tail -20 "$out.tmp" >&2
    exit 1
fi
if [ -n "$ae_file" ]; then
    mv "$ae_file.tmp" "$ae_file"
fi
mv "$out.tmp" "$out"
