#!/bin/bash
# Checks the TOML input for one test input: the program must give the same
# output for INPUT.toml (read with -t) as for INPUT.dat (read from stdin).
# tests/data/toml holds the TOML form of every tests/data input, made with
# scripts/input_text_to_toml.py. The outputs are compared with fldiff.pl
# at its default tolerance (1e-10, absolute or relative), not byte for byte:
# with -ffast-math (the study build) a sum can round differently depending
# on the alignment of its arrays, which differs between the two runs.
#
# usage: toml_test.sh EXE INPUT.dat INPUT.toml WORKDIR
set -euo pipefail
exe=$1; dat=$2; toml=$3; work=$4
mkdir -p "$work"
scripts=$(realpath "$(dirname "$0")")
cd "$work"

"$exe" < "$dat" > dat.out
"$exe" -t "$toml" > toml.out
# fldiff.pl exits non-zero when lines differ, even within the tolerance
perl "$scripts/fldiff.pl" dat.out toml.out > toml.diff || true
if ! grep -q 'Summary  : no significant difference' toml.diff; then
    cat toml.diff
    exit 1
fi
echo "TOML output the same for $(basename "$toml")"
