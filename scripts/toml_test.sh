#!/bin/bash
# Checks the TOML input for one test input: the program must give the same
# output for INPUT.toml (read with -t) as for INPUT.dat (read from stdin).
# tests/data/toml holds the TOML form of every tests/data input, made with
# scripts/input_text_to_toml.py.
#
# usage: toml_test.sh EXE INPUT.dat INPUT.toml WORKDIR
set -euo pipefail
exe=$1; dat=$2; toml=$3; work=$4
mkdir -p "$work"
cd "$work"

"$exe" < "$dat" > dat.out
"$exe" -t "$toml" > toml.out
cmp dat.out toml.out
echo "TOML output identical for $(basename "$toml")"
