#!/bin/bash
# Checks the TOML input for one test input: the program must give the same
# output for INPUT.toml (read with -t) as for INPUT.dat (read from stdin).
# tests/data/toml holds the TOML form of every tests/data input, made with
# scripts/input_text_to_toml.py. The outputs are compared with
# compare_output.py --identical (every number to 1e-10, absolute or
# relative), not byte for byte: with -ffast-math (the study build) a sum can
# round differently depending on the alignment of its arrays, which differs
# between the two runs.
#
# usage: toml_test.sh EXE INPUT.dat INPUT.toml WORKDIR
# The python, with pyoncvpsp and upf-tools, is $ONCVPSP_PYTHON, or python3.
# $ONCVPSP_RUN, if set, is the output of EXE for INPUT.dat (run_test.sh),
# read in place of running it.
set -euo pipefail
exe=$1; dat=$2; toml=$3; work=$4
mkdir -p "$work"
scripts=$(realpath "$(dirname "$0")")
python=${ONCVPSP_PYTHON:-python3}
cd "$work"

if [ -n "${ONCVPSP_RUN:-}" ]; then
    cp "$ONCVPSP_RUN" dat.out
else
    "$exe" < "$dat" > dat.out
fi
"$exe" -t "$toml" > toml.out
if ! "$python" "$scripts/compare_output.py" --identical dat.out toml.out > toml.diff; then
    cat toml.diff
    exit 1
fi
echo "TOML output the same for $(basename "$toml")"
