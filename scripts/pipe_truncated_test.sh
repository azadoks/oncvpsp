#!/bin/bash
# Checks that an input piped in that ends before its data do is reported with
# the line the reader stopped at, rather than a run-time error.
#
# usage: pipe_truncated_test.sh EXE INPUT.dat
set -euo pipefail
exe=$1; dat=$2
out=$(head -n 20 "$dat" | "$exe" 2>&1 || true)
grep -q 'Read ERROR, input data file line  21' <<< "$out" || { echo "$out" | tail -5 >&2; exit 1; }
echo "$(basename "$exe") reports the input ending at line 21"
