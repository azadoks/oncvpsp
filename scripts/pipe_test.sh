#!/bin/bash
# Checks that the text input can be read from a pipe (as from a Python
# subprocess), which the program cannot backspace, as well as from a file:
# both must give the same output (compare_output.py --identical: every
# number to 1e-10, the text exactly, apart from the date).
#
# usage: pipe_test.sh EXE INPUT.dat WORKDIR
# The python, with pyoncvpsp and upf-tools, is $ONCVPSP_PYTHON, or python3.
set -euo pipefail
exe=$(realpath "$1"); dat=$(realpath "$2"); work=$3
scripts=$(realpath "$(dirname "$0")")
python=${ONCVPSP_PYTHON:-python3}
rm -rf "$work"
mkdir -p "$work"
cd "$work"
"$exe" < "$dat" > file.out
# the input must come through a pipe, not a file
# shellcheck disable=SC2002
cat "$dat" | "$exe" > pipe.out
if ! "$python" "$scripts/compare_output.py" --identical file.out pipe.out > pipe.diff; then
    cat pipe.diff
    exit 1
fi
grep -q '^# atsym.*psfile' pipe.out
echo "$(basename "$exe") gives the same output for $(basename "$dat") from a pipe as from a file"
