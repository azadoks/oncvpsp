#!/bin/bash
# Checks that the text input can be read from a pipe (as from a Python
# subprocess), which the program cannot backspace, as well as from a file:
# both must give the same output (compare_output.py --identical: every
# number to 1e-10, the text exactly, apart from the date).
#
# usage: pipe_test.sh EXE INPUT.dat WORKDIR
# The python, with pyoncvpsp and upf-tools, is $ONCVPSP_PYTHON, or python3.
# $ONCVPSP_RUN, if set, is the output of EXE for INPUT.dat (run_test.sh),
# read in place of running it.
# shellcheck source=scripts/test_lib.sh
source "$(dirname "$0")/test_lib.sh"
exe=$(realpath "$1"); dat=$(realpath "$2"); work=$3
rm -rf "$work"
mkdir -p "$work"
cd "$work"
if [ -n "${ONCVPSP_RUN:-}" ]; then
    cp "$ONCVPSP_RUN" file.out
else
    "$exe" < "$dat" > file.out
fi
# the input must come through a pipe, not a file
# shellcheck disable=SC2002
cat "$dat" | "$exe" > pipe.out
same file.out pipe.out
grep -q '^# atsym.*psfile' pipe.out
echo "$(basename "$exe") gives the same output for $(basename "$dat") from a pipe as from a file"
