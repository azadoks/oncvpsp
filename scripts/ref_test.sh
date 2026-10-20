#!/bin/bash
# Checks the output for one test input against its reference in tests/refs
# with compare_output.py, which parses both and compares them quantity by
# quantity (see its usage): the numbers to their tolerances, the
# pseudopotential by its nonlocal kernels |beta>D<beta| (to 5e-5) and its
# local potential, never by the projectors and their energies, which are not
# unique. The PSML file is left out: its values near r = 0 depend on the
# compiler flags (see psml_test.sh).
#
# usage: ref_test.sh EXE INPUT.dat REF.out WORKDIR
# The output and the comparison are left in WORKDIR as REF.out and REF.diff
# (with REF the basename of the reference). The python, with pyoncvpsp and
# upf-tools, is $ONCVPSP_PYTHON, or python3.
set -euo pipefail
exe=$(realpath "$1"); dat=$(realpath "$2"); ref=$(realpath "$3"); work=$4
scripts=$(realpath "$(dirname "$0")")
python=${ONCVPSP_PYTHON:-python3}
mkdir -p "$work"
cd "$work"
name=$(basename "$ref" .out)

if ! "$exe" < "$dat" > "$name.out"; then
    echo "$(basename "$exe") failed on $(basename "$dat")" >&2
    tail -20 "$name.out" >&2
    exit 1
fi

status=0
"$python" "$scripts/compare_output.py" "$ref" "$name.out" > "$name.diff" || status=$?
cat "$name.diff"
exit $status
