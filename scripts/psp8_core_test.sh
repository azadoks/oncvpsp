#!/bin/bash
# Checks that the psp8 file carries the model core for an input changed by
# SED: its header must give fchrg > 0, and the file must have the model core
# block, as the one for the unchanged input (whose model core is written)
# has. Both runs write psp8, read with pyoncvpsp.
#
# usage: psp8_core_test.sh EXE INPUT.dat SED WORKDIR
# The python, with pyoncvpsp, is $ONCVPSP_PYTHON, or python3.
# shellcheck source=scripts/test_lib.sh
source "$(dirname "$0")/test_lib.sh"
exe=$1; dat=$2; edit=$3; work=$4
rm -rf "$work"
mkdir -p "$work"
cd "$work"

# psfile is the last field of the first data line
awk '!/^#/ && !done {sub(/[^ ]+[ ]*$/, "psp8"); done=1} {print}' "$dat" > base.dat
sed -e "$edit" base.dat > in.dat
if cmp -s base.dat in.dat; then
    echo "psp8_core_test.sh: the sed script did not change the input" >&2; exit 1
fi
"$exe" < base.dat > base.out
"$exe" < in.dat > in.out

"$python" "$scripts/psp8_core_check.py" base.out in.out
