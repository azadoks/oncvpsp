#!/bin/bash
# Checks that the psp8 file carries the model core for an input changed by
# SED: its header must give fchrg > 0, and the file must have as many lines
# as the one for the unchanged input (whose model core is written). Both
# runs write psp8.
#
# usage: psp8_core_test.sh EXE INPUT.dat SED WORKDIR
set -euo pipefail
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

# psp8 OUT: the psp8 file in OUT, without the echo of the input
psp8() { awk '/^Begin PSPCODE8$/ {on=1} /^<INPUT>$/ {on=0} on' "$1"; }

fchrg=$(psp8 in.out | awk '/rchrg fchrg qchrg/ {print $2}')
if ! awk -v f="$fchrg" 'BEGIN {exit !(f > 0)}'; then
    echo "psp8_core_test.sh: fchrg=$fchrg, no model core" >&2; exit 1
fi
if [ "$(psp8 base.out | wc -l)" -ne "$(psp8 in.out | wc -l)" ]; then
    echo "psp8_core_test.sh: the psp8 file has $(psp8 in.out | wc -l) lines," \
         "$(psp8 base.out | wc -l) for $(basename "$dat")" >&2
    exit 1
fi
echo "the psp8 file has the model core, fchrg=$fchrg"
