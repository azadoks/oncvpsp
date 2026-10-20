#!/bin/bash
# Checks that the psp8 file carries the model core for an input changed by
# SED: its header must give fchrg > 0, and the file must have the model core
# block, as the one for the unchanged input (whose model core is written)
# has. Both runs write psp8, read with pyoncvpsp.
#
# usage: psp8_core_test.sh EXE INPUT.dat SED WORKDIR
# The python, with pyoncvpsp, is $ONCVPSP_PYTHON, or python3.
set -euo pipefail
exe=$1; dat=$2; edit=$3; work=$4
python=${ONCVPSP_PYTHON:-python3}
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

"$python" - base.out in.out <<'EOF'
import sys

from pyoncvpsp.io import OncvpspTextParser, parse_psp8

base, new = (
    parse_psp8(OncvpspTextParser(path=path, check=False).psp8_string)
    for path in sys.argv[1:]
)
if base["model_core"] is None:
    sys.exit("psp8_core_test.sh: the psp8 file of the unchanged input has no model core")
if new["fchrg"] <= 0:
    sys.exit(f"psp8_core_test.sh: fchrg={new['fchrg']}, no model core")
if new["model_core"] is None or new["model_core"].shape != base["model_core"].shape:
    sys.exit("psp8_core_test.sh: the psp8 file has no model core block")
print(f"the psp8 file has the model core, fchrg={new['fchrg']}")
EOF
