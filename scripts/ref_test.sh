#!/bin/bash
# Checks the output for one test input against its reference in tests/refs.
# Everything after "ATOM AND REFERENCE CONFIGURATION" is compared with
# fldiff.pl -easy, apart from the lines naming pspd or a date, the PSML
# block (whose values near r = 0 depend on the compiler flags; see
# psml_test.sh) and the projectors and their energies, which are not unique
# (kernel_diff.py --mask). The pseudopotential is compared instead by its
# nonlocal kernels |beta>D<beta|, which must agree to 5e-5 (kernel_diff.py).
# The test passes if no number differs and the kernels agree.
#
# usage: ref_test.sh EXE INPUT.dat REF.out WORKDIR
# The output and the fldiff report are left in WORKDIR as REF.out and
# REF.diff (with REF the basename of the reference). The python is
# $ONCVPSP_PYTHON, or python3.
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

compared() {
    "$python" "$scripts/kernel_diff.py" --mask "$1" |
        awk 'out {print} /ATOM AND REFERENCE CONFIGURATION/ {out=1}' |
        sed -e '/pspd/s/^/-/' -e '/date/s/^/-/' -e '/^Begin PSML/,/^END_PSP/s/^/-/'
}
compared "$ref" > ref.cmp
compared "$name.out" > out.cmp
# fldiff.pl exits non-zero when lines differ, which they may only in text
perl "$scripts/fldiff.pl" -easy ref.cmp out.cmp > "$name.diff" || true
rm -f ref.cmp out.cmp
cat "$name.diff"
status=0
grep -qE 'Summary  : no significant difference|Summary  : different lines=\s*[0-9]+ , max abs_diff= -?0\.000e[+-]00 \(.*\), max rel_diff= -?0\.000e[+-]00 \(.*\)' "$name.diff" || status=1
printf "\nNonlocal kernels and local potential against the reference:\n"
"$python" "$scripts/kernel_diff.py" "$ref" "$name.out" --tol 5e-5 || status=1
exit $status
