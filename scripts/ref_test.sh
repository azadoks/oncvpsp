#!/bin/bash
# Checks the output for one test input against its reference in tests/refs.
# Everything after "ATOM AND REFERENCE CONFIGURATION" is compared with
# fldiff.pl -easy, apart from the lines naming pspd or a date and the PSML
# block (whose values near r = 0 depend on the compiler flags; see
# psml_test.sh), and the test passes if no number differs.
#
# usage: ref_test.sh EXE INPUT.dat REF.out WORKDIR
# The output and the fldiff report are left in WORKDIR as REF.out and
# REF.diff (with REF the basename of the reference).
set -euo pipefail
exe=$(realpath "$1"); dat=$(realpath "$2"); ref=$(realpath "$3"); work=$4
fldiff=$(realpath "$(dirname "$0")")/fldiff.pl
mkdir -p "$work"
cd "$work"
name=$(basename "$ref" .out)

if ! "$exe" < "$dat" > "$name.out"; then
    echo "$(basename "$exe") failed on $(basename "$dat")" >&2
    tail -20 "$name.out" >&2
    exit 1
fi

compared() {
    awk 'out {print} /ATOM AND REFERENCE CONFIGURATION/ {out=1}' "$1" |
        sed -e '/pspd/s/^/-/' -e '/date/s/^/-/' -e '/^Begin PSML/,/^END_PSP/s/^/-/'
}
compared "$ref" > ref.cmp
compared "$name.out" > out.cmp
# fldiff.pl exits non-zero when lines differ, which they may only in text
perl "$fldiff" -easy ref.cmp out.cmp > "$name.diff" || true
rm -f ref.cmp out.cmp
cat "$name.diff"
grep -qE 'Summary  : no significant difference|Summary  : different lines=\s*[0-9]+ , max abs_diff= -?0\.000e[+-]00 \(.*\), max rel_diff= -?0\.000e[+-]00 \(.*\)' "$name.diff"
