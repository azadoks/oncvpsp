#!/bin/bash
# Checks the PSML output for one input: with psfile = psml, the program must
# write the document into its output between "Begin PSML" and "END_PSP",
# as it does the psp8 and upf files, and leave no files behind.
#
# usage: psml_test.sh EXE INPUT.dat WORKDIR
# shellcheck source=scripts/test_lib.sh
source "$(dirname "$0")/test_lib.sh"
exe=$1; dat=$2; work=$3
rm -rf "$work"
mkdir -p "$work/run"
# psfile is the last field of the first data line
awk '!/^#/ && !done {sub(/[^ ]+[ ]*$/, "psml"); done=1} {print}' "$dat" > "$work/psml.dat"
cd "$work/run"
"$exe" < ../psml.dat > ../psml.out

left=$(ls -A)
if [ -n "$left" ]; then
    echo "psml_test.sh: files left in the working directory: $left" >&2
    exit 1
fi
grep -q '^# atsym.*psfile' ../psml.out
[ "$(grep -c '^Begin PSML$' ../psml.out)" -eq 1 ]
[ "$(grep -c '^Begin PSP' ../psml.out || true)" -eq 0 ]
# the document: the XML declaration, a psml root element closed before END_PSP
awk '/^Begin PSML$/ {f=1; next} f && /^END_PSP$/ {exit} f' ../psml.out > ../doc.psml
head -1 ../doc.psml | grep -q '^<?xml version="1.0" encoding="UTF-8" ?>'
grep -q '<psml ' ../doc.psml
grep -v '^$' ../doc.psml | tail -1 | grep -q '</psml>$'
grep -q '<input-file name="oncvpsp-input">' ../doc.psml
echo "PSML written to the output of $(basename "$exe") for $(basename "$dat"), no files left"
