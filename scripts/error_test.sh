#!/bin/bash
# Checks that an input error stops the program with exit status 1 and the
# expected message (see "Error messages" in doc/coding_standards.md).
#
# usage: error_test.sh EXE MESSAGE WORKDIR [INPUT.dat SED_EXPR]
# With INPUT.dat, EXE reads INPUT.dat changed by the sed expression
# SED_EXPR on stdin; without it, EXE is given an input file that does not
# exist.
set -uo pipefail
exe=$(realpath "$1"); message=$2; work=$3
mkdir -p "$work"
cd "$work"

if [ $# -ge 5 ]; then
    sed -e "$5" "$4" > input.dat
    if cmp -s "$4" input.dat; then
        echo "error_test: the sed expression did not change $(basename "$4")" >&2
        exit 1
    fi
    "$exe" < input.dat > out.txt 2>&1
else
    "$exe" -i no_such_file.dat > out.txt 2>&1
fi
status=$?
tail -20 out.txt
if [ $status -ne 1 ]; then
    echo "error_test: exit status $status, expected 1" >&2
    exit 1
fi
grep -qF -- "$message" out.txt
