#!/bin/bash
# Checks that a test configuration whose occupation the pseudo-atom SCF
# ramps to from the reference one is reached exactly: with 3s going from
# 2 to 0.1 in INPUT.dat, 2 + (0.1 - 2) rounds to 0.10000000000000009, and
# psatom's convergence test, which compares the occupations exactly, once
# never passed, so the configuration was dropped ("no output for
# configuration").
#
# usage: occupation_ramp_test.sh EXE INPUT.dat WORKDIR
set -euo pipefail
exe=$(realpath "$1"); dat=$(realpath "$2"); work=$3
mkdir -p "$work"
cd "$work"

# the test configurations replaced by 3s 0.1, 3p 2
awk '/^# ncnf/ {print; print "    1"; print "#"; print "    2";
                print "    3    0    0.1"; print "    3    1    2.0"; exit}
     {print}' "$dat" > ramp.dat
"$exe" < ramp.dat > ramp.out 2>/dev/null
if grep -q "no output for configuration\|failed to converge" ramp.out; then
    echo "occupation_ramp_test: ERROR the configuration 3s 0.1 was dropped"
    exit 1
fi
grep -q "^Test configuration 1" ramp.out || {
    echo "occupation_ramp_test: ERROR no output for configuration 1"
    exit 1
}
echo "configuration 3s 0.1, 3p 2 converged"
