#!/bin/bash
# Checks the onset of the well that confines the well states (wellstate,
# well_onset in [pseudopotentials]; doc/smoothness_plan.md, R-12): it is
# rc times a factor that does not depend on the log mesh, so the well
# states found on two meshes agree, and the input sets the factor.
#
# usage: well_onset_test.sh EXE BASE.toml WORKDIR
# BASE.toml must have well states and no [log_mesh] table. Runs it on the
# default mesh (amesh 1.006) and on amesh 1.012, with the default factor
# and with well_onset = 1, and checks that
#   - the well states (their n and half-point radius) on the two meshes
#     agree to 5e-4 bohr for each factor (with the onset rc amesh**5 they
#     differed by 1e-2 and more);
#   - the half-point radii with well_onset = 1 differ from those of the
#     default by more than 1e-3 (the input is used);
#   - the echo of the input carries well_onset (on the lmax line) only
#     when it is set.
# shellcheck source=scripts/test_lib.sh
source "$(dirname "$0")/test_lib.sh"
exe=$1; base=$2; work=$3
mkdir -p "$work"
cd "$work"

# wells OUT: one line "n radius" per well state of OUT
wells() {
    awk '/Wellstate for l =/ {n = $NF}
         /half-point radius =/ {print n, $NF}' "$1"
}

run() {
    local name=$1 onset=$2 amesh=$3
    if [ -n "$onset" ]; then
        sed "/^\[pseudopotentials\]\$/a well_onset = $onset" "$base" > "$name.toml"
    else
        cp "$base" "$name.toml"
    fi
    if [ -n "$amesh" ]; then
        printf '\n[log_mesh]\namesh = %s\n' "$amesh" >> "$name.toml"
    fi
    "$exe" -t "$name.toml" > "$name.out"
    wells "$name.out" > "$name.wells"
    if [ ! -s "$name.wells" ]; then
        echo "$name: no well states in the output" >&2
        exit 1
    fi
}

# same_wells A B TOL: the same n and radii within TOL
same_wells() {
    paste "$1" "$2" | awk -v tol="$3" -v a="$1" -v b="$2" '
        $1 != $3 {print a ", " b ": n " $1 " against " $3; bad = 1}
        $1 == $3 && ($2 - $4 > tol || $4 - $2 > tol) {
            print a ", " b ": half-point radius " $2 " against " $4; bad = 1}
        END {exit bad}' || return 1
    [ "$(wc -l < "$1")" -eq "$(wc -l < "$2")" ] || {
        echo "$1, $2: different numbers of well states"; return 1; }
}

run default "" ""
run default_coarse "" 1.012
run one 1.0 ""
run one_coarse 1.0 1.012

# the half-point radii of the well states on two meshes (bohr): the onset at
# a fixed radius, not a number of mesh points (with rc amesh**5 they differed
# by 1e-2 and more)
mesh_tol=5e-4
# the least change of a radius from the default onset to 1 bohr
onset_change=1e-3

failed=0
same_wells default.wells default_coarse.wells "$mesh_tol" || failed=1
same_wells one.wells one_coarse.wells "$mesh_tol" || failed=1
if same_wells default.wells one.wells "$onset_change" > /dev/null; then
    echo "well_onset = 1 gives the well states of the default" >&2
    failed=1
fi
if grep -q 'well_onset' default.out; then
    echo "the default echoes well_onset" >&2
    failed=1
fi
if ! grep -A1 '^# lmax, ep_as_given, well_onset, well_tail$' one.out |
        grep -Eq '^ +[0-9]+ +0 +1\.0+E\+00 +0$'; then
    echo "well_onset = 1 is not echoed on the lmax line" >&2
    failed=1
fi
if [ $failed -ne 0 ]; then
    exit 1
fi
echo "well_onset: the well states agree on both meshes, and the input sets the onset"
