#!/bin/bash
# Checks the all-electron / pseudopotential split for one input:
# oncvpsp_ae.x writes the all-electron file, and oncvpsp_ps.x must then
# reproduce the output of the single program exactly, both for the input
# itself and for a copy with every rc scaled by 1.1, which reuses the
# same all-electron file. The build lines, which name the executable,
# are checked separately.
#
# usage: split_test.sh FULL_EXE AE_EXE PS_EXE INPUT.dat WORKDIR
set -euo pipefail
full=$1; ae=$2; ps=$3; dat=$4; work=$5
mkdir -p "$work"
cd "$work"

# same_output A.out A_EXE B.out B_EXE: the outputs are identical apart from
# their build lines, which name their executables
same_output() {
    grep -q "^build: $(basename "$2") (" "$1"
    grep -q "^build: $(basename "$4") (" "$3"
    cmp <(grep -v '^build: ' "$1") <(grep -v '^build: ' "$3")
}

"$ae" -a atom.ae < "$dat" > ae.out
"$full" < "$dat" > full.out
"$ps" -a atom.ae < "$dat" > ps.out
same_output full.out "$full" ps.out "$ps"

# rc is the second field of the lmax+1 lines after "# l, rc, ..."
awk '/^# l, rc/ {inrc=1; print; next}
     inrc && /^#/ {inrc=0}
     inrc {$2=sprintf("%.2f", $2*1.1)}
     {print}' "$dat" > rc.dat
if cmp -s rc.dat "$dat"; then
    echo "split_test.sh: could not change rc in $dat" >&2
    exit 1
fi
"$full" < rc.dat > rc_full.out
"$ps" -a atom.ae < rc.dat > rc_ps.out
same_output rc_full.out "$full" rc_ps.out "$ps"
echo "split output identical for $(basename "$dat") and its rc*1.1 variant"
