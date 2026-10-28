#!/bin/bash
# The Wigner-Seitz sphere of oncvpsp_confine.x, from the all-electron and
# pseudopotential files (the checks are those of confine_check.py):
# - the band edges of the valence bands agree with the all-electron atom's;
# - each radius gives the same as alone;
# - the edges depend neither on the step of the scan that brackets them nor
#   on a window that takes in core bands and energies where the
#   pseudopotential no longer follows the atom (which once paired each band
#   with the next one down);
# - the cutoffs rise as the error falls;
# - in the walls, the free atoms are the run's reference atoms, and the
#   pseudo-atom's energy rises as the all-electron atom's does;
# - the summary line is that of the other lines.
# With MODE ghost (60_Nd_GHOST), only that the band of the core 4f is
# ws_core, that the window has no pseudopotential-only edge (its f wave
# function gains a node inside the projectors, which once made one at
# +0.051 Ha), and that the f ghost of the run's bound-ghost test at -50 Ha
# is ws_extra in a window around it.
#
# usage: confine_test.sh FULL_EXE AE_EXE CONFINE_EXE INPUT.dat WORKDIR [MODE]
# shellcheck source=scripts/test_lib.sh
source "$(dirname "$0")/test_lib.sh"
full=$1; ae=$2; confine=$3; dat=$4; work=$5
mode=${6:-}
mkdir -p "$work"
cd "$work"

if [ -n "${ONCVPSP_AE_FILE:-}" ]; then
    cp "$ONCVPSP_AE_FILE" atom.ae
else
    "$ae" -a atom.ae < "$dat" > ae.out
fi
"$full" --no-log-derivatives --ps-file atom.ps < "$dat" > full.out

"$confine" -a atom.ae -p atom.ps > conf.out
if [ "$mode" = ghost ]; then
    "$confine" -a atom.ae -p atom.ps -l 3 -e -55,-45,0.1 -w none > ghost.out
else
    read -r r1 r2 < <(awk '$1 == "ws_r" {print $4}' conf.out | sort -g -u | head -2 | xargs)
    "$confine" -a atom.ae -p atom.ps -r "$r1" -w none > r1.out
    "$confine" -a atom.ae -p atom.ps -r "$r2" -w none > r2.out
    "$confine" -a atom.ae -p atom.ps -r "$r1,$r2" -w none > r12.out
    "$confine" -a atom.ae -p atom.ps -r "$r1,$r2" -l 1 -w none > l1.out
    "$confine" -a atom.ae -p atom.ps -r "$r1,$r2" -e -2.0,2.0,0.137 -w none > step.out
    "$confine" -a atom.ae -p atom.ps -r "$r1,$r2" -e -12.0,13.0,0.1 -w none > wide.out
fi

"$python" "$scripts/confine_check.py" "$mode"
