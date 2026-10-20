#!/bin/bash
# Checks the all-electron atom with a 4f shell (lschfb, ldiracfb):
#
#   scf   neutral Yb ([Xe] 4f14 6s2): in an early SCF potential its 4f has
#         no bound state, which once stopped the run; it now gets the box
#         state, and the SCF must converge
#   mesh  neutral Tm ([Xe] 4f13 6s2): its shallow 4f decays slowly beyond
#         the turning point, and the inward integration once started at 10
#         times it, so that the total energy converged to first order in
#         the mesh spacing (1e-3 Ha between amesh 1.006 and 1.0015); it must
#         now agree to 5e-6 Ha (it does to 1.8e-6)
#
# usage: ae_4f_test.sh AE_EXE INPUT.dat CASE WORKDIR
set -euo pipefail
ae=$(realpath "$1"); dat=$(realpath "$2"); case=$3; work=$4
mkdir -p "$work"
cd "$work"

energy() { awk '$1 == "total" {print $2}' "$1.out"; }

case $case in
    scf)
        "$ae" -a ae.bin < "$dat" > scf.out || true
        e=$(energy scf)
        if [ -z "$e" ]; then
            echo "ae_4f_test: ERROR the SCF did not converge:"; tail -5 scf.out
            exit 1
        fi
        echo "total energy $e Ha" ;;
    mesh)
        for a in 1.006 1.0015; do
            { cat "$dat"; echo "# amesh, rr1, rmax"; echo "$a 1.45e-10 40.0"; } |
                "$ae" -a ae.bin > "m$a.out"
        done
        e1=$(energy m1.006); e2=$(energy m1.0015)
        echo "amesh 1.006: $e1 Ha; amesh 1.0015: $e2 Ha"
        awk -v a="$e1" -v b="$e2" \
            'BEGIN {d = a - b; if (d < 0) d = -d; exit !(d < 5e-6)}' || {
            echo "ae_4f_test: ERROR the total energy is not converged in the mesh"
            exit 1
        } ;;
    *)
        echo "unknown case $case"; exit 1 ;;
esac
