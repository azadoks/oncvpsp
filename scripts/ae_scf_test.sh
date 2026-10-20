#!/bin/bash
# Checks that the all-electron SCF (sratom) does not stop before it is
# self-consistent: on a log mesh starting at rr1 = 1e-8 bohr, the H atom of
# INPUT.dat once stopped at iteration 1, 0.29 Ha above the converged total
# energy, because its eigenvalue alone did not move. The total energy must
# be that of the built-in mesh, to 1e-8 Ha, after more than one iteration.
#
# usage: ae_scf_test.sh AE_EXE INPUT.dat WORKDIR
set -euo pipefail
ae=$(realpath "$1"); dat=$(realpath "$2"); work=$3
mkdir -p "$work"
cd "$work"

run() { "$ae" -a ae.bin > "$1.out"; }
run def < "$dat"
{ cat "$dat"; echo "# amesh, rr1, rmax"; echo "1.006 1e-08 -1"; } | run rr1

energy() { awk '$1 == "total" {print $2}' "$1.out"; }
iterations() { awk '$1 == "iterations" {print $2}' "$1.out"; }
e_def=$(energy def); e_rr1=$(energy rr1); it_rr1=$(iterations rr1)
echo "built-in mesh: $e_def Ha; rr1 = 1e-8: $e_rr1 Ha in $it_rr1 iterations"
awk -v a="$e_def" -v b="$e_rr1" -v it="$it_rr1" \
    'BEGIN {d = a - b; if (d < 0) d = -d; exit !(it > 1 && d < 1e-8)}' || {
    echo "ae_scf_test: ERROR the SCF on the rr1 = 1e-8 mesh is not converged"
    exit 1
}
