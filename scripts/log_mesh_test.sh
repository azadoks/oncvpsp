#!/bin/bash
# Checks the optional log mesh input (amesh, rr1, rmax), the last line of the
# text input or the [log_mesh] table of the TOML input, for one input:
#
#   defaults      0 or less for every value: the output of the input without it
#   explicit      the built-in amesh and rmax given: the output of the input
#                 without them, apart from the echo of the mesh
#   rerun         another mesh: an output that differs from the built-in mesh's,
#                 and whose echo of the input, and the input in its psp8 file
#                 (if it writes one), run again, give the same output
#   upf           another mesh, with psfile upf: the input in the UPF file, run
#                 again, gives the same output
#   split         another mesh for oncvpsp_ae.x and oncvpsp_ps.x: the output of
#                 the single program
#   split_differ  another mesh for oncvpsp_ae.x only: oncvpsp_ps.x must stop
#   toml          another mesh as [log_mesh] in INPUT.toml: the output of the
#                 text input
#   amesh, rr1, rmax, read  a bad value or line: the program must stop with a
#                 message naming it
#
# usage: log_mesh_test.sh FULL_EXE AE_EXE PS_EXE INPUT.dat CASE WORKDIR [INPUT.toml]
set -euo pipefail
full=$1; ae=$2; ps=$3; dat=$4; case=$5; work=$6; toml=${7:-}
mkdir -p "$work"
cd "$work"

# another mesh: finer, and further out, than the built-in one
mesh='1.004 -1 50.0'

# with_mesh LINE < IN > OUT: the input with the log mesh line LINE after it
with_mesh() { cat; echo "# amesh, rr1, rmax"; echo "$1"; }

# drop_echo < OUT: the output without the echoes of the log mesh
drop_echo() {
    python3 -c 'import re, sys
sys.stdout.write(re.sub(r"#\n# LOG MESH\n[^\n]*\n[^\n]*\n", "", sys.stdin.read()))'
}

# rerun_block OUT START END: run the input between the lines START and END of
# OUT again, which must set the mesh and give OUT
rerun_block() {
    awk -v start="$2" -v end="$3" '$0 ~ end {on=0} on; $0 ~ start {on=1}' "$1" > block.dat
    grep -q '^# LOG MESH$' block.dat
    "$full" < block.dat > block.out
    cmp "$1" block.out
}

# expect_stop OUT MESSAGE: the output ends with MESSAGE, and has no results
expect_stop() {
    if ! grep -qF -- "$2" "$1"; then
        echo "expected '$2' in $1:"; tail -5 "$1"; exit 1
    fi
    if grep -q '^Reference configufation results' "$1" && [[ $2 != *'radial mesh'* ]]; then
        echo "$1 went on to the calculation"; exit 1
    fi
}

"$full" < "$dat" > base.out
case $case in
    defaults)
        with_mesh '0 -1 0' < "$dat" > in.dat
        "$full" < in.dat > mesh.out
        cmp base.out mesh.out ;;
    explicit)
        with_mesh '1.006 -1 45' < "$dat" > in.dat
        "$full" < in.dat > mesh.out
        grep -q '^# LOG MESH$' mesh.out
        cmp base.out <(drop_echo < mesh.out) ;;
    rerun)
        with_mesh "$mesh" < "$dat" > in.dat
        "$full" < in.dat > mesh.out
        if cmp -s <(grep -v '^build: ' base.out) <(drop_echo < mesh.out | grep -v '^build: '); then
            echo "the mesh input did not change the output"; exit 1
        fi
        awk '/^# ATOM AND REFERENCE CONFIGURATION$/ {on=1}
             /^Reference configufation results/ {on=0}
             on' mesh.out > echo.dat
        grep -q '^# LOG MESH$' echo.dat
        "$full" < echo.dat > echo.out
        cmp mesh.out echo.out
        if grep -q '^<INPUT>$' mesh.out; then
            rerun_block mesh.out '^<INPUT>$' '^</INPUT>$'
        fi ;;
    upf)
        with_mesh "$mesh" < "$dat" | sed -E '0,/psp8/s//upf/' > in.dat
        "$full" < in.dat > mesh.out
        rerun_block mesh.out '<PP_INPUTFILE>' '</PP_INPUTFILE>' ;;
    split)
        with_mesh "$mesh" < "$dat" > in.dat
        "$ae" -a atom.ae < in.dat > ae.out
        "$full" < in.dat > mesh.out
        "$ps" -a atom.ae < in.dat > ps.out
        cmp <(grep -v '^build: ' mesh.out) <(grep -v '^build: ' ps.out) ;;
    split_differ)
        with_mesh "$mesh" < "$dat" > in.dat
        "$ae" -a atom.ae < in.dat > ae.out
        "$ps" -a atom.ae < "$dat" > ps.out || true
        expect_stop ps.out 'read_ae_file: ERROR radial mesh' ;;
    toml)
        with_mesh "$mesh" < "$dat" > in.dat
        { cat "$toml"; printf '\n[log_mesh]\namesh = 1.004\nrmax = 50\n'; } > in.toml
        "$full" < in.dat > mesh.out
        "$full" -t in.toml > toml.out
        cmp mesh.out toml.out ;;
    amesh)
        with_mesh '0.999 -1 -1' < "$dat" > in.dat
        "$full" < in.dat > mesh.out || true
        expect_stop mesh.out 'must have amesh > 1.0' ;;
    rr1)
        with_mesh '-1 10.0 -1' < "$dat" > in.dat
        "$full" < in.dat > mesh.out || true
        expect_stop mesh.out 'must have rr1 < rcmax' ;;
    rmax)
        with_mesh '-1 -1 1.0' < "$dat" > in.dat
        "$full" < in.dat > mesh.out || true
        expect_stop mesh.out 'must have rmax > rlmax' ;;
    read)
        with_mesh '1.004 50.0' < "$dat" > in.dat
        "$full" < in.dat > mesh.out || true
        expect_stop mesh.out 'Read ERROR, input data file line' ;;
    *)
        echo "unknown case $case"; exit 1 ;;
esac
echo "log mesh $case: OK"
