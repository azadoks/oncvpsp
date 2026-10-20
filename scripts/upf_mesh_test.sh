#!/bin/bash
# Checks the optional upf_mesh input, the third value of the output grid line
# (rlmax, drl, upf_mesh) of the text input or the upf_mesh key of [pp_output]
# in the TOML input, for one input with psfile upf or both:
#
#   defaults  upf_mesh 0: the output of the input without it
#   log       upf_mesh 3500: the output of the input without it apart from the
#             UPF file and the echo of the grid line; a UPF file on every point
#             of the log mesh (ratio amesh, rab = log(amesh) r) out to the
#             radius of the linear mesh, whose input, run again, gives the same
#             output
#   stride    upf_mesh 500: at most 500 points, every istr-th of the log mesh
#   toml      upf_mesh in [pp_output] of INPUT.toml: the output of the text
#             input
#   negative, read  a bad value: the program must stop with a message naming it
#
# The same output is one that compare_output.py --identical accepts (every
# number to 1e-10, the text exactly); the UPF mesh is checked by
# upf_mesh_check.py.
#
# usage: upf_mesh_test.sh FULL_EXE INPUT.dat CASE WORKDIR [INPUT.toml]
# The python, with pyoncvpsp and upf-tools, is $ONCVPSP_PYTHON, or python3.
set -euo pipefail
full=$1; dat=$(realpath "$2"); case=$3; work=$4; toml=${5:+$(realpath "$5")}
scripts=$(realpath "$(dirname "$0")")
python=${ONCVPSP_PYTHON:-python3}
mkdir -p "$work"
cd "$work"

# with_upf_mesh VALUE < IN > OUT: the input with VALUE after rlmax, drl
with_upf_mesh() {
    awk -v v="$1" 'done != 1 && prev ~ /^#.*rlmax/ {$0 = $0 "  " v; done = 1}
                   {print; prev = $0}'
}

# same A.out B.out [ARGS...]: B is the same output as A (compare_output.py
# --identical, with ARGS)
same() {
    if ! "$python" "$scripts/compare_output.py" --identical "${@:3}" "$1" "$2" \
            > "${2%.out}.diff"; then
        cat "${2%.out}.diff"
        exit 1
    fi
}

# check_mesh OUT MAXPTS [--linear BASE.out]: the UPF file of OUT is on the
# log mesh (upf_mesh_check.py)
check_mesh() {
    "$python" "$scripts/upf_mesh_check.py" "$1" 1.006 "${@:2}"
}

# rerun_block OUT START END: run the input between the lines START and END of
# OUT again, which must give OUT
rerun_block() {
    awk -v start="$2" -v end="$3" '$0 ~ end {on=0} on; $0 ~ start {on=1}' "$1" > block.dat
    grep -q 'rlmax, drl, upf_mesh' block.dat
    "$full" < block.dat > block.out
    same "$1" block.out
}

# expect_stop OUT MESSAGE: the output ends with MESSAGE, and has no results
expect_stop() {
    if ! grep -qF -- "$2" "$1"; then
        echo "expected '$2' in $1:"; tail -5 "$1"; exit 1
    fi
    if grep -q '^Reference configufation results' "$1"; then
        echo "$1 went on to the calculation"; exit 1
    fi
}

"$full" < "$dat" > base.out
case $case in
    defaults)
        with_upf_mesh 0 < "$dat" > in.dat
        "$full" < in.dat > mesh.out
        same base.out mesh.out ;;
    log)
        with_upf_mesh 3500 < "$dat" > in.dat
        "$full" < in.dat > mesh.out
        grep -q '^# rlmax, drl, upf_mesh$' mesh.out
        same base.out mesh.out --no-lines --ignore upf --ignore input.output_grid.upf_mesh
        # out to the radius of the linear mesh
        check_mesh mesh.out 3500 --linear base.out
        rerun_block mesh.out '<PP_INPUTFILE>' '</PP_INPUTFILE>' ;;
    stride)
        with_upf_mesh 500 < "$dat" > in.dat
        "$full" < in.dat > mesh.out
        check_mesh mesh.out 500 | grep -qv 'stride 1$' ;;
    toml)
        with_upf_mesh 500 < "$dat" > in.dat
        awk '{print} /^\[pp_output\]/ {print "upf_mesh = 500"}' "$toml" > in.toml
        "$full" < in.dat > mesh.out
        "$full" -t in.toml > toml.out
        same mesh.out toml.out ;;
    negative)
        with_upf_mesh -1 < "$dat" > in.dat
        "$full" < in.dat > mesh.out || true
        expect_stop mesh.out 'must have upf_mesh >= 0' ;;
    read)
        with_upf_mesh 2.5 < "$dat" > in.dat
        "$full" < in.dat > mesh.out || true
        expect_stop mesh.out 'Read ERROR, input data file line' ;;
    *)
        echo "unknown case $case"; exit 1 ;;
esac
echo "upf mesh $case: OK"
