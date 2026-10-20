#!/bin/bash
# Checks the optional upf_mesh input, the third value of the output grid line
# (rlmax, drl, upf_mesh) of the text input or the upf_mesh key of [pp_output]
# in the TOML input, for one input with psfile upf or both:
#
#   defaults  upf_mesh 0: the output of the input without it
#   log       upf_mesh 3500: the output of the input without it up to the UPF
#             file, apart from the echo of the grid line; a UPF file on every
#             point of the log mesh (ratio amesh, rab = log(amesh) r) out to the
#             radius of the linear mesh, whose input, run again, gives the same
#             output
#   stride    upf_mesh 500: at most 500 points, every istr-th of the log mesh
#   toml      upf_mesh in [pp_output] of INPUT.toml: the output of the text
#             input
#   negative, read  a bad value: the program must stop with a message naming it
#
# usage: upf_mesh_test.sh FULL_EXE INPUT.dat CASE WORKDIR [INPUT.toml]
set -euo pipefail
full=$1; dat=$(realpath "$2"); case=$3; work=$4; toml=${5:+$(realpath "$5")}
mkdir -p "$work"
cd "$work"

# with_upf_mesh VALUE < IN > OUT: the input with VALUE after rlmax, drl
with_upf_mesh() {
    awk -v v="$1" 'done != 1 && prev ~ /^#.*rlmax/ {$0 = $0 "  " v; done = 1}
                   {print; prev = $0}'
}

# before_upf < OUT: the output up to the UPF file, without the build line and
# the echo of the grid line
before_upf() {
    awk '/^Begin PSP_UPF/ {exit} /^build: / {next}
         /^# rlmax, drl/ {skip = 2} skip > 0 {skip--; next} {print}'
}

# upf_block OUT TAG: the data of the UPF element TAG in OUT, one per line
upf_block() {
    awk -v tag="$2" '$0 ~ "</" tag ">" {on = 0} on {print} $0 ~ "<" tag " " {on = 1}' "$1" |
        tr -s ' ' '\n' | grep -v '^$'
}

# check_mesh OUT MAXPTS: PP_R is every istr-th point of a log mesh of ratio
# amesh, PP_RAB is istr log(amesh) r, mesh_size is at most MAXPTS
check_mesh() {
    upf_block "$1" PP_R > r.txt
    upf_block "$1" PP_RAB > rab.txt
    paste r.txt rab.txt | awk -v amesh=1.006 -v maxpts="$2" '
        {r[NR] = $1; rab[NR] = $2; n = NR}
        END {
            al = log(amesh); s = int(log(r[2] / r[1]) / al + 0.5)
            if (s < 1) bad("ratio " r[2] / r[1])
            for (i = 2; i <= n; i++)
                if (abs(log(r[i] / r[i - 1]) / al - s) > 1e-6) bad("ratio at " i)
            for (i = 1; i <= n; i++)
                if (abs(rab[i] / (r[i] * s * al) - 1) > 1e-9) bad("rab at " i)
            if (n > maxpts) bad(n " points")
            print n " points, stride " s
        }
        function abs(x) {return x < 0 ? -x : x}
        function bad(m) {print "bad mesh: " m; exit 1}'
    grep -q "mesh_size=\" *$(wc -l < r.txt)\"" "$1"
}

# rerun_block OUT START END: run the input between the lines START and END of
# OUT again, which must give OUT
rerun_block() {
    awk -v start="$2" -v end="$3" '$0 ~ end {on=0} on; $0 ~ start {on=1}' "$1" > block.dat
    grep -q 'rlmax, drl, upf_mesh' block.dat
    "$full" < block.dat > block.out
    cmp "$1" block.out
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
        cmp base.out mesh.out ;;
    log)
        with_upf_mesh 3500 < "$dat" > in.dat
        "$full" < in.dat > mesh.out
        grep -q '^# rlmax, drl, upf_mesh$' mesh.out
        cmp <(before_upf < base.out) <(before_upf < mesh.out)
        check_mesh mesh.out 3500
        # out to the radius of the linear mesh
        lin_end=$(upf_block base.out PP_R | tail -1)
        log_end=$(tail -1 r.txt)
        awk -v a="$lin_end" -v b="$log_end" 'BEGIN {exit !(b <= a && b > 0.99 * a)}'
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
        cmp mesh.out toml.out ;;
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
