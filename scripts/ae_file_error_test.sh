#!/bin/bash
# Checks how oncvpsp_ps.x handles one change to a valid all-electron data
# file: it must either stop with a non-zero status and the expected message
# on stdout, or, for a change within round-off, give the output of the
# single program.
#
# The same output is one that compare_output.py --identical accepts (every
# number to 1e-10, the text exactly, apart from the build line).
#
# usage: ae_file_error_test.sh FULL_EXE AE_EXE PS_EXE INPUT.dat CASE WORKDIR
# The python, with pyoncvpsp and upf-tools, is $ONCVPSP_PYTHON, or python3.
# $ONCVPSP_RUN and $ONCVPSP_AE_FILE, if set, are the output of FULL_EXE and
# the all-electron file of AE_EXE for INPUT.dat (run_test.sh), read in place
# of running them.
# shellcheck source=scripts/test_lib.sh
source "$(dirname "$0")/test_lib.sh"
full=$1; ae=$2; ps=$3; dat=$4; case=$5; work=$6
mkdir -p "$work"
cd "$work"

if [ -n "${ONCVPSP_AE_FILE:-}" ]; then
    cp "$ONCVPSP_AE_FILE" atom.ae
else
    "$ae" -a atom.ae < "$dat" > ae.out
fi

# patch FIELD FACTOR: multiply a real in the header of atom.ae by FACTOR, or
# with FIELD=version reverse the bytes of the format version. The header is
# magic(8 chars), version, irel, atsym(2 chars), zz, nc, nv, iexc,
# na(nc+nv), la(nc+nv), fa(nc+nv), mmax, rr(mmax), c_light in native
# byte order
patch() {
    "$python" "$scripts/ae_file_patch.py" atom.ae "$@"
}

expect=same
case $case in
    byte_order)  patch version;          expect='other byte order' ;;
    zz)          patch zz 1.001;         expect='atsym, z, nc, nv or iexc differ' ;;
    fa)          patch fa 1.001;         expect='reference configuration differs' ;;
    mesh)        patch rr 1.000001;      expect='radial mesh differs' ;;
    c_light)     patch c_light 1.001;    expect='c_light differs' ;;
    zz_roundoff) patch zz 1.00000000000001 ;;
    fa_roundoff) patch fa 1.00000000000001 ;;
    mesh_roundoff) patch rr 1.00000000000001 ;;
    *) echo "ae_file_error_test.sh: unknown case $case" >&2; exit 1 ;;
esac

status=0
"$ps" -a atom.ae < "$dat" > ps.out || status=$?

if [ "$expect" = same ]; then
    if [ $status -ne 0 ]; then
        echo "case $case: oncvpsp_ps.x stopped with status $status" >&2
        tail -3 ps.out >&2
        exit 1
    fi
    if [ -n "${ONCVPSP_RUN:-}" ]; then
        cp "$ONCVPSP_RUN" full.out
    else
        "$full" < "$dat" > full.out
    fi
    same full.out ps.out
    echo "case $case: the same output as the single program"
else
    if [ $status -eq 0 ]; then
        echo "case $case: oncvpsp_ps.x did not stop" >&2
        exit 1
    fi
    if ! grep -q "read_ae_file: ERROR.*$expect" ps.out; then
        echo "case $case: no message matching '$expect'; got:" >&2
        tail -3 ps.out >&2
        exit 1
    fi
    echo "case $case: stopped with '$expect'"
fi
