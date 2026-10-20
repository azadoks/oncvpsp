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
set -euo pipefail
full=$1; ae=$2; ps=$3; dat=$4; case=$5; work=$6
scripts=$(realpath "$(dirname "$0")")
python=${ONCVPSP_PYTHON:-python3}
mkdir -p "$work"
cd "$work"

"$ae" -a atom.ae < "$dat" > ae.out

# patch FIELD FACTOR: multiply a real in the header of atom.ae by FACTOR, or
# with FIELD=version reverse the bytes of the format version. The header is
# magic(8 chars), version, irel, atsym(2 chars), zz, nc, nv, iexc,
# na(nc+nv), la(nc+nv), fa(nc+nv), mmax, rr(mmax), c_light in native
# byte order
patch() {
    python3 - atom.ae "$@" <<'PY'
import struct
import sys

path, field = sys.argv[1:3]
with open(path, "rb") as f:
    buf = bytearray(f.read())
zz = 8 + 4 + 4 + 2
nc, nv = struct.unpack_from("=ii", buf, zz + 8)
nt = nc + nv
fa = zz + 8 + 12 + 8 * nt
rr = fa + 8 * nt + 4
(mmax,) = struct.unpack_from("=i", buf, rr - 4)
at = {"zz": zz, "fa": fa + 8 * (nt - 1), "rr": rr + 8, "c_light": rr + 8 * mmax}
if field == "version":
    buf[8:12] = buf[8:12][::-1]
else:
    (x,) = struct.unpack_from("=d", buf, at[field])
    struct.pack_into("=d", buf, at[field], x * float(sys.argv[3]))
with open(path, "wb") as f:
    f.write(buf)
PY
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
    "$full" < "$dat" > full.out
    if ! "$python" "$scripts/compare_output.py" --identical full.out ps.out > ps.diff; then
        cat ps.diff
        exit 1
    fi
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
