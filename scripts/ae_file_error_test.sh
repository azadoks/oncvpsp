#!/bin/bash
# Checks how oncvpsp_ps.x handles one change to a valid all-electron data
# file: it must either stop with a non-zero status and the expected message
# on stdout, or, for a change within round-off, give the output of the
# single program.
#
# usage: ae_file_error_test.sh FULL_EXE AE_EXE PS_EXE INPUT.dat CASE WORKDIR
set -euo pipefail
full=$1; ae=$2; ps=$3; dat=$4; case=$5; work=$6
mkdir -p "$work"
cd "$work"

"$ae" -a atom.ae < "$dat" > ae.out

# patch FIELD FACTOR: multiply a real in the header of atom.ae by FACTOR, or
# with FIELD=version reverse the bytes of the format version. The header is
# magic(8 chars), version, irel, atsym(2 chars), zz, nc, nv, iexc,
# na(nc+nv), la(nc+nv), fa(nc+nv), mmax, rr(mmax), c_light in native
# byte order
patch() {
    perl -e '
        my ($field, $factor) = @ARGV[1, 2];
        open(my $fh, "+<:raw", $ARGV[0]) or die "patch: $!";
        local $/; my $buf = <$fh>;
        my $zz = 8 + 4 + 4 + 2;
        my ($nc, $nv) = unpack("l l", substr($buf, $zz + 8, 8));
        my $nt = $nc + $nv;
        my $fa = $zz + 8 + 12 + 8 * $nt;
        my $rr = $fa + 8 * $nt + 4;
        my $mmax = unpack("l", substr($buf, $rr - 4, 4));
        my %at = (zz => $zz, fa => $fa + 8 * ($nt - 1), rr => $rr + 8,
                  c_light => $rr + 8 * $mmax);
        if ($field eq "version") {
            substr($buf, 8, 4) = reverse substr($buf, 8, 4);
        } else {
            my $off = $at{$field} // die "patch: no field $field";
            my $x = unpack("d", substr($buf, $off, 8));
            substr($buf, $off, 8) = pack("d", $x * $factor);
        }
        seek($fh, 0, 0); print $fh $buf; close($fh);
    ' atom.ae "$@"
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
    cmp <(grep -v '^build: ' full.out) <(grep -v '^build: ' ps.out)
    echo "case $case: output identical to the single program"
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
