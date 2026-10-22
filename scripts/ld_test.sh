#!/bin/bash
# Checks oncvpsp_ld.x for one input: FULL_EXE --no-log-derivatives --ps-file
# writes the pseudopotential file, and must give the output of FULL_EXE but
# the log derivatives; oncvpsp_ld.x with that file and the all-electron file
# of AE_EXE must then give the log derivatives of FULL_EXE's DIAGNOSTICS
# block (every number to 1e-10), and the phases of each radius of a list of
# radii (one inside the core radii, one beyond) must be those of that radius
# alone (to 1e-12), and the channels of -l those of every channel.
#
# usage: ld_test.sh FULL_EXE AE_EXE LD_EXE INPUT.dat WORKDIR
# The python is $ONCVPSP_PYTHON, or python3. $ONCVPSP_RUN and
# $ONCVPSP_AE_FILE, if set, are the output of FULL_EXE and the all-electron
# file of AE_EXE for INPUT.dat (run_test.sh), read in place of running them.
set -euo pipefail
full=$1; ae=$2; ld=$3; dat=$4; work=$5
python=${ONCVPSP_PYTHON:-python3}
mkdir -p "$work"
cd "$work"

if [ -n "${ONCVPSP_AE_FILE:-}" ]; then
    cp "$ONCVPSP_AE_FILE" atom.ae
else
    "$ae" -a atom.ae < "$dat" > ae.out
fi
if [ -n "${ONCVPSP_RUN:-}" ]; then
    cp "$ONCVPSP_RUN" full.out
else
    "$full" < "$dat" > full.out
fi
"$full" --no-log-derivatives --ps-file atom.ps < "$dat" > skip.out

# the output without the log derivatives: their DIAGNOSTICS lines and the
# blocks for plotting, and the blank lines between those
no_ld() {
    grep -vE "^(ld|ld_r) |^# ld|^! |^log derivativve|^atan\(r|^l, energy|^build: " "$1" |
        grep -v '^[[:space:]]*$'
}
if ! diff <(no_ld full.out) <(no_ld skip.out) > skip.diff; then
    echo "ld_test.sh: --no-log-derivatives changed more than the log derivatives" >&2
    cat skip.diff >&2
    exit 1
fi
if grep -qE "^(ld|ld_r) " skip.out; then
    echo "ld_test.sh: --no-log-derivatives left log derivatives" >&2
    exit 1
fi

"$ld" -a atom.ae -p atom.ps > ld.out
"$ld" -a atom.ae -p atom.ps -r 1.0 > r1.out
"$ld" -a atom.ae -p atom.ps -r 3.5 > r2.out
"$ld" -a atom.ae -p atom.ps -r 1.0,3.5 > r12.out
"$ld" -a atom.ae -p atom.ps -r 1.0,3.5 -l 1 > l1.out

"$python" - <<'EOF'
import sys

def blocks(path):
    """the phases of each (l, ikap, radius): rows of (e, ae, ps)"""
    out, key = {}, None
    for line in open(path):
        w = line.split()
        if not w:
            continue
        if w[0] == "ld_r":
            key = (int(w[1]), int(w[2]), float(w[3]))
            out[key] = []
        elif w[0] == "ld":
            out[key].append([float(x) for x in w[3:6]])
    return out

def compare(a, b, tol, what):
    if a.keys() != b.keys():
        sys.exit(f"ld_test.sh: {what}: different channels or radii")
    worst = max((abs(x - y) for k in a for ra, rb in zip(a[k], b[k])
                 for x, y in zip(ra, rb)), default=0.0)
    if any(len(a[k]) != len(b[k]) for k in a) or worst > tol:
        sys.exit(f"ld_test.sh: {what}: differ by {worst:.2e}")
    print(f"{what}: {sum(map(len, a.values()))} energies, max difference {worst:.1e}")

full = blocks("full.out")
if not full:
    sys.exit("ld_test.sh: no log derivatives in the output of the run")
compare(blocks("ld.out"), full, 1e-10, "oncvpsp_ld.x against the run")
r12 = blocks("r12.out")
alone = {**blocks("r1.out"), **blocks("r2.out")}
compare(r12, alone, 1e-12, "two radii against each alone")
compare(blocks("l1.out"), {k: v for k, v in r12.items() if k[0] == 1}, 0.0,
        "l = 1 against every l")
EOF
