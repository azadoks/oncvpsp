#!/usr/bin/env bash
# Wigner-Seitz sphere of oncvpsp_confine.x from the all-electron and pseudopotential files: the band
# edges of the valence bands agree with the all-electron atom's, each radius gives the same as alone,
# the edges do not depend on the step of the scan that brackets them, and the cutoffs rise as the
# error falls.
#
# Usage: confine_test.sh FULL_EXE AE_EXE CONFINE_EXE INPUT.dat WORKDIR
set -euo pipefail
full=$1; ae=$2; confine=$3; dat=$4; work=$5
python=${ONCVPSP_PYTHON:-python3}
mkdir -p "$work"
cd "$work"

if [ -n "${ONCVPSP_AE_FILE:-}" ]; then
    cp "$ONCVPSP_AE_FILE" atom.ae
else
    "$ae" -a atom.ae < "$dat" > ae.out
fi
"$full" --no-log-derivatives --ps-file atom.ps < "$dat" > full.out

"$confine" -a atom.ae -p atom.ps > conf.out
read -r r1 r2 < <(awk '$1 == "ws_r" {print $4}' conf.out | sort -g -u | head -2 | xargs)
"$confine" -a atom.ae -p atom.ps -r "$r1" > r1.out
"$confine" -a atom.ae -p atom.ps -r "$r2" > r2.out
"$confine" -a atom.ae -p atom.ps -r "$r1,$r2" > r12.out
"$confine" -a atom.ae -p atom.ps -r "$r1,$r2" -l 1 > l1.out
"$confine" -a atom.ae -p atom.ps -r "$r1,$r2" -e -2.0,2.0,0.137 > step.out

"$python" - <<'EOF'
import sys

def blocks(path):
    """per (l, ikap, rs): the edges {(edge, n): (e_ae, e_ps)}, the cutoffs
    {(edge, n): [(tol, ecut)]}, the floors and the missing edges"""
    out, key = {}, None
    for line in open(path):
        w = line.split()
        if not w or w[0].startswith("#"):
            continue
        if w[0] == "ws_r":
            key = (int(w[1]), int(w[2]), float(w[3]))
            out[key] = {"ws": {}, "ecut": {}, "floor": {}, "missing": []}
        elif w[0] in ("ws", "ws_ecut", "ws_floor"):
            edge = (int(w[3]), int(w[4]))
            vals = tuple(float(x) for x in w[5:7])
            kind = w[0][3:] or "ws"
            if kind == "ecut":
                out[key]["ecut"].setdefault(edge, []).append(vals)
            else:
                out[key][kind][edge] = vals
        elif w[0] == "ws_missing":
            out[key]["missing"].append(line.strip())
    return out

def compare(a, b, tol, what, window=None):
    if a.keys() != b.keys():
        sys.exit(f"confine_test.sh: {what}: different channels or radii")
    worst, count = 0.0, 0
    for k in a:
        ea, eb = a[k]["ws"], b[k]["ws"]
        if window:
            ea = {e: v for e, v in ea.items() if window[0] < v[0] < window[1]}
            eb = {e: v for e, v in eb.items() if window[0] < v[0] < window[1]}
        if ea.keys() != eb.keys():
            sys.exit(f"confine_test.sh: {what}: different edges at {k}")
        for e in ea:
            worst = max(worst, abs(ea[e][0] - eb[e][0]), abs(ea[e][1] - eb[e][1]))
            count += 1
    if worst > tol:
        sys.exit(f"confine_test.sh: {what}: differ by {worst:.2e}")
    print(f"{what}: {count} edges, max difference {worst:.1e}")

conf = blocks("conf.out")
if not conf or not any(b["ws"] for b in conf.values()):
    sys.exit("confine_test.sh: no band edges in the output")
for k, b in conf.items():
    if b["missing"]:
        sys.exit(f"confine_test.sh: edges with no pseudopotential partner at {k}: {b['missing']}")
    for e, (eae, eps) in b["ws"].items():
        # the bands below the vacuum level are those the pseudopotential is made for
        if eae < 0.0 and abs(eps - eae) > 1e-4:
            sys.exit(f"confine_test.sh: edge {e} at {k}: e_ps - e_ae = {eps - eae:.2e}")
    for e, rows in b["ecut"].items():
        cuts = [c for _, c in sorted(rows, reverse=True)]
        if any(c2 < c1 for c1, c2 in zip(cuts, cuts[1:])):
            sys.exit(f"confine_test.sh: edge {e} at {k}: cutoffs fall as the error falls: {cuts}")
    for e, (_, err) in b["floor"].items():
        if b["ws"][e][0] < 0.0 and abs(err) > 1e-7:
            sys.exit(f"confine_test.sh: edge {e} at {k}: error {err:.2e} at the largest cutoff")
print(f"{sum(len(b['ws']) for b in conf.values())} edges, those below 0 within 1e-4 Ha of the"
      " all-electron atom's and converged to 1e-7 Ha")

r12 = blocks("r12.out")
compare(r12, {**blocks("r1.out"), **blocks("r2.out")}, 1e-12, "two radii against each alone")
compare(blocks("l1.out"), {k: v for k, v in r12.items() if k[0] == 1}, 0.0, "l = 1 against every l")
compare(blocks("step.out"), r12, 1e-11, "a scan step of 0.137 Ha against the input's", (-1.9, 1.9))
EOF
