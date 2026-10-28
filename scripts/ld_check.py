"""The checks of ld_test.sh on the outputs in the current directory: the log
derivatives of oncvpsp_ld.x (ld.out) against those of the run's DIAGNOSTICS
block (full.out), two radii (r12.out) against each alone (r1.out, r2.out),
and l = 1 alone (l1.out) against every l.

usage: ld_check.py
"""

import sys

# oncvpsp_ld.x against the run, which compute the phases from the same
# potentials read back from the files: their 17-digit printing
LD_TOL = 1e-10
# two radii in one run against each alone: round-off
RADII_TOL = 1e-12


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
        sys.exit(f"ld_check.py: {what}: different channels or radii")
    worst = max(
        (abs(x - y) for k in a for ra, rb in zip(a[k], b[k]) for x, y in zip(ra, rb)),
        default=0.0,
    )
    if any(len(a[k]) != len(b[k]) for k in a) or worst > tol:
        sys.exit(f"ld_check.py: {what}: differ by {worst:.2e}")
    print(f"{what}: {sum(map(len, a.values()))} energies, max difference {worst:.1e}")


full = blocks("full.out")
if not full:
    sys.exit("ld_check.py: no log derivatives in the output of the run")
compare(blocks("ld.out"), full, LD_TOL, "oncvpsp_ld.x against the run")
r12 = blocks("r12.out")
alone = {**blocks("r1.out"), **blocks("r2.out")}
compare(r12, alone, RADII_TOL, "two radii against each alone")
compare(
    blocks("l1.out"),
    {k: v for k, v in r12.items() if k[0] == 1},
    0.0,
    "l = 1 against every l",
)
