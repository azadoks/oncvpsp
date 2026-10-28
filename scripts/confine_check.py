"""The checks of confine_test.sh on the outputs of oncvpsp_confine.x in the
current directory: conf.out and full.out, then r1.out, r2.out, r12.out,
l1.out, step.out and wide.out, or ghost.out with MODE ghost.

usage: confine_check.py [MODE]

Relative differences are |a - b| / |b|, b the expected value.
"""

import math
import sys

from scripts_common import fortran_float

# the summary line against the lines it summarizes: both printed to 1e-12
# relative (max_dde absolute, from changes of total energies near 1e-16)
SUMMARY_TOL = 1e-12
# a band edge below the vacuum level, pseudopotential against all-electron
# (Ha): the accuracy the pseudopotential is made for
EDGE_TOL = 1e-4
# the error at the largest cutoff of an edge below the vacuum level (Ha)
FLOOR_TOL = 1e-7
# the f ghost of 60_Nd_GHOST, from the run's bound-ghost test (Ha)
GHOST_E, GHOST_TOL = -50.06, 0.05
# the same edges from two radii and from each alone (Ha): round-off
RADII_TOL = 1e-12
# the same edges from another scan step or window (Ha): the bisection's
# convergence
SCAN_TOL = 1e-11
# the free atoms against the run's reference energies, relative
FREE_TOL = 1e-8
# the energy changes in a wall, pseudo-atom against all-electron atom (Ha)
WALL_TOL = 1e-4


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
        sys.exit(f"confine_check.py: {what}: different channels or radii")
    worst, count = 0.0, 0
    for k in a:
        ea, eb = a[k]["ws"], b[k]["ws"]
        if window:
            ea = {e: v for e, v in ea.items() if window[0] < v[0] < window[1]}
            eb = {e: v for e, v in eb.items() if window[0] < v[0] < window[1]}
        if ea.keys() != eb.keys():
            sys.exit(f"confine_check.py: {what}: different edges at {k}")
        for e in ea:
            worst = max(worst, abs(ea[e][0] - eb[e][0]), abs(ea[e][1] - eb[e][1]))
            count += 1
    if worst > tol:
        sys.exit(f"confine_check.py: {what}: differ by {worst:.2e}")
    print(f"{what}: {count} edges, max difference {worst:.1e}")


conf = blocks("conf.out")
if not conf or not any(b["ws"] for b in conf.values()):
    sys.exit("confine_check.py: no band edges in the output")

# the summary from the other lines
lines = [l.split() for l in open("conf.out") if l.strip() and not l.startswith("#")]
summary = [w for w in lines if w[0] == "confine_summary"]
if len(summary) != 1 or lines[-1][0] != "confine_summary":
    sys.exit("confine_check.py: not one summary line, last")
w = summary[0]
got = (
    [float(x) for x in w[1:4]] + [int(x) for x in w[4:8]] + [float(x) for x in w[8:11]]
)


def count(kind):
    return sum(1 for l in lines if l[0] == kind)


neg = {(k, e) for k, b in conf.items() for e, v in b["ws"].items() if v[0] < 0.0}


def worst_cut(tol, top=False):
    cuts = []
    for k, e in neg:
        if top and e[0] != 3:
            continue
        c = [c for t, c in conf[k]["ecut"].get(e, []) if t == tol]
        cuts.append(c[0] if c else math.inf)
    return max(cuts, default=math.nan)


kwalls = {float(l[1]): (float(l[4]), float(l[5])) for l in lines if l[0] == "wall"}
want = [
    max(
        (abs(conf[k]["ws"][e][1] - conf[k]["ws"][e][0]) for k, e in neg),
        default=math.nan,
    ),
    worst_cut(1e-4),
    worst_cut(1e-5),
    count("ws_missing"),
    count("ws_extra"),
    count("ws_paired_off"),
    count("wall_failed"),
    max((abs(p - a) for k, (a, p) in kwalls.items() if k > 0), default=math.nan),
    worst_cut(1e-4, True),
    worst_cut(1e-5, True),
]
for name, g, x in zip(
    (
        "max_de",
        "ecut4",
        "ecut5",
        "n_missing",
        "n_extra",
        "n_paired_off",
        "n_wall_failed",
        "max_dde",
        "ecut4_top",
        "ecut5_top",
    ),
    got,
    want,
):
    # max_dde from the printed changes, each rounded at the 1e-16 of a total energy
    tol = SUMMARY_TOL if name == "max_dde" else SUMMARY_TOL * abs(x)
    if not (g == x or (math.isnan(g) and math.isnan(x)) or abs(g - x) <= tol):
        sys.exit(f"confine_check.py: summary {name} = {g}, from the lines {x}")
print(f"summary: {' '.join(w[1:])}")

if sys.argv[1:] == ["ghost"]:
    core = [l for l in lines if l[0] == "ws_core"]
    extra = [l for l in lines if l[0] == "ws_extra"]
    if not core or any(l[1] != "3" or l[4] != "4" for l in core):
        sys.exit(f"confine_check.py: the core 4f band is not ws_core: {core}")
    if extra or count("ws_missing"):
        sys.exit(
            f"confine_check.py: pseudopotential-only or missing edges in the window: {extra}"
        )
    ghost = [float(l.split()[4]) for l in open("ghost.out") if l.startswith("ws_extra")]
    if not ghost or any(abs(e - GHOST_E) > GHOST_TOL for e in ghost):
        sys.exit(f"confine_check.py: the f ghost near -50 Ha is not ws_extra: {ghost}")
    print(
        f"{len(core)} core-state edges, none pseudopotential-only; the f ghost at {ghost[0]:.4f} Ha"
    )
    sys.exit(0)
for k, b in conf.items():
    if b["missing"]:
        sys.exit(
            f"confine_check.py: edges with no pseudopotential partner at {k}: {b['missing']}"
        )
    for e, (eae, eps) in b["ws"].items():
        # the bands below the vacuum level are those the pseudopotential is made for
        if eae < 0.0 and abs(eps - eae) > EDGE_TOL:
            sys.exit(
                f"confine_check.py: edge {e} at {k}: e_ps - e_ae = {eps - eae:.2e}"
            )
    for e, rows in b["ecut"].items():
        cuts = [c for _, c in sorted(rows, reverse=True)]
        if any(c2 < c1 for c1, c2 in zip(cuts, cuts[1:])):
            sys.exit(
                f"confine_check.py: edge {e} at {k}: cutoffs fall as the error falls: {cuts}"
            )
    for e, (_, err) in b["floor"].items():
        if b["ws"][e][0] < 0.0 and abs(err) > FLOOR_TOL:
            sys.exit(
                f"confine_check.py: edge {e} at {k}: error {err:.2e} at the largest cutoff"
            )
print(
    f"{sum(len(b['ws']) for b in conf.values())} edges, those below 0 within {EDGE_TOL:.0e} Ha"
    f" of the all-electron atom's and converged to {FLOOR_TOL:.0e} Ha"
)

r12 = blocks("r12.out")
compare(
    r12,
    {**blocks("r1.out"), **blocks("r2.out")},
    RADII_TOL,
    "two radii against each alone",
)
compare(
    blocks("l1.out"),
    {k: v for k, v in r12.items() if k[0] == 1},
    0.0,
    "l = 1 against every l",
)
compare(
    blocks("step.out"),
    r12,
    SCAN_TOL,
    "a scan step of 0.137 Ha against the input's",
    (-1.9, 1.9),
)
compare(
    blocks("wide.out"),
    r12,
    SCAN_TOL,
    "a window of -12..13 Ha against the input's",
    (-1.9, 1.9),
)

# the atom in walls: the free atoms are the run's reference, the energy rises with the wall,
# and the pseudo-atom's rise is the all-electron atom's
walls = {}
eigs = {}
for line in open("conf.out"):
    w = line.split()
    if w and w[0] == "wall":
        walls[float(w[1])] = [float(x) for x in w[2:6]]
    elif w and w[0] == "wall_eig":
        eigs.setdefault(float(w[1]), []).append((float(w[5]), float(w[6])))
    elif w and w[0] == "wall_failed":
        sys.exit(f"confine_check.py: no self-consistent atom in a wall: {line.strip()}")
if sorted(walls) != [0.0, 0.01, 0.03, 0.1] or sorted(eigs) != sorted(walls):
    sys.exit(f"confine_check.py: walls {sorted(walls)}, eigenvalues for {sorted(eigs)}")
ref = {}
for line in open("full.out"):
    for key in ("AE_ref=", "PS_ref="):
        if key in line and key not in ref:
            ref[key] = fortran_float(line.split(key)[1].split()[0])
etae, etps = walls[0.0][:2]
if abs(etae - ref["AE_ref="]) > FREE_TOL * abs(ref["AE_ref="]) or abs(
    etps - ref["PS_ref="]
) > FREE_TOL * abs(ref["PS_ref="]):
    sys.exit(f"confine_check.py: free atoms {etae} {etps} against the run's {ref}")
ks = sorted(walls)
for k1, k2 in zip(ks, ks[1:]):
    if not (walls[k2][2] > walls[k1][2] and walls[k2][3] > walls[k1][3]):
        sys.exit(f"confine_check.py: the energy does not rise from wall {k1} to {k2}")
worst = max(abs(v[3] - v[2]) for v in walls.values())
if worst > WALL_TOL:
    sys.exit(
        f"confine_check.py: the PS and AE energy changes in the walls differ by {worst:.2e}"
    )
print(
    f"walls {ks[1:]}: free atoms the run's, energy changes AE - PS within {worst:.1e} Ha"
)
