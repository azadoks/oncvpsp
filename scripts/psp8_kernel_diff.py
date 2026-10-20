#!/usr/bin/env python3
"""Compares the pseudopotentials of two psp8 files by what a calculation uses.

The projectors and their energies of a psp8 file are not unique (any rotation
of the projectors of one l that keeps sum_i |beta_i> e_i <beta_i| leaves the
pseudopotential the same), so a change in them can be harmless.  This prints,
for each projector block (scalar-relativistic, then spin-orbit), the largest
difference of the nonlocal kernel

    K(r, r') = sum_i beta_i(r) e_i beta_i(r')

relative to its largest value, and the same for the local potential: the
invariant to look at when the reference outputs change.  The files can be
psp8 files or outputs that contain one (between "Begin PSPCODE8" and
"END_PSP").  With --tol, the exit status is 1 if any difference exceeds it.

usage: psp8_kernel_diff.py A B [--tol TOL]
"""

import argparse
import sys


def fortran_float(text):
    return float(text.replace("D", "E").replace("d", "e"))


def read_psp8(path):
    """[(l, energies, projectors)] of the file, projectors[k] the values at the
    k-th radial point; energies empty for the local potential (one column)."""
    lines = open(path).read().split("\n")
    start = next((k + 1 for k, line in enumerate(lines) if "Begin PSPCODE8" in line), 0)
    head = lines[start : start + 8]
    lmax, lloc, mmax = (int(x) for x in head[2].split()[2:5])
    extension = int(head[5].split()[0])
    pos = start + 6 + (1 if extension >= 2 else 0)
    nblock = (lmax + 1 - (1 if lloc <= lmax else 0)) + 1
    if extension >= 2:
        nblock += lmax
    blocks = []
    for _ in range(nblock):
        tokens = lines[pos].split()
        energies = [fortran_float(x) for x in tokens[1:]]
        values = [
            [fortran_float(x) for x in lines[pos + 1 + k].split()[2:]]
            for k in range(mmax)
        ]
        blocks.append((int(tokens[0]), energies, values))
        pos += mmax + 1
    return blocks


def kernel(energies, values):
    return [
        [
            sum(e * bi * bj for e, bi, bj in zip(energies, row_i, row_j))
            for row_j in values
        ]
        for row_i in values
    ]


def max_relative(a, b):
    diff = max(abs(x - y) for row_a, row_b in zip(a, b) for x, y in zip(row_a, row_b))
    scale = max(abs(x) for row in a for x in row)
    return diff / scale if scale > 0 else diff


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("a")
    parser.add_argument("b")
    parser.add_argument("--tol", type=float, help="fail above this difference")
    args = parser.parse_args()

    worst = 0.0
    for (l, ea, va), (_, eb, vb) in zip(read_psp8(args.a), read_psp8(args.b)):
        if not ea:
            diff = max_relative([[r[0]] for r in va], [[r[0]] for r in vb])
            print(f"  V_loc (l={l})              max rel diff {diff:.1e}")
        else:
            diff = max_relative(kernel(ea, va), kernel(eb, vb))
            proj = max_relative(va, vb)
            print(
                f"  l={l} projectors {len(ea)}  V_NL max rel diff {diff:.1e}"
                f"   (projectors themselves {proj:.1e})"
            )
        worst = max(worst, diff)
    if args.tol is not None and worst > args.tol:
        sys.exit(1)


if __name__ == "__main__":
    main()
