#!/usr/bin/env python3
"""Compares the pseudopotentials of two outputs by what a calculation uses.

The projectors and their energies are not unique (any rotation of the
projectors of one l, or one l and j, that keeps sum_ij |beta_i> D_ij <beta_j|
leaves the pseudopotential the same), so a change in them can be harmless, and
they are never compared directly.  This prints, for each projector block, the
difference of the nonlocal kernel

    K(r, r') = sum_ij beta_i(r) D_ij beta_j(r')

in the Frobenius norm on the file's radial grid, relative to the norm of the
first file's kernel (for a spin-orbit kernel of the psp8 file, of the
scalar-relativistic kernel of its l, to which it is a correction, as it can be
small; computed from the overlaps of the projectors, it resolves differences
down to about 1e-8), and the largest difference of the local potential relative
to its largest value.  The blocks are those of the psp8 file (scalar-relativistic,
then spin-orbit) and of the UPF file (one per l, or per l and j) that the
outputs contain (between "Begin PSPCODE8" and "END_PSP", and "<UPF" and
"</UPF>"); the files themselves can be given too.  With --tol, the exit status
is 1 if any difference exceeds it, or the blocks do not match.

It needs only the standard library, and reads bare psp8 and UPF files too;
scripts/compare_output.py, which the tests use, compares the kernels in the
same way along with everything else in two outputs.

usage: kernel_diff.py A B [--tol TOL]
"""

import argparse
import math
import re
import sys

from scripts_common import fortran_float


def find(lines, marker, start=0):
    return next((k for k in range(start, len(lines)) if marker in lines[k]), None)


def psp8_blocks(lines):
    """[(name, D, projectors, first line, last line)] of the psp8 file in
    lines (projectors[i] the values of beta_i on the grid; D None for the local
    potential), [] if there is none."""
    begin = find(lines, "Begin PSPCODE8")
    if begin is None:
        return []
    head = lines[begin + 1 : begin + 8]
    lmax, lloc, mmax = (int(x) for x in head[2].split()[2:5])
    nproj = [int(x) for x in head[4].split()[: lmax + 1]]
    extension = int(head[5].split()[0])
    pos = begin + 7
    nprojso = [0] * (lmax + 1)
    if extension >= 2:
        nprojso[1:] = [int(x) for x in lines[pos].split()[:lmax]]
        pos += 1
    blocks = []

    def block(name, nvalue, local):
        nonlocal pos
        tokens = lines[pos].split()
        if int(tokens[0]) != int(name.split("=")[1]):
            raise ValueError(f"psp8 block {name} not at line {pos + 1}")
        rows = [lines[pos + 1 + k].split()[2 : 2 + nvalue] for k in range(mmax)]
        values = [[fortran_float(row[i]) for row in rows] for i in range(nvalue)]
        if local:
            matrix = None
        else:
            energies = [fortran_float(x) for x in tokens[1 : 1 + nvalue]]
            matrix = [
                [e if i == j else 0.0 for j in range(nvalue)]
                for i, e in enumerate(energies)
            ]
        blocks.append((name, matrix, values, pos, pos + mmax))
        pos += mmax + 1

    for ll in range(lmax + 1):
        if ll == lloc:
            block(f"psp8 V_loc l={ll}", 1, True)
        elif nproj[ll] > 0:
            block(f"psp8 l={ll}", nproj[ll], False)
    if lloc > lmax:
        block(f"psp8 V_loc l={lloc}", 1, True)
    for ll in range(1, lmax + 1):
        if nprojso[ll] > 0:
            block(f"psp8 spin-orbit l={ll}", nprojso[ll], False)
    return blocks


UPF_NUMBERS = re.compile(r'(\w+)="\s*([^"]*?)\s*"')


def upf_array(lines, start):
    """The numbers of the UPF element starting at line start, and the line of
    its end tag."""
    tag = lines[start].split()[0][1:]
    end = find(lines, f"</{tag}", start)
    text = " ".join(lines[start : end + 1])
    body = text[text.index(">") + 1 : text.rindex("<")]
    return [fortran_float(x) for x in body.split()], end


def upf_blocks(lines):
    """[(name, D, projectors, first line, last line)] of the UPF file in lines,
    the projectors grouped by l (and j), D the block of PP_DIJ of each group,
    and the local potential; [] if there is none."""
    begin = find(lines, "<UPF")
    if begin is None:
        return []
    end = find(lines, "</UPF>", begin)
    betas = {}
    dij = None
    relbeta = {}
    local = None
    k = begin
    while k < end:
        line = lines[k].lstrip()
        if line.startswith("<PP_BETA."):
            head_end = find(lines, ">", k)
            attrs = dict(UPF_NUMBERS.findall(" ".join(lines[k : head_end + 1])))
            values, last = upf_array(lines, k)
            betas[int(attrs["index"])] = (
                int(attrs["angular_momentum"]),
                values,
                k,
                last,
            )
            k = last
        elif line.startswith("<PP_DIJ"):
            dij, last = upf_array(lines, k)
            dij_lines = (k, last)
            k = last
        elif line.startswith("<PP_LOCAL"):
            local, last = upf_array(lines, k)
            local_lines = (k, last)
            k = last
        elif line.startswith("<PP_RELBETA."):
            attrs = dict(UPF_NUMBERS.findall(line))
            relbeta[int(attrs["index"])] = attrs["jjj"]
        k += 1
    blocks = []
    if local is not None:
        blocks.append(("UPF V_loc", None, [local], *local_lines))
    nbeta = len(betas)
    groups = {}
    for i in sorted(betas):
        key = (betas[i][0], relbeta.get(i))
        groups.setdefault(key, []).append(i)
    for (ll, jj), members in groups.items():
        name = f"UPF l={ll}" + (f" j={jj}" if jj is not None else "")
        matrix = [[dij[(i - 1) * nbeta + (j - 1)] for j in members] for i in members]
        first = min(betas[i][2] for i in members)
        last = max(betas[i][3] for i in members)
        blocks.append((name, matrix, [betas[i][1] for i in members], first, last))
    if nbeta:
        blocks.append(("UPF D_ij", None, None, *dij_lines))
    return blocks


def kernel_dot(da, pa, db, pb):
    """tr(K_a K_b) = sum_ijkl Da_ij Db_kl (a_i . b_k) (a_j . b_l), on the grid."""
    n = min(len(pa[0]), len(pb[0]))
    overlap = [
        [math.fsum(x * y for x, y in zip(ai[:n], bk[:n])) for bk in pb] for ai in pa
    ]
    return math.fsum(
        da[i][j] * db[k][l] * overlap[i][k] * overlap[j][l]
        for i in range(len(pa))
        for j in range(len(pa))
        for k in range(len(pb))
        for l in range(len(pb))
    )


def kernel_diff(da, pa, db, pb, scale=None):
    """||K_a - K_b|| / ||K_a|| in the Frobenius norm, or relative to
    sqrt(scale) (the squared norm of another kernel) when given."""
    aa = kernel_dot(da, pa, da, pa)
    bb = kernel_dot(db, pb, db, pb)
    ab = kernel_dot(da, pa, db, pb)
    norm = aa if scale is None else scale
    return math.sqrt(max(aa + bb - 2.0 * ab, 0.0) / norm) if norm > 0 else math.sqrt(bb)


def local_diff(va, vb):
    scale = max(abs(x) for x in va)
    diff = max(abs(x - y) for x, y in zip(va, vb))
    return diff / scale if scale > 0 else diff


def blocks_of(path):
    lines = open(path).read().split("\n")
    return [b for b in psp8_blocks(lines) + upf_blocks(lines) if b[2] is not None]


def differences(a, b):
    """[(name, kind, difference)] of the blocks of the files a and b, kind
    "V_NL" for a kernel and "V_loc" for the local potential; None if the
    files have different blocks."""
    blocks_a, blocks_b = blocks_of(a), blocks_of(b)
    names_a = [(x[0], len(x[2])) for x in blocks_a]
    names_b = [(x[0], len(x[2])) for x in blocks_b]
    if names_a != names_b:
        print(f"  the blocks differ: {names_a} and {names_b}")
        return None
    # the norms of the scalar-relativistic kernels, the scale of the spin-orbit ones
    scalar = {
        x[0]: kernel_dot(x[1], x[2], x[1], x[2]) for x in blocks_a if x[1] is not None
    }
    out = []
    for (name, da, pa, *_), (_, db, pb, *_) in zip(blocks_a, blocks_b):
        if da is None:
            out.append((name, "V_loc", local_diff(pa[0], pb[0])))
        else:
            scale = scalar.get(name.replace("spin-orbit ", ""))
            out.append((name, "V_NL", kernel_diff(da, pa, db, pb, scale)))
    return out


def compare(a, b, tol):
    diffs = differences(a, b)
    if diffs is None:
        return False
    if not diffs:
        print("  no psp8 or UPF file")
    ok = True
    for name, kind, diff in diffs:
        label = f"{name:24s} " + (
            "max rel diff" if kind == "V_loc" else "V_NL rel diff"
        )
        bad = tol is not None and diff > tol
        ok = ok and not bad
        print(f"  {label} {diff:.1e}" + ("   > tolerance" if bad else ""))
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("files", nargs="+")
    parser.add_argument("--tol", type=float, help="fail above this difference")
    args = parser.parse_args()
    if len(args.files) != 2:
        parser.error("give two files to compare")
    if not compare(*args.files, args.tol):
        sys.exit(1)


if __name__ == "__main__":
    main()
