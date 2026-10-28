"""Checks that the two kernel comparisons agree: kernel_diff.py, which needs
only the standard library and works from the overlaps of the projectors,
and compare_output.py, which reads the files with pyoncvpsp and forms the
kernels on the grid. B is the output A with the projectors of its psp8
file changed smoothly, each by its own factor (1 + EPS k r for projector k
at radius r), written to WORKDIR/B.out. The relative difference of the
nonlocal kernel of each block (relative to the norm of A's kernel, or of
the scalar-relativistic kernel of its l for a spin-orbit one) must be the
same from both, to TOL relative.

usage: kernel_agreement_test.py A.out WORKDIR
"""

import pathlib
import re
import sys

import compare_output
import kernel_diff

# the two forms of the same difference, for kernels that differ well above
# round-off (the overlap form loses precision only as they come to agree)
TOL = 1e-6
# the smallest difference compared: the overlap form resolves about 1e-8
MIN_DIFF = 1e-8
# the change of the projectors, per unit radius and projector
EPS = 1e-3
NAME = re.compile(r"(psp8|UPF)( spin-orbit)? l=(\d+)(?: j=(\S+))?$")


def key_of(name):
    """compare_output.py's key of the block kernel_diff.py calls name."""
    form, so, ll, jj = NAME.match(name).groups()
    if form == "psp8":
        return ("psp8", "spin_orbit" if so else "scalar", int(ll))
    return ("upf", int(ll), None if jj is None else float(jj))


def perturbed(path_a, path_b):
    """Writes the output path_a to path_b with the projectors of its psp8
    file changed."""
    lines = pathlib.Path(path_a).read_text().split("\n")
    blocks = [b for b in kernel_diff.psp8_blocks(lines) if b[1] is not None]
    if not blocks:
        sys.exit(f"kernel_agreement_test.py: no psp8 projectors in {path_a}")
    for _, _, betas, first, last in blocks:
        for row in range(first + 1, last + 1):
            words = lines[row].split()
            rr = float(words[1])
            values = [
                float(x) * (1.0 + EPS * k * rr) for k, x in enumerate(words[2:], 1)
            ]
            lines[row] = f"{words[0]:>6s}" + "".join(
                f" {x: .13E}" for x in [rr, *values]
            )
    pathlib.Path(path_b).write_text("\n".join(lines))


def main():
    path_a, work = sys.argv[1:3]
    pathlib.Path(work).mkdir(parents=True, exist_ok=True)
    path_b = str(pathlib.Path(work) / "B.out")
    perturbed(path_a, path_b)
    blocks_a = compare_output.load(path_a, with_psml=False)[1]
    blocks_b = compare_output.load(path_b, with_psml=False)[1]
    diffs = kernel_diff.differences(path_a, path_b)
    if diffs is None:
        sys.exit("kernel_agreement_test.py: kernel_diff.py found different blocks")
    kernels = [(name, diff) for name, kind, diff in diffs if kind == "V_NL"]
    if {key_of(name) for name, _ in kernels} != set(blocks_a):
        sys.exit("kernel_agreement_test.py: the two tools find different blocks")
    worst = 0.0
    for name, diff in kernels:
        key = key_of(name)
        a = blocks_a[key]
        scale = compare_output.kernel_norm2(
            blocks_a[a.scale_key] if a.scale_key is not None else a
        )
        other = compare_output.kernel_diff(a, blocks_b[key], scale)
        rel = abs(diff - other) / other
        worst = max(worst, rel)
        print(f"  {name:24s} {diff:.6e} {other:.6e}  {rel:.1e}")
        if not other > MIN_DIFF or rel > TOL:
            sys.exit(
                f"kernel_agreement_test.py: {name}: {diff:.6e} against {other:.6e}"
            )
    print(f"{len(kernels)} kernels, the two tools within {worst:.1e} relative")


if __name__ == "__main__":
    main()
