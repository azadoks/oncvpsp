#!/usr/bin/env python3
"""Updates a reference output minimally to a new output.

The reference tests compare with scripts/compare_output.py, which allows a
relative (or absolute) difference of 1e-5 in every number (the base test, that
of fldiff.pl -easy, which they used before), tighter for some quantities, and
compares the projectors and their energies by the nonlocal kernels instead.
Regenerating a reference from a new output would move every number that
changed in its last digits; this keeps the reference line wherever the base
test accepts it, and takes the new line only where a number moves beyond it or
the text differs.  Lines are aligned by their shape (the text with the numbers
masked), so runs of lines that appear or disappear come from the new output.
The header above the input echo (citation and build lines) is the reference's.
Prints how many lines changed.

The lines of the quantities whose tolerance in compare_output.py is close to
their printed precision (the pseudo-atom energies, the eigenvalues of the test
configurations, the UPF file's total energy, the rows of the d2Exc matrices)
take the new line whenever its text differs.  Look at what changed before
patching (compare_output.py REF.out NEW.out), and run the reference test
after: another quantity with a tighter tolerance can still differ.

The kernels are compared to 1e-6, tighter than the 1e-5 by which a projector
line is kept, so a kernel can fail although every line passes. With --files,
the pseudopotential files (Begin PSPCODE8 or PSP_UPF to END_PSP) take every
line whose text differs too.

usage: ref_patch.py NEW.out REF.out [-o OUT] [--files]
       (OUT defaults to REF.out)
"""

import argparse
import difflib
import re

from scripts_common import fortran_float

TOL = 1.01e-5
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
ECHO = "# ATOM AND REFERENCE"
# lines taken from the new output whenever their text differs
EXACT = [
    re.compile(r"^Pseudoatom total energy"),
    re.compile(r"PS_ref=.*PS_tst="),
    re.compile(r"total_psenergy="),
    # n, l, (kap,) f, eae, eps, diff of a test configuration
    re.compile(
        r"^\s*\d+\s+\d+(\s+-?\d+)?\s+\d+\.\d{4}(\s+-?\d+\.\d{8}){2}"
        r"\s+\S+[DE][-+]\d+$"
    ),
    # a row of a d2Exc matrix
    re.compile(r"^(\s+-?\d\.\d{6}D[-+]\d\d)+$"),
]


def shape(line):
    return " ".join(NUMBER.sub("#", line).split())


def close(a, b):
    """Whether the base test accepts line b for reference line a."""
    if a == b or "pspd" in a or "date" in a:
        return True
    if any(pattern.search(a) for pattern in EXACT):
        return False
    na, nb = NUMBER.findall(a), NUMBER.findall(b)
    if shape(a) != shape(b) or len(na) != len(nb):
        return False
    for x, y in zip(na, nb):
        diff = abs(fortran_float(x) - fortran_float(y))
        if diff > TOL and diff > TOL * (abs(fortran_float(x)) + abs(fortran_float(y))):
            return False
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("new")
    parser.add_argument("ref")
    parser.add_argument("-o", "--out")
    parser.add_argument("--files", action="store_true")
    args = parser.parse_args()

    new = open(args.new).read().split("\n")
    ref = open(args.ref).read().split("\n")
    i_new = next(k for k, line in enumerate(new) if line.startswith(ECHO))
    i_ref = next(k for k, line in enumerate(ref) if line.startswith(ECHO))
    head, ref, new = ref[:i_ref], ref[i_ref:], new[i_new:]

    matcher = difflib.SequenceMatcher(
        None, [shape(x) for x in ref], [shape(x) for x in new], autojunk=False
    )
    out = []
    patched = replaced = 0
    in_file = False
    for op, i1, i2, j1, j2 in matcher.get_opcodes():
        if op == "equal":
            for a, b in zip(ref[i1:i2], new[j1:j2]):
                if a.startswith(("Begin PSPCODE8", "Begin PSP_UPF")):
                    in_file = True
                elif a.startswith("END_PSP"):
                    in_file = False
                if close(a, b) and not (args.files and in_file and a != b):
                    out.append(a)
                else:
                    out.append(b)
                    patched += 1
        else:
            out += new[j1:j2]
            replaced += (i2 - i1) + (j2 - j1)
    open(args.out or args.ref, "w").write("\n".join(head + out))
    print(f"{args.ref}: {patched} lines patched, {replaced} lines in replaced runs")


if __name__ == "__main__":
    main()
