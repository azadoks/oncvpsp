#!/usr/bin/env python3
"""Updates a reference output minimally to a new output.

The reference tests compare with fldiff, which allows a relative (or absolute)
difference of 1e-5 in every number, apart from the projectors and their
energies, which they compare by the nonlocal kernels instead
(scripts/kernel_diff.py).  Regenerating a reference from a new output
would move every number that changed in its last digits; this keeps the
reference line wherever fldiff would accept it, and takes the new line only
where a number moves beyond the tolerance or the text differs.  Lines are
aligned by their shape (the text with the numbers masked), so runs of lines
that appear or disappear come from the new output.  The header above the input
echo (citation and build lines) is the reference's.  Prints how many lines
changed.

Look at what changed before patching: for the pseudopotential itself,
scripts/kernel_diff.py compares the nonlocal kernel and the local
potential, which a change of the projectors alone leaves the same.

usage: ref_patch.py NEW.out REF.out [-o OUT]   (OUT defaults to REF.out)
"""

import argparse
import difflib
import re

TOL = 1.01e-5
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
ECHO = "# ATOM AND REFERENCE"


def value(text):
    return float(text.replace("D", "E").replace("d", "e"))


def shape(line):
    return " ".join(NUMBER.sub("#", line).split())


def close(a, b):
    """Whether fldiff accepts line b for reference line a."""
    if a == b or "pspd" in a or "date" in a:
        return True
    na, nb = NUMBER.findall(a), NUMBER.findall(b)
    if shape(a) != shape(b) or len(na) != len(nb):
        return False
    for x, y in zip(na, nb):
        diff = abs(value(x) - value(y))
        if diff > TOL and diff > TOL * (abs(value(x)) + abs(value(y))):
            return False
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("new")
    parser.add_argument("ref")
    parser.add_argument("-o", "--out")
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
    for op, i1, i2, j1, j2 in matcher.get_opcodes():
        if op == "equal":
            for a, b in zip(ref[i1:i2], new[j1:j2]):
                if close(a, b):
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
