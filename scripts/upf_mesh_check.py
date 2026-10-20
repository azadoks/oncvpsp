#!/usr/bin/env python3
"""Checks the mesh of the UPF file in an output, read with upf-tools.

PP_R must be every istr-th point of a log mesh of ratio AMESH (for some
stride istr >= 1), PP_RAB must be istr log(AMESH) r, and the number of points,
mesh_size, at most MAXPTS.  With --linear BASE.out, the mesh must also end at
the end of the linear mesh of the UPF file in BASE.out, within 1 %, and not
beyond it.  Prints the number of points and the stride.

usage: upf_mesh_check.py OUT AMESH MAXPTS [--linear BASE.out]
"""

import argparse
import math
import sys

import numpy as np
from pyoncvpsp.io import OncvpspTextParser
from upf_tools import UPFDict


def upf(path):
    text = OncvpspTextParser(path=path, check=False).upf_string
    if not text:
        sys.exit(f"upf_mesh_check.py: no UPF file in {path}")
    return UPFDict.from_str(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("out")
    parser.add_argument("amesh", type=float)
    parser.add_argument("maxpts", type=int)
    parser.add_argument("--linear", metavar="BASE.out")
    args = parser.parse_args()

    pseudo = upf(args.out)
    r = np.asarray(pseudo["mesh"]["r"], dtype=float)
    rab = np.asarray(pseudo["mesh"]["rab"], dtype=float)
    n = len(r)
    log_amesh = math.log(args.amesh)
    stride = round(math.log(r[1] / r[0]) / log_amesh)
    problems = []
    if stride < 1:
        problems.append(f"ratio {r[1] / r[0]}")
    steps = np.log(r[1:] / r[:-1]) / log_amesh
    if np.abs(steps - stride).max() > 1e-6:
        problems.append(f"ratio at {int(np.argmax(np.abs(steps - stride))) + 2}")
    if np.abs(rab / (r * stride * log_amesh) - 1).max() > 1e-9:
        problems.append("rab")
    if n > args.maxpts:
        problems.append(f"{n} points")
    if pseudo["header"]["mesh_size"] != n:
        problems.append(f"mesh_size {pseudo['header']['mesh_size']}, {n} points")
    if args.linear:
        end = float(upf(args.linear)["mesh"]["r"][-1])
        if not 0.99 * end < r[-1] <= end:
            problems.append(f"ends at {r[-1]}, the linear mesh at {end}")
    if problems:
        sys.exit("bad mesh: " + "; ".join(problems))
    print(f"{n} points, stride {stride}")


if __name__ == "__main__":
    main()
