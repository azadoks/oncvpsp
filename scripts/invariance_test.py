#!/usr/bin/env python3
"""Checks that the pseudopotential does not depend on one input parameter.

Runs the program on a TOML input with KEY set to each of VALUES in turn, takes
from each output the features of smooth_scan.py (the kernels, the local
potential and the model core at the sample radii, the log derivatives, the
convergence tables, the ghost test) and prints, for each group, the largest
difference from the first value's features, relative to the group's scale as
smooth_scan.py forms it. With --tol (one value, or group=value pairs), the
exit status is 1 if a difference exceeds it.

The well states' tail in E_resid is the use (doc/smoothness_plan.md, phase 8):
with well_tail = "free", nothing depends on the well's onset (well_onset).

usage: invariance_test.py EXE BASE.toml KEY VALUE VALUE [VALUE ...]
                          [--set KEY=VALUE ...] [--groups G,G]
                          [--tol TOL | --tol G=TOL,G=TOL] [--sample S]
                          [--workdir DIR]
"""

import argparse
import concurrent.futures
import os
import pathlib
import sys
import tempfile
import tomllib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from smooth_scan import (  # noqa: E402
    GROUPS,
    apply_sets,
    features,
    file_mesh,
    group_scales,
    parse_tol,
    run_one,
)


def group_of(name):
    """The group of a feature, kernel for kernel.log:..."""
    return name.split(":", 1)[0].split(".")[0]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("exe")
    parser.add_argument("toml")
    parser.add_argument("key")
    parser.add_argument("values", nargs="+", help="values of KEY, in TOML")
    parser.add_argument("--set", action="append", metavar="KEY=VALUE")
    parser.add_argument("--groups", default="kernel,vloc,core")
    parser.add_argument("--tol")
    parser.add_argument("--sample", type=float, default=0.25)
    parser.add_argument("--workdir")
    args = parser.parse_args()
    if len(args.values) < 2:
        parser.error("give at least two values")
    exe = os.path.realpath(args.exe)
    data = apply_sets(tomllib.loads(pathlib.Path(args.toml).read_text()), args.set)
    values = [tomllib.loads(f"v = {v}")["v"] for v in args.values]

    workdir = args.workdir or tempfile.mkdtemp(prefix="invariance_")
    with concurrent.futures.ThreadPoolExecutor(len(values)) as pool:
        outputs = list(
            pool.map(
                lambda kv: run_one(exe, data, args.key, kv[1], workdir, kv[0], ["log"]),
                enumerate(values),
            )
        )

    mesh = file_mesh(outputs[0], "log")
    samples = sorted(
        {
            min(range(len(mesh)), key=lambda k: abs(mesh[k] - r))
            for r in (
                args.sample * j for j in range(1, int(mesh[-1] / args.sample) + 1)
            )
        }
    )
    groups = set(args.groups.split(","))
    runs = [
        {
            n: v
            for n, v in features(out, {"log": samples}).items()
            if group_of(n) in groups
        }
        for out in outputs
    ]
    names = sorted(set.intersection(*(set(r) for r in runs)))
    scales = group_scales({n: runs[0][n] for n in names})
    tols = parse_tol(args.tol)

    print(f"{args.key}: {', '.join(args.values)}, outputs in {workdir}")
    print(f"{'group':8s} {'value':>12s} {'difference':>11s}  worst feature")
    status = 0
    for group in GROUPS:
        if group not in groups:
            continue
        mine = [n for n in names if group_of(n) == group]
        if not mine:
            print(f"{group:8s} no features  FAIL")
            status = 1
            continue
        tol = tols.get(group, tols.get(None))
        for text, run in zip(args.values[1:], runs[1:]):
            diff, name = max((abs(run[n] - runs[0][n]) / scales[n], n) for n in mine)
            flag = ""
            if tol is not None and diff > tol:
                flag = "  FAIL"
                status = 1
            print(f"{group:8s} {text:>12s} {diff:11.2e}  {name}{flag}")
    return status


if __name__ == "__main__":
    sys.exit(main())
