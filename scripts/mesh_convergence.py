#!/usr/bin/env python3
"""Measures how the nonlocal kernels converge as the log mesh is refined.

Runs the program on a TOML input at each --amesh and at the finer --ref, and
prints for each kernel block (each l, or l and kind) the relative difference
of its kernel from the --ref run's, as compare_output.py forms it, and the
observed order between successive meshes, log(d1/d2)/log(al1/al2), al =
ln(amesh). The kernels are those of the psp8 file, on its linear mesh, which
does not change with amesh (doc/smoothness_plan.md, T4: the rates must not
get worse in a phase that changes numerics).

With --set KEY=VALUE (repeatable, VALUE in TOML) the input is changed first,
as in smooth_scan.py.

usage: mesh_convergence.py EXE BASE.toml [--amesh A A ...] [--ref A]
                           [--set KEY=VALUE ...] [--workdir DIR]
"""

import argparse
import math
import pathlib
import subprocess
import sys
import tempfile
import tomllib

from compare_output import kernel_diff, kernel_norm2, load
from scripts_common import apply_sets, dump_toml


def run(exe, data, amesh, workdir):
    data = {k: (dict(v) if isinstance(v, dict) else v) for k, v in data.items()}
    data.setdefault("log_mesh", {})
    data["log_mesh"] = dict(data["log_mesh"], amesh=amesh)
    data["pp_output"] = dict(data.get("pp_output", {}), psfile="psp8")
    path = pathlib.Path(workdir) / f"amesh{amesh}"
    path.mkdir(parents=True, exist_ok=True)
    (path / "in.toml").write_text(dump_toml(data))
    result = subprocess.run(
        [exe, "-t", "in.toml"], cwd=path, capture_output=True, text=True
    )
    (path / "out").write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"{exe} failed at amesh={amesh} in {path}")
    return load(path / "out", False)[1]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("exe")
    parser.add_argument("toml")
    parser.add_argument("--amesh", type=float, nargs="+", default=[1.012, 1.006, 1.003])
    parser.add_argument("--ref", type=float, default=1.0015)
    parser.add_argument("--set", action="append", metavar="KEY=VALUE")
    parser.add_argument("--workdir")
    args = parser.parse_args()

    data = apply_sets(tomllib.loads(pathlib.Path(args.toml).read_text()), args.set)
    exe = str(pathlib.Path(args.exe).resolve())
    workdir = args.workdir or tempfile.mkdtemp(prefix="mesh_convergence_")
    ref = run(exe, data, args.ref, workdir)
    diffs = {}
    for amesh in args.amesh:
        blocks = run(exe, data, amesh, workdir)
        for key, block in blocks.items():
            if key not in ref or block.beta.shape[0] != ref[key].beta.shape[0]:
                continue
            scale_block = ref.get(ref[key].scale_key) if ref[key].scale_key else None
            scale = kernel_norm2(scale_block if scale_block is not None else ref[key])
            diffs.setdefault(key, []).append(kernel_diff(block, ref[key], scale))

    print(f"{args.toml}: kernels against amesh={args.ref}, outputs in {workdir}")
    head = "".join(f"{a:>11.4f}" for a in args.amesh)
    print(f"{'block':28s}{head}   orders")
    als = [math.log(a) for a in args.amesh]
    for key, ds in sorted(diffs.items(), key=str):
        orders = [
            math.log(ds[i] / ds[i + 1]) / math.log(als[i] / als[i + 1])
            if ds[i] > 0 and ds[i + 1] > 0
            else float("nan")
            for i in range(len(ds) - 1)
        ]
        cells = "".join(f"{d:11.2e}" for d in ds)
        print(f"{str(key[1:]):28s}{cells}   " + " ".join(f"{o:5.2f}" for o in orders))
    return 0


if __name__ == "__main__":
    sys.exit(main())
