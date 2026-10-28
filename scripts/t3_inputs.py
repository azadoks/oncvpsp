#!/usr/bin/env python3
"""Writes the inputs of the on-node equivalence test (T3 of
doc/smoothness_plan.md): each reference input with its core radii, its
local-potential radius and its qcut already where the code before phase 3
rounds them to.

That code truncates each rc to 6 decimals and takes the first log-mesh point
beyond it, and puts qcut on the nearest point of eresid's q grid, a multiple
of dq = 2 pi / (150 rc) at the rounded rc. Given those values, it computes
what it did for the original input, and code that no longer rounds computes
at the same radii and wave vectors, so the two differ by their numerics
alone. The values are written to 17 digits.

usage: t3_inputs.py OUTDIR   (OUTDIR/NAME.toml for each tests/refs/NAME.out)
"""

import math
import pathlib
import sys

from input_text_to_toml import OncvpspInput
from scripts_common import dump_toml

ROOT = pathlib.Path(__file__).resolve().parent.parent


def log_mesh(data, irel):
    """The log mesh of oncvpsp_main for the input data and irel."""
    mesh = data.get("log_mesh", {})
    zz = float(data["oncvpsp"]["z"])
    amesh = mesh.get("amesh", 0.0) if mesh.get("amesh", 0.0) > 0 else 1.006
    rr1 = mesh.get("rr1", 0.0)
    if rr1 <= 0.0:
        rr1 = 0.0005 / zz
        if irel != 2:
            rr1 = min(rr1, 0.0005 / 10)
    rmax = mesh.get("rmax", 0.0) if mesh.get("rmax", 0.0) > 0 else 45.0
    al = math.log(amesh)
    mmax = int(math.log(rmax / rr1) / al)
    return [rr1 * math.exp(al * ii) for ii in range(mmax)]


def truncated(rc):
    """rc truncated to 6 decimals, as oncvpsp_main does."""
    return int(rc * 10.0e5) / 10.0e5


def snapped(rr, rct):
    """The first mesh point beyond the truncated radius rct."""
    return next(r for r in rr[1:] if r > rct)


def main():
    outdir = pathlib.Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)
    for ref in sorted((ROOT / "tests/refs").glob("*.out")):
        name = ref.stem
        irel = 0 if name.endswith("_nr") else 2 if name.endswith("_r") else 1
        stem = name.removesuffix("_nr").removesuffix("_r")
        dat = ROOT / "tests/data" / f"{stem}.dat"
        if not dat.is_file():
            dat = ROOT / "tests/data" / f"{name}.dat"
        data = OncvpspInput.from_text_file(dat).to_dict()
        rr = log_mesh(data, irel)
        pp = data["pseudopotentials"]
        rct = [truncated(rc) for rc in pp["rc"]]
        rcmax = max(rct)
        rcs = [snapped(rr, rc if rc != 0.0 else rcmax) for rc in rct]
        pp["rc"] = rcs
        pp["qcut"] = [
            (2.0 * math.pi / (150.0 * rc)) * round(q / (2.0 * math.pi / (150.0 * rc)))
            for q, rc in zip(pp["qcut"], rcs)
        ]
        loc = data["local_potential"]
        if loc["lloc"] == 4:
            loc["rcloc"] = snapped(rr, truncated(loc["rcloc"]))
        (outdir / f"{name}.toml").write_text(dump_toml(data))
    return 0


if __name__ == "__main__":
    sys.exit(main())
