#!/usr/bin/env python3
"""Checks the projectors of the psp8 and UPF files on the linear output mesh
(doc/smoothness_plan.md, phase 4).

Runs EXE on BASE.toml with psfile "both" and upf_mesh 0 (the linear mesh),
in WORKDIR, and checks in the output that

  - every projector of the psp8 file (scalar-relativistic and spin-orbit) and
    of the UPF file is exactly 0 at every radius beyond the core radius of its
    l (O-1, R-15: the writers evaluate each projector from its pieces, which
    end at rc; an interpolant of the projector on the log mesh rings past it);
  - in the fully-relativistic psp8 file, every term of the decomposition into
    scalar-relativistic and spin-orbit projectors is written: nproj of them
    for l = 0 and 2 nproj of each kind for l > 0 (O-2, D-2; terms below 1e-3
    of the first were dropped).

usage: projector_output_test.py EXE BASE.toml WORKDIR
"""

import argparse
import pathlib
import subprocess
import sys
import tomllib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_output import psp8_part, upf_part  # noqa: E402
from pyoncvpsp.io import OncvpspTextParser  # noqa: E402
from smooth_scan import dump_toml  # noqa: E402


def beyond_rc(name, r, blocks, rc):
    """Failures for the projectors of blocks {(l, label): Block} not 0 at the
    radii r beyond rc[l]."""
    failures = []
    for (ll, label), block in blocks.items():
        outside = r > rc[ll]
        values = block.beta[:, outside]
        nonzero = int((values != 0.0).sum())
        if nonzero:
            failures.append(
                f"{name} {label}: {nonzero} values beyond rc = {rc[ll]:g},"
                f" the largest {abs(values).max():.2e}"
            )
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("exe")
    parser.add_argument("toml")
    parser.add_argument("workdir")
    args = parser.parse_args()
    exe = str(pathlib.Path(args.exe).resolve())

    data = tomllib.loads(pathlib.Path(args.toml).read_text())
    data["pp_output"]["psfile"] = "both"
    data["pp_output"]["upf_mesh"] = 0
    work = pathlib.Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    (work / "in.toml").write_text(dump_toml(data))
    result = subprocess.run(
        [exe, "-t", "in.toml"], cwd=work, capture_output=True, text=True
    )
    (work / "out").write_text(result.stdout)
    if result.returncode != 0:
        sys.exit(f"{exe} failed in {work}")

    pseudo = data["pseudopotentials"]
    rc = dict(zip(pseudo["l"], pseudo["rc"]))
    vkb = data["vkb_projectors"]
    nproj = dict(zip(vkb["l"], vkb["nproj"]))
    lloc = data["local_potential"]["lloc"]

    output = OncvpspTextParser(path=work / "out", check=False)
    psp8, blocks8 = psp8_part(output.psp8_string)
    upf, blocks_upf = upf_part(output.upf_string)
    failures = beyond_rc(
        "psp8",
        psp8["r"],
        {(ll, f"{kind} l={ll}"): b for (_, kind, ll), b in blocks8.items()},
        rc,
    )
    failures += beyond_rc(
        "UPF",
        upf["mesh"]["r"],
        {
            (ll, f"l={ll}" if jj is None else f"l={ll} j={jj:g}"): b
            for (_, ll, jj), b in blocks_upf.items()
        },
        rc,
    )

    if any(kind == "spin_orbit" for _, kind, _ in blocks8):
        for ll, count in nproj.items():
            if ll == lloc or count == 0:
                continue
            expected = {"scalar": count if ll == 0 else 2 * count}
            if ll > 0:
                expected["spin_orbit"] = 2 * count
            for kind, want in expected.items():
                block = blocks8.get(("psp8", kind, ll))
                got = 0 if block is None else block.beta.shape[0]
                if got != want:
                    failures.append(f"psp8 {kind} l={ll}: {got} projectors, not {want}")

    for failure in failures:
        print("FAILED:", failure)
    if failures:
        return 1
    print("projectors: 0 beyond rc in psp8 and UPF; every fr term in psp8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
