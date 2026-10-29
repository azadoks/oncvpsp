"""The confined atom of liboncvpsp (oncv_confine) through ctypes against
oncvpsp_confine.x: the run's all-electron and pseudopotential data files
from the executables, the confined atom of each, and the same run through
the library, whose confined atom, in its text output's buffer, must be the
executable's output line for line, byte for byte (but its build line), and
whose summary fields must be its confine_summary line; with the default
options and with radii, a channel, a window, a wall, its onset and a
cutoff given.

usage: c_api_confine.py LIBONCVPSP FULL_EXE AE_EXE CONFINE_EXE IREL
                        INPUT.toml WORKDIR
"""

import ctypes as C
import math
import os
import subprocess
import sys

from c_api_ctypes import get, library, output_text


def run(args, out):
    """runs an executable, its stdout into the file out"""
    with open(out, "w") as fh:
        subprocess.run(args, stdout=fh, check=True)
    return open(out).read()


def summary_fields(lib, handle):
    """the confine_summary line's values from the fields: NaN for a part
    left out, Infinity for a cutoff never reached"""

    def value(xx, inf=0):
        if inf:
            return math.inf
        return math.nan if xx < 0 else xx

    cut = get(lib, handle, "r", "confine.cut")[1]
    inf = get(lib, handle, "i", "confine.inf")[1]
    return (
        [
            value(get(lib, handle, "r", "confine.max_de")[1][0]),
            value(cut[0], inf[0]),
            value(cut[1], inf[1]),
        ]
        + [
            get(lib, handle, "i", f"confine.{name}")[1][0]
            for name in ("nmiss", "nextra", "noff", "nwfail")
        ]
        + [
            value(get(lib, handle, "r", "confine.dde")[1][0]),
            value(cut[2], inf[2]),
            value(cut[3], inf[3]),
        ]
    )


def same(aa, bb):
    """the values the same, NaN the same as NaN"""
    return len(aa) == len(bb) and all(
        (math.isnan(x) and math.isnan(y))
        if isinstance(x, float) and math.isnan(x)
        else x == y
        for x, y in zip(aa, bb)
    )


def main():
    lib_path, full, ae, confine, irel, toml, work = sys.argv[1:8]
    irel = int(irel)
    os.makedirs(work, exist_ok=True)
    os.chdir(work)
    nfail = 0

    def check(label, ok):
        nonlocal nfail
        if not ok:
            print("FAILED:", label, file=sys.stderr)
            nfail += 1

    run([ae, "-a", "atom.ae", "--toml-input", toml], "ae.out")
    run(
        [full, "--no-log-derivatives", "--ps-file", "atom.ps", "--toml-input", toml],
        "full.out",
    )

    lib = library(lib_path)
    text = open(toml, "rb").read()
    hrun, hc = C.c_void_p(lib.oncv_new()), C.c_void_p(lib.oncv_new())
    check("run", lib.oncv_run_toml(hrun, irel, text, len(text), 1) == 0)
    lib.oncv_set_output_buffer(hc, 1)

    out = run([confine, "-a", "atom.ae", "-p", "atom.ps"], "conf.out")
    rs = sorted(
        {
            float(line.split()[3])
            for line in out.splitlines()
            if line.startswith("ws_r ")
        }
    )
    options = [
        ("defaults", [], (None, 0, None, 0, None, None, -1, 0.0, 0.0)),
        (
            "options",
            [
                "-r",
                f"{rs[0]!r},{rs[-1]!r}",
                "-l",
                "1",
                "-e",
                "-2.0,2.0,0.137",
                "-w",
                "0.05",
                "-R",
                f"{rs[0]!r}",
                "-c",
                "100",
            ],
            (
                (C.c_double * 2)(rs[0], rs[-1]),
                2,
                (C.c_int * 1)(1),
                1,
                (C.c_double * 3)(-2.0, 2.0, 0.137),
                (C.c_double * 1)(0.05),
                1,
                rs[0],
                100.0,
            ),
        ),
    ]
    for label, args, call in options:
        out = run([confine, "-a", "atom.ae", "-p", "atom.ps"] + args, f"{label}.out")
        check(f"oncv_confine, {label}", lib.oncv_confine(hc, hrun, *call) == 0)
        lines = [
            line
            for line in out.splitlines(keepends=True)
            if not line.startswith("build: ")
        ]
        buffer = output_text(lib, hc)
        open(f"{label}.lib.out", "w").write(buffer)
        check(f"the text output, {label}", buffer == "".join(lines))
        last = [line for line in lines if line.startswith("confine_summary ")]
        check(f"a summary, {label}", len(last) == 1)
        if last:
            words = last[0].split()[1:]
            values = (
                [float(w) for w in words[:3]]
                + [int(w) for w in words[3:7]]
                + [float(w) for w in words[7:]]
            )
            check(
                f"the summary's fields, {label}", same(summary_fields(lib, hc), values)
            )
    lib.oncv_free(hrun)
    lib.oncv_free(hc)

    if nfail:
        print(f"c_api_confine: {nfail} check(s) failed", file=sys.stderr)
        return 1
    print("c_api_confine: all checks passed", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
