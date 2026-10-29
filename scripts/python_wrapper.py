"""The Python package oncvpsp (python/oncvpsp) against the executable: the
run's text output and psp8 file the executable's byte for byte (its output
but for the banner); the input as the TOML file, its dictionary and tables
given as keyword arguments, the same; the arrays in Fortran order against
the dimensions; the library's build the executable's; the stages (the all-electron atom, a run from it, each
channel with its projectors, the node test, a confined atom) against the
run; and a failure raised, or left in the handle, and a handle closed.

usage: python_wrapper.py FULL_EXE IREL INPUT.toml WORKDIR
(ONCVPSP_LIBRARY the library, the package on PYTHONPATH)
"""

import gc
import os
import re
import subprocess
import sys

import numpy as np

import oncvpsp


def main():
    full, irel, toml, work = sys.argv[1:5]
    irel = int(irel)
    os.makedirs(work, exist_ok=True)
    os.chdir(work)
    nfail = 0

    def check(label, ok):
        nonlocal nfail
        if not ok:
            print("FAILED:", label, file=sys.stderr)
            nfail += 1

    with open("exe.out", "w") as fh:
        subprocess.run([full, "--toml-input", toml], stdout=fh, check=True)
    exe = open("exe.out").read()

    # the run: its text output and psp8 file, the executable's
    run = oncvpsp.run(toml, irel)
    text = run.text
    run.write("lib.psp8")
    psp8 = open("lib.psp8").read()
    at = exe.find(text)
    check("the text output, the executable's", at > 0 and "build: " in exe[:at])
    # between them, the fully relativistic executable's lines of the
    # spin-orbit projectors' conversion, which the file's writing makes
    between, _, after = exe[at + len(text) :].partition("Begin PSPCODE8\n")
    check(
        "the psp8 file, the executable's",
        after.startswith(psp8)
        and all(
            line.startswith("l1-1,npr_sr,npr_so") or not line
            for line in between.splitlines()
        ),
    )

    # the library's build, the executable's (one build)
    info = oncvpsp.build_info()
    check(
        "the build, the executable's",
        all(f"{kk} {vv}" in exe.replace(",", "") for kk, vv in info.items()),
    )

    dims = run.dims
    check("dims", dims["version"] == oncvpsp.API_VERSION and dims["irel"] == irel)
    vkb = run["pot.vkb"]
    check("pot.rr", run["pot.rr"].shape == (dims["mmax"],))
    check(
        "pot.vkb in Fortran order",
        vkb.shape[0] == dims["mmax"]
        and vkb.shape[3] == dims["nkap"]
        and vkb.flags["F_CONTIGUOUS"],
    )
    check("a scalar", np.ndim(run["pot.zion"]) == 0)
    check("no channel in a run", "channel.pswf" not in run)
    fields = {name: run[name] for name in run.keys()}
    check("every field of a run", len(fields) > 50)

    def same(handle, names):
        return all(np.array_equal(handle[nn], fields[nn]) for nn in names)

    # the input as a dictionary, and as tables in keyword arguments
    data = oncvpsp.read_toml(toml)
    check("the dictionary", same(oncvpsp.run(data, irel), fields))
    nones = {
        kk: dict(vv, nothing=None) if isinstance(vv, dict) else vv
        for kk, vv in data.items()
    }
    check("None values left out", same(oncvpsp.run(nones, irel), fields))
    check(
        "the tables as keyword arguments", same(oncvpsp.run({}, irel, **data), fields)
    )

    # the stages against the run
    ae = oncvpsp.all_electron(data, irel)
    check("ae's fields", all(nn.startswith("ae.") for nn in ae.keys()))
    check("a run from ae", same(oncvpsp.run(data, ae=ae), fields))
    for ll in range(dims["lmax"] + 1):
        ch = oncvpsp.channel(data, ae, ll, "projectors")
        if ll != dims["lloc"]:
            check(
                f"channel {ll}'s projectors",
                np.array_equal(ch["channel.evkb"][:, 0], fields["pot.evkb"][:, ll, 0]),
            )
        check(f"no node in channel {ll}", not oncvpsp.has_node(data, ae, ll))
    # the node test, over qcuts of channel 0 that reach a node
    qcut = data["pseudopotentials"]["qcut"]
    for factor in (2.0, 3.0):
        pp = dict(data["pseudopotentials"], qcut=[qcut[0] * factor] + qcut[1:])
        ch = oncvpsp.channel(data, ae, 0, check=False, pseudopotentials=pp)
        check(
            f"the node test at qcut {pp['qcut'][0]}",
            oncvpsp.has_node(data, ae, 0, pseudopotentials=pp)
            == (ch.ifail != 0 and ch.failure.msg == 122),
        )
    conf = oncvpsp.confine(run, radii=(), walls=())
    check("a confined atom of nothing", conf["confine.nmiss"] == -1)
    conf = oncvpsp.confine(run, walls=())
    check("the sphere", conf["confine.nmiss"] >= 0 and "\nws_r " in conf.text)

    # failures
    try:
        oncvpsp.run({"oncvpsp": {"atsym": "O"}}, irel)
        check("a broken input raises", False)
    except oncvpsp.OncvpspError as err:
        check("its kind and stage", err.kind_name == "input" == err.stage_name)
        check("its ifail", err.ifail % 100 == err.kind and err.msg > 0)
    broken = oncvpsp.run({"oncvpsp": {"atsym": "O"}}, irel, check=False)
    check("or is left in the handle", broken.ifail != 0 and broken.keys() == [])
    try:
        oncvpsp.confine(run, radii=0.1)
        check("a radius inside rc raises", False)
    except oncvpsp.OncvpspError as err:
        check("at the confined atom", err.stage_name == "confine")
        check("with its message", re.search("confine_atom: ERROR", str(err)))

    # handles freed
    run.close()
    try:
        run["pot.rr"]
        check("a closed handle refuses", False)
    except ValueError:
        pass
    for _ in range(3):
        oncvpsp.all_electron(data, irel, text=False)
    gc.collect()

    if nfail:
        print(f"python_wrapper: {nfail} check(s) failed", file=sys.stderr)
        return 1
    print(f"python_wrapper: {len(fields)} fields, all checks passed", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
