"""The Python package oncvpsp (python/oncvpsp) against the executable: the
run's text output and psp8 file the executable's byte for byte (its output
but for the banner); the input as the TOML file, its dictionary and tables
given as keyword arguments, the same; the arrays in Fortran order against
the dimensions and each field's axes; the library's build the
executable's; the accessors (projectors, kernel, densities, convergence)
against the fields; the stages (the all-electron atom, a run from it, each
channel with its projectors, the node test, the model core again, a
confined atom) against the run; and a failure raised as its class, or left
in the handle, and a handle closed.

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
    check("dims", dims.version == oncvpsp.API_VERSION and dims.irel == irel)
    check("a run", isinstance(run, oncvpsp.Run) and isinstance(dims.irel, oncvpsp.Irel))
    vkb = run["pot.vkb"]
    check("pot.rr", run["pot.rr"].shape == (dims.mmax,))
    check(
        "pot.vkb in Fortran order",
        vkb.shape[0] == dims.mmax
        and vkb.shape[3] == dims.nkap
        and vkb.flags["F_CONTIGUOUS"],
    )
    check("a scalar", np.ndim(run["pot.zion"]) == 0)
    check("no channel in a run", "channel.pswf" not in run)
    fields = dict(run)
    check("every field of a run", len(fields) > 50 and len(run) == len(fields))

    def axes_ranks(handle):
        """each field's axes (field_info), as many as its rank"""
        return all(
            len(oncvpsp.field_info(nn).axes) == np.ndim(handle[nn])
            and oncvpsp.field_info(nn).doc
            for nn in handle
        )

    check("the run's axes", axes_ranks(run))

    # the accessors against the fields
    check("mesh", np.array_equal(run.mesh, fields["pot.rr"]))
    check(
        "local potential",
        np.array_equal(run.local_potential, fields["pot.vpuns"][:, dims.lloc]),
    )
    for ll in range(dims.lmax + 1):
        for ikap in range(1, dims.nkap + 1):
            ee, beta = run.projectors(ll, ikap)
            nn = dims.nproj[ll]
            kk = (
                sum(ee[ii] * np.outer(beta[:, ii], beta[:, ii]) for ii in range(nn))
                if nn
                else np.zeros((dims.mmax, dims.mmax))
            )
            check(f"l={ll} ikap={ikap}'s kernel", np.allclose(run.kernel(ll, ikap), kk))
            if nn:
                table = run.convergence(ll, ikap)
                check(
                    f"l={ll} ikap={ikap}'s convergence",
                    table.shape[1] == 2 and np.all(np.diff(table[:, 1]) < 0),
                )
    rr = run.mesh
    charge = np.trapezoid(4 * np.pi * rr**2 * run.valence_density, rr)
    check("the valence charge", abs(charge - fields["pot.zion"]) < 1e-3)

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
    check("ae's fields", all(nn.startswith("ae.") for nn in ae) and axes_ranks(ae))
    check("a run from ae", same(oncvpsp.run(data, ae=ae), fields))
    for ll in range(dims.lmax + 1):
        ch = oncvpsp.channel(data, ae, ll, "projectors")
        check(
            f"channel {ll}'s axes", isinstance(ch, oncvpsp.Channel) and axes_ranks(ch)
        )
        if ll != dims.lloc:
            check(
                f"channel {ll}'s projectors",
                np.array_equal(ch["channel.evkb"][:, 0], fields["pot.evkb"][:, ll, 0])
                and np.allclose(ch.kernel(), run.kernel(ll)),
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
        check(
            "a node raised as NodeError",
            ch.ifail == 0
            or isinstance(ch.error(), oncvpsp.NodeError) == (ch.failure.msg == 122),
        )
    # the run assembled from channels made on the local potential alone
    loc = oncvpsp.local_potential(data, ae, text=False)
    check("the local potential", isinstance(loc, oncvpsp.Local) and axes_ranks(loc))
    tiles = []
    for ll in range(dims.lmax + 1):
        ch = oncvpsp.channel(
            data, ae, ll, "projectors", local=loc, tests=True, text=False
        )
        check(f"channel {ll}'s tests' axes", axes_ranks(ch))
        mk = ch["channel.logder.ae"].shape[1]
        check(
            f"channel {ll}'s log derivatives",
            np.array_equal(
                ch["channel.logder.ae"], fields["tests.logder.ae"][:, ll, :mk]
            )
            and (
                not ch.exact
                or np.array_equal(
                    ch["channel.logder.ps"], fields["tests.logder.ps"][:, ll, :mk]
                )
            ),
        )
        check(
            f"channel {ll}'s ghosts near the run's",
            np.allclose(
                ch["channel.ghosts.ee"],
                fields["tests.ghosts.ee"][:, :, ll, :mk],
                atol=1e-5,
            ),
        )
        if ll != dims.lloc:
            tiles.append(ch)
    check(
        "the assembled run, the run's",
        same(oncvpsp.assemble(data, ae, tiles, local=loc, text=False), fields),
    )
    try:
        oncvpsp.assemble(data, ae, tiles[1:], text=False)
        check("an l missing raises", False)
    except oncvpsp.InputError as err:
        check("an l missing", err.msg == 190)
    pp = dict(data["pseudopotentials"], qcut=[qq * 1.01 for qq in qcut])
    try:
        oncvpsp.assemble(data, ae, tiles, local=loc, text=False, pseudopotentials=pp)
        check("a channel of another input raises", False)
    except oncvpsp.InputError as err:
        check("a channel of another input", err.msg in (186, 187))

    # the icmod 6 search from a guess of its radius: the same results
    six = oncvpsp.coupled(run, 6, fcfact=-1.0, rcfact=1e-4, check=False, text=False)
    if six.ifail == 0:
        check("the search's passes", six["core.budget.passes"] >= 1)
        for factor in (0.5, 1.05, 1.5):
            guessed = oncvpsp.coupled(
                run,
                6,
                fcfact=-1.0,
                rcfact=1e-4,
                text=False,
                rs_guess=float(six["core.budget.r_s"]) * factor,
            )
            check(
                f"the search from {factor} r_s",
                all(np.array_equal(guessed[nn], six[nn]) for nn in six),
            )

    with oncvpsp.coupled(run, 0, text=False) as again:
        check("the model core again", np.array_equal(again["pot.vkb"], vkb))
    check("closed at the end of the with block", repr(again).endswith("closed>"))
    conf = oncvpsp.confine(run, radii=(), walls=())
    check("a confined atom of nothing", conf["confine.nmiss"] == -1)
    conf = oncvpsp.confine(run, walls=())
    check("the sphere", conf["confine.nmiss"] >= 0 and "\nws_r " in conf.text)
    check("the confined atom's axes", axes_ranks(conf))

    # failures
    try:
        oncvpsp.run({"oncvpsp": {"atsym": "O"}}, irel)
        check("a broken input raises", False)
    except oncvpsp.InputError as err:
        check("its kind and stage", err.kind_name == "input" == err.stage_name)
        check(
            "as names",
            err.kind is oncvpsp.FailureKind.INPUT and err.stage is oncvpsp.Stage.INPUT,
        )
        check("its ifail", err.ifail % 100 == err.kind and err.msg > 0)
    broken = oncvpsp.run({"oncvpsp": {"atsym": "O"}}, irel, check=False)
    check("or is left in the handle", broken.ifail != 0 and len(broken) == 0)
    check("as its class", isinstance(broken.error(), oncvpsp.InputError))
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
