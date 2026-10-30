"""The C interface of liboncvpsp (include/oncvpsp.h) through ctypes, as a
Python binding will call it: the structs' layouts against the getters (the
sizes of oncv_dims against the fields' shapes), every field copied out and
the same after a second run in the same handle, bit for bit, the UPF file
written with its mesh and densities (what a check of the model core reads),
a broken input's failure record, and the stage handles: the all-electron
atom, each channel from it with its fields' shapes, and the model core made
again on the run's stages with the input's, the same as the run.

usage: c_api_ctypes.py LIBONCVPSP IREL INPUT.toml
"""

import ctypes as C
import os
import re
import sys
import tempfile

MXRANK = 5
MXLIN = 5
FAIL_INPUT = 1
STAGE_INPUT = 1
STAGE_AE = 2
MSG_AE_GIVEN_ATOM = 176
MSG_PSPOT_NODE = 122


class Failure(C.Structure):
    """oncv_failure_t"""

    _fields_ = [
        (name, C.c_int)
        for name in ("kind", "msg", "ierr", "stage", "ll", "ikap", "iproj")
    ]


class Dims(C.Structure):
    """oncv_dims_t"""

    _fields_ = [
        ("version", C.c_int),
        ("irel", C.c_int),
        ("nkap", C.c_int),
        ("iexc", C.c_int),
        ("mmax", C.c_int),
        ("lmax", C.c_int),
        ("lloc", C.c_int),
        ("nproj", C.c_int * (MXLIN + 1)),
        ("nrl", C.c_int),
    ]


def library(path):
    """the library, with the functions' signatures; loaded global, as MKL's
    dispatch libraries need its symbols (RTLD_LOCAL aborts the process at
    the first LAPACK call)"""
    lib = C.CDLL(path, mode=C.RTLD_GLOBAL)
    vp, sz, cp = C.c_void_p, C.c_size_t, C.c_char_p
    lib.oncv_new.restype = vp
    lib.oncv_free.argtypes = [vp]
    lib.oncv_run_toml.argtypes = [vp, C.c_int, cp, sz, C.c_int]
    lib.oncv_failure.argtypes = [vp, C.POINTER(Failure)]
    lib.oncv_dims.argtypes = [vp, C.POINTER(Dims)]
    lib.oncv_field_name.argtypes = [C.c_int, cp, sz]
    lib.oncv_shape.argtypes = [vp, cp, sz, C.POINTER(C.c_int), C.c_int * MXRANK]
    lib.oncv_get_double.argtypes = [vp, cp, sz, C.POINTER(C.c_double), sz]
    lib.oncv_get_int.argtypes = [vp, cp, sz, C.POINTER(C.c_int), sz]
    lib.oncv_write_files.argtypes = [vp, cp, sz, cp, sz]
    lib.oncv_ae_toml.argtypes = [vp, C.c_int, cp, sz]
    lib.oncv_channel_toml.argtypes = [vp, vp, cp, sz, C.c_int, C.c_int]
    lib.oncv_coupled.argtypes = [vp, vp, C.c_int, C.c_double, C.c_double, C.c_int]
    lib.oncv_set_output_buffer.argtypes = [vp, C.c_int]
    lib.oncv_build_info.argtypes = [cp, sz, C.POINTER(sz)]
    lib.oncv_output_text.argtypes = [vp, cp, sz, C.POINTER(sz)]
    pd, pi = C.POINTER(C.c_double), C.POINTER(C.c_int)
    lib.oncv_confine.argtypes = [
        vp,
        vp,
        pd,
        C.c_int,
        pi,
        C.c_int,
        pd,
        pd,
        C.c_int,
        C.c_double,
        C.c_double,
    ]
    return lib


def output_text(lib, handle):
    """the buffered text output of the handle, as a string"""
    nout = C.c_size_t()
    lib.oncv_output_text(handle, None, 0, C.byref(nout))
    buf = C.create_string_buffer(nout.value)
    lib.oncv_output_text(handle, buf, nout.value, C.byref(nout))
    return buf.raw[: nout.value].decode()


def fields(lib):
    """(kind, name) of every field"""
    out, buf = [], C.create_string_buffer(64)
    while True:
        nn = lib.oncv_field_name(len(out) + 1, buf, 64)
        if nn == 0:
            return out
        kind, name = buf.raw[:nn].decode().split()
        out.append((kind, name))


def get(lib, handle, kind, name):
    """the shape and the flat values of a field"""
    key = name.encode()
    rank, shp = C.c_int(), (C.c_int * MXRANK)()
    assert lib.oncv_shape(handle, key, len(key), C.byref(rank), shp) == 0, name
    shape = tuple(shp[: rank.value])
    size = 1
    for nn in shape:
        size *= nn
    if kind == "r":
        buf = (C.c_double * max(size, 1))()
        ifail = lib.oncv_get_double(handle, key, len(key), buf, size)
    else:
        buf = (C.c_int * max(size, 1))()
        ifail = lib.oncv_get_int(handle, key, len(key), buf, size)
    assert ifail == 0, (name, ifail)
    return shape, list(buf[:size])


def slice_l(field, lpos, ll):
    """the values of a field (shape, flat values in Fortran order) at
    index ll of its axis lpos, in Fortran order"""
    shape, values = field
    inner = 1
    for nn in shape[:lpos]:
        inner *= nn
    return [vv for ii, vv in enumerate(values) if ii // inner % shape[lpos] == ll]


def stages(lib, handle, irel, text, dims, first, check):
    """the stage handles of the input text, against the run of handle: the
    all-electron atom, each channel from it, and the model core again with
    the input's (core.teter and the run's model core in first)"""
    ae, chan, again = (C.c_void_p(lib.oncv_new()) for _ in range(3))
    check("ae", lib.oncv_ae_toml(ae, irel, text, len(text)) == 0)
    for kind, name in (kn for kn in fields(lib) if kn[1].startswith("ae.")):
        check(f"{name} of the ae handle", get(lib, ae, kind, name) == first[name])
    nkap = 2 if irel == 2 else 1
    for ll in range(dims.lmax + 1):
        check(
            f"channel {ll}",
            lib.oncv_channel_toml(chan, ae, text, len(text), ll, 1) == 0,
        )
        shape, curve = get(lib, chan, "r", "channel.ecurve")
        check(
            "channel.ecurve's shape",
            len(shape) == 4 and shape[0] == 2 and shape[3] == nkap,
        )
        check("a curve", max(curve) > 0)
        check(
            "channel.pswf's shape",
            get(lib, chan, "r", "channel.pswf")[0][0] == dims.mmax,
        )
        check("channel.vloc", get(lib, chan, "r", "channel.vloc")[0] == (dims.mmax,))
        if ll != dims.lloc:
            check(
                "channel.evkb against the run's",
                get(lib, chan, "r", "channel.evkb")[1]
                == slice_l(get(lib, handle, "r", "pot.evkb"), 1, ll),
            )
    # mode 2, the node test alone, fails as the channel does at a node
    # and at the steps before it, over qcuts of channel 0 that reach one
    qcut = re.search(rb"qcut = \[ *([0-9.]+),", text)
    for factor in (1.0, 2.0, 2.5, 3.0):
        qq = b"%r" % (float(qcut.group(1)) * factor)
        other = text[: qcut.start(1)] + qq + text[qcut.end(1) :]
        full = lib.oncv_channel_toml(chan, ae, other, len(other), 0, 0)
        node = lib.oncv_channel_toml(chan, ae, other, len(other), 0, 2)
        at_node = full // 100 % 1000 == MSG_PSPOT_NODE
        check(
            f"mode 2 at qcut {qq.decode()}", node == full or node == 0 and not at_node
        )
    # another functional is another atom (check_ae)
    other = re.sub(
        rb"iexc = (\d+)", lambda m: b"iexc = %d" % (7 - int(m.group(1))), text
    )
    ifail = lib.oncv_channel_toml(chan, ae, other, len(other), 0, 0)
    failure = Failure()
    lib.oncv_failure(chan, C.byref(failure))
    check("another atom", ifail % 100 == FAIL_INPUT and failure.stage == STAGE_AE)
    check("its message", failure.msg == MSG_AE_GIVEN_ATOM)
    # the input's model core, whatever it is, from the run's stages
    lines = text.decode().split("[model_core_charge]")[-1].split("\n[")[0].splitlines()
    core = dict(line.replace(" ", "").split("=", 1) for line in lines if "=" in line)
    check(
        "coupled",
        lib.oncv_coupled(
            again,
            handle,
            int(core.get("icmod", 0)),
            float(core.get("fcfact", 0)),
            float(core.get("rcfact", 0)),
            0,
        )
        == 0,
    )
    for name in ("pot.vkb", "pot.rhomod", "tests.convergence.cross", "core.teter"):
        check(f"{name} of coupled", get(lib, again, "r", name) == first[name])
    for hh in (ae, chan, again):
        lib.oncv_free(hh)


def main():
    lib = library(sys.argv[1])
    irel = int(sys.argv[2])
    text = open(sys.argv[3], "rb").read()
    nfail = 0

    def check(label, ok):
        nonlocal nfail
        if not ok:
            print("FAILED:", label, file=sys.stderr)
            nfail += 1

    check("api version", lib.oncv_api_version() == 8)
    nout = C.c_size_t()
    lib.oncv_build_info(None, 0, C.byref(nout))
    buf = C.create_string_buffer(nout.value)
    lib.oncv_build_info(buf, nout.value, C.byref(nout))
    keys = [line.split(" ", 1)[0] for line in buf.raw.decode().splitlines()]
    check("build info", keys == ["commit", "source", "compiled"])
    handle = C.c_void_p(lib.oncv_new())
    check("run", lib.oncv_run_toml(handle, irel, text, len(text), 0) == 0)
    dims = Dims()
    lib.oncv_dims(handle, C.byref(dims))
    check("dims", dims.version == 8 and dims.irel == irel and dims.mmax > 0)
    names = [
        (kind, name)
        for kind, name in fields(lib)
        if not name.startswith(("channel.", "confine.", "local."))
    ]
    first = {name: get(lib, handle, kind, name) for kind, name in names}
    check("rr against dims", first["pot.rr"][0] == (dims.mmax,))
    check("nproj against dims", first["pot.nproj"][1] == list(dims.nproj))
    check("vkb against dims", first["pot.vkb"][0][0] == dims.mmax)

    lib.oncv_set_output_buffer(handle, 1)
    check("rerun", lib.oncv_run_toml(handle, irel, text, len(text), 0) == 0)
    out = output_text(lib, handle)
    check("the text output in the buffer", "DIAGNOSTICS" in out)
    lib.oncv_set_output_buffer(handle, 0)
    for kind, name in names:
        check(f"{name} after a rerun", get(lib, handle, kind, name) == first[name])

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "psp.upf").encode()
        check("UPF", lib.oncv_write_files(handle, b"upf", 3, path, len(path)) == 0)
        upf = open(path).read()
    sizes = {
        tag: re.search(rf'<{tag} type="real"\s+size="\s*(\d+)"', upf)
        for tag in ("PP_R", "PP_RAB", "PP_RHOATOM")
    }
    check("UPF mesh and densities", all(sizes.values()))
    check(
        "UPF sizes",
        len({int(m.group(1)) for m in sizes.values() if m}) == 1,
    )

    stages(lib, handle, irel, text, dims, first, check)

    broken = b"[oncvpsp]\natsym = 'O'\n"
    ifail = lib.oncv_run_toml(handle, irel, broken, len(broken), 0)
    failure = Failure()
    lib.oncv_failure(handle, C.byref(failure))
    check("broken input", ifail % 100 == FAIL_INPUT == failure.kind)
    check("its message", failure.msg == ifail // 100 % 1000 > 0)
    check("at the input", failure.stage == STAGE_INPUT and failure.ll == -1)
    lib.oncv_free(handle)

    if nfail:
        print(f"c_api_ctypes: {nfail} check(s) failed", file=sys.stderr)
        return 1
    print(f"c_api_ctypes: {len(names)} fields, all checks passed", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
