"""The C interface of liboncvpsp (include/oncvpsp.h) through ctypes, as a
Python binding will call it: the structs' layouts against the getters (the
sizes of oncv_dims against the fields' shapes), every field copied out and
the same after a second run in the same handle, bit for bit, the UPF file
written with its mesh and densities (what a check of the model core reads),
and a broken input's failure record.

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
    """the library, with the functions' signatures"""
    lib = C.CDLL(path)
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
    return lib


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

    check("api version", lib.oncv_api_version() == 1)
    handle = C.c_void_p(lib.oncv_new())
    check("run", lib.oncv_run_toml(handle, irel, text, len(text), 0) == 0)
    dims = Dims()
    lib.oncv_dims(handle, C.byref(dims))
    check("dims", dims.version == 1 and dims.irel == irel and dims.mmax > 0)
    names = fields(lib)
    first = {name: get(lib, handle, kind, name) for kind, name in names}
    check("rr against dims", first["pot.rr"][0] == (dims.mmax,))
    check("nproj against dims", first["pot.nproj"][1] == list(dims.nproj))
    check("vkb against dims", first["pot.vkb"][0][0] == dims.mmax)

    check("rerun", lib.oncv_run_toml(handle, irel, text, len(text), 0) == 0)
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
