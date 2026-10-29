"""A handle of liboncvpsp: one calculation, a run or a stage of one, its
fields as NumPy arrays, its failure and its text output"""

import ctypes as C
import os
import weakref

import numpy as np

from ._clib import API_VERSION, MXRANK, Dims, Failure, lib

FAILURES = (
    "none",
    "input",
    "convergence",
    "not_found",
    "construction",
    "numerical",
    "io",
    "internal",
)
STAGES = (
    "none",
    "input",
    "ae",
    "channel",
    "vkb",
    "valence",
    "model_core",
    "screening",
    "diag",
    "ghosts",
    "unscreen",
    "configs",
    "plot",
    "phsft",
    "tests",
    "output",
    "confine",
)
MSG_PSPOT_NODE = 122


def _field_kinds():
    """{name: 'r' or 'i'} of every field the library has"""
    out, buf = {}, C.create_string_buffer(64)
    while True:
        nn = lib.oncv_field_name(len(out) + 1, buf, 64)
        if nn == 0:
            return out
        kind, name = buf.raw[:nn].decode().split()
        out[name] = kind


FIELDS = _field_kinds()
# the fields of each kind of handle, by prefix ("" for the run's, every
# one not of another kind); the library writes an error for each field
# asked of a handle that holds none, so only these are asked
STAGE_PREFIXES = ("ae.", "channel.", "confine.")
HOLDS = {
    "run": ("",) + STAGE_PREFIXES[:1],
    "ae": STAGE_PREFIXES[:1],
    "channel": STAGE_PREFIXES[1:2],
    "confine": STAGE_PREFIXES[2:],
    "none": (),
}


def _prefix(name):
    return next((pp for pp in STAGE_PREFIXES if name.startswith(pp)), "")


class OncvpspError(RuntimeError):
    """a failed run or stage: ifail and its parts (oncv_failure_t), the
    names of its kind and stage, and the text output (the message is
    there) if the handle kept it"""

    def __init__(self, ifail, failure, text):
        self.ifail = ifail
        for name, _ in Failure._fields_:
            setattr(self, name, getattr(failure, name))
        self.kind_name = FAILURES[self.kind] if self.kind < len(FAILURES) else "?"
        self.stage_name = STAGES[self.stage] if self.stage < len(STAGES) else "?"
        self.text = text
        errors = [line.strip() for line in (text or "").splitlines() if "ERROR" in line]
        where = f"stage {self.stage_name}"
        if self.ll >= 0:
            where += f", l={self.ll}"
        super().__init__(
            f"{self.kind_name} failure (ifail {ifail}, msg {self.msg}) at {where}"
            + (f": {errors[-1]}" if errors else "")
        )


class Handle:
    """one calculation of liboncvpsp, of a kind: a run, an all-electron
    atom (ae), a channel or a confined atom (confine). handle[name] is a
    field (Fortran order, as the Fortran arrays), handle.keys() those it
    holds; the handle is freed when it is no longer referenced, or by
    close()"""

    def __init__(self, kind, text=True):
        ptr = lib.oncv_new()
        if not ptr:
            raise MemoryError("oncv_new")
        self._ptr = C.c_void_p(ptr)
        self._free = weakref.finalize(self, lib.oncv_free, self._ptr)
        self._keep_text = text is True
        if self._keep_text:
            lib.oncv_set_output_buffer(self._ptr, 1)
        elif text:
            path = os.fsencode(text)
            lib.oncv_set_output(self._ptr, path, len(path))
        self.ifail = 0
        self._holds = HOLDS[kind]
        self.kind = kind

    def close(self):
        """frees the handle now"""
        self._free()

    @property
    def ptr(self):
        if not self._free.alive:
            raise ValueError("the handle is closed")
        return self._ptr

    def _done(self, ifail, check):
        self.ifail = ifail
        if ifail != 0:
            self._holds = ()
        if ifail != 0 and check:
            raise self.error()
        return self

    def error(self):
        """the OncvpspError of the handle's last run or stage"""
        return OncvpspError(self.ifail, self.failure, self.text)

    @property
    def failure(self):
        """the failure record (oncv_failure_t) of the last run or stage"""
        ff = Failure()
        lib.oncv_failure(self.ptr, C.byref(ff))
        return ff

    @property
    def dims(self):
        """the dimensions (oncv_dims_t), as a dictionary"""
        dd = Dims()
        lib.oncv_dims(self.ptr, C.byref(dd))
        if dd.version != API_VERSION:
            raise RuntimeError(f"oncv_dims has version {dd.version}")
        out = {name: getattr(dd, name) for name, _ in Dims._fields_}
        out["nproj"] = list(dd.nproj)
        return out

    @property
    def text(self):
        """the text output of the last run or stage, None if not kept"""
        if not self._keep_text:
            return None
        nout = C.c_size_t()
        lib.oncv_output_text(self.ptr, None, 0, C.byref(nout))
        buf = C.create_string_buffer(nout.value)
        lib.oncv_output_text(self.ptr, buf, nout.value, C.byref(nout))
        return buf.raw[: nout.value].decode()

    def shape(self, name):
        """the shape of a field (Fortran order), None if the handle holds
        none"""
        if name not in FIELDS or _prefix(name) not in self._holds:
            return None
        key = name.encode()
        rank, shp = C.c_int(), (C.c_int * MXRANK)()
        if lib.oncv_shape(self.ptr, key, len(key), C.byref(rank), shp) != 0:
            return None
        return tuple(shp[: rank.value])

    def __getitem__(self, name):
        if name not in FIELDS:
            raise KeyError(f"no field {name}")
        shape = self.shape(name)
        if shape is None:
            raise KeyError(f"the handle holds no {name}")
        key = name.encode()
        if FIELDS[name] == "r":
            out = np.empty(shape, dtype=np.float64, order="F")
            ptr = out.ctypes.data_as(C.POINTER(C.c_double))
            ifail = lib.oncv_get_double(self.ptr, key, len(key), ptr, out.size)
        else:
            out = np.empty(shape, dtype=np.intc, order="F")
            ptr = out.ctypes.data_as(C.POINTER(C.c_int))
            ifail = lib.oncv_get_int(self.ptr, key, len(key), ptr, out.size)
        if ifail != 0:
            raise RuntimeError(f"oncv_get of {name}: ifail {ifail}")
        return out[()] if shape == () else out

    def __contains__(self, name):
        return self.shape(name) is not None

    def keys(self):
        """the names of the fields the handle holds"""
        return [name for name in FIELDS if self.shape(name) is not None]

    def get(self, name, default=None):
        return self[name] if name in self else default

    def write(self, path, format=None):
        """writes the pseudopotential file of a run: psp8, upf or psml, by
        default from path's extension"""
        if format is None:
            format = os.path.splitext(os.fspath(path))[1].lstrip(".")
        fmt, pp = format.encode(), os.fsencode(path)
        ifail = lib.oncv_write_files(self.ptr, fmt, len(fmt), pp, len(pp))
        if ifail != 0:
            self.ifail = ifail
            raise self.error()

    def __repr__(self):
        state = "closed" if not self._free.alive else f"ifail={self.ifail}"
        return f"<oncvpsp.Handle {state}>"
