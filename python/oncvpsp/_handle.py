"""A handle of liboncvpsp: one calculation, a run or a stage of one, its
fields as NumPy arrays, its failure and its text output"""

import ctypes as C
import enum
import math
import os
import weakref
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, fields
from typing import Any, NamedTuple

import numpy as np

from ._clib import API_VERSION, MXRANK, lib
from ._clib import Dims as _CDims
from ._clib import Failure as _CFailure


class Irel(enum.IntEnum):
    """the relativistic treatment"""

    NR = 0
    SR = 1
    FR = 2


class FailureKind(enum.IntEnum):
    """the kind of a failure (failure_m's fail_*, ONCV_FAIL_*)"""

    NONE = 0
    INPUT = 1
    CONVERGENCE = 2
    NOT_FOUND = 3
    CONSTRUCTION = 4
    NUMERICAL = 5
    IO = 6
    INTERNAL = 7


class Stage(enum.IntEnum):
    """the stage of a run a failure is of (failure_m's stage_*,
    ONCV_STAGE_*)"""

    NONE = 0
    INPUT = 1
    AE = 2
    CHANNEL = 3
    VKB = 4
    VALENCE = 5
    MODEL_CORE = 6
    SCREENING = 7
    DIAG = 8
    GHOSTS = 9
    UNSCREEN = 10
    CONFIGS = 11
    PLOT = 12
    PHSFT = 13
    TESTS = 14
    OUTPUT = 15
    CONFINE = 16


MSG_PSPOT_NODE = 122


def _enum(cls, value):
    """the member of value, or value itself if a newer library has more"""
    try:
        return cls(value)
    except ValueError:
        return value


@dataclass(frozen=True)
class FieldInfo:
    """a field of the results (oncv_field_info): its name, its kind ('r'
    real, 'i' integer), its axes, the first index first (an axis of
    different quantities as their names in brackets, "[ecut,err]"), its
    units ('' for none) and what it is"""

    name: str
    kind: str
    axes: tuple[str, ...]
    units: str
    doc: str

    def __str__(self):
        axes = ", ".join(self.axes) or "scalar"
        units = f" [{self.units}]" if self.units else ""
        return f"{self.name} ({axes}){units}: {self.doc}"


def _info_text(ifield, what):
    nout = C.c_size_t()
    lib.oncv_field_info(ifield, what, None, 0, C.byref(nout))
    buf = C.create_string_buffer(max(nout.value, 1))
    lib.oncv_field_info(ifield, what, buf, nout.value, C.byref(nout))
    return buf.raw[: nout.value].decode()


def _field_infos():
    """{name: FieldInfo} of every field the library has"""
    out, buf = {}, C.create_string_buffer(64)
    while True:
        ifield = len(out) + 1
        nn = lib.oncv_field_name(ifield, buf, 64)
        if nn == 0:
            return out
        kind, name = buf.raw[:nn].decode().split()
        axes, units, doc = (_info_text(ifield, what) for what in (1, 2, 3))
        out[name] = FieldInfo(name, kind, tuple(axes.split()), units, doc)


FIELDS: dict[str, FieldInfo] = _field_infos()
# the fields of each kind of handle, by prefix ("" for the run's, every
# one not of another kind); the library writes an error for each field
# asked of a handle that holds none, so only these are asked
STAGE_PREFIXES = ("ae.", "channel.", "confine.", "local.")


def field_info(name: str) -> FieldInfo:
    """what the field name is: its axes, units and a sentence"""
    try:
        return FIELDS[name]
    except KeyError:
        raise KeyError(f"no field {name}") from None


def _prefix(name):
    return next((pp for pp in STAGE_PREFIXES if name.startswith(pp)), "")


@dataclass(frozen=True)
class Failure:
    """the failure record of a run or stage (oncv_failure_t): its kind,
    its message's number, the solver status, the stage, and where it was:
    ll (-1 none), ikap and iproj (0 none)"""

    kind: FailureKind
    msg: int
    ierr: int
    stage: Stage
    ll: int
    ikap: int
    iproj: int


@dataclass(frozen=True)
class Dims:
    """the sizes of a run's results (oncv_dims_t)"""

    version: int
    irel: Irel
    nkap: int
    iexc: int
    mmax: int
    lmax: int
    lloc: int
    nproj: tuple[int, ...]
    nrl: int


class OncvpspError(RuntimeError):
    """a failed run or stage: ifail, the parts of its record (failure, and
    each as an attribute) and the text output (the message is there) if
    the handle kept it. It is raised as the class of its kind
    (InputError, ConvergenceError, ...), a node in a channel's first
    pseudo-wave as NodeError"""

    def __init__(self, ifail: int, failure: Failure, text: str | None):
        self.ifail = ifail
        self.failure = failure
        for ff in fields(Failure):
            setattr(self, ff.name, getattr(failure, ff.name))
        self.text = text
        errors = [line.strip() for line in (text or "").splitlines() if "ERROR" in line]
        where = f"stage {self.stage_name}"
        if self.ll >= 0:
            where += f", l={self.ll}"
        super().__init__(
            f"{self.kind_name} failure (ifail {ifail}, msg {self.msg}) at {where}"
            + (f": {errors[-1]}" if errors else "")
        )

    @property
    def kind_name(self) -> str:
        return self.kind.name.lower() if isinstance(self.kind, FailureKind) else "?"

    @property
    def stage_name(self) -> str:
        return self.stage.name.lower() if isinstance(self.stage, Stage) else "?"

    @classmethod
    def of(cls, ifail: int, failure: Failure, text: str | None) -> "OncvpspError":
        """the error of the failure, of the class of its kind"""
        if failure.msg == MSG_PSPOT_NODE:
            return NodeError(ifail, failure, text)
        return _ERRORS.get(failure.kind, OncvpspError)(ifail, failure, text)


class InputError(OncvpspError):
    """the input is invalid, inconsistent or asks for something not built
    in"""


class ConvergenceError(OncvpspError):
    """an iteration reached its limit"""


class NotFoundError(OncvpspError):
    """something the construction needs does not exist for these
    parameters: a bound state, a root in a bracket, a crossing"""


class ConstructionError(OncvpspError):
    """the pseudopotential cannot be constructed as asked"""


class NodeError(ConstructionError):
    """a channel's first pseudo-wave has a node (has_node)"""


class NumericalError(OncvpspError):
    """a numerical failure"""


class FileIOError(OncvpspError):
    """a file that could not be read or written"""


class InternalError(OncvpspError):
    """a failure that should not happen"""


_ERRORS = {
    FailureKind.INPUT: InputError,
    FailureKind.CONVERGENCE: ConvergenceError,
    FailureKind.NOT_FOUND: NotFoundError,
    FailureKind.CONSTRUCTION: ConstructionError,
    FailureKind.NUMERICAL: NumericalError,
    FailureKind.IO: FileIOError,
    FailureKind.INTERNAL: InternalError,
}


class Projectors(NamedTuple):
    """a channel's projectors, nproj of them: their coefficients (Ha) and
    the projectors r*beta on the mesh (mmax, nproj), normalized. Only
    their kernel, sum_i vkb[:, i] energies[i] vkb[:, i]^T, is independent
    of their basis"""

    energies: np.ndarray
    vkb: np.ndarray


def _kernel(pp):
    return (pp.vkb * pp.energies) @ pp.vkb.T


class Handle(Mapping[str, Any]):
    """one calculation of liboncvpsp: a Mapping of the fields it holds to
    NumPy arrays (Fortran order, as the Fortran arrays; a scalar as a
    NumPy scalar). The handle is freed when it is no longer referenced, at
    the end of a with block, or by close()"""

    # the prefixes of the fields it holds when its run or stage succeeds
    _HOLDS: tuple[str, ...] = ()

    def __init__(self, text: bool | str | os.PathLike = True):
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
        self._holds = self._HOLDS

    def close(self) -> None:
        """frees the handle now"""
        self._free()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def ptr(self) -> C.c_void_p:
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

    def error(self) -> OncvpspError:
        """the OncvpspError of the handle's last run or stage"""
        return OncvpspError.of(self.ifail, self.failure, self.text)

    @property
    def failure(self) -> Failure:
        """the failure record of the last run or stage"""
        ff = _CFailure()
        lib.oncv_failure(self.ptr, C.byref(ff))
        values = {name: getattr(ff, name) for name, _ in _CFailure._fields_}
        values["kind"] = _enum(FailureKind, values["kind"])
        values["stage"] = _enum(Stage, values["stage"])
        return Failure(**values)

    @property
    def dims(self) -> Dims:
        """the sizes of the results"""
        dd = _CDims()
        lib.oncv_dims(self.ptr, C.byref(dd))
        if dd.version != API_VERSION:
            raise RuntimeError(f"oncv_dims has version {dd.version}")
        values = {name: getattr(dd, name) for name, _ in _CDims._fields_}
        values["nproj"] = tuple(dd.nproj)
        values["irel"] = _enum(Irel, dd.irel)
        return Dims(**values)

    @property
    def text(self) -> str | None:
        """the text output of the last run or stage, None if not kept"""
        if not self._keep_text:
            return None
        nout = C.c_size_t()
        lib.oncv_output_text(self.ptr, None, 0, C.byref(nout))
        buf = C.create_string_buffer(nout.value)
        lib.oncv_output_text(self.ptr, buf, nout.value, C.byref(nout))
        return buf.raw[: nout.value].decode()

    def shape(self, name: str) -> tuple[int, ...] | None:
        """the shape of a field (Fortran order), None if the handle holds
        none"""
        if name not in FIELDS or _prefix(name) not in self._holds:
            return None
        key = name.encode()
        rank, shp = C.c_int(), (C.c_int * MXRANK)()
        if lib.oncv_shape(self.ptr, key, len(key), C.byref(rank), shp) != 0:
            return None
        return tuple(shp[: rank.value])

    def __getitem__(self, name: str) -> Any:
        if name not in FIELDS:
            raise KeyError(f"no field {name}")
        shape = self.shape(name)
        if shape is None:
            raise KeyError(f"the handle holds no {name}")
        key = name.encode()
        if FIELDS[name].kind == "r":
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

    def __contains__(self, name: object) -> bool:
        return isinstance(name, str) and self.shape(name) is not None

    def __iter__(self) -> Iterator[str]:
        return (name for name in FIELDS if self.shape(name) is not None)

    def __len__(self) -> int:
        return sum(1 for _ in self)

    # one calculation: equal only to itself (not Mapping's equal fields)
    __eq__ = object.__eq__
    __hash__ = object.__hash__

    def _ipython_key_completions_(self):
        return list(self)

    @staticmethod
    def info(name: str) -> FieldInfo:
        """what the field name is (field_info)"""
        return field_info(name)

    def __repr__(self):
        name = f"oncvpsp.{type(self).__name__}"
        if not self._free.alive:
            return f"<{name} closed>"
        if self.ifail:
            return f"<{name} ifail={self.ifail}>"
        return f"<{name} {len(self)} fields>"


def _kap(ikap, nkap):
    if not 1 <= ikap <= nkap:
        raise ValueError(f"ikap is 1 to {nkap}, not {ikap}")
    return ikap - 1


def _l(l, lmax):  # noqa: E741
    if not 0 <= l <= lmax:
        raise ValueError(f"l is 0 to {lmax}, not {l}")
    return l


class Run(Handle):
    """a whole run (run, coupled): the pseudopotential pot.*, its tests
    tests.*, the model core core.*, the warnings, and its all-electron
    atom ae.*; write writes its pseudopotential file. ikap is 1, or 2 for
    j=l-1/2 when fully relativistic"""

    _HOLDS = ("", "ae.")

    def write(self, path: str | os.PathLike, format: str | None = None) -> None:
        """writes the pseudopotential file: psp8, upf or psml, by default
        from path's extension"""
        if format is None:
            format = os.path.splitext(os.fspath(path))[1].lstrip(".")
        fmt, pp = format.encode(), os.fsencode(path)
        ifail = lib.oncv_write_files(self.ptr, fmt, len(fmt), pp, len(pp))
        if ifail != 0:
            self.ifail = ifail
            raise self.error()

    @property
    def mesh(self) -> np.ndarray:
        """the log radial mesh (bohr)"""
        return self["pot.rr"]

    @property
    def local_potential(self) -> np.ndarray:
        """the local potential (Ha) on the mesh"""
        return self["pot.vpuns"][:, self.dims.lloc]

    def semilocal(self, l: int) -> np.ndarray:  # noqa: E741
        """the unscreened semi-local potential of l (Ha) on the mesh"""
        return self["pot.vpuns"][:, _l(l, self.dims.lmax)]

    def projectors(self, l: int, ikap: int = 1) -> Projectors:  # noqa: E741
        """the projectors of l, nproj of them"""
        dims = self.dims
        kk, nn = _kap(ikap, dims.nkap), dims.nproj[_l(l, dims.lmax)]
        return Projectors(
            self["pot.evkb"][:nn, l, kk].copy(), self["pot.vkb"][:, :nn, l, kk].copy()
        )

    def kernel(self, l: int, ikap: int = 1) -> np.ndarray:  # noqa: E741
        """the nonlocal operator of l on the mesh, sum_i vkb_i(r) e_i
        vkb_i(r') (Ha/bohr), (mmax, mmax): of the projectors, what does
        not depend on their basis"""
        return _kernel(self.projectors(l, ikap))

    @property
    def valence_density(self) -> np.ndarray:
        """the pseudo valence density (e/bohr^3) of the reference
        configuration: integral 4 pi r^2 rho dr is the valence charge (the
        field pot.rho is 4 pi times it)"""
        return self["pot.rho"] / (4 * math.pi)

    @property
    def model_core(self) -> np.ndarray:
        """the model core density (e/bohr^3), as valence_density"""
        return self["pot.rhomod"][:, 0] / (4 * math.pi)

    @property
    def core_density(self) -> np.ndarray:
        """the all-electron core density (e/bohr^3), as valence_density"""
        return self["rhoc"] / (4 * math.pi)

    def convergence(self, l: int, ikap: int = 1, iproj: int = 1) -> np.ndarray:  # noqa: E741
        """the convergence of the pseudo-wave iproj of l with the cutoff:
        rows (Ecut, err) (Ha) of the levels it reaches, err from 1e-2 down
        by half decades"""
        dims = self.dims
        kk = _kap(ikap, dims.nkap)
        if not 1 <= iproj <= dims.nproj[_l(l, dims.lmax)]:
            raise ValueError(f"iproj is 1 to {dims.nproj[l]}, not {iproj}")
        table = self["tests.convergence.cross"][:, :, iproj - 1, l, kk].T
        return table[table[:, 0] > 0].copy()


class AllElectron(Handle):
    """the all-electron atom alone (all_electron): fields ae.*, from which
    runs and channels start"""

    _HOLDS = ("ae.",)

    @property
    def mesh(self) -> np.ndarray:
        """the log radial mesh (bohr)"""
        return self["ae.rr"]


class Channel(Handle):
    """one channel (channel): fields channel.*"""

    _HOLDS = ("channel.",)

    @property
    def l(self) -> int:  # noqa: E743
        """its angular momentum"""
        return int(self["channel.l"])

    def projectors(self, ikap: int = 1) -> Projectors:
        """its projectors, with mode 'projectors'"""
        kk = _kap(ikap, int(self["channel.mkap"]))
        evkb = self["channel.evkb"]
        if evkb.size == 0:
            raise ValueError("the channel has no projectors without mode 'projectors'")
        nn = int(self["channel.nproj"])
        return Projectors(evkb[:nn, kk].copy(), self["channel.vkb"][:, :nn, kk].copy())

    def kernel(self, ikap: int = 1) -> np.ndarray:
        """its nonlocal operator, as Run.kernel, with mode 'projectors'"""
        return _kernel(self.projectors(ikap))

    @property
    def exact(self) -> bool:
        """with tests, whether its log derivatives are those of the
        input's whole run bit for bit (they read the local potential no
        further than rcmax, where the run's ionic tail starts); the ghost
        tests always reach beyond it"""
        return bool(self["channel.logder.exact"])


class Local(Handle):
    """the local potential alone (local_potential): fields local.*, from
    which channels make their projectors without making it again, and
    for lloc < 4 the local channel a run assembled from channels takes"""

    _HOLDS = ("local.",)

    @property
    def vloc(self) -> np.ndarray:
        """the screened local potential (Ha) on the mesh, before a whole
        run's ionic tail beyond rcmax"""
        return self["local.vloc"]


class Confinement(Handle):
    """the confined atom of a run (confine): fields confine.*, the summary
    line of oncvpsp_confine.x"""

    _HOLDS = ("confine.",)
