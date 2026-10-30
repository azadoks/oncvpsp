"""ONCVPSP's runs and stages in Python, over liboncvpsp (include/oncvpsp.h).

The input is a dictionary in the layout of the TOML input, the TOML text,
or the path of a TOML file; keyword arguments are tables merged over it,
key by key, and None values are left out (TOML has no null). Each function
returns a handle (Run, AllElectron, Channel, Confinement), a Mapping of
its fields to NumPy arrays in Fortran order, which field_info describes; a
failure raises an OncvpspError of its kind, or with check=False is left in
the handle (handle.ifail, handle.failure). text=True keeps the text output
(handle.text), False discards it, a path writes it there.

    import oncvpsp
    with oncvpsp.run("tests/data/toml/14_Si.toml", irel="sr") as run:
        run["pot.vkb"].shape   # (mmax, mxprj, mxl+1, nkap)
        run.kernel(0)          # the nonlocal operator of l=0, (mmax, mmax)
        run.write("Si.upf")
"""

import ctypes as C
import os
from typing import Any

import numpy as np

from ._clib import API_VERSION, lib
from ._handle import (
    FIELDS,
    MSG_PSPOT_NODE,
    AllElectron,
    Channel,
    Confinement,
    ConstructionError,
    ConvergenceError,
    Dims,
    Failure,
    FailureKind,
    FieldInfo,
    FileIOError,
    Handle,
    InputError,
    InternalError,
    Irel,
    Local,
    NodeError,
    NotFoundError,
    NumericalError,
    OncvpspError,
    Projectors,
    Run,
    Stage,
    field_info,
)
from ._toml import dumps as toml_dumps
from ._toml import input_text

__all__ = [
    "API_VERSION",
    "FIELDS",
    "AllElectron",
    "Channel",
    "Confinement",
    "ConstructionError",
    "ConvergenceError",
    "Dims",
    "Failure",
    "FailureKind",
    "FieldInfo",
    "FileIOError",
    "Handle",
    "InputError",
    "InternalError",
    "Irel",
    "Local",
    "NodeError",
    "NotFoundError",
    "NumericalError",
    "OncvpspError",
    "Projectors",
    "Run",
    "Stage",
    "all_electron",
    "assemble",
    "build_info",
    "channel",
    "confine",
    "coupled",
    "field_info",
    "has_node",
    "local_potential",
    "read_toml",
    "run",
    "toml_dumps",
]

_RELS = {"nr": Irel.NR, "sr": Irel.SR, "fr": Irel.FR, "r": Irel.FR}
_MODES = {"channel": 0, "projectors": 1, "node": 2}
Input = dict[str, Any] | str | os.PathLike | None
Text = bool | str | os.PathLike


def build_info() -> dict[str, str]:
    """the library's build: {"commit", "source", "compiled"}, its git
    commit, the state of the source against it and the build time"""
    nout = C.c_size_t()
    lib.oncv_build_info(None, 0, C.byref(nout))
    buf = C.create_string_buffer(nout.value)
    lib.oncv_build_info(buf, nout.value, C.byref(nout))
    return dict(line.split(" ", 1) for line in buf.raw.decode().splitlines())


def _irel(irel):
    """0 non-, 1 scalar- and 2 fully-relativistic (Irel), or 'nr', 'sr',
    'fr'"""
    if isinstance(irel, str):
        return int(_RELS[irel.lower()])
    if irel not in (0, 1, 2):
        raise ValueError(f"irel is 0, 1 or 2, not {irel}")
    return int(irel)


def run(
    inp: Input = None,
    irel: int | str = 1,
    *,
    ae: "AllElectron | Run | None" = None,
    skip_ld: bool = False,
    check: bool = True,
    text: Text = True,
    **tables: Any,
) -> Run:
    """the whole run of the input (oncv_run_toml), or from the all-electron
    atom of ae (all_electron, or a run), bit for bit the same
    (oncv_run_ae_toml); skip_ld leaves out the log derivatives"""
    data = input_text(inp, tables)
    hh = Run(text)
    if ae is None:
        ifail = lib.oncv_run_toml(hh.ptr, _irel(irel), data, len(data), int(skip_ld))
    else:
        ifail = lib.oncv_run_ae_toml(hh.ptr, ae.ptr, data, len(data), int(skip_ld))
    return hh._done(ifail, check)


def all_electron(
    inp: Input = None,
    irel: int | str = 1,
    *,
    check: bool = True,
    text: Text = True,
    **tables: Any,
) -> AllElectron:
    """the all-electron atom of the input (oncv_ae_toml): fields ae.*"""
    data = input_text(inp, tables)
    hh = AllElectron(text)
    return hh._done(lib.oncv_ae_toml(hh.ptr, _irel(irel), data, len(data)), check)


def channel(
    inp: Input,
    ae: "AllElectron | Run",
    l: int,  # noqa: E741
    mode: str = "channel",
    *,
    local: Local | None = None,
    tests: bool = False,
    energies: Any = None,
    check: bool = True,
    text: Text = True,
    **tables: Any,
) -> Channel:
    """the channel of angular momentum l of the input, from the
    all-electron atom of ae (oncv_channel_toml): fields channel.*. mode
    'projectors' also makes the local potential (or takes local's,
    local_potential) and the channel's projectors; 'node' stops at the
    node test of the first projector (has_node). tests=True, with mode
    'projectors', also makes its log derivatives (channel.logder.*, on
    energies (emin, emax, de) in Ha, None for the input's window) and its
    ghost tests (channel.ghosts.*), on the local potential before a whole
    run's ionic tail beyond rcmax: Channel.exact says whether the log
    derivatives are the run's bit for bit"""
    data = input_text(inp, tables)
    hh = Channel(text)
    if mode == "node":
        hh._holds = ()
    if local is None and not tests:
        ifail = lib.oncv_channel_toml(hh.ptr, ae.ptr, data, len(data), l, _MODES[mode])
        return hh._done(ifail, check)
    if energies is not None and len(energies) != 3:
        raise ValueError("energies is (emin, emax, de)")
    pe, ke = (None, None) if energies is None else _array(list(energies), C.c_double)
    ifail = lib.oncv_channel_local_toml(
        hh.ptr,
        ae.ptr,
        None if local is None else local.ptr,
        data,
        len(data),
        l,
        _MODES[mode],
        int(tests),
        pe,
    )
    del ke
    return hh._done(ifail, check)


def local_potential(
    inp: Input,
    ae: "AllElectron | Run",
    *,
    check: bool = True,
    text: Text = True,
    **tables: Any,
) -> Local:
    """the local potential alone of the input, from the all-electron atom
    of ae (oncv_local_toml): fields local.*, for channel(..., local=) and
    assemble"""
    data = input_text(inp, tables)
    hh = Local(text)
    return hh._done(lib.oncv_local_toml(hh.ptr, ae.ptr, data, len(data)), check)


def assemble(
    inp: Input,
    ae: "AllElectron | Run",
    channels: "list[Channel]",
    *,
    local: Local | None = None,
    skip_ld: bool = False,
    check: bool = True,
    text: Text = True,
    **tables: Any,
) -> Run:
    """the run of the input from the all-electron atom of ae with the
    channels made before (channel, modes 'channel' or 'projectors', in any
    order) in place of optimizing them (oncv_assemble_toml): one of each l
    from 0 to lmax, the local channel (lloc < 4) local's unless one of the
    channels is of l = lloc. Each must be of this input and atom; the run
    is then run(inp, ae=ae)'s bit for bit, its text output without the
    channels' lines"""
    data = input_text(inp, tables)
    ptrs = (C.c_void_p * max(len(channels), 1))(*[ch.ptr for ch in channels])
    hh = Run(text)
    ifail = lib.oncv_assemble_toml(
        hh.ptr,
        ae.ptr,
        None if local is None else local.ptr,
        C.cast(ptrs, C.POINTER(C.c_void_p)),
        len(channels),
        data,
        len(data),
        int(skip_ld),
    )
    return hh._done(ifail, check)


def has_node(
    inp: Input,
    ae: "AllElectron | Run",
    l: int,  # noqa: E741
    *,
    text: Text = False,
    **tables: Any,
) -> bool:
    """whether the first pseudo-wave of channel l has a node (pspot's
    test), from the all-electron atom of ae; a failure before the test
    raises OncvpspError"""
    hh = channel(inp, ae, l, "node", check=False, text=text, **tables)
    if hh.ifail == 0:
        return False
    if hh.failure.msg == MSG_PSPOT_NODE:
        return True
    raise hh.error()


def coupled(
    base: Run,
    icmod: int,
    *,
    fcfact: float = 0.0,
    rcfact: float = 0.0,
    r_s: float = 0.0,
    skip_ld: bool = False,
    check: bool = True,
    text: Text = True,
) -> Run:
    """the model core, screening and tests again on the stages of base's
    run with the model core icmod, fcfact, rcfact and r_s (icmod 8, 9; 0 for
    the input's) (oncv_coupled): a run"""
    hh = Run(text)
    ifail = lib.oncv_coupled(hh.ptr, base.ptr, icmod, fcfact, rcfact, r_s, int(skip_ld))
    return hh._done(ifail, check)


def _array(values, ctype):
    arr = (ctype * max(len(values), 1))(*values)
    return C.cast(arr, C.POINTER(ctype)), arr


def confine(
    base: Run,
    radii: Any = None,
    channels: Any = None,
    energies: Any = None,
    walls: Any = None,
    r0: float = 0.0,
    ecutmax: float = 0.0,
    *,
    check: bool = True,
    text: Text = True,
) -> Confinement:
    """the confined atom of base's run, as oncvpsp_confine.x (oncv_confine):
    fields confine.*. radii of the sphere (bohr), None for the default
    (1.25, 1.5, 2 times the largest core radius), () for no sphere;
    channels' l, None for all; energies (emin, emax, de) in Ha, None for
    the input's window and 0.1; walls (Ha/bohr**2), None for the default,
    () for none; r0 the walls' onset, ecutmax the basis's largest cutoff
    (Ha), 0 for their defaults"""
    radii = None if radii is None else np.atleast_1d(np.asarray(radii, float))
    nradii = 0 if radii is None else (len(radii) if len(radii) else -1)
    pr, kr = _array([] if radii is None else radii.tolist(), C.c_double)
    chans = [] if channels is None else [int(ll) for ll in np.atleast_1d(channels)]
    pc, kc = _array(chans, C.c_int)
    pe, ke = (None, None) if energies is None else _array(list(energies), C.c_double)
    if energies is not None and len(energies) != 3:
        raise ValueError("energies is (emin, emax, de)")
    ww = None if walls is None else np.atleast_1d(np.asarray(walls, float))
    nwalls = -1 if ww is None else len(ww)
    pw, kw = _array([] if ww is None else ww.tolist(), C.c_double)
    hh = Confinement(text)
    ifail = lib.oncv_confine(
        hh.ptr, base.ptr, pr, nradii, pc, len(chans), pe, pw, nwalls, r0, ecutmax
    )
    del kr, kc, ke, kw
    return hh._done(ifail, check)


def read_toml(path: str | os.PathLike) -> dict[str, Any]:
    """the dictionary of a TOML input file"""
    import tomllib

    with open(os.fspath(path), "rb") as ff:
        return tomllib.load(ff)
