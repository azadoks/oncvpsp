"""ONCVPSP's runs and stages in Python, over liboncvpsp (include/oncvpsp.h).

The input is a dictionary in the layout of the TOML input, the TOML text,
or the path of a TOML file; keyword arguments are tables merged over it,
key by key, and None values are left out (TOML has no null). Each function
returns a Handle, whose fields are NumPy arrays in Fortran order; a failure raises OncvpspError, or with check=False is
left in the handle (handle.ifail, handle.failure). text=True keeps the text
output (handle.text), False discards it, a path writes it there.

    import oncvpsp
    run = oncvpsp.run("tests/data/toml/14_Si.toml", irel=1)
    run["pot.vkb"].shape   # (mmax, mxprj, lmax+1, nkap)
    run.write("Si.upf")
"""

import ctypes as C
import os

import numpy as np

from ._clib import API_VERSION, lib
from ._handle import FIELDS, MSG_PSPOT_NODE, Handle, OncvpspError
from ._toml import dumps as toml_dumps
from ._toml import input_text

__all__ = [
    "API_VERSION",
    "FIELDS",
    "Handle",
    "OncvpspError",
    "all_electron",
    "build_info",
    "channel",
    "confine",
    "coupled",
    "has_node",
    "read_toml",
    "run",
    "toml_dumps",
]

_RELS = {"nr": 0, "sr": 1, "fr": 2, "r": 2}
_MODES = {"channel": 0, "projectors": 1, "node": 2}


def build_info():
    """the library's build: {"commit", "source", "compiled"}, its git
    commit, the state of the source against it and the build time"""
    nout = C.c_size_t()
    lib.oncv_build_info(None, 0, C.byref(nout))
    buf = C.create_string_buffer(nout.value)
    lib.oncv_build_info(buf, nout.value, C.byref(nout))
    return dict(line.split(" ", 1) for line in buf.raw.decode().splitlines())


def _irel(irel):
    """0 non-, 1 scalar- and 2 fully-relativistic, or 'nr', 'sr', 'fr'"""
    if isinstance(irel, str):
        return _RELS[irel]
    if irel not in (0, 1, 2):
        raise ValueError(f"irel is 0, 1 or 2, not {irel}")
    return irel


def run(inp=None, irel=1, *, ae=None, skip_ld=False, check=True, text=True, **tables):
    """the whole run of the input (oncv_run_toml), or from the all-electron
    atom of ae (all_electron, or a run), bit for bit the same
    (oncv_run_ae_toml); skip_ld leaves out the log derivatives"""
    data = input_text(inp, tables)
    hh = Handle("run", text)
    if ae is None:
        ifail = lib.oncv_run_toml(hh.ptr, _irel(irel), data, len(data), int(skip_ld))
    else:
        ifail = lib.oncv_run_ae_toml(hh.ptr, ae.ptr, data, len(data), int(skip_ld))
    return hh._done(ifail, check)


def all_electron(inp=None, irel=1, *, check=True, text=True, **tables):
    """the all-electron atom of the input (oncv_ae_toml): fields ae.*"""
    data = input_text(inp, tables)
    hh = Handle("ae", text)
    return hh._done(lib.oncv_ae_toml(hh.ptr, _irel(irel), data, len(data)), check)


def channel(inp, ae, l, mode="channel", *, check=True, text=True, **tables):
    """the channel of angular momentum l of the input, from the
    all-electron atom of ae (oncv_channel_toml): fields channel.*. mode
    'projectors' also makes the local potential and the channel's
    projectors; 'node' stops at the node test of the first projector
    (has_node)"""
    data = input_text(inp, tables)
    hh = Handle("channel" if mode != "node" else "none", text)
    ifail = lib.oncv_channel_toml(hh.ptr, ae.ptr, data, len(data), l, _MODES[mode])
    return hh._done(ifail, check)


def has_node(inp, ae, l, *, text=False, **tables):
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
    base, icmod, fcfact=0.0, rcfact=0.0, *, skip_ld=False, check=True, text=True
):
    """the model core, screening and tests again on the stages of base's
    run with the model core icmod, fcfact, rcfact (oncv_coupled): a run"""
    hh = Handle("run", text)
    ifail = lib.oncv_coupled(hh.ptr, base.ptr, icmod, fcfact, rcfact, int(skip_ld))
    return hh._done(ifail, check)


def _array(values, ctype):
    arr = (ctype * max(len(values), 1))(*values)
    return C.cast(arr, C.POINTER(ctype)), arr


def confine(
    base,
    radii=None,
    channels=None,
    energies=None,
    walls=None,
    r0=0.0,
    ecutmax=0.0,
    *,
    check=True,
    text=True,
):
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
    hh = Handle("confine", text)
    ifail = lib.oncv_confine(
        hh.ptr, base.ptr, pr, nradii, pc, len(chans), pe, pw, nwalls, r0, ecutmax
    )
    del kr, kc, ke, kw
    return hh._done(ifail, check)


def read_toml(path):
    """the dictionary of a TOML input file"""
    import tomllib

    with open(os.fspath(path), "rb") as ff:
        return tomllib.load(ff)
