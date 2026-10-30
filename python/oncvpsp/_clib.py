"""liboncvpsp through ctypes: the library, its structs and the functions'
signatures (include/oncvpsp.h)"""

import ctypes as C
import ctypes.util
import glob
import os

API_VERSION = 5
MXLIN = 5
MXRANK = 5


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


def _candidates():
    """where the library may be: ONCVPSP_LIBRARY, beside the package (a
    wheel's), then the system's search path"""
    if os.environ.get("ONCVPSP_LIBRARY"):
        yield os.environ["ONCVPSP_LIBRARY"]
        return
    here = os.path.dirname(os.path.abspath(__file__))
    for pattern in ("liboncvpsp.so*", "liboncvpsp*.dylib", "oncvpsp*.dll"):
        yield from sorted(glob.glob(os.path.join(here, pattern)))
    found = ctypes.util.find_library("oncvpsp")
    if found:
        yield found


def load():
    """the library, with the functions' signatures. It is loaded global:
    MKL's dispatch libraries need the symbols of the MKL it links, and
    with RTLD_LOCAL abort the process at the first LAPACK call"""
    errors = []
    for path in _candidates():
        try:
            lib = C.CDLL(path, mode=C.RTLD_GLOBAL)
            break
        except OSError as err:
            errors.append(f"{path}: {err}")
    else:
        raise ImportError(
            "liboncvpsp not found; set ONCVPSP_LIBRARY to its path"
            + "".join(f"\n  {ee}" for ee in errors)
        )
    vp, sz, cp, ci, cd = C.c_void_p, C.c_size_t, C.c_char_p, C.c_int, C.c_double
    pd, pi = C.POINTER(C.c_double), C.POINTER(C.c_int)
    signatures = {
        "oncv_api_version": ([], ci),
        "oncv_build_info": ([cp, sz, C.POINTER(sz)], ci),
        "oncv_new": ([], vp),
        "oncv_free": ([vp], None),
        "oncv_set_output": ([vp, cp, sz], ci),
        "oncv_set_output_buffer": ([vp, ci], ci),
        "oncv_output_text": ([vp, cp, sz, C.POINTER(sz)], ci),
        "oncv_run_toml": ([vp, ci, cp, sz, ci], ci),
        "oncv_ae_toml": ([vp, ci, cp, sz], ci),
        "oncv_run_ae_toml": ([vp, vp, cp, sz, ci], ci),
        "oncv_channel_toml": ([vp, vp, cp, sz, ci, ci], ci),
        "oncv_coupled": ([vp, vp, ci, cd, cd, ci], ci),
        "oncv_confine": ([vp, vp, pd, ci, pi, ci, pd, pd, ci, cd, cd], ci),
        "oncv_write_files": ([vp, cp, sz, cp, sz], ci),
        "oncv_failure": ([vp, C.POINTER(Failure)], ci),
        "oncv_dims": ([vp, C.POINTER(Dims)], ci),
        "oncv_field_name": ([ci, cp, sz], ci),
        "oncv_field_info": ([ci, ci, cp, sz, C.POINTER(sz)], ci),
        "oncv_shape": ([vp, cp, sz, pi, C.c_int * MXRANK], ci),
        "oncv_get_double": ([vp, cp, sz, pd, sz], ci),
        "oncv_get_int": ([vp, cp, sz, pi, sz], ci),
    }
    for name, (args, res) in signatures.items():
        func = getattr(lib, name)
        func.argtypes = args
        func.restype = res
    version = lib.oncv_api_version()
    if version != API_VERSION:
        raise ImportError(
            f"liboncvpsp {path} has interface version {version}, "
            f"this package {API_VERSION}"
        )
    return lib


lib = load()
