"""The text output of liboncvpsp with the executable's text
(oncv_set_output_exe) through ctypes against the executables: a whole run
(oncv_run_toml) against oncvpsp.x, and a run from oncvpsp_ae.x's
all-electron data file (oncv_run_ps_toml) against oncvpsp_ps.x, each byte
for byte, the banner, the build line and the pseudopotential file
included; and without oncv_set_output_exe, the same text less those.

usage: c_api_exe_text.py LIBONCVPSP FULL_EXE AE_EXE PS_EXE IREL INPUT.toml
                         WORKDIR
"""

import ctypes as C
import os
import subprocess
import sys

from c_api_ctypes import library, output_text


def run(args, out):
    """runs an executable, its stdout into the file out"""
    with open(out, "w") as fh:
        subprocess.run(args, stdout=fh, check=True)
    return open(out).read()


def main():
    lib_path, full, ae, ps, irel, toml, work = sys.argv[1:8]
    irel = int(irel)
    os.makedirs(work, exist_ok=True)
    os.chdir(work)
    nfail = 0

    def check(label, ok):
        nonlocal nfail
        if not ok:
            print("FAILED:", label, file=sys.stderr)
            nfail += 1

    want_full = run([full, "--toml-input", toml], "full.out")
    run([ae, "-a", "atom.ae", "--toml-input", toml], "ae.out")
    want_ps = run([ps, "-a", "atom.ae", "--toml-input", toml], "ps.out")

    lib = library(lib_path)
    vp, sz, cp = C.c_void_p, C.c_size_t, C.c_char_p
    lib.oncv_set_output_exe.argtypes = [vp, C.c_int]
    lib.oncv_run_ps_toml.argtypes = [vp, C.c_int, cp, sz, cp, sz, C.c_int]
    text = open(toml, "rb").read()
    path = os.fsencode(os.path.abspath("atom.ae"))
    hh = C.c_void_p(lib.oncv_new())
    lib.oncv_set_output_buffer(hh, 1)

    lib.oncv_set_output_exe(hh, 1)
    check("run", lib.oncv_run_toml(hh, irel, text, len(text), 0) == 0)
    got = output_text(lib, hh)
    open("full.lib", "w").write(got)
    check("oncvpsp.x's text", got == want_full)
    check(
        "run from the file",
        lib.oncv_run_ps_toml(hh, irel, text, len(text), path, len(path), 0) == 0,
    )
    got = output_text(lib, hh)
    open("ps.lib", "w").write(got)
    check("oncvpsp_ps.x's text", got == want_ps)

    # without the executable's text: its run's text, less the banner and
    # the file at the end
    lib.oncv_set_output_exe(hh, 0)
    check("run, plain", lib.oncv_run_toml(hh, irel, text, len(text), 0) == 0)
    plain = output_text(lib, hh)
    check("the plain text inside the executable's", plain in want_full)
    check("the plain text has no banner", "ONCVPSP  (" not in plain)
    check("the plain text is shorter", len(plain) < len(want_full))

    lib.oncv_free(hh)
    print(f"c_api_exe_text: {nfail} failure(s)")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
