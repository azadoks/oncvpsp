#!/usr/bin/env python3
"""Compares the all-electron atom of ONCVPSP with the NIST atomic reference data.

The data (tests/data/nist, see its README.md) are total energies, their terms and
eigenvalues of the neutral atoms H-U and their first cations, for three sets:
LDA (non-relativistic), ScRLDA (scalar-relativistic) and RLDA (Dirac).  NIST
uses Slater exchange and VWN correlation, with the relativistic exchange of
MacDonald and Vosko in ScRLDA and RLDA, and the speed of light 137.0359895.
The inputs written here set exactly that: iexc=6 for LDA, iexc=7 and c_light
for ScRLDA and RLDA, on a log mesh fine enough for the stated accuracy of the
data.  Only the all-electron program is run (oncvpspnr_ae.x for LDA,
oncvpsp_ae.x for ScRLDA, oncvpspr_ae.x for RLDA); the pseudopotential part of
the input is a placeholder that it does not use.

A case passes when it converges and its total energy and every eigenvalue NIST
gives are within NIST's stated accuracy (TOL_* below), against the data
exactly as NIST gives them.  NIST states no accuracy for the energy terms;
their largest difference is printed but not tested.  The cases in
EXPECTED_FAILURES miss the stated accuracy by a little, as other codes do; they
are expected to fail, and one that passes is reported (it should leave the
list).  The script prints one line per case and exits non-zero if a case does
not do what is expected of it.

usage: nist_test.py EXE SET DATA_DIR WORK_DIR [--z 1-18,29,...] [--neutrals]
       [--jobs N] [--per-decade N]
"""

import argparse
import concurrent.futures
import pathlib
import re
import struct
import subprocess
import sys

# The accuracy NIST states for the data (Ha): "numerical convergence to the
# target accuracy of 1 microHartree in total energy" (Procedure), and "all
# individual orbital eigenvalues agreed to within 2 microHartrees"
# (Reliability)
TOL_ETOT = 1.0e-6
TOL_EIGENVALUE = 2.0e-6

# the cases that miss the stated accuracy, at the default mesh: (set, Z,
# charge).  ScRLDA lanthanides and Fr: Etot +1.0e-6 to +1.3e-6, where the
# solvers of RadialAtomicSolvers.jl miss the lanthanides by 1.1e-6 to 1.3e-6
# too, and Lu 4s -2.2e-6.  Which of them pass at a few 1e-7 under the limit
# moves with the eigenvalue tolerance and the mixing of the SCF.
EXPECTED_FAILURES = {
    ("ScRLDA", 62, "neutral"),
    ("ScRLDA", 62, "cation"),
    ("ScRLDA", 65, "neutral"),
    ("ScRLDA", 68, "neutral"),
    ("ScRLDA", 69, "neutral"),
    ("ScRLDA", 70, "cation"),
    ("ScRLDA", 71, "neutral"),
    ("ScRLDA", 87, "neutral"),
}

# speed of light of the NIST data (CODATA 1986)
C_LIGHT_NIST = 137.0359895

# iexc of each set, and whether it needs c_light
SETS = {"LDA": (6, False), "ScRLDA": (7, True), "RLDA": (7, True)}

# labels that NIST's files get wrong: (set, Z, charge, label) -> label.  The
# ScRLDA file of neutral Pm labels its 3d eigenvalue "32".
LABEL_ERRATA = {("ScRLDA", 61, "neutral", "32"): "3d"}

# the NIST names of the energy terms
TERMS = ("Etot", "Ekin", "Eenuc", "Ecoul", "Exc")

NOBLE = {
    "He": "1s2",
    "Ne": "[He] 2s2 2p6",
    "Ar": "[Ne] 3s2 3p6",
    "Kr": "[Ar] 3d10 4s2 4p6",
    "Xe": "[Kr] 4d10 5s2 5p6",
    "Rn": "[Xe] 4f14 5d10 6s2 6p6",
}
LETTERS = "spdf"

# the input: the configuration, every state but the last as core, and a
# placeholder pseudopotential
TEMPLATE = """# ATOM AND REFERENCE CONFIGURATION
# atsym, z, nc, nv, iexc, psfile[, c_light]
  {atsym} {z}.0 {nc} 1 {iexc} none{c_light}
# n, l, f
{states}
# PSEUDOPOTENTIAL AND OPTIMIZATION
 0
 0 2.0 0.0 4 8 6.0
# LOCAL POTENTIAL
 4 5 2.0 0.0
# VANDERBILT-KLEINMAN-BYLANDER PROJECTORS
 0 2 1.0
# MODEL CORE CHARGE
 0 0.0 0.0
# LOG DERIVATIVE ANALYSIS
 -12.0 12.0 0.02
# OUTPUT GRID
 6.0 0.01
# TEST CONFIGURATIONS
 0
# LOG MESH: amesh, rr1, rmax
 {amesh!r} {rr1!r} {rmax!r}
"""


def log_mesh(z, per_decade):
    """amesh, rr1, rmax of a mesh fine enough for the stated accuracy.

    rr1 = 1e-8/z resolves the nuclear cusp (the built-in rr1 does not, for the
    non-relativistic equation), and rmax = 40 bohr holds the most diffuse
    cation and neutral states.  The error of the total energy falls as the
    fourth power of the spacing: for U, 3e-5, 1e-6 and 8e-8 Ha at 300, 600 and
    1200 points per decade of r.
    """
    rr1 = 1.0e-8 / z
    rmax = 40.0
    return 10.0 ** (1.0 / per_decade), rr1, rmax


def expand(configuration):
    """[(n, l, f)] of a configuration such as '[Ar] 3d10 4s1'."""
    states = []
    for token in configuration.split():
        noble = re.fullmatch(r"\[(\w+)\]", token)
        if noble:
            states += expand(NOBLE[noble[1]])
        else:
            state = re.fullmatch(r"(\d)([spdf])(\d+)", token)
            states.append((int(state[1]), LETTERS.index(state[2]), int(state[3])))
    return states


def read_data(path, set_name):
    """{(Z, charge): case} of a data file, a case being a dict with symbol,
    configuration, terms ({name: value}) and eigenvalues ({label: value})."""
    cases = {}
    case = None
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        header = re.fullmatch(r"(\d+) (\w+) (neutral|cation) (.+)", line)
        if header:
            case = {
                "symbol": header[2],
                "configuration": header[4],
                "terms": {},
                "eigenvalues": {},
            }
            key = (int(header[1]), header[3])
            cases[key] = case
            continue
        tokens = line.split()
        if tokens[1] == "=":
            case["terms"][tokens[0]] = float(tokens[2])
        else:
            label = LABEL_ERRATA.get((set_name, *key, tokens[0]), tokens[0])
            case["eigenvalues"][label] = float(tokens[1])
    return cases


def read_ae_file(path):
    """terms ({NIST name: value}) and eigenvalues ({NIST label: value}) of the
    all-electron data file, at full precision (the output prints eigenvalues
    to 8 digits): the layout of src/ae_file_m.f90, native byte order."""
    buf = path.read_bytes()
    offset = 0

    def take(fmt):
        nonlocal offset
        values = struct.unpack_from("=" + fmt, buf, offset)
        offset += struct.calcsize("=" + fmt)
        return values

    _, version, irel = take("8sii")
    if version != 3:
        raise ValueError(f"{path}: all-electron file version {version}, not 3")
    _, _, nc, nv, _ = take("2sdiii")
    nt = nc + nv
    na, la = take(f"{nt}i"), take(f"{nt}i")
    take(f"{nt}d")
    (mmax,) = take("i")
    take(f"{mmax}d")
    take("d")  # c_light
    _, etot, ekin, enuc, ecoul, exc, _ = take("i6d")
    nkap = 2 if irel == 2 else 1
    ea = take(f"{nt * nkap}d")
    terms = {"Etot": etot, "Ekin": ekin, "Eenuc": enuc, "Ecoul": ecoul, "Exc": exc}
    eigenvalues = {}
    for ii in range(nt):
        label = f"{na[ii]}{LETTERS[la[ii]]}"
        if irel != 2:
            eigenvalues[label] = ea[ii]
        else:
            # kap = -(l+1), j = l + 1/2; kap = l, j = l - 1/2 for l > 0
            eigenvalues[label + "P"] = ea[ii]
            if la[ii] > 0:
                eigenvalues[label + "M"] = ea[nt + ii]
    return {"terms": terms, "eigenvalues": eigenvalues}


def run_case(exe, set_name, z, charge, case, work, per_decade):
    """Runs one case; returns (passed, message)."""
    iexc, needs_c = SETS[set_name]
    states = expand(case["configuration"])
    amesh, rr1, rmax = log_mesh(z, per_decade)
    text = TEMPLATE.format(
        atsym=case["symbol"],
        z=z,
        nc=len(states) - 1,
        iexc=iexc,
        c_light=f" {C_LIGHT_NIST}" if needs_c else "",
        states="\n".join(f" {n} {l} {f}.0" for n, l, f in states),
        amesh=amesh,
        rr1=rr1,
        rmax=rmax,
    )
    name = f"{z:02d}{case['symbol']}_{charge}"
    directory = work / name
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "ae.dat").write_text(text)
    with open(directory / "ae.dat") as stdin, open(directory / "ae.out", "w") as stdout:
        subprocess.run(
            [exe, "-a", "ae.bin"],
            stdin=stdin,
            stdout=stdout,
            cwd=directory,
            check=False,
        )
    output = (directory / "ae.out").read_text(errors="replace")
    if "All-electron data written" not in output or "ERROR" in output:
        return False, f"{name:14s} FAILED: no converged result (see {directory}/ae.out)"
    result = read_ae_file(directory / "ae.bin")

    failures = []
    d_etot = result["terms"]["Etot"] - case["terms"]["Etot"]
    if abs(d_etot) > TOL_ETOT:
        failures.append(f"Etot {d_etot:+.1e}")
    d_term = max(abs(result["terms"][t] - case["terms"][t]) for t in TERMS[1:])
    d_eig, d_label = 0.0, ""
    for label, ref in case["eigenvalues"].items():
        if label not in result["eigenvalues"]:
            failures.append(f"{label} not computed")
            continue
        delta = result["eigenvalues"][label] - ref
        if abs(delta) > abs(d_eig):
            d_eig, d_label = delta, label
        if abs(delta) > TOL_EIGENVALUE:
            failures.append(f"{label} {delta:+.1e}")
    summary = (
        f"{name:14s} dEtot {d_etot:+.1e}  max|dterm| {d_term:.1e}"
        f"  max|deig| {abs(d_eig):.1e} ({d_label})"
    )
    if failures:
        return False, summary + "  FAILED: " + ", ".join(failures)
    return True, summary


def parse_z(text):
    zs = []
    for part in text.split(","):
        first, _, last = part.partition("-")
        zs += range(int(first), int(last or first) + 1)
    return zs


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("exe", help="the all-electron program of the set")
    parser.add_argument("set", choices=SETS)
    parser.add_argument("data", type=pathlib.Path, help="tests/data/nist")
    parser.add_argument("work", type=pathlib.Path, help="directory for the runs")
    parser.add_argument("--z", default="1-92", help="atomic numbers, e.g. 1-18,29")
    parser.add_argument("--neutrals", action="store_true", help="no cations")
    parser.add_argument("--jobs", type=int, default=1, help="runs at a time")
    parser.add_argument(
        "--per-decade", type=int, default=2000, help="log mesh points per decade of r"
    )
    args = parser.parse_args()

    data = read_data(args.data / f"{args.set}.dat", args.set)
    charges = ("neutral",) if args.neutrals else ("neutral", "cation")
    cases = [(z, c) for z in parse_z(args.z) for c in charges if (z, c) in data]
    exe = str(pathlib.Path(args.exe).resolve())
    work = args.work.resolve()
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        results = list(
            pool.map(
                lambda zc: run_case(
                    exe, args.set, zc[0], zc[1], data[zc], work, args.per_decade
                ),
                cases,
            )
        )
    unexpected = 0
    for (z, charge), (passed, message) in zip(cases, results):
        expected = (args.set, z, charge) not in EXPECTED_FAILURES
        if passed != expected:
            unexpected += 1
        if not expected:
            message += "  (unexpected pass)" if passed else "  (expected)"
        print(message)
    npass = sum(passed for passed, _ in results)
    print(
        f"{args.set}: {npass} of {len(results)} cases within"
        f" {TOL_ETOT:.0e} (Etot) and {TOL_EIGENVALUE:.0e} (eigenvalues) Ha;"
        f" {unexpected} not as expected"
    )
    sys.exit(1 if unexpected else 0)


if __name__ == "__main__":
    main()
