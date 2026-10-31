"""The arrays of the library (python/oncvpsp) against a reference output of
tests/refs, number by number: each field against the numbers the
executables print of it, so that the arrays are checked by themselves, not
only through the text the library writes.

  - the DIAGNOSTICS block, the tests to 17 digits (diagout_m): the
    convergence of the pseudo-waves (tests.convergence: ecut, opt,
    well_share, ekin_tail, ekin_check), of the local energy
    (tests.vloc_cut), the log derivatives (tests.logder), the ghosts
    (tests.ghosts), the eigenvalues of the B matrices (tests.bbeig) and
    the test configurations (tests.configs);
  - the all-electron atom: its eigenvalues (the input echo), energy terms
    and iterations (ae.*), and the pseudo-atom's total energy (epstot);
  - the radial functions of the plots, on the points plotted: the valence
    charge, the semi-local and local potentials (!p, !L: pot.rho,
    pot.vpuns) and the core charges (!r: rhoc, pot.rhomod).

The projectors and their energies are not compared: they are not unique
(compare_output.py); the kernel they make is, through the psp8 file of
python_wrapper.py. The tolerances are those of compare_output.py, with
which the reference tests compare the outputs: the base test (fldiff -easy),
the radial functions relative to their largest value, and its rules for the
energies and eigenvalues, the noise of the compiler flags and of the changes
since the references were made.

usage: library_arrays.py IREL INPUT.toml REF.out
(ONCVPSP_LIBRARY the library, the package on PYTHONPATH)
"""

import re
import sys

import numpy as np

import oncvpsp

# the base test of compare_output.py (fldiff -easy)
BASE = 1.01e-5


class Checks:
    """the comparisons, the largest difference of each group, the failures"""

    def __init__(self):
        self.worst = {}
        self.failed = []

    def close(self, group, got, ref, atol=0.0, rtol=0.0, scale=False):
        """got against ref within the base test and, if given, the rule
        (atol, rtol), as compare_output.py: |a - b| <= atol or <= rtol
        (|a| + |b|), with the largest |a| + |b| of the array if scale"""
        got = np.atleast_1d(np.asarray(got, dtype=float))
        ref = np.atleast_1d(np.asarray(ref, dtype=float))
        if got.shape != ref.shape:
            self.failed.append(f"{group}: shape {got.shape}, ref {ref.shape}")
            return
        diff = np.abs(got - ref)
        size = np.abs(got) + np.abs(ref)
        if scale:
            size = np.full_like(size, size.max(initial=0.0))
        ok = (diff <= BASE) | (diff <= BASE * size)
        if atol or rtol:
            ok &= (diff <= atol) | (diff <= rtol * size)
        worst = self.worst.setdefault(group, [0, 0.0])
        worst[0] += diff.size
        worst[1] = max(worst[1], diff.max(initial=0.0))
        for ii in np.flatnonzero(~ok):
            self.failed.append(f"{group}[{ii}]: {got.flat[ii]!r}, ref {ref.flat[ii]!r}")

    def same(self, group, got, ref):
        self.close(group, got, ref, atol=0.0, rtol=1e-15)


def records(text):
    """the DIAGNOSTICS block: its records by tag, each a list of numbers"""
    block = text.partition("Begin DIAGNOSTICS\n")[2].partition("End DIAGNOSTICS")[0]
    recs = {}
    for line in block.splitlines():
        if line and not line.startswith("#"):
            tag, *values = line.split()
            recs.setdefault(tag, []).append([float(vv) for vv in values])
    return recs


def plotted(text, tag):
    """the lines of a plot tag (!p, !L, !r) as rows of numbers"""
    rows = [
        line.split()[1:]
        for line in text.splitlines()
        if line.lstrip().startswith(tag + " ")
    ]
    return np.array(rows, dtype=float)


def points(chk, tag, rr, radii):
    """the mesh points of the radii plotted (f12.7), a failure for those
    off the mesh"""
    idx = np.searchsorted(rr, radii - 1e-7)
    idx = np.minimum(idx, len(rr) - 1)
    lower = np.maximum(idx - 1, 0)
    idx = np.where(np.abs(rr[lower] - radii) < np.abs(rr[idx] - radii), lower, idx)
    for ii in np.flatnonzero(np.abs(rr[idx] - radii) > 6e-8):
        chk.failed.append(f"{tag}: the radius {radii[ii]} off the mesh")
    return idx


def main():
    irel, toml, ref = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    text = open(ref).read()
    run = oncvpsp.run(toml, irel, text=False)
    dims = run.dims
    chk = Checks()
    recs = records(text)
    lmax, lloc, nproj = dims.lmax, dims.lloc, dims.nproj
    # the entries diagout prints: l to lmax, ikap 1 for l = 0, the
    # projectors of nproj (none for the local potential's l)
    nkap = dims.nkap
    kaps = np.zeros((6, 2), bool)
    kaps[: lmax + 1, :nkap] = True
    kaps[0, 1] = False
    prjs = np.zeros((5, 6, 2), bool)
    for ll in range(lmax + 1):
        prjs[: nproj[ll], ll] = kaps[ll]

    # the convergence of the pseudo-waves: (Ecut, err) where E_resid
    # crosses each level, at qcut (E_resid, the regularization's term, e1,
    # |f1|), the share of a well's tail, the tail and the kinetic check
    cross = run["tests.convergence.cross"]
    got, want = [], []
    for ll, ikap, iprj, err, ecut in recs.get("ecut", []):
        ll, ikap, iprj = int(ll), int(ikap), int(iprj)
        levels = cross[1, :, iprj - 1, ll, ikap - 1]
        jj = int(np.argmin(np.abs(levels - err)))
        got.append(cross[:, jj, iprj - 1, ll, ikap - 1])
        want.append([ecut, err])
    nlevel = np.count_nonzero(
        cross[0] * prjs[None, :, : cross.shape[3], : cross.shape[4]]
    )
    chk.same("tests.convergence.cross levels", nlevel, len(want))
    chk.close("tests.convergence.cross", got, want)
    opt = run["tests.convergence.opt"]
    for ll, ikap, iprj, *values in recs.get("opt", []):
        at = (slice(None), int(iprj) - 1, int(ll), int(ikap) - 1)
        chk.close("tests.convergence.opt", opt[at], values)
    wshare = run["tests.convergence.wshare"]
    for ll, ikap, iprj, *values in recs.get("well_share", []):
        at = (slice(None), int(iprj) - 1, int(ll), int(ikap) - 1)
        chk.same("tests.convergence.wshare flag", wshare[2][at[1:]], 1.0)
        chk.close("tests.convergence.wshare", wshare[:2][at], values)
    tail = run["tests.convergence.tail"]
    for kk, tag in enumerate(("ekin_tail", "ekin_check")):
        for ll, ikap, iprj, value in recs.get(tag, []):
            chk.close(
                f"tests.convergence.tail {tag}",
                tail[kk, int(iprj) - 1, int(ll), int(ikap) - 1],
                value,
            )

    # the local energy: (Ecut, err) per level and the transforms' error
    vcut = run["tests.vloc_cut.cut"]
    want = [[ecut, err] for err, ecut in recs.get("ecut_vloc", [])]
    chk.close("tests.vloc_cut.cut", vcut[:, : len(want)].T, want)
    chk.same("tests.vloc_cut.cut levels", np.count_nonzero(vcut[0]), len(want))
    chk.close(
        "tests.vloc_cut.check", run["tests.vloc_cut.check"], recs["vloc_check"][0]
    )

    # the log derivatives: the radii, the energies and atan(r u'/u)
    ne = run["tests.logder.ne"]
    if ne > 0:
        rl = run["tests.logder.rl"]
        ee = run["tests.logder.ee"]
        ae, ps = run["tests.logder.ae"], run["tests.logder.ps"]
        for ll, ikap, radius in recs["ld_r"]:
            chk.close("tests.logder.rl", rl[int(ll)], radius)
        rows = np.array(recs["ld"])
        for ll, ikap in sorted({(int(a), int(b)) for a, b in rows[:, :2]}):
            sel = rows[(rows[:, 0] == ll) & (rows[:, 1] == ikap)]
            # the grid epsh2 - depsh*(i - 1) rounds differently where the
            # compiler fuses it (-march=native: -6.2e-17 for 0, a few ulps
            # near it), so to round-off of its largest energy
            chk.close("tests.logder.ee", ee, sel[:, 2], rtol=1e-15, scale=True)
            chk.close("tests.logder.ae", ae[:, ll, ikap - 1], sel[:, 3])
            chk.close("tests.logder.ps", ps[:, ll, ikap - 1], sel[:, 4])
    else:
        chk.same("tests.logder (none)", len(recs.get("ld", [])), 0)

    # the ghosts: the bound states (n, e_ref, e_basis) and the basis state
    # of least <r>/rc
    nrow, nn = run["tests.ghosts.nrow"], run["tests.ghosts.nn"]
    gee, least = run["tests.ghosts.ee"], run["tests.ghosts.least"]
    seen = {}
    for ll, ikap, n, *values in recs.get("ghost", []):
        ll, ikap = int(ll), int(ikap)
        jj = seen.get((ll, ikap), 0)
        seen[(ll, ikap)] = jj + 1
        chk.same("tests.ghosts.nn", nn[jj, ll, ikap - 1], n)
        chk.close("tests.ghosts.ee", gee[:, jj, ll, ikap - 1], values)
    for (ll, ikap), count in seen.items():
        chk.same("tests.ghosts.nrow", nrow[ll, ikap - 1], count)
    shown = kaps[: nrow.shape[0], : nrow.shape[1]]
    chk.same("tests.ghosts.nrow total", nrow[shown].sum(), sum(seen.values()))
    nleast = 0
    for ll, ikap, *values in recs.get("ghost_r", []):
        chk.close("tests.ghosts.least", least[:, int(ll), int(ikap) - 1], values)
        nleast += 1
    chk.same("tests.ghosts.least count", np.count_nonzero(least[1][shown] > 0), nleast)

    # the eigenvalues of the B matrices and their condition
    bbeig = run["tests.bbeig"]
    for ll, ikap, iprj, value in recs.get("vkb_b", []):
        chk.close("tests.bbeig", bbeig[int(iprj) - 1, int(ll), int(ikap) - 1], value)
    for ll, ikap, value in recs.get("vkb_cond", []):
        bb = np.abs(bbeig[: nproj[int(ll)], int(ll), int(ikap) - 1])
        chk.close("tests.bbeig cond", bb.max() / bb.min(), value)

    # the test configurations: the energy differences and the excitation
    # error, each row's n, l, ikap, occupation and eigenvalues
    configs = recs.get("test_config", [])
    chk.same("tests.configs.ncfg", run["tests.configs.ncfg"], len(configs))
    de, crow = run["tests.configs.de"], run["tests.configs.nrow"]
    state, cee = run["tests.configs.state"], run["tests.configs.ee"]
    for icfg, *values in configs:
        chk.close("tests.configs.de", de[:, int(icfg)], values)
    rows = {}
    for icfg, n, ll, ikap, *values in recs.get("test_eig", []):
        icfg = int(icfg)
        ii = rows.get(icfg, 0)
        rows[icfg] = ii + 1
        chk.same("tests.configs.state", state[:, ii, icfg], [n, ll, ikap])
        chk.same("tests.configs.ee occupation", cee[0, ii, icfg], values[0])
        chk.close("tests.configs.ee", cee[1:, ii, icfg], values[1:], atol=1e-7)
    for icfg, count in rows.items():
        chk.same("tests.configs.nrow", crow[icfg], count)

    # the all-electron atom: the eigenvalues of the input echo (8 digits),
    # the energy terms (10 decimals) and the SCF's iterations (one step)
    echo = text.partition("#   n    l    f        energy (Ha)\n")[2]
    eig = []
    for line in echo.splitlines():
        if line.startswith("#"):
            break
        eig.append(float(line.split()[3].replace("D", "E")))
    ea = run["ae.ea"]
    chk.close("ae.ea", ea[: len(eig), 0], eig, rtol=1e-7)
    terms = text.partition("energy terms (Ha)\n")[2].splitlines()[:5]
    terms = [float(line.split()[-1]) for line in terms]
    chk.close("ae.etot", run["ae.etot"], terms[0], atol=1e-8)
    chk.close("ae.eterm", run["ae.eterm"], terms[1:], atol=1e-8)
    its = int(re.search(r"iterations +(\d+)", text).group(1))
    # (the base test does not allow a step, compare_output.py's STEP_COUNTS)
    chk.same("ae.it", min(abs(run["ae.it"] - its), 1), abs(run["ae.it"] - its))
    eps = float(re.search(r"Pseudoatom total energy *(\S+)", text).group(1))
    chk.close("epstot", run["epstot"], eps, atol=2e-6)

    # the radial functions plotted (f12.7), on the points plotted, each
    # relative to its largest value
    rr = run["pot.rr"]
    rho, vpuns = run["pot.rho"], run["pot.vpuns"]
    plot = plotted(text, "!p")
    at = points(chk, "!p", rr, plot[:, 0])
    chk.close("pot.rho (!p)", rho[at], plot[:, 1], scale=True)
    chk.close("pot.vpuns (!p)", vpuns[at, : lmax + 1], plot[:, 2:], scale=True)
    if lloc == 4:
        plot = plotted(text, "!L")
        at = points(chk, "!L", rr, plot[:, 0])
        chk.close("pot.vpuns local (!L)", vpuns[at, lloc], plot[:, 1], scale=True)
    plot = plotted(text, "!r")
    at = points(chk, "!r", rr, plot[:, 0])
    rhoc, rhomod = run["rhoc"], run["pot.rhomod"]
    chk.close("pot.rho (!r)", rho[at], plot[:, 1], scale=True)
    # the core charge, capped where it is above the plot's range
    cut = plot[:, 2] < rhoc[at] - 1e-6
    chk.close("rhoc (!r)", rhoc[at][~cut], plot[~cut, 2], scale=True)
    chk.close("pot.rhomod (!r)", rhomod[at, 0], plot[:, 3], scale=True)

    for group, (count, worst) in sorted(chk.worst.items()):
        print(f"{group:40s} {count:6d}  max |diff| {worst:.3e}")
    for failure in chk.failed:
        print("FAILED:", failure, file=sys.stderr)
    if not chk.worst:
        print("FAILED: nothing compared", file=sys.stderr)
        return 1
    return 1 if chk.failed else 0


if __name__ == "__main__":
    sys.exit(main())
