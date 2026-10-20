#!/usr/bin/env python3
"""Compares two ONCVPSP outputs quantity by quantity.

Both outputs are parsed with pyoncvpsp (and the UPF files they contain with
upf-tools), so the comparison does not depend on how the quantities are laid
out in lines: every parsed quantity of A must be in B, with the same structure,
the same text and numbers that agree within the tolerance of that quantity
(TOLERANCES below).  A number passes if it passes both

  - the base test, that of fldiff.pl -easy, with which the tests compared
    the outputs before: |a - b| <= 1.01e-5 or |a - b| <= 1.01e-5 (|a| + |b|)
    (1.01e-10, fldiff's default, with --identical), and
  - the rule for its quantity, (atol, rtol), the first in TOLERANCES that
    matches its path: |a - b| <= atol or |a - b| <= rtol (|a| + |b|), with
    the largest |a| + |b| of the whole array for radial functions; quantities
    without a rule take the base test alone,

so a rule can only tighten the comparison.  The relative differences are
fldiff's, |a - b| / (|a| + |b|), about half the usual ones.  Text is compared
with its numbers masked, and the numbers in it as numbers.  Integers are
numbers too (fldiff compared them as text, and let them differ).

The projectors and their energies are not unique (any change of basis among
the projectors of one l, or one l and j, that keeps the nonlocal kernel
sum_ij |beta_i> D_ij <beta_j| leaves the pseudopotential the same), so they are
never compared: neither the projector coefficients and the plotted projectors
of the output, nor the projectors of the psp8, UPF and PSML files.  The
pseudopotential is compared by the nonlocal kernel of each block of the files
(psp8: scalar-relativistic and spin-orbit per l; UPF: per l, or l and j; PSML:
per set, l and j),

    K(r, r') = sum_ij beta_i(r) D_ij beta_j(r'),

whose difference in the Frobenius norm on the file's grid, relative to the
norm of A's kernel (for a spin-orbit kernel, of the scalar-relativistic kernel
of its l, to which it is a correction), must be at most --kernel-tol (5e-5;
1e-10 with --identical), and by the local potential, which is compared as any
other radial function.

Never compared: the date, the build line (it names the executable) and what
follows from it (pyoncvpsp's stage), the lines of the messages (only their
names), and whatever --ignore names; without --identical, nor the numbers of
steps of the solvers (REFERENCE_IGNORE).  The PSML file is left out unless
--psml is given: its values near r = 0 depend on the compiler flags.

--identical is for outputs that must be the same (the same input read in two
ways, or run in two ways): the base test and the kernels to 1e-10 (with
-ffast-math a sum can round differently between two runs), the PSML files
too, and, unless --no-lines is given, the lines: the outputs must have the
same lines, apart from the build line and the lines naming a date, with the
same text and numbers that pass the base test, which covers what pyoncvpsp
does not parse.

Prints, for each group of quantities, the largest differences found, and each
quantity beyond its tolerance; the exit status is 1 if any quantity differs
beyond its tolerance or is in only one output, 2 if they cannot be compared.

usage: compare_output.py A.out B.out [--identical [--no-lines]] [--psml]
                         [--ignore PATH]... [--kernel-tol TOL]

PATH is a quantity as printed, with or without its indices, e.g.
input.log_mesh, upf, plot_data.wavefunctions[].wfn_ps, or a pattern of them
(fnmatch); leaving out a quantity leaves out what is in it too (--ignore upf,
the UPF file and its kernels).
"""

import argparse
import fnmatch
import math
import re
import sys
import traceback
from dataclasses import dataclass, field

import numpy as np
from pyoncvpsp.io import OncvpspTextParser, parse_psml, parse_psp8
from upf_tools import UPFDict

# (pattern, atol, rtol): the first pattern matching a quantity's path, its
# indices written [], sets the tolerance of its numbers (on top of the base
# test, so that only tightens it).  Radial functions (arrays) are compared
# relative to their largest |a| + |b|.  The rules tighten the base test where
# its absolute 1e-5, or its relative 1e-5 of large numbers, says little; each
# is three to ten times the largest difference between the references in
# tests/refs and the outputs of the code they were last checked with.
TOLERANCES = [
    # the input echoes, written by the same code from the same input
    ("input.*", 0.0, 1e-12),
    ("psp8.input", 0.0, 1e-12),
    ("upf.info.*", 0.0, 1e-12),
    ("psml.input", 0.0, 1e-12),
    # total energies (Ha; the UPF's in Ry): 1e-5 relative is 0.4 Ha for Bi;
    # they agree to 1e-3 Ha, and the terms of the all-electron energy to 4e-3 Ha
    ("ae_reference_energies.total_energy", 3e-3, 0.0),
    ("ae_reference_energies.terms.total", 3e-3, 0.0),
    ("ae_reference_energies.terms.*", 1e-2, 0.0),
    ("test_configuration_energies[].etot_*", 3e-3, 0.0),
    ("pseudoatom_total_energy", 3e-3, 0.0),
    ("upf.header.total_psenergy", 3e-3, 0.0),
    # d2Exc/dn_i dn_j (about 1e-2 Ha, 7 digits) and its rms error (1e-4 to
    # 1e-3 Ha): they agree to 3e-8 and 4e-9 Ha
    ("d2exc.*", 3e-7, 0.0),
    ("d2exc_rmse", 5e-8, 0.0),
]

# Left out of the comparison with a reference (not with --identical): the
# steps of the solvers, which depend on the way they converge (the mixing of
# the SCF), not on what they converge to, which is compared.
REFERENCE_IGNORE = ["ae_reference_energies.iterations", "ionic_tail.steps"]

NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
# record fields that name a record in the report
NAME_FIELDS = ("set", "kind", "n", "l", "j", "iproj", "jproj", "index", "seq")


@dataclass
class Stats:
    """Largest differences of the numbers of one quantity."""

    count: int = 0
    abs_diff: float = 0.0
    rel_diff: float = 0.0


@dataclass
class Report:
    base: float
    ignore: list
    failures: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)

    def ignored(self, path):
        generic = re.sub(r"\[\d+\]", "[]", path)
        return any(
            fnmatch.fnmatchcase(p, pattern + suffix)
            for pattern in self.ignore
            for suffix in ("", ".*", "[[]*")
            for p in (path, generic)
        )

    def fail(self, path, message):
        self.failures.append(f"{path}: {message}")

    def record(self, path, abs_diff, rel_diff, count):
        group = re.sub(r"\[\d+\]", "[]", path)
        stats = self.stats.setdefault(group, Stats())
        stats.count += count
        stats.abs_diff = max(stats.abs_diff, abs_diff)
        stats.rel_diff = max(stats.rel_diff, rel_diff)


def rule(path):
    """(atol, rtol) of the first rule of TOLERANCES matching path, or None."""
    generic = re.sub(r"\[\d+\]", "[]", path)
    for pattern, atol, rtol in TOLERANCES:
        if fnmatch.fnmatchcase(generic, pattern):
            return atol, rtol
    return None


def label(value):
    """The fields of a record that name it, e.g. ' (l=1, iproj=2)'."""
    if not isinstance(value, dict):
        return ""
    names = [
        f"{key}={value[key]}"
        for key in NAME_FIELDS
        if key in value and not isinstance(value[key], (dict, list, np.ndarray))
    ]
    return f" ({', '.join(names)})" if names else ""


def compare_numbers(report, path, a, b, scale=False):
    """Compare arrays (or scalars) of numbers a and b at path; with scale, a
    rule is relative to the largest |a| + |b| of the array."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        report.fail(path, f"shape {a.shape} and {b.shape}")
        return
    if a.size == 0:
        return
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    if (nan_a != nan_b).any():
        report.fail(path, "NaN in only one of them")
        return
    a, b = np.where(nan_a, 0.0, a), np.where(nan_b, 0.0, b)
    diff = np.abs(a - b)
    pointwise = np.abs(a) + np.abs(b)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(diff > 0, diff / pointwise, 0.0)
    report.record(path, float(diff.max()), float(rel.max()), a.size)
    bad = (diff > report.base) & (diff > report.base * pointwise)
    tolerance = f"{report.base:.3g} abs or rel"
    specific = rule(path)
    if specific is not None:
        atol, rtol = specific
        size = float(pointwise.max()) if scale else pointwise
        bad |= (diff > atol) & (diff > rtol * size)
        tolerance += f", and {atol:.0e} abs or {rtol:.0e} rel"
        if scale:
            tolerance += " to the largest value"
    if bad.any():
        k = int(np.argmax(np.where(bad, diff, -1.0)))
        where = f"[{k}]" if a.ndim else ""
        report.fail(
            path + where,
            f"{a.flat[k]:.10g} and {b.flat[k]:.10g} differ by {diff.flat[k]:.3e}"
            f" (tolerance {tolerance}; {int(bad.sum())} of {a.size} beyond)",
        )


def compare_text(report, path, a, b):
    if NUMBER.sub("#", a).split() != NUMBER.sub("#", b).split():
        report.fail(path, f"text differs:\n    {a[:200]!r}\n    {b[:200]!r}")
        return
    na = [fortran_float(x) for x in NUMBER.findall(a)]
    nb = [fortran_float(x) for x in NUMBER.findall(b)]
    compare_numbers(report, path, na, nb)


def is_number(value):
    return isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(
        value, (bool, np.bool_)
    )


def compare(report, path, a, b):
    """Compare the parsed values a and b at path."""
    if report.ignored(path):
        return
    if isinstance(a, dict) and isinstance(b, dict):
        for key in list(a) + [k for k in b if k not in a]:
            sub = f"{path}.{key}" if path else str(key)
            if report.ignored(sub):
                continue
            if key not in b:
                report.fail(sub, "only in A")
            elif key not in a:
                report.fail(sub, "only in B")
            else:
                compare(report, sub, a[key], b[key])
    elif isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            report.fail(path, f"{len(a)} and {len(b)} items")
        elif a and all(is_number(x) for x in a) and all(is_number(x) for x in b):
            compare_numbers(report, path, a, b)
        else:
            for k, (x, y) in enumerate(zip(a, b)):
                name = label(x)
                if name and name != label(y):
                    report.fail(f"{path}[{k}]", f"records{name} and{label(y)}")
                else:
                    compare(report, f"{path}[{k}]", x, y)
    elif isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        if a.dtype.kind in "fiu" and b.dtype.kind in "fiu":
            compare_numbers(report, path, a, b, scale=True)
        elif not np.array_equal(a, b):
            report.fail(path, "differs")
    elif is_number(a) and is_number(b):
        compare_numbers(report, path, a, b)
    elif isinstance(a, str) and isinstance(b, str):
        compare_text(report, path, a, b)
    elif a != b:
        report.fail(path, f"{a!r} and {b!r}")


def compare_lines(report, path_a, path_b, limit=10):
    """With --identical: the outputs have the same lines, apart from the build
    line and the lines naming a date, with the same text and numbers that pass
    the base test.  This covers what pyoncvpsp does not parse."""

    def lines(path):
        with open(path, encoding="latin-1") as f:
            return [
                (k + 1, line.rstrip("\n"))
                for k, line in enumerate(f)
                if not line.startswith("build: ") and "date" not in line
            ]

    lines_a, lines_b = lines(path_a), lines(path_b)
    if len(lines_a) != len(lines_b):
        report.fail("lines", f"{len(lines_a)} and {len(lines_b)} lines")
    failed = 0
    for (ka, a), (kb, b) in zip(lines_a, lines_b):
        if a == b:
            continue
        where = f"line {ka} of A, {kb} of B"
        if NUMBER.sub("#", a).split() != NUMBER.sub("#", b).split():
            report.fail(where, f"text differs:\n    {a!r}\n    {b!r}")
            failed += 1
        else:
            na = np.array([fortran_float(x) for x in NUMBER.findall(a)])
            nb = np.array([fortran_float(x) for x in NUMBER.findall(b)])
            diff = np.abs(na - nb)
            sizes = np.abs(na) + np.abs(nb)
            rel = np.where(diff > 0, diff / np.where(sizes > 0, sizes, 1.0), 0.0)
            report.record("(lines)", float(diff.max()), float(rel.max()), 1)
            if ((diff > report.base) & (diff > report.base * sizes)).any():
                report.fail(
                    where,
                    f"numbers differ beyond {report.base:.3g}:\n    {a!r}\n    {b!r}",
                )
                failed += 1
        if failed >= limit:
            report.fail("lines", f"stopped after {limit} differing lines")
            break


def fortran_float(text):
    return float(text.replace("D", "E").replace("d", "e"))


# --- the pseudopotential files, with the projectors replaced by kernels ------


@dataclass
class Block:
    """A block of projectors: their D matrix and their values on the grid."""

    dij: np.ndarray
    beta: np.ndarray
    scale_key: tuple | None = None  # the block whose norm scales this one's


def kernel_norm2(block):
    """||K||^2 = tr(D S D S), S the overlaps of the projectors on the grid."""
    s = block.beta @ block.beta.T
    return float(np.trace(block.dij @ s @ block.dij @ s))


def kernel_diff(a, b, scale_norm2):
    """||K_a - K_b|| / sqrt(scale_norm2), the kernels formed on the grid (the
    difference of the formed kernels does not lose the precision that the
    overlap form does when they agree)."""
    n = min(a.beta.shape[1], b.beta.shape[1])
    # the projectors vanish beyond their cutoff: form the kernels up to it
    support = np.flatnonzero(
        (np.abs(a.beta[:, :n]) > 0).any(axis=0)
        | (np.abs(b.beta[:, :n]) > 0).any(axis=0)
    )
    last = int(support[-1]) + 1 if support.size else 0
    ka = a.beta[:, :last].T @ a.dij @ a.beta[:, :last]
    kb = b.beta[:, :last].T @ b.dij @ b.beta[:, :last]
    diff2 = float(np.sum((ka - kb) ** 2))
    if scale_norm2 <= 0:
        return math.sqrt(diff2)
    return math.sqrt(diff2 / scale_norm2)


def psp8_part(text):
    psp8 = parse_psp8(text)
    blocks = {}
    for block in psp8.pop("projectors"):
        key = ("psp8", block["kind"], block["l"])
        scale = (
            ("psp8", "scalar", block["l"]) if block["kind"] == "spin_orbit" else None
        )
        blocks[key] = Block(np.diag(block["ekb"]), block["beta"], scale)
    psp8.pop("pspd")  # the date
    return psp8, blocks


def plain(value):
    """upf-tools' OrderedDicts as dicts."""
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, list):
        return [plain(item) for item in value]
    return value


def upf_part(text):
    upf = plain(UPFDict.from_str(text))
    nonlocal_ = upf.pop("nonlocal", {})
    betas = nonlocal_.get("beta", [])
    if isinstance(betas, dict):
        betas = [betas]
    nbeta = len(betas)
    dij = np.asarray(nonlocal_.get("dij", []), dtype=float).reshape(nbeta, nbeta)
    jjj = {}
    spin_orb = upf.get("spin_orb", {})
    relbeta = spin_orb.get("relbeta", []) if isinstance(spin_orb, dict) else []
    for rel in relbeta if isinstance(relbeta, list) else [relbeta]:
        jjj[int(rel["index"])] = float(rel["jjj"])
    groups = {}
    for k, beta in enumerate(betas):
        index = int(beta.get("index", k + 1))
        key = ("upf", int(beta["angular_momentum"]), jjj.get(index))
        groups.setdefault(key, []).append(k)
    blocks = {}
    for key, members in groups.items():
        beta = np.array([betas[k]["content"] for k in members])
        blocks[key] = Block(dij[np.ix_(members, members)], beta)
    upf.get("header", {}).pop("date", None)
    return upf, blocks


def psml_part(text):
    psml = parse_psml(text)
    groups = {}
    for proj in psml.pop("projectors"):
        key = ("psml", proj["set"], proj["l"], proj["j"])
        groups.setdefault(key, []).append(proj)
    blocks = {}
    for key, projs in groups.items():
        scale = None
        if key[1] == "spin_orbit":
            scale = ("psml", "scalar_relativistic", key[2], None)
        blocks[key] = Block(
            np.diag([p["ekb"] for p in projs]),
            np.array([p["values"] for p in projs]),
            scale,
        )
    for provenance in psml["provenance"]:
        provenance.pop("date", None)
    return psml, blocks


def load(path, with_psml):
    """The parsed output of path, ready to compare, and its kernel blocks."""
    parser = OncvpspTextParser(path=path, check=False)
    tree = dict(parser.output)
    # the build line names the executable; pyoncvpsp's stage follows from it
    info = dict(tree.pop("program_information"))
    info.pop("date", None)
    info.pop("build", None)
    # the header of the fully-relativistic program, without a build line
    if info.get("relativistic_treatment") == "relativistic":
        info["relativistic_treatment"] = "fully-relativistic"
    tree["program_information"] = info
    tree.pop("stage", None)
    # messages by name: their lines and line numbers depend on the layout
    tree["errors"] = [message["name"] for message in tree["errors"]]
    tree["warnings"] = [message["name"] for message in tree["warnings"]]
    tree["test_configuration_status"] = [
        {**status, "messages": [m["name"] for m in status["messages"]]}
        for status in tree["test_configuration_status"]
    ]
    # projectors and their energies: compared by the kernels instead
    tree.pop("vkb_scalar_projector_coefficients", None)
    tree.pop("vkb_spin_orbit_projector_coefficients", None)
    tree["plot_data"] = {
        key: value
        for key, value in tree["plot_data"].items()
        if key != "vkb_projectors"
    }
    blocks = {}
    psp8_text, upf_text, psml_text = (tree.pop(key) for key in ("psp8", "upf", "psml"))
    if psp8_text:
        tree["psp8"], part = psp8_part(psp8_text)
        blocks.update(part)
    if upf_text:
        tree["upf"], part = upf_part(upf_text)
        blocks.update(part)
    if psml_text and with_psml:
        tree["psml"], part = psml_part(psml_text)
        blocks.update(part)
    return tree, blocks


def compare_kernels(report, blocks_a, blocks_b, tol):
    """Compare the kernels of the blocks; returns the lines to print."""
    lines = []
    keys_a = [k for k in blocks_a if not report.ignored(kernel_path(k))]
    keys_b = [k for k in blocks_b if not report.ignored(kernel_path(k))]
    for key in sorted(set(keys_a) ^ set(keys_b), key=str):
        report.fail(kernel_path(key), "only in " + ("A" if key in blocks_a else "B"))
    for key in [k for k in keys_a if k in blocks_b]:
        a, b = blocks_a[key], blocks_b[key]
        path = kernel_path(key)
        if a.beta.shape[0] != b.beta.shape[0]:
            report.fail(path, f"{a.beta.shape[0]} and {b.beta.shape[0]} projectors")
            continue
        scale_block = blocks_a.get(a.scale_key) if a.scale_key else None
        scale = kernel_norm2(scale_block if scale_block is not None else a)
        diff = kernel_diff(a, b, scale)
        bad = diff > tol
        lines.append(f"  {path:42s} {diff:9.2e}" + ("   > tolerance" if bad else ""))
        if bad:
            report.fail(path, f"relative difference {diff:.3e} > {tol:.0e}")
    return lines


def kernel_path(key):
    names = {"psp8": ("kind", "l"), "upf": ("l", "j"), "psml": ("set", "l", "j")}[
        key[0]
    ]
    fields = ", ".join(f"{n}={v}" for n, v in zip(names, key[1:]) if v is not None)
    return f"{key[0]}.kernel({fields})"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n", 1)[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("\n\n", 1)[1],
    )
    parser.add_argument("a", help="the reference output")
    parser.add_argument("b", help="the output to check")
    parser.add_argument(
        "--identical",
        action="store_true",
        help="outputs that must be the same: base test 1e-10, kernels 1e-10, the"
        " PSML files and the lines too",
    )
    parser.add_argument(
        "--psml", action="store_true", help="compare the PSML files too"
    )
    parser.add_argument(
        "--no-lines",
        action="store_true",
        help="with --identical, compare the parsed quantities only, not the lines",
    )
    parser.add_argument(
        "--ignore", action="append", default=[], metavar="PATH", help="leave out PATH"
    )
    parser.add_argument("--kernel-tol", type=float, help="tolerance of the kernels")
    args = parser.parse_args()

    base = 1.01e-10 if args.identical else 1.01e-5
    kernel_tol = args.kernel_tol or (1e-10 if args.identical else 5e-5)
    ignore = args.ignore + ([] if args.identical else REFERENCE_IGNORE)
    report = Report(base, ignore)
    with_psml = args.psml or args.identical
    tree_a, blocks_a = load(args.a, with_psml)
    tree_b, blocks_b = load(args.b, with_psml)
    compare(report, "", tree_a, tree_b)
    kernel_lines = compare_kernels(report, blocks_a, blocks_b, kernel_tol)
    if args.identical and not args.no_lines:
        compare_lines(report, args.a, args.b)

    print(f"A: {args.a}\nB: {args.b}")
    print(
        f"base test: {base:.3g} abs or rel; kernels {kernel_tol:.0e} rel"
        + (f"; ignored: {', '.join(ignore)}" if ignore else "")
    )
    print(
        "\nBy quantity: how many numbers, their largest |a - b| and"
        " |a - b| / (|a| + |b|):"
    )
    for group in sorted(report.stats):
        stats = report.stats[group]
        if stats.count:
            print(
                f"  {group:60s} {stats.count:7d} {stats.abs_diff:9.2e} {stats.rel_diff:9.2e}"
            )
    print("\nNonlocal kernels, relative difference:")
    print("\n".join(kernel_lines) if kernel_lines else "  no psp8, UPF or PSML file")
    if report.failures:
        print(f"\n{len(report.failures)} differences beyond the tolerances:")
        for failure in report.failures:
            print(f"  {failure}")
        sys.exit(1)
    print("\nNo difference beyond the tolerances.")


if __name__ == "__main__":
    try:
        main()
    except Exception:  # a failure to compare, not a difference: exit status 2
        traceback.print_exc()
        sys.exit(2)
