#!/usr/bin/env python3
"""Measures how smoothly the pseudopotential depends on one input parameter.

Runs the program on a TOML input with one parameter (KEY, e.g.
pseudopotentials.rc[0]) stepped through NPTS equally spaced values over a
small window, takes from each output the features below, fits each feature
with a polynomial of degree DEGREE in the parameter by least squares, and
prints the largest residual of each group relative to its scale, and the
ratio of the largest step between neighbouring values to the median step
(about 1 for a smooth function, large for a jump, infinite for a staircase).
A feature that is a smooth function of the parameter over the window is
fitted to round-off (or to the window's truncation error, kept small by the
window's width); a jump or a kink in it, such as one from rounding the
parameter to the log mesh, leaves a residual of the order of its size. With
--tol (one value, or group=value pairs), the exit status is 1 if a residual
exceeds it or a jump ratio exceeds --max-jump. A step no larger than the
group's --tol is not a jump: a feature that does not depend on the parameter
at all (an optimum the search lands on from wherever it starts) has only
noise for steps, and their ratio is ~3 for noise alone.

The features come from the UPF file, written on the log mesh itself
(upf_mesh, every few points of it) and read with upf-tools, and from the
output read with pyoncvpsp (compare_output.py), so they are the generator's own values
with no interpolation onto another grid in between; the UPF writer gives 17
digits. With --output psp8 or upf they come instead from the psp8 file (14
digits) or the UPF file on the linear output mesh, onto which the writers
evaluate the functions (doc/smoothness_plan.md, O-1). They are what a
calculation uses, never the projectors themselves (they are not unique; see
kernel_diff.py):

  kernel  K(r, r') = sum_ij beta_i(r) D_ij beta_j(r') of each projector block
          (each l, or l and j; in psp8, scalar-relativistic and spin-orbit
          per l), at all pairs of the sample radii, relative to the largest
          |K| over the pairs
  vloc    the local potential at the sample radii, relative to its largest
          value there
  core    the model core charge at the sample radii, relative to its largest
          value there; from the psp8 file also its four derivatives, each
          relative to its own largest value
  tests   the pseudo-atom eigenvalues of the test configurations (printed to
          8 decimals, so only a coarser --tol is meaningful for them)
  ld      the pseudopotential's log-derivative phase at rxpsh, each l and
          energy, as cos and sin of twice it: the phase is defined only up
          to pi, and the doubled angle has no branch
  ecut    the cutoffs of each projector's convergence table, at its
          energy errors 1e-2..1e-5 Ha, relative to the largest in the table
          (both from the output's DIAGNOSTICS block, 17 digits; without one,
          from the tables printed to 6 decimals and 0.01 Ha)
  ghost   the energies of the bound-ghost test, the radial equation's and
          the basis's (printed to 6 decimals)

The last three are what the optimizer reads (doc/smoothness_plan.md, order
of the remaining work); until O-6 prints them in full, only a --tol at
their print resolution is meaningful.

The sample radii are the points of the file's mesh nearest to every --sample
bohr. For a parameter that is a radius (a key named rc, rcloc or rxpsh), the
radii inside the scanned window (widened by two steps of the file's log mesh,
or of the generator's for a linear one), and inside 0.95 times it (where the
projectors' taper starts), are dropped: there the functions are only as
smooth as the matching at the radius makes them, which is a property of the
construction, not of its numerics. So are the radii swept by a breakpoint
the output reports (the r_s of icmod 6 and 7), whatever the parameter. The 8
points of the file's mesh on each side of each window are sampled too: an
interpolant across the radius would reach them.

usage: smooth_scan.py EXE BASE.toml KEY [--center X] [--cells N | --span S]
                      [--npts N] [--degree D] [--groups G,G]
                      [--output log | psp8 | upf]
                      [--tol TOL | --tol G=TOL,G=TOL] [--max-jump R]
                      [--jobs J] [--workdir DIR]
       smooth_scan.py --self-test
"""

import argparse
import concurrent.futures
import math
import os
import pathlib
import re
import subprocess
import sys
import tempfile
import tomllib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_output import psp8_part, upf_part  # noqa: E402
from pyoncvpsp.io import OncvpspTextParser  # noqa: E402

RADIUS_KEYS = {"rc", "rcloc", "rxpsh"}
GROUPS = ("kernel", "vloc", "core", "tests", "ld", "ecut", "ghost")
# points of the UPF file's log mesh: every few of the generator's out to the
# linear mesh's extent
UPF_MESH = 600
# points of the file's mesh sampled on each side of the window of a radius
NEAR = 8
# breakpoints of the model core that move with the inputs, as the output
# reports them: the radius r_s of icmod 6 and 7, a C2 joint (C4 with
# match_order = 4)
BREAKPOINTS = (re.compile(r"^\s*model radius r_s\s+([-+.0-9eE]+)"),)
# a core radius rc moves the start of the projectors' taper, TAPER rc
# (run_vkb), as well as rc
TAPER = 0.95
# --output: the file the features come from, as (psfile, upf_mesh)
OUTPUTS = {"log": ("upf", UPF_MESH), "psp8": ("psp8", 0), "upf": ("upf", 0)}


# TOML input: read with tomllib, one value changed, written back


def parse_key(key):
    """('table', 'name', index or None) of KEY, e.g. pseudopotentials.rc[0]."""
    match = re.fullmatch(r"(\w+)\.(\w+)(?:\[(\d+)\])?", key)
    if not match:
        raise ValueError(f"KEY must be table.name or table.name[i], not {key!r}")
    table, name, index = match.groups()
    return table, name, None if index is None else int(index)


def get_value(data, key):
    table, name, index = parse_key(key)
    value = data[table][name]
    return value if index is None else value[index]


def set_value(data, key, value):
    table, name, index = parse_key(key)
    if index is None:
        data[table][name] = value
    else:
        data[table][name][index] = value


def toml_value(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(value, list):
        return "[ " + ", ".join(toml_value(v) for v in value) + ",]"
    raise TypeError(f"no TOML form for {value!r}")


def dump_toml(data):
    """TOML text for the input tables (tables of scalars and arrays, and arrays
    of such tables), with floats written so that they read back exactly."""
    out = []
    for name, value in data.items():
        if isinstance(value, list) and value and isinstance(value[0], dict):
            for entry in value:
                out.append(f"[[{name}]]")
                out += [f"{k} = {toml_value(v)}" for k, v in entry.items()]
                out.append("")
    for name, value in data.items():
        if isinstance(value, dict):
            out.append(f"[{name}]")
            out += [f"{k} = {toml_value(v)}" for k, v in value.items()]
            out.append("")
    return "\n".join(out)


# Features of one output


def file_part(parser, output):
    """(mesh, {name: Block}, local potential, model core or None) of the file
    of --output in the parsed output; the psp8 model core with its four
    derivatives, shape (5, mesh)."""
    if output == "psp8":
        psp8, blocks = psp8_part(parser.psp8_string)
        named = {
            (f"l{ll}" if kind == "scalar" else f"so_l{ll}"): block
            for (_, kind, ll), block in blocks.items()
        }
        core = psp8["model_core"] if psp8["fchrg"] > 0.0 else None
        return psp8["r"], named, psp8["local_potential"], core
    upf, blocks = upf_part(parser.upf_string)
    named = {
        (f"l{ll}" if jj is None else f"l{ll}j{jj:g}"): block
        for (_, ll, jj), block in blocks.items()
    }
    return upf["mesh"]["r"], named, upf["local"], upf.get("nlcc")


def features(path, samples, output):
    """{name: value} of the output at path, the names of the form group:label."""
    parser = OncvpspTextParser(path=path, check=False)
    _, blocks, local, core = file_part(parser, output)
    out = {}
    for name, block in blocks.items():
        beta = block.beta[:, samples]
        kern = beta.T @ block.dij @ beta
        for a, ka in enumerate(samples):
            for b in range(a, len(samples)):
                out[f"kernel:{name}:r{ka},r{samples[b]}"] = float(kern[a, b])
    for k in samples:
        out[f"vloc:r{k}"] = float(local[k])
    if core is not None:
        rows = core if hasattr(core[0], "__len__") else [core]
        for d, row in enumerate(rows):
            block = "core" if d == 0 else f"core:d{d}"
            for k in samples:
                out[f"{block}:r{k}"] = float(row[k])
    eigs = (
        state["eig_ps"]
        for config in parser.output["test_configurations"]
        for state in config
        if state["eig_ps"] is not None
    )
    for i, eps in enumerate(eigs):
        out[f"tests:eps{i}"] = eps
    plots = parser.output["plot_data"]
    diagnostics = parser.output.get("diagnostics") or {}
    curves = diagnostics.get("log_derivatives", plots.get("log_derivatives"))
    tables = diagnostics.get("ecut_crossings", plots.get("convergence_profiles"))
    for curve in curves or []:
        tag = f"l{curve['l']}" + ("" if curve["j"] is None else f"j{curve['j']}")
        for i, phase in enumerate(curve["log_deriv_ps"]):
            out[f"ld:{tag}:cos:e{i}"] = math.cos(2.0 * phase)
            out[f"ld:{tag}:sin:e{i}"] = math.sin(2.0 * phase)
    for table in tables or []:
        tag = f"l{table['l']}p{table['iproj']}"
        if table["j"] is not None:
            tag += f"j{table['j']}"
        for err, cut in zip(table["eresid"], table["ecut"]):
            out[f"ecut:{tag}:{err:.0e}"] = float(cut)
    for test in parser.output.get("bound_ghost_tests") or []:
        tag = f"n{test['n']}l{test['l']}"
        out[f"ghost:{tag}:eig"] = test["eig"]
        out[f"ghost:{tag}:basis"] = test["eig_basis"]
    return out


def file_mesh(path, output):
    """The radii of the file of one output."""
    parser = OncvpspTextParser(path=path, check=False)
    return [float(r) for r in file_part(parser, output)[0]]


# Least-squares polynomial fit (Legendre basis, modified Gram-Schmidt)


def legendre(degree, t):
    values = [1.0, t]
    for n in range(1, degree):
        values.append(((2 * n + 1) * t * values[n] - n * values[n - 1]) / (n + 1))
    return values[: degree + 1]


def fit_residuals(ts, ys, degree):
    """y - p(t) at each t, p the least-squares polynomial of the degree."""
    columns = [list(c) for c in zip(*(legendre(degree, t) for t in ts))]
    resid = list(ys)
    basis = []
    for col in columns:
        for q in basis:
            dot = sum(a * b for a, b in zip(q, col))
            col = [a - dot * b for a, b in zip(col, q)]
        norm = math.sqrt(sum(a * a for a in col))
        q = [a / norm for a in col]
        basis.append(q)
        dot = sum(a * b for a, b in zip(q, resid))
        resid = [a - dot * b for a, b in zip(resid, q)]
    return resid


def group_scales(center):
    """Scale of each feature: the largest |value| at the window's centre over
    the features of its block (kernel block, vloc, core derivative order)."""
    blocks = {}
    for name, value in center.items():
        block = name.rsplit(":", 1)[0]
        blocks[block] = max(blocks.get(block, 0.0), abs(value))
    return {name: blocks[name.rsplit(":", 1)[0]] or 1.0 for name in center}


def analyse(ts, runs, degree):
    """{group: (worst relative residual, feature, jump ratio, largest step
    relative to the scale)}."""
    names = sorted(set.intersection(*(set(r) for r in runs)))
    center = runs[len(runs) // 2]
    scales = group_scales({n: center[n] for n in names})
    worst = {}
    for name in names:
        ys = [r[name] for r in runs]
        resid = fit_residuals(ts, ys, degree)
        rel = max(abs(x) for x in resid) / scales[name]
        steps = sorted(abs(b - a) for a, b in zip(ys, ys[1:]))
        median = steps[len(steps) // 2]
        jump = (
            steps[-1] / median if median > 0.0 else math.inf if steps[-1] > 0 else 1.0
        )
        group = name.split(":", 1)[0]
        if group not in worst or rel > worst[group][0]:
            worst[group] = (rel, name, jump, steps[-1] / scales[name])
    return worst


def fails(rel, jump, step, tol, max_jump):
    """A residual over tol, or a jump ratio over max_jump in steps over tol."""
    if tol is not None and rel > tol:
        return True
    return max_jump is not None and jump > max_jump and (tol is None or step > tol)


# Running the program


def run_one(exe, data, key, value, workdir, index, output):
    data = {k: (dict(v) if isinstance(v, dict) else v) for k, v in data.items()}
    for k, v in data.items():
        if isinstance(v, dict):
            data[k] = {
                kk: (list(vv) if isinstance(vv, list) else vv) for kk, vv in v.items()
            }
    set_value(data, key, value)
    # the file of --output, whatever the input asks for
    data["pp_output"]["psfile"], data["pp_output"]["upf_mesh"] = OUTPUTS[output]
    path = pathlib.Path(workdir) / f"scan{index:03d}"
    path.mkdir(parents=True, exist_ok=True)
    (path / "in.toml").write_text(dump_toml(data))
    result = subprocess.run(
        [exe, "-t", "in.toml"], cwd=path, capture_output=True, text=True
    )
    (path / "out").write_text(result.stdout)
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-3:])
        raise RuntimeError(
            f"{os.path.basename(exe)} failed for {key}={value!r} in {path}:\n{tail}"
        )
    return path / "out"


def breakpoints(outputs):
    """The breakpoint radii the outputs report (BREAKPOINTS), over all runs:
    the functions are only as smooth as the joint there, a property of the
    construction, so the radii it sweeps are dropped as for a radius key."""
    radii = []
    for out in outputs:
        for line in pathlib.Path(out).read_text().splitlines():
            for pattern in BREAKPOINTS:
                m = pattern.match(line)
                if m:
                    radii.append(float(m.group(1)))
    return radii


def scan(args):
    data = tomllib.loads(pathlib.Path(args.toml).read_text())
    center = float(get_value(data, args.key)) if args.center is None else args.center
    if args.span is not None:
        span = args.span
    else:
        amesh = data.get("log_mesh", {}).get("amesh", 1.006)
        span = center * (amesh**args.cells - 1.0)
    values = [center + span * (k / (args.npts - 1) - 0.5) for k in range(args.npts)]
    ts = [2.0 * k / (args.npts - 1) - 1.0 for k in range(args.npts)]

    workdir = args.workdir or tempfile.mkdtemp(prefix="smooth_scan_")
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        outputs = list(
            pool.map(
                lambda kv: run_one(
                    args.exe, data, args.key, kv[1], workdir, kv[0], args.output
                ),
                enumerate(values),
            )
        )

    mesh = file_mesh(outputs[0], args.output)
    samples = sorted(
        {
            min(range(len(mesh)), key=lambda k: abs(mesh[k] - r))
            for r in (
                args.sample * j for j in range(1, int(mesh[-1] / args.sample) + 1)
            )
        }
    )
    if args.output == "log":
        step = max(b / a for a, b in zip(mesh[1:], mesh[2:]))
    else:
        step = data.get("log_mesh", {}).get("amesh", 1.006)
    windows = []
    if parse_key(args.key)[1] in RADIUS_KEYS:
        lo, hi = min(values) / step**2, max(values) * step**2
        windows += [(lo, hi), (TAPER * lo, TAPER * hi)]
    moved = breakpoints(outputs)
    if moved:
        windows += [(min(moved) / step**2, max(moved) * step**2)]
    if windows:
        near = set()
        for a, b in windows:
            # the NEAR points of the mesh on each side of the window, which an
            # interpolant across the radius would reach
            below = [k for k in range(len(mesh)) if mesh[k] < a][-NEAR:]
            above = [k for k in range(len(mesh)) if mesh[k] > b][:NEAR]
            near |= set(below) | set(above)
        samples = [
            k
            for k in sorted(set(samples) | near)
            if not any(a <= mesh[k] <= b for a, b in windows)
        ]

    groups = set(args.groups.split(","))
    runs = [
        {
            n: v
            for n, v in features(out, samples, args.output).items()
            if n.split(":", 1)[0] in groups
        }
        for out in outputs
    ]
    worst = analyse(ts, runs, args.degree)

    print(
        f"{args.key}: {args.npts} values from {values[0]:.10g} to {values[-1]:.10g}"
        f" (span {span:.4g}), degree {args.degree} fit, {args.output} file,"
        f" outputs in {workdir}"
    )
    print(f"{'group':8s} {'residual':>10s} {'jump':>10s}  worst feature")
    tols = parse_tol(args.tol)
    status = 0
    for group in GROUPS:
        if group not in worst:
            continue
        rel, name, jump, step = worst[group]
        flag = ""
        tol = tols.get(group, tols.get(None))
        if fails(rel, jump, step, tol, args.max_jump):
            flag = "  FAIL"
            status = 1
        print(f"{group:8s} {rel:10.2e} {jump:10.2e}  {name}{flag}")
    missing = groups - set(worst)
    if missing:
        print(f"no features of group(s) {', '.join(sorted(missing))}  FAIL")
        status = 1
    return status


def parse_tol(text):
    """{group: tol} from 'TOL' (key None: every group) or 'G=TOL,G=TOL'."""
    if text is None:
        return {}
    if "=" not in text:
        return {None: float(text)}
    pairs = (item.split("=") for item in text.split(","))
    return {group: float(value) for group, value in pairs}


# Checks of the analysis on functions whose smoothness is known


def self_test():
    """The fit must pass a smooth function to round-off and catch a jump and a
    kink well below the scale of the function."""
    npts, degree = 41, 6
    ts = [2.0 * k / (npts - 1) - 1.0 for k in range(npts)]
    span = 0.03
    smooth = [{"kernel:a:x": math.exp(math.sin(1.0 + span * t))} for t in ts]
    jump = [
        {"kernel:a:x": math.exp(math.sin(1.0 + span * t)) + (1e-7 if t > 0.13 else 0.0)}
        for t in ts
    ]
    kink = [
        {"kernel:a:x": math.exp(math.sin(1.0 + span * t)) + 1e-6 * max(0.0, t - 0.13)}
        for t in ts
    ]
    failures = []
    for degree in (degree, 10):
        if analyse(ts, smooth, degree)["kernel"][0] > 1e-12:
            failures.append(f"smooth function not fitted to 1e-12 at degree {degree}")
        if analyse(ts, jump, degree)["kernel"][0] < 1e-8:
            failures.append(f"jump of 1e-7 not seen at 1e-8 at degree {degree}")
        if analyse(ts, kink, degree)["kernel"][0] < 1e-9:
            failures.append(f"kink of slope 1e-6 not seen at 1e-9 at degree {degree}")
    # a constant with noise of ~1e-14: steps under a tol of 1e-12 are no jump,
    # however their ratio; a staircase (a parameter rounded) is one
    noise = [
        {"kernel:a:x": 1.0 + 1e-15 * (-1) ** k + (6e-15 if k == 20 else 0.0)}
        for k in range(npts)
    ]
    rel, _, ratio, step = analyse(ts, noise, degree)["kernel"]
    if ratio < 3.0 or fails(rel, ratio, step, 1e-12, 3.0):
        failures.append("noise of ~1e-14 failed at tol 1e-12")
    stairs = [
        {"kernel:a:x": math.exp(math.sin(1.0 + span * round(4.0 * t) / 4.0))}
        for t in ts
    ]
    _, _, ratio, step = analyse(ts, stairs, degree)["kernel"]
    if not fails(0.0, ratio, step, 1e-9, 3.0):
        failures.append("staircase passed the jump ratio at tol 1e-9")
    data = {
        "t": [{"n": [1, 2]}],
        "a": {"x": 0.1 + 0.2, "s": 'p"q', "v": [2.7, 1e-5, 3]},
    }
    if tomllib.loads(dump_toml(data)) != data:
        failures.append("TOML written does not read back")
    for f in failures:
        print("FAILED:", f, file=sys.stderr)
    if not failures:
        print("smooth_scan self-test: all checks passed")
    return 1 if failures else 0


def main():
    if sys.argv[1:] == ["--self-test"]:
        return self_test()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("exe")
    parser.add_argument("toml")
    parser.add_argument("key")
    parser.add_argument(
        "--center", type=float, help="centre of the window (default: the input's value)"
    )
    parser.add_argument(
        "--cells",
        type=float,
        default=2.0,
        help="width of the window in log-mesh cells at the centre (default 2)",
    )
    parser.add_argument(
        "--span", type=float, help="width of the window, instead of --cells"
    )
    parser.add_argument("--npts", type=int, default=41)
    parser.add_argument("--degree", type=int, default=6)
    parser.add_argument(
        "--sample",
        type=float,
        default=0.25,
        help="spacing of the sample radii (bohr, default 0.25)",
    )
    parser.add_argument(
        "--groups",
        default="kernel,vloc,core",
        help=f"comma-separated feature groups of {','.join(GROUPS)}",
    )
    parser.add_argument(
        "--tol", help="largest residual: one value, or group=value pairs"
    )
    parser.add_argument(
        "--max-jump",
        type=float,
        help="largest ratio of the largest step to the median step",
    )
    parser.add_argument(
        "--output",
        choices=sorted(OUTPUTS),
        default="log",
        help="file the features come from: the UPF file on the log mesh (default),"
        " the psp8 file, or the UPF file on the linear mesh",
    )
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--workdir")
    args = parser.parse_args()
    args.exe = os.path.realpath(args.exe)
    return scan(args)


if __name__ == "__main__":
    sys.exit(main())
