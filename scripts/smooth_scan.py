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
digits. With --output psp8, upf or psml they come instead from the psp8 file
(14 digits), the UPF file on the linear output mesh, onto which the writers
evaluate the functions (doc/smoothness_plan.md, O-1), or the PSML file (on the
generator's log mesh, 13 digits). A comma-separated list (e.g. log,psp8,psml)
takes them from several files of each run (psfile = "both" writes all three;
log and upf exclude each other, as one run writes one UPF file), and the
groups below drawn from a file are reported for each, e.g. kernel.psp8. They
are what a
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
          (the block's levels to 1e-7 are left out: their last crossings
          jump where a shoulder of the curve passes them)
          (both from the output's DIAGNOSTICS block, 17 digits; without one,
          from the tables printed to 6 decimals and 0.01 Ha)
  ghost   the energies of the bound-ghost test, the radial equation's and
          the basis's, and <r>/rc of the basis state the positive-energy
          test comes closest to flagging (from the DIAGNOSTICS block; without
          one, the energies of the table printed to 6 decimals)

The last three are what the optimizer reads (doc/smoothness_plan.md, order
of the remaining work); from an output without the DIAGNOSTICS block only a
--tol at their print resolution is meaningful.

The sample radii are the points of the file's mesh nearest to every --sample
bohr. For a parameter that is a radius (a key named rc, rcloc or rxpsh), the
radii inside the scanned window (widened by two steps of the file's log mesh,
or of the generator's for a linear one), and inside 0.95 times it (where the
projectors' taper starts), are dropped: there the functions are only as
smooth as the matching at the radius makes them, which is a property of the
construction, not of its numerics. So are the radii swept by a breakpoint
the output reports (r* of icmod 1, r_s of icmod 6 and 7), whatever the
parameter. The 8 points of the file's mesh on each side of each window are
sampled too: an interpolant across the radius would reach them.

KEY[*] moves every element of an array input by the same offset (e.g.
pseudopotentials.rc[*], each core radius at once): --center is then the
offset's centre (default 0), --cells counts cells at the first element, and
each element's radius window is dropped as above. Each channel's features
stay their own, so one scan does the work of one per channel for the
outputs that depend on their own channel's inputs only.

With --set KEY=VALUE (repeated as needed; VALUE in TOML, e.g.
pseudopotentials.well_tail="free") the input is changed before the scan.

The groups kernel, vloc and core come from each file of --output: a group
named alone is taken from all of them, one named with its file (core.log)
from that file only. --tol and --degree take a value for a group of one file
the same way (--tol core.psp8=1e-9, --degree 18,psp8=24: 24 for the groups of
the psp8 file, 18 for the rest).

usage: smooth_scan.py EXE BASE.toml KEY [--center X] [--cells N | --span S]
                      [--set KEY=VALUE ...]
                      [--npts N] [--degree D[,OUT=D]] [--groups G,G]
                      [--output OUT[,OUT]] (OUT: log, psp8, upf, psml)
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
from compare_output import psml_part, psp8_part, upf_part  # noqa: E402
from pyoncvpsp.io import OncvpspTextParser  # noqa: E402

RADIUS_KEYS = {"rc", "rcloc", "rxpsh"}
GROUPS = ("kernel", "vloc", "core", "tests", "ld", "ecut", "ghost")
# the groups drawn from the pseudopotential file, one of each per file
FILE_GROUPS = ("kernel", "vloc", "core")
# points of the UPF file's log mesh: every few of the generator's out to the
# linear mesh's extent
UPF_MESH = 600
# points of the file's mesh sampled on each side of the window of a radius
NEAR = 8
# breakpoints of the model core that move with the inputs, as the output
# reports them: r* of icmod 1, a C4 joint, and r_s of icmod 6 and 7, a C2
# joint (C4 with match_order = 4)
BREAKPOINTS = (re.compile(r"^\s*model radius r(?:\*|_s)\s+([-+.0-9eE]+)"),)
# a core radius rc moves the start of the projectors' taper, TAPER rc
# (run_vkb), as well as rc
TAPER = 0.95
# --output: the file the features come from, as (psfile, upf_mesh)
OUTPUTS = {
    "log": ("upf", UPF_MESH),
    "psp8": ("psp8", 0),
    "upf": ("upf", 0),
    "psml": ("psml", 0),
}


def parse_outputs(text):
    """The files of --output, a comma-separated list of OUTPUTS."""
    outputs = list(dict.fromkeys(text.split(",")))
    unknown = [o for o in outputs if o not in OUTPUTS]
    if unknown:
        raise ValueError(f"--output: no file {', '.join(unknown)}")
    if "log" in outputs and "upf" in outputs:
        raise ValueError("--output: log and upf are both the UPF file of a run")
    return outputs


def psfile_of(outputs):
    """(psfile, upf_mesh) of the input that writes the files of outputs."""
    if len(outputs) == 1:
        return OUTPUTS[outputs[0]]
    return "both", UPF_MESH if "log" in outputs else 0


# TOML input: read with tomllib, one value changed, written back


def parse_key(key):
    """('table', 'name', index, '*' or None) of KEY, e.g. pseudopotentials.rc[0]
    or pseudopotentials.rc[*] (every element)."""
    match = re.fullmatch(r"(\w+)\.(\w+)(?:\[(\d+|\*)\])?", key)
    if not match:
        raise ValueError(
            f"KEY must be table.name, table.name[i] or table.name[*], not {key!r}"
        )
    table, name, index = match.groups()
    return table, name, index if index in (None, "*") else int(index)


def get_value(data, key):
    table, name, index = parse_key(key)
    value = data[table][name]
    return value if index in (None, "*") else value[index]


def set_value(data, key, value):
    table, name, index = parse_key(key)
    if index in (None, "*"):
        data[table][name] = value
    else:
        data[table][name][index] = value


def apply_sets(data, sets):
    """data with each KEY=VALUE of sets applied, VALUE read as TOML."""
    for item in sets or []:
        key, _, text = item.partition("=")
        set_value(data, key.strip(), tomllib.loads(f"v = {text}")["v"])
    return data


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
    if output == "psml":
        psml, blocks = psml_part(parser.psml_string)
        named = {
            f"{kind}_l{ll}" + ("" if jj is None else f"j{jj:g}"): block
            for (_, kind, ll, jj), block in blocks.items()
        }
        core = psml["core_charge"]
        return (
            psml["grid"],
            named,
            psml["local_potential"]["values"],
            None if core is None else core["values"],
        )
    upf, blocks = upf_part(parser.upf_string)
    named = {
        (f"l{ll}" if jj is None else f"l{ll}j{jj:g}"): block
        for (_, ll, jj), block in blocks.items()
    }
    return upf["mesh"]["r"], named, upf["local"], upf.get("nlcc")


def features(path, samples, groups=GROUPS):
    """{name: value} of the output at path, the names of the form group:label,
    those of the files group.file:label, samples {file: indices on its mesh}.
    Only the parts of the output that the groups read are parsed."""
    parser = OncvpspTextParser(path=path, check=False)
    groups = {g.split(".")[0] for g in groups}
    out = {}
    for output, points in samples.items() if groups & set(FILE_GROUPS) else ():
        _, blocks, local, core = file_part(parser, output)
        for name, block in blocks.items():
            beta = block.beta[:, points]
            kern = beta.T @ block.dij @ beta
            for a, ka in enumerate(points):
                for b in range(a, len(points)):
                    out[f"kernel.{output}:{name}:r{ka},r{points[b]}"] = float(
                        kern[a, b]
                    )
        for k in points:
            out[f"vloc.{output}:r{k}"] = float(local[k])
        if core is not None:
            rows = core if hasattr(core[0], "__len__") else [core]
            for d, row in enumerate(rows):
                block = f"core.{output}" if d == 0 else f"core.{output}:d{d}"
                for k in points:
                    out[f"{block}:r{k}"] = float(row[k])
    if "tests" in groups:
        eigs = (
            state["eig_ps"]
            for config in parser.test_configurations
            for state in config
            if state["eig_ps"] is not None
        )
        for i, eps in enumerate(eigs):
            out[f"tests:eps{i}"] = eps
    if not groups & {"ld", "ecut", "ghost"}:
        return out
    diagnostics = parser.diagnostics
    if "ld" in groups:
        curves = (
            diagnostics["log_derivatives"]
            if "log_derivatives" in diagnostics
            else parser.log_derivatives
        )
        for curve in curves or []:
            tag = f"l{curve['l']}" + ("" if curve["j"] is None else f"j{curve['j']}")
            for i, phase in enumerate(curve["log_deriv_ps"]):
                out[f"ld:{tag}:cos:e{i}"] = math.cos(2.0 * phase)
                out[f"ld:{tag}:sin:e{i}"] = math.sin(2.0 * phase)
    if "ecut" in groups:
        tables = (
            diagnostics["ecut_crossings"]
            if "ecut_crossings" in diagnostics
            else parser.convergence_profiles
        )
        for table in tables or []:
            tag = f"l{table['l']}p{table['iproj']}"
            if table["j"] is not None:
                tag += f"j{table['j']}"
            for err, cut in zip(table["eresid"], table["ecut"]):
                # the table's levels: below 1e-5 the last crossing jumps where a
                # shoulder of the curve passes the level (G5)
                if err > 0.99e-5:
                    out[f"ecut:{tag}:{err:.0e}"] = float(cut)
        weighted = diagnostics.get("ecut_weighted")
        if weighted:
            out[f"ecut:weighted:{weighted['error']:.0e}"] = float(weighted["ecut"])
    if "ghost" in groups:
        margins = diagnostics.get("ghost_margins") or parser.bound_ghost_tests
        for test in margins or []:
            tag = f"n{test['n']}l{test['l']}"
            if test.get("j") is not None:
                tag += f"j{test['j']}"
            elif test.get("set") is not None:
                tag += f"s{test['set']}"
            out[f"ghost:{tag}:eig"] = test["eig"]
            out[f"ghost:{tag}:basis"] = test["eig_basis"]
        for test in diagnostics.get("ghost_localization") or []:
            tag = f"l{test['l']}" + ("" if test["j"] is None else f"j{test['j']}")
            out[f"ghost:{tag}:r_over_rc"] = test["r_over_rc"]
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
    relative to the scale)}; degree a number, or {group: degree}."""
    names = sorted(set.intersection(*(set(r) for r in runs)))
    center = runs[len(runs) // 2]
    scales = group_scales({n: center[n] for n in names})
    worst = {}
    for name in names:
        ys = [r[name] for r in runs]
        group = name.split(":", 1)[0]
        resid = fit_residuals(
            ts, ys, degree if isinstance(degree, int) else degree[group]
        )
        rel = max(abs(x) for x in resid) / scales[name]
        steps = sorted(abs(b - a) for a, b in zip(ys, ys[1:]))
        median = steps[len(steps) // 2]
        jump = (
            steps[-1] / median if median > 0.0 else math.inf if steps[-1] > 0 else 1.0
        )
        if group not in worst or rel > worst[group][0]:
            worst[group] = (rel, name, jump, steps[-1] / scales[name])
    return worst


def fails(rel, jump, step, tol, max_jump):
    """A residual over tol, or a jump ratio over max_jump in steps over tol."""
    if tol is not None and rel > tol:
        return True
    return max_jump is not None and jump > max_jump and (tol is None or step > tol)


# Running the program


def run_one(exe, data, key, value, workdir, index, outputs):
    data = {k: (dict(v) if isinstance(v, dict) else v) for k, v in data.items()}
    for k, v in data.items():
        if isinstance(v, dict):
            data[k] = {
                kk: (list(vv) if isinstance(vv, list) else vv) for kk, vv in v.items()
            }
    set_value(data, key, value)
    # the files of --output, whatever the input asks for
    data["pp_output"]["psfile"], data["pp_output"]["upf_mesh"] = psfile_of(outputs)
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


def sample_points(mesh, output, data, key, values, moved, spacing):
    """The indices of the sample radii on the mesh of one file: every spacing
    bohr, less the windows swept by a radius key and by the breakpoints moved,
    plus the NEAR points on each side of each window."""
    samples = sorted(
        {
            min(range(len(mesh)), key=lambda k: abs(mesh[k] - r))
            for r in (spacing * j for j in range(1, int(mesh[-1] / spacing) + 1))
        }
    )
    if output in ("log", "psml"):
        step = max(b / a for a, b in zip(mesh[1:], mesh[2:]))
    else:
        step = data.get("log_mesh", {}).get("amesh", 1.006)
    windows = []
    if parse_key(key)[1] in RADIUS_KEYS:
        # each element's values for KEY[*]
        for column in zip(*values) if parse_key(key)[2] == "*" else [values]:
            lo, hi = min(column) / step**2, max(column) * step**2
            windows += [(lo, hi), (TAPER * lo, TAPER * hi)]
    if moved:
        windows += [(min(moved) / step**2, max(moved) * step**2)]
    if not windows:
        return samples
    near = set()
    for a, b in windows:
        # the NEAR points of the mesh on each side of the window, which an
        # interpolant across the radius would reach
        below = [k for k in range(len(mesh)) if mesh[k] < a][-NEAR:]
        above = [k for k in range(len(mesh)) if mesh[k] > b][:NEAR]
        near |= set(below) | set(above)
    return [
        k
        for k in sorted(set(samples) | near)
        if not any(a <= mesh[k] <= b for a, b in windows)
    ]


def scan(args):
    outputs = parse_outputs(args.output)
    data = apply_sets(tomllib.loads(pathlib.Path(args.toml).read_text()), args.set)
    every = parse_key(args.key)[2] == "*"
    if every:
        # an offset to every element
        base = [float(v) for v in get_value(data, args.key)]
        center = 0.0 if args.center is None else args.center
        start = base[0] + center
    else:
        center = (
            float(get_value(data, args.key)) if args.center is None else args.center
        )
        start = center
    if args.span is not None:
        span = args.span
    else:
        amesh = data.get("log_mesh", {}).get("amesh", 1.006)
        span = start * (amesh**args.cells - 1.0)
    offsets = [center + span * (k / (args.npts - 1) - 0.5) for k in range(args.npts)]
    values = [[b + o for b in base] for o in offsets] if every else offsets
    ts = [2.0 * k / (args.npts - 1) - 1.0 for k in range(args.npts)]

    workdir = args.workdir or tempfile.mkdtemp(prefix="smooth_scan_")
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        runs_out = list(
            pool.map(
                lambda kv: run_one(
                    args.exe, data, args.key, kv[1], workdir, kv[0], outputs
                ),
                enumerate(values),
            )
        )

    moved = breakpoints(runs_out)
    samples = {
        output: sample_points(
            file_mesh(runs_out[0], output),
            output,
            data,
            args.key,
            values,
            moved,
            args.sample,
        )
        for output in outputs
    }

    groups = set(args.groups.split(","))
    # parsing an output is Python-bound (~0.1 s, of the groups read): in processes,
    # as many as the runs
    with concurrent.futures.ProcessPoolExecutor(args.jobs) as pool:
        parsed = list(
            pool.map(
                features,
                runs_out,
                [samples] * len(runs_out),
                [groups] * len(runs_out),
            )
        )

    def wanted(name):
        group = name.split(":", 1)[0]
        return group in groups or group.split(".")[0] in groups

    runs = [{n: v for n, v in run.items() if wanted(n)} for run in parsed]
    degrees = parse_tol(args.degree, int)
    degree = {
        g: degrees.get(g.partition(".")[2], degrees.get(None, 6))
        for g in {n.split(":", 1)[0] for run in runs for n in run}
    }
    worst = analyse(ts, runs, degree)

    print(
        f"{args.key}: {args.npts} {'offsets' if every else 'values'} from"
        f" {offsets[0]:.10g} to {offsets[-1]:.10g}"
        f" (span {span:.4g}), degree {args.degree} fit, {args.output} file,"
        f" outputs in {workdir}"
    )
    print(f"{'group':11s} {'residual':>10s} {'jump':>10s}  worst feature")
    tols = parse_tol(args.tol)
    status = 0
    order = {g: i for i, g in enumerate(GROUPS)}
    for group in sorted(worst, key=lambda g: (order[g.split(".")[0]], g)):
        rel, name, jump, step = worst[group]
        flag = ""
        tol = tols.get(group, tols.get(group.split(".")[0], tols.get(None)))
        if fails(rel, jump, step, tol, args.max_jump):
            flag = "  FAIL"
            status = 1
        print(f"{group:11s} {rel:10.2e} {jump:10.2e}  {name}{flag}")
    missing = {
        g for g in groups if not any(w == g or w.split(".")[0] == g for w in worst)
    }
    if missing:
        print(f"no features of group(s) {', '.join(sorted(missing))}  FAIL")
        status = 1
    return status


def parse_tol(text, kind=float):
    """{group: value} of --tol or --degree, 'V,G=V,G=V': V for every group
    (key None), G=V for group G (kernel, or kernel.psp8 for that file's
    group; for --degree a file, psp8, for its groups)."""
    out = {}
    for item in (text or "").split(",") if text else []:
        name, _, value = item.rpartition("=")
        out[name or None] = kind(value)
    return out


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
        "--center",
        type=float,
        help="centre of the window (default: the input's value; for KEY[*], of"
        " the offset, default 0)",
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
    parser.add_argument(
        "--set",
        action="append",
        metavar="KEY=VALUE",
        help="change the input before the scan (VALUE in TOML; repeatable)",
    )
    parser.add_argument("--npts", type=int, default=41)
    parser.add_argument(
        "--degree",
        default="6",
        help="degree of the fit: D, or D,OUT=D for the groups of file OUT",
    )
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
        default="log",
        help="files the features come from, comma-separated: log (the UPF file on"
        " the log mesh, the default), psp8, upf (the UPF file on the linear mesh),"
        " psml",
    )
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--workdir")
    args = parser.parse_args()
    args.exe = os.path.realpath(args.exe)
    return scan(args)


if __name__ == "__main__":
    sys.exit(main())
