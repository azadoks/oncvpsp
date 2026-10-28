#!/usr/bin/env python3
"""Coverage of src/ by each test, and the fewest inputs and tests that keep it.

usage:
  coverage.py run BUILD OUTDIR [-j N] [-L REGEX] [-LE REGEX] [-R REGEX]
  coverage.py report OUTDIR [--compare OLDDIR] [--uncovered] [--cover]

run: runs the tests of BUILD (the coverage preset, build-cov/) that ctest selects with -L, -LE
and -R, the fixtures first, each with its own GCOV_PREFIX, so that the .gcda files of each test
land in OUTDIR/gcda/<test>; then reads them with gcov into OUTDIR/cov.json:
{test: {"lines": [[file, line], ...], "branches": [[file, line, i], ...]}, "_all": {...}},
"_all" every instrumented line and branch of src/. The logs are in OUTDIR/log, the times
(-O0 wall seconds) in OUTDIR/times.json.

report: the lines and branches covered, in total and per file; with --compare, those OLDDIR
covered and OUTDIR does not (lost) and the reverse; with --uncovered, the lines no test runs, as
ranges per file; with --cover, the fewest inputs, then the cheapest tests on them, that keep the
coverage of each kind of test (the tests of one kind assert the same thing: ref, toml, echo,
split, ...), by greedy choice. The scan label is its own measure and the nist tests check
accuracy, not coverage: run them with -LE "scan|nist".
"""

import argparse
import collections
import gzip
import json
import os
import re
import subprocess
import tempfile
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

KINDS = (
    "ref",
    "toml_input",
    "toml",
    "echo",
    "log_mesh",
    "upf_mesh",
    "split",
    "ld",
    "confine",
    "ae_file",
    "ae_scf",
    "ae_4f",
    "occupation_ramp",
    "error",
    "pipe",
    "psml",
    "psp8_core",
    "well_onset",
    "well_tail",
    "well_rmax",
    "well_norm",
    "unit",
    "lint",
    "projector_output",
    "run_ae",
    "run",
    "scan",
    "nist",
)


def props(test):
    return {p["name"]: p["value"] for p in test.get("properties", [])}


def run_tests(args):
    sel = ["ctest", "--test-dir", args.build, "--show-only=json-v1"]
    for flag in ("L", "LE", "R"):
        if getattr(args, flag):
            sel += [f"-{flag}", getattr(args, flag)]
    tests = json.loads(
        subprocess.run(sel, capture_output=True, check=True, text=True).stdout
    )
    tests = tests["tests"]
    os.makedirs(f"{args.outdir}/log", exist_ok=True)
    cond = threading.Condition()
    free = [args.jobs]
    times = {}

    def run(t):
        p = props(t)
        name = t["name"]
        nproc = min(int(p.get("PROCESSORS", 1)), args.jobs)
        with cond:
            cond.wait_for(lambda: free[0] >= nproc)
            free[0] -= nproc
        try:
            env = dict(os.environ)
            for kv in p.get("ENVIRONMENT", []):
                k, v = kv.split("=", 1)
                env[k] = v
            env["GCOV_PREFIX"] = os.path.abspath(f"{args.outdir}/gcda/{name}")
            env["GCOV_PREFIX_STRIP"] = "0"
            t0 = time.time()
            with open(f"{args.outdir}/log/{name}", "w") as log:
                r = subprocess.run(
                    t["command"],
                    cwd=p.get("WORKING_DIRECTORY"),
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
            dt = time.time() - t0
            times[name] = (r.returncode, dt)
            print(f"{name} rc={r.returncode} {dt:.1f}s", flush=True)
        finally:
            with cond:
                free[0] += nproc
                cond.notify_all()

    # the fixtures first, as their tests read what they write; in each phase the long ones first
    setup = [t for t in tests if "FIXTURES_SETUP" in props(t)]
    rest = [t for t in tests if "FIXTURES_SETUP" not in props(t)]
    for phase in (setup, rest):
        phase.sort(key=lambda t: not t["name"].startswith(("scan_", "nist_")))
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            list(ex.map(run, phase))
    with open(f"{args.outdir}/times.json", "w") as f:
        json.dump(times, f, indent=1)
    with open(f"{args.outdir}/tests.json", "w") as f:
        json.dump(tests, f)
    failed = sorted(n for n, (rc, _) in times.items() if rc)
    if failed:
        print("failed (their coverage is kept):", " ".join(failed))
    collect(args)


def gcov_json(gcno, gcda_root):
    """gcov's JSON of one object; its .gcda under gcda_root (None: no counts, every line)."""
    objdir = os.path.dirname(gcno)
    stem = os.path.basename(gcno)[:-5]
    with tempfile.TemporaryDirectory() as tmp:
        os.symlink(gcno, os.path.join(tmp, stem + ".gcno"))
        if gcda_root:
            g = os.path.join(gcda_root + objdir, stem + ".gcda")
            if not os.path.exists(g):
                return []
            os.symlink(g, os.path.join(tmp, stem + ".gcda"))
        subprocess.run(
            ["gcov", "-j", "-b", "-o", tmp, os.path.join(tmp, stem + ".gcda")],
            cwd=tmp,
            capture_output=True,
        )
        return [
            json.load(gzip.open(os.path.join(tmp, f)))
            for f in os.listdir(tmp)
            if f.endswith(".gcov.json.gz")
        ]


def summarize(jsons, src, executed):
    lines, branches = set(), set()
    for j in jsons:
        for f in j["files"]:
            path = f["file"]
            if not os.path.isabs(path) or not path.startswith(src):
                continue
            rel = os.path.relpath(path, os.path.dirname(src))
            for ln in f["lines"]:
                if not executed or ln["count"] > 0:
                    lines.add((rel, ln["line_number"]))
                for i, b in enumerate(ln.get("branches", [])):
                    if not executed or b["count"] > 0:
                        branches.add((rel, ln["line_number"], i))
    return lines, branches


def all_items(task):
    gcno, src = task
    return summarize(gcov_json(gcno, None), src, False)


def test_items(task):
    name, root, gcnos, src = task
    L, B = set(), set()
    for g in gcnos:
        lines, branches = summarize(gcov_json(g, root), src, True)
        L |= lines
        B |= branches
    return name, sorted(L), sorted(B)


def collect(args):
    build = os.path.abspath(args.build)
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src")
    src = os.path.realpath(src)
    gcnos = [
        os.path.join(d, f)
        for d, _, fs in os.walk(build)
        if "/_deps" not in d
        for f in fs
        if f.endswith(".gcno")
    ]
    names = sorted(os.listdir(os.path.join(args.outdir, "gcda")))
    out = {}
    with ProcessPoolExecutor(args.jobs) as ex:
        AL, AB = set(), set()
        for lines, branches in ex.map(all_items, [(g, src) for g in gcnos]):
            AL |= lines
            AB |= branches
        out["_all"] = {"lines": sorted(AL), "branches": sorted(AB)}
        roots = [
            (n, os.path.abspath(os.path.join(args.outdir, "gcda", n)), gcnos, src)
            for n in names
        ]
        for name, lines, branches in ex.map(test_items, roots):
            out[name] = {"lines": lines, "branches": branches}
    with open(os.path.join(args.outdir, "cov.json"), "w") as f:
        json.dump(out, f)
    print(len(names), "tests;", len(AL), "lines;", len(AB), "branches")


def load(outdir):
    with open(os.path.join(outdir, "cov.json")) as f:
        cov = json.load(f)
    al = cov.pop("_all")
    tot = {k: {tuple(x) for x in al[k]} for k in ("lines", "branches")}
    tests = {
        n: {k: {tuple(x) for x in c[k]} for k in ("lines", "branches")}
        for n, c in cov.items()
    }
    return tot, tests


def union(tests, kind):
    return set().union(set(), *(c[kind] for c in tests.values()))


def ranges(numbers):
    out = []
    for n in sorted(numbers):
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in out)


def kind_of(name):
    for k in KINDS:
        if name == k or name.startswith(k + "_"):
            return k
    return name


def inputs_of(test):
    """The inputs of tests/data a test runs (several for the tests that combine them)."""
    return frozenset(
        m.group(1)
        for a in test["command"]
        for m in [re.search(r"/tests/data/([^/]+)\.(dat|toml)$", a)]
        if m
    )


def cover(outdir, tests):
    """Fewest inputs, then cheapest tests on them, keeping the coverage of each kind of test."""
    with open(os.path.join(outdir, "tests.json")) as f:
        info = {t["name"]: t for t in json.load(f)}
    with open(os.path.join(outdir, "times.json")) as f:
        cost = {n: v[1] for n, v in json.load(f).items()}
    idx = {}
    S, I, K = {}, {}, {}
    none = {"lines": set(), "branches": set()}
    for n in info:
        if kind_of(n) in ("run", "run_ae"):
            continue  # a fixture's coverage counts for the tests that read its run
        K[n] = kind_of(n)
        I[n] = inputs_of(info[n])
        c = tests.get(n, none)
        for f in props(info[n]).get("FIXTURES_REQUIRED", []):
            c = {k: c[k] | tests.get(f, none)[k] for k in c}
        S[n] = frozenset(
            idx.setdefault((K[n], k) + x, len(idx))
            for k in ("lines", "branches")
            for x in c[k]
        )
    U = frozenset().union(*S.values())
    by_in = collections.defaultdict(set)
    for n in S:
        by_in[I[n]].add(n)

    def covered(chosen):
        return set().union(
            set(), *(S[n] for ins, ns in by_in.items() if ins <= chosen for n in ns)
        )

    allin = sorted(set().union(*I.values()))
    chosen = set()
    cur = covered(chosen)
    while len(cur) < len(U):
        best = max(
            (i for i in allin if i not in chosen),
            key=lambda i: len(covered(chosen | {i})),
        )
        new = covered(chosen | {best})
        if len(new) == len(cur):
            break
        chosen.add(best)
        cur = new
    for i in sorted(chosen):
        if len(covered(chosen - {i})) == len(U):
            chosen.discard(i)
    print(
        f"\n{len(chosen)} of {len(allin)} inputs keep the coverage of each kind of test:"
    )
    print("  " + " ".join(sorted(chosen)))
    print("not needed: " + " ".join(sorted(set(allin) - chosen)))

    left = set(U)
    sel = []
    cand = [n for n in S if I[n] <= chosen]
    while left:
        b = max(cand, key=lambda n: len(S[n] & left) / max(cost.get(n, 1.0), 0.05))
        if not S[b] & left:
            break
        sel.append(b)
        left -= S[b]
    print(f"\n{len(sel)} of {len(S)} tests on them, by kind (tests, -O0 seconds):")
    for k in sorted(set(K.values())):
        ns = [n for n in S if K[n] == k]
        ks = [n for n in sel if K[n] == k]
        print(
            f"  {k:18s} {len(ks):3d} / {len(ns):3d}  {sum(cost.get(n, 0) for n in ks):7.1f}"
            f" / {sum(cost.get(n, 0) for n in ns):7.1f}"
        )
    print("not needed:")
    for n in sorted(set(S) - set(sel)):
        print(f"  {n}")


def report(args):
    tot, tests = load(args.outdir)
    for kind in ("lines", "branches"):
        u = union(tests, kind)
        print(
            f"{kind}: {len(u)} of {len(tot[kind])} ({100 * len(u) / len(tot[kind]):.1f} %)"
        )
    u = union(tests, "lines")
    print(f"\n{'file':40s} {'lines':>13s}")
    for f in sorted({x[0] for x in tot["lines"]}):
        a = sum(1 for x in tot["lines"] if x[0] == f)
        c = sum(1 for x in u if x[0] == f)
        print(f"{f:40s} {c:5d} / {a:5d}  {100 * c / a:5.1f} %")
    if args.compare:
        _, old = load(args.compare)
        for kind in ("lines", "branches"):
            a, b = union(old, kind), union(tests, kind)
            print(
                f"\n{kind}: {len(a)} to {len(b)}; lost {len(a - b)}, gained {len(b - a)}"
            )
            for x in sorted(a - b):
                print(
                    "  lost",
                    ":".join(map(str, x)),
                    "(covered by "
                    + " ".join(sorted(n for n, c in old.items() if x in c[kind])[:5])
                    + ")",
                )
    if args.uncovered:
        print("\nlines no test runs:")
        miss = collections.defaultdict(set)
        for f, n in tot["lines"] - u:
            miss[f].add(n)
        for f in sorted(miss):
            print(f"  {f}: {ranges(miss[f])}")
    if args.cover:
        cover(args.outdir, tests)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("build")
    r.add_argument("outdir")
    r.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    r.add_argument("-L")
    r.add_argument("-LE")
    r.add_argument("-R")
    p = sub.add_parser("report")
    p.add_argument("outdir")
    p.add_argument("--compare")
    p.add_argument("--uncovered", action="store_true")
    p.add_argument("--cover", action="store_true")
    args = ap.parse_args()
    run_tests(args) if args.cmd == "run" else report(args)


if __name__ == "__main__":
    main()
