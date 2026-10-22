# Test suite speed

The plan for a faster test suite with at least today's coverage (the
user's request, 2026-10-22; the item "Optimize the test suite" of
[smoothness_plan](smoothness_plan.md) holds the measurements). The aim: the
fewest inputs for the broadest code coverage, first matching today's
coverage, then filling its gaps. The `nist` tests are outside it: they
check the all-electron atom against reference data, not coverage, and are
rerun when the all-electron or solver numerics change.

## Done

- [x] Coverage per test: a `-O0 --coverage` build, each test run with its
      own `GCOV_PREFIX`, read with gcov. The non-smooth suite without nist
      runs 88 % of the lines of `src/` and 62 % of the branches gcov
      counts.
- [x] The non-smooth suite cut to the fewest inputs, and the cheapest tests
      on them, that keep the coverage of each kind of test (reference,
      TOML, echo, split, log and UPF mesh, pipe, PSML, psp8 core, well
      onset and tail, projector output, unit): 30 of 37 inputs (29 for
      the coverage, and 73_Ta, the example of the users guide), 179 of 264
      tests (350fe06; `ctest -LE smooth -j 16`: 185 s to 112 s wall, 2729 to 1743
      s·proc; without nist, 150 s to 88 s wall; the same parser for
      both).
- [x] `smooth_scan.py` reads several files of each run (`--output
      log,psp8,psml`), the PSML file among them (c2a6ce2).

## Next

1. [ ] Share one run among the tests that check it: CTest fixtures
       (`FIXTURES_SETUP`/`FIXTURES_REQUIRED`) run each input once, and the
       reference, TOML, echo, PSML and pipe checks read that output.
2. [ ] A cheaper fully-relativistic input for the mode tests (log mesh,
       UPF mesh, echo, split, pipe), which run 52_Te or 74_W at 30–50 s a
       run but test the handling, not the atom.
3. [ ] Merge the smooth rows that differ only in the file they read
       (14_Si rc0, 52_Te rc1, 08_O rc1, and the fcfact and rcfact psp8
       pairs), now that one scan reads several files; and scan every
       channel of a key at once (`rc[*]`, `qcut[*]`, `debl[*]`) in place
       of one row per channel, each channel's kernel still its own
       feature. The tooling is done (c2a6ce2); the rows are not merged.
4. [ ] Commit the coverage tooling: a `coverage` preset in
       `CMakePresets.json` and a script that writes the coverage map and
       the covering-set report, so that both can be redone after any
       change.
5. [ ] Run the suite once in the Debug build (`-fcheck=all`, floating-point
       traps), which it has not been, to find out-of-bounds accesses and
       invalid operations that `-O2` hides.
6. [ ] Turn the gaps in coverage into a list: code that wants a test (the
       last step of the aim) or code that cannot be reached and can go.
7. [ ] Measure how many scan points the smooth rows need: rerun the
       self-test's jumps and kinks at 21 and 31 points; if 21 still catch
       them, the smooth label halves.
8. [ ] Parse only what each scan reads: `smooth_scan.py` parses the whole
       output of each run, and with three files a run that grows
       (pyoncvpsp's parser became 4.7 times faster on 2026-10-22,
       2b17106, which helps either way).
9. [ ] Move `regen.sh`, which regenerates the references in parallel,
       into `scripts/`; with fixtures it reuses the same runs.
