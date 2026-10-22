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

1. [x] Share one run among the tests that check it: the CTest fixtures
       `run_<rel>_<input>` and `run_ae_<rel>_<input>` (`scripts/run_test.sh`)
       run each input, and its all-electron file, once; the reference, TOML,
       pipe, split and `ae_file_*` tests read them. The log and UPF mesh
       scripts run the input itself only in the cases that compare with it
       (none of today's). Echo and PSML tests change the input, so run their
       own. `ctest -LE "smooth|nist" -j 16`: 88 s to 74 s wall, 1354 to
       1128 s·proc; with nist, 1737 to 1492 s·proc, the wall time (112 s)
       set by `nist_RLDA_heavy`.
2. [x] The fully-relativistic mode tests (log and UPF mesh, echo, split,
       pipe, PSML, well onset, psp8 core, TOML) on inputs that run in
       seconds in place of 52_Te, 74_W and 80_Hg: 14_Si, 14_Si_UPF,
       14_Si_lder, 40_Zr_grid; 03_Li for the split (its test configurations
       differ only in l, as 80_Hg's, which the all-electron file's matching
       needs) and 60_Nd_GHOST for the psp8 core (it reaches the edges of
       `modcore_budget` that only 74_W did). The coverage of `src/`
       without nist and smooth loses nothing and gains 11 lines and 8
       branches (8702 to 8713, 9795 to 9803). 52_Te_r and 80_Hg_r stay as
       references, the inputs the docs measure. `ctest -LE "smooth|nist"
       -j 16`: 75 s to 70 s wall, 1152 to 959 s·proc; with nist 1499 to
       1336 s·proc (the wall time, 115 s, still `nist_RLDA_heavy`'s).
3. [x] Merge the smooth rows that differ only in the file they read
       (14_Si rc0, 52_Te rc1, 08_O rc1, and the fcfact and rcfact psp8
       pairs), now that one scan reads several files; and scan every
       channel of a key at once (`rc[*]`, `qcut[*]`, `debl[*]`: every
       element moved by the same offset, `--center` the offset's) in place
       of one row per channel, each channel's kernel still its own
       feature: 41 rows to 30 (`ctest -L smooth -j 16`: 599 s to 469 s
       wall, 9006 to 6866 s·proc). Every channel at once on the linear meshes
       found PZ81's jump at r_s = 1 spread by their interpolants (F-1 in
       [smoothness_plan](smoothness_plan.md)), which the rows with one
       channel, or the core on the psp8 file's core only, had not moved:
       52_Te keeps rc[1] for the linear meshes, and the psp8 file's V_loc
       is left out of the core rows.
4. [x] Commit the coverage tooling: the `coverage` preset (`build-cov/`)
       and `scripts/coverage.py` (`run`: each test with its own counters,
       the fixtures first; `report`: totals per file, `--compare` with an
       earlier run, `--uncovered`, and `--cover`, the fewest inputs and tests
       that keep the coverage of each kind of test). It reproduces the
       measure of item 2 (8713 lines, 9803 branches; 2.5 min on 28 cores)
       and its cover the 28 inputs of the suite (73_Ta is the users guide's
       example).
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
