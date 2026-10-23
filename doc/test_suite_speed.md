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
5. [x] Run the suite once in the Debug build (`-fcheck=all`, floating-point
       traps), which it has not been, to find out-of-bounds accesses and
       invalid operations that `-O2` hides. Done 2026-10-23 (the `debug`
       preset, gfortran 13.3, on a Zen 4 node): one defect, an overflow in
       `qjacobi`, whose sweeps chased off-diagonal elements into underflow
       until theta**2 overflowed (tt = 0 at -O2, the right limit; 88 runs
       trapped in `sbf_basis_con`). Fixed in 98f6850 (nothing printed
       moves); then 290 of 290 pass, 246 s wall at -j48, with no bounds,
       invalid-operation or division trap.
6. [ ] Turn the gaps in coverage into a list: code that wants a test (the
       last step of the aim) or code that cannot be reached and can go.
7. [x] Measure how many scan points the smooth rows need. The self-test's
       jump of 1e-7 and kink of slope 1e-6, moved through the window, are
       caught (residual above 1e-8 and 1e-9) with 21 points to degree 10
       (the jump at 1.09e-8), 31 to degree 14 and 41 to degree 18, so 21
       points do not do for every row. NPTS now defaults by the fit's
       highest degree: 21 to degree 6 (16 rows), 31 to 14 (13), 41 above
       (14_Si_bc_ecut), 780 runs in place of 1230; the self-test checks the
       rule. Every row passes at its new count, the largest residual 0.83
       of its tolerance (34_Se rcfact, psp8 core; nr_08_O_rc's kernel
       7.3e-11 against 7.0e-11 at 41), and the four rows that see F-1 (the
       psp8 or linear-mesh V_loc of 52_Te rc[*], 14_Si fcfact, 34_Se
       rcfact) still fail at 21, by 6e-9 to 8e-8. `ctest -L smooth -j 16`:
       469 to 372 s wall, 6866 to 5901 s·proc, measured under a load of
       about 20 on 28 cores (other sessions' runs), so less than the 37 %
       fewer runs; the wall time is set by 14_Si_bc_ecut (116 s).
8. [x] Parse only what each scan reads: `smooth_scan.py` asked for the
       parser's whole output; it now asks only for the parts of its groups
       (the files for kernel, vloc and core, the test configurations, the
       DIAGNOSTICS block). The features of the 16 scans of item 3 are the
       same; a parse takes 0.034 s in place of 0.092 s (pyoncvpsp 5deacd6,
       whose parser became 4.7 times faster on 2026-10-22), about 1 % of
       the smooth label, whose time is the oncvpsp runs.
9. [ ] Move `regen.sh`, which regenerates the references in parallel,
       into `scripts/`; with fixtures it reuses the same runs.
