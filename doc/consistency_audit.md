# Consistency audit of ONCVPSP, branch `smooooth` (2026-10-28)

The audit covers commit `cca5f9a`. It looks for one kind of problem:
places where the code does the same thing in two or more ways. That
includes names, numerical methods, errors and messages, interfaces and
data layout, and the scripts, tests and build. For each finding it
proposes a standard and a way to enforce it.

It does not repeat what [Coding Standards](coding_standards.md),
[Fortran best practices](fortran_best_practices.md) and the
[Numerics audit](numerics_audit.md) already cover. Their items appear
here only where the code at `cca5f9a` still breaks them. Files are
relative to `src/` unless they start with another directory.

**How it was made.** Five reviews ran in parallel, one per dimension.
Each parsed the sources with small scripts: 513 procedures, 3333 dummy
arguments, 2424 `@param` lines. Each then read the code to confirm the
counts.
- **Re-checked by hand** for this page: the `-ffast-math` finding
  (M-1), the unchecked `ierr` (E-1), the copies of the mesh step (M-2),
  the literal counts (I-1), the unguarded tests (S-1) and four of the
  wrong `@param` texts (N-1).
- **Approximate:** the other counts may be off by a few, because they
  come from the reviews' parsers.
- **Not part of the audit:** the imported files (`exc_off_pbe_m.f`,
  `functionals.F90`, `m_libxc_list.F90`, `xmlf90-wxml/`,
  `tests/unit/gg1cc|gp1cc|gpp1cc.f90`).

**Status.** Open, all of it. Nothing has been fixed. The proposed
standards are not in [Coding Standards](coding_standards.md) yet.

**What holds already.**
- Every call of a routine with `ifail` checks it afterwards (317 calls).
- `iout,ifail` come last in 108 of 109 routines.
- Every `use` has `only`, and there are no `save` or `data` statements
  and no `d`-specific intrinsics.
- l is indexed the same way everywhere (`l1`, 1-based).
- The GPL block is byte-identical in all 114 files that have it.
- `iout`, `ifail` and `c_light` each have one `@param` text throughout.
- The 55 TOML keys match the users guide.
- The shared smooth layer of logmesh_m and pieces_m is uniform.

The rest is of three kinds:
- helpers kept private and re-derived elsewhere;
- legacy code next to the shared code;
- tolerances, literals and descriptions written differently at each
  site.

## Summary

Ranked by payoff against effort. "Output" says whether a fix changes
what the program writes. A fix that moves the references needs them
regenerated with `scripts/ref_patch.py`.

| # | ID | finding | fix | output |
|---|---|---|---|---|
| 1 | E-1 | `scatterstate`'s `ierr` ignored at 4 call sites; `ug_brent`'s at 1 | check it, or `error stop` in `ldiracfs` | no |
| 2 | S-1 | 2 fr tests registered without the FR guard; 9 Python tests without the pyoncvpsp guard | 5 lines in `tests/CMakeLists.txt` | no |
| 3 | M-1 | `-ffast-math` (Release, `study`) removes the Kahan compensation | `-fno-associative-math` on 3 files, plus a unit test | round-off |
| 4 | N-1 | 5 `@param` texts state the physics wrongly | 5 lines | no |
| 5 | M-3 | C-12's wrong sign in `derivs` still there; dead code in pspot | 3 lines | no |
| 6 | S-2 | ruff without an explicit `select`; dead shellcheck excludes; CI runs the scans | config | no |
| 7 | I-1 | literals `1.d0` (193) and `_dp` (47) beside `1.0d0`; `qp` defined in 7 files; `dp` from `iso_fortran_env` in 2 readers | fortitude S291, C021 with `--fix` | no |
| 8 | I-2 | 14 declarations use a bound before declaring it; `-std=f2008` rejects them | reorder; `-std=f2008` in Debug | no |
| 9 | E-2 | `iout` under `!Output variables` (101) and 19 other label mismatches | scripted move | no |
| 10 | M-2 | the mesh step re-derived 5 ways at 30 sites; index from radius at 8 | public `lm_step`, `lm_index` | round-off |
| 11 | E-3 | error-message form: ~47 wrong prefixes, ~90 without `/`, 49 other formats | message checker, then fix | error text only |
| 12 | M-4 | Brent twice (`lm_root` is a copy of `ug_brent`); 6 hand-written polish loops, 5 stopping rules | one Brent in a leaf module | no |
| 13 | M-5 | tolerances and floors restated per routine (`eps=1.0d-12` in 5 solvers, the 0.2 floor ×6) | named constants | no |
| 14 | N-2 | `@param` texts drift (310 lines against a 29-name glossary); units stated for a minority | glossary + checker | no |
| 15 | I-3 | kappa extent `nkap`/`mxkap`/`2`; 4 idioms for the kappa count; `mxprj` both parameter and argument | one idiom | no |
| 16 | E-4 | 64 module procedures without `@param`, 56 without section labels, 12 descriptions fortls cannot show | write them file by file | no |
| 17 | N-3 | one-letter type components (`config_t%n`, `%l`, `%f`, ...) against rule 5 | rename (~70 sites) | no |
| 18 | S-3 | test names in 6 patterns; labels mix test kind and treatment | one registration helper | no |
| 19 | I-4 | `use` at module level in 53 modules; `implicit none` repeated in 101 module procedures | move; fortitude S201 | no |
| 20 | M-6 | three Fourier–Bessel conventions; four j_l routines; sbfq still normalised by j_0 (C-13) | document; fix C-13 | DIAGNOSTICS |
| 21 | E-5 | `confine_*` stops (`error stop 1`) on a callee's `ifail` instead of passing it on (9) | pass `ifail` on | no |
| 22 | M-7 | relatom and psatom mix in one step; sratom uses Anderson-8 | `ae_mix` everywhere | fr refs ~1e-11 |
| 23 | E-6 | sr/nr energies printed with D exponents, fr with E | with the next regeneration | every sr/nr ref |
| 24 | M-8 | the all-electron energies still use order-4 Gregory sums, copied 4 times | merge now; retire later | retire: every ref |

The sections below give the details by dimension. The last section
collects the proposed standards and how each would be checked.

## Naming (N)

**N-1. Wrong `@param` texts.**
- `ae_scf.f90:38`: the density's normalisation has a spurious 4π.
- `der2exc_m.f90:42`: `zion` is called "ion potential"; it is the ion
  charge.
- `modcore_m.f90:782`: `rhotae` is called "all-electron total charge";
  it is the total valence charge.
- `ea` is documented as all-electron eigenvalues in `results_m.f90:46`
  and as pseudo-orbital eigenvalues in `upfout_m.f90:66`.
- `linout_m.f90:342`: `ee` is called "energies of the projectors"; it
  receives `evkb`, their coefficients.

There are also about 15 typos (`psuedo` ×6).

**N-2. One argument, many descriptions; units missing.**
- `mmax` has 16 different `@param` texts, and `ee` has 14.
- The Header rule asks for units "where they apply". They are stated
  for 5 of 81 potentials, 38 of 133 energies and 49 of 83 radii.
- 874 dummies in 153 mostly private procedures have no `@param` at all
  (see E-4).

**N-3. One-letter names (rule 5).**
- 13 type components in 7 types, e.g. `config_t%n`, `%l`, `%f`.
- The dummy `e` in two internal functions.
- About 70 sites; the compiler catches every one.

**N-4. One quantity, several names.**
- `epgiven` and `ep_as_given`.
- `iunit` and `iout`; `main_setup_m.F90:168` passes a literal 6 where
  `iout` is in scope.
- `ik` and `ikap`; `wmesh` and `wm`; `nrlin` and `nrl`; `vae` and
  `vfull`.
- `psp%echo%fcfact0` is `fcfact0` in linout but `fcfact` in upfout and
  psmlout.
- The "as read" suffix is written three ways: `0`, `_in` and `_f`.
- The configuration triple has four sets of names: `na/la/fa`,
  `cnf%n/l/f`, `nat/lat/fat` and `nacnf_f`.
- Three `input_t` components are abbreviations, unlike the other 52:
  `wreg`, `ecut_wt` and `budget_cut`.

**N-5. One name, several meanings.**
- `uu` is ln r_s at `modcore_budget_m.f90:903`.
- `vv` holds eigenvectors or projector values in four routines.
- `pw` is a real power in two places.
- `kk`, `jj`, `hh` and `fl` also change meaning.
- `mch` is the matching point the bound solvers return, but the point
  integrated to that other routines are given.

**N-6. Argument order.**
- Pairs (`irel,c_light`, `uu,up`) always come in the same order, but
  their place in the list varies.
- `boundstate`, `ldiracfb` and `lschfb` take the same arguments in
  three orders. `zz` and `c_light` are both real scalars, so swapping
  them would compile.
- Proposed: a rule for new code now, and reordering only when other
  work touches these routines.

**N-7. Module, file and procedure names.**
- Public procedures carry a prefix from their module (`lm_`, `ug_`,
  `bs_`, `pw_`, `sbf_`). Nothing writes that convention down, and 9 of
  the prefixes are abbreviations.
- `free_tail_m`, `quad_m` and `kappa_m` break their own scheme.
- `ae_scf.f90` holds `ae_scf_m`. It is the one file whose name breaks
  the `<name>_m.f90` rule.
- 14 public names are used nowhere outside their module.

## Numerical methods (M)

**M-1. `-ffast-math` removes the Kahan compensation.**
- **Where:** `cmake/OncvpspCompilerOptions.cmake:28` gives Release
  `-ffast-math`, and the `study` preset is built the same way: it is the
  build the optimizer runs on. `-fassociative-math` folds
  `(tmp-sum_hi)-yy` to 0.
- **What it breaks:** the compensated sums in `cumint_to`
  (`logmesh_m.f90:232`), `ug_cumint_down` (`ugrid_m.f90:93`) and
  outwf's 4-lane sum.
- **Shown in a test program, not on a full run:** 1e7 additions of
  about 1e-16 to 1. At `-O2` the compensated sum gives 1.0e-9; at
  `-O3 -ffast-math`, with or without `-march=znver3`, it gives 0.
- **Effect, estimated:** about 1e-14 relative, below the noise floor.
  But it makes false the claim at `logmesh_m.f90:233` that the
  round-off does not grow with mmax, and the same claim in the
  numerics audit (§1) and `developer_setup.md`.
- **Fix:** `set_source_files_properties(... COMPILE_OPTIONS
  -fno-associative-math)` for GNU on those three files, plus a unit test
  that fails when the compensation is lost (about 20 lines). The test
  catches it on any compiler.

**M-2. The mesh step, and the index from a radius.**
- **The mesh step:** `al` is derived 5 ways at 30 sites. 22 of them are
  `0.01d0*log(rr(101)/rr(1))`, in `ae_scf.f90`, `lschfs_m.f90` and
  others. It is also passed as an argument to 15 procedures.
- **The index from a radius:** found from a closed form at 8 sites. Two
  of them correct the guess upward only, so they could be off by one at
  a radius exactly on a node. That is possible from the code, but not
  observed.
- **Fix:** `mesh_step`, `mesh_position` and `node_index` already exist
  in logmesh_m but are private. Made public (`lm_step`, `lm_index`),
  they replace all 38 sites.
- **Output:** moves at 1e-16. The reviews disagree on whether that
  reaches a reference digit; the suite will say.

**M-3. Ad-hoc derivative stencils.**
- C-12 of the numerics audit, the sign of the ρ″ end formulas in
  `derivs` (`derivs_m.f90:93,98`), is still wrong.
- pspot's second stencil pass for u″, and the values it writes beyond
  rc, are dead code.
- pspot takes u′ from a 7-point stencil, where the convention elsewhere
  is the solver's `up`. Changing it moves the kinetic-energy test line.

**M-4. Root finders.**
- `lm_root` is a line-for-line copy of `ug_brent`, because ugrid_m uses
  logmesh_m and neither can call the other. One Brent in a leaf module
  would remove about 70 lines and change no output.
- Six hand-written bracket-and-polish loops (qroots and others) stop by
  five rules: 4 eps, 8 eps, `4*spacing`, eps·(2|x|+1e-3), and fixed
  absolute values.
- Proposed: one relative-plus-absolute rule, written once.

**M-5. Tolerances, floors and limits restated.**
- `eps=1.0d-12` is a local in 5 eigen-solvers.
- The solvers' iteration limits are 60 or 100.
- The 0.2 Ha energy floor is written 6 times; the SCF floors at 1.
- `0.95d0*rc` (the relativistic window) appears 9 times.
- The `1.0d-17` series cut-off appears 3 times.
- `eps` means five different things.
- Proposed: named parameters in constants_m (or in the module that
  owns them), such as `eig_tol`, `eig_floor`, `rel_window`, `qp`,
  `pi_qp` and `ha_ev`. `Ha_eV` is local to run_optimize now. `qp` is
  defined in 7 files.

**M-6. Bessel functions and Fourier–Bessel transforms.**
- There are four routines for j_l values, plus a private `j0`.
- `sbfq` is still normalised by j_0 (C-13 of the numerics audit).
- `sbf_rc_der` copies an uninitialised value that is never used
  (`sb_out(llin+6)`).
- The three transform engines (eresid, qroots and vloc_cutoff) use
  different prefactor conventions.
- vloc_cutoff's sine recurrence is never re-seeded. Fixing that moves
  DIAGNOSTICS.
- Proposed: document one convention for the transform (prefactor, sign,
  q units) and keep one source of j_l per precision.

**M-7. SCF mixing.**
- sratom mixes with Anderson-8 (`ae_mix`); relatom and psatom still mix
  in one step.
- The SCF tolerances are spelled three ways.
- Moving relatom to `ae_mix` moves the fr references by about 1e-11.

**M-8. Legacy quadrature.**
- The all-electron energies still use the order-4 Gregory sums: in
  `radint` and in three inline copies in vout.
- Merging the copies into one function changes nothing.
- Retiring them in favour of `lm_cumint` moves every reference.
  Moving lschpsbar to outwf would too. Both were deferred earlier.

**M-9. Differences that stay, to be documented only.**
- The `adams.inc` include: inlining it is worth about 9 % of run time.
- The three RK4 loops. Two of them use a "2e-3" step, in different
  units.
- The finite-difference steps kept as parameters.
- eresid's 16-digit Gauss–Legendre table. The review assumed it is
  within an ulp of `gauss_legendre`, but did not compare them.
- The mesh-end conventions (`mch`, `mtp`, `n2`, zeroing beyond the last
  point) and the sign conventions of wave functions and projectors.
  These should be one table in program_overview.md.
- The legacy sr/fr difference in upfout.

## Errors, output and in-source documentation (E)

**E-1. Ignored `ierr`.**
- `scatterstate` prints `ldiracfs: ERROR` for a bad argument and
  returns `ierr` 2 or 3.
- Four callers ignore it and carry on with zero wave functions:
  `phsft_m.f90:491` and `:614`, `confine_ws_m.f90:603` and
  `run_diag_m.f90:195`.
- `ws_basis_m.f90:133` ignores `ug_brent`'s `ierr`.
- The arguments are checked upstream, so this is a bug if it happens.
  Rule 5 of the error messages calls for `error stop 1` on `error_unit`
  in `ldiracfs`. That covers 4 of the 5 sites, and the fifth needs its
  own check.

**E-2. Section labels.**
- 120 declarations sit under a label that doesn't match their intent.
- 101 of them are `iout`, `intent(in)` under `!Output variables`, in 51
  files. A script can move them.
- 7 routines have the labels out of order, and 3 have declarations
  before the first label.

**E-3. Message form (Error messages, rules 2–4).**
- About 47 new messages name another routine or none:
  - helpers using their public routine's name (`read_ae_file`,
    `modcore_budget`);
  - `oncvpsp` in code that has moved out of the main program;
  - module names such as `quad_m:`.

  Internal procedures use their host's name in 22 of 23 cases. The
  rules don't settle helpers. Proposed: the routine's own name, its
  host's, or a public procedure of the same module.
- About 90 new messages have no leading `/`.
- 49 use edit descriptors other than `i0`, `f12.6` and `es12.4`. `rc`
  alone is printed as `f10.6`, `es12.4` and `f12.6`.
- 13 messages in modcore_budget say "ERROR stop -" although the routine
  returns `ifail`.
- 9 messages continue a string across lines, which copies the
  indentation into the text.
- Fixing them changes error text only, except `ae_floor_warning`
  (6 references) and three regexes in `tests/CMakeLists.txt`.
- The 59 messages from 4.0.1 keep their text (a generated allowlist).

**E-4. Header documentation.**
- 64 module procedures have no `@param`, mostly in new or recently split
  files: `confine_ws_m`, `main_setup_m`, `ae_projectors_m`, `cli_m`,
  `ps_valence_m`, `kappa_m`.
- 56 have no section labels.
- 7 `@param` lists name the wrong dummies, and 10 are in the wrong
  order.
- 12 descriptions are plain `!` comments, which fortls doesn't show.
- Only 2 of 85 functions have `@return`. Proposed: drop `@return` from
  the rule rather than add 83 of them.
- 7 types document their components only above the `type` statement.
  14 declarations put one `!>` over several names, so about 75
  components show nothing on hover.

**E-5. Stops in library routines.**
- 9 `error stop 1` in the `confine_*` modules stop on a callee's
  `ifail` instead of passing it on.
- `psmlout_m.f90:744,749` use `error stop "<string>"`.
- `core_radius_points` writes an ERROR to unit 6 outside the command
  line.
- The 10 bare `stop`s are all success exits in the programs or the
  command line.

**E-6. Output formats.**
- `run_config_m.f90:254` (sr and nr) prints the total energies with D
  exponents; `:296` (fr) prints the same lines with E. `1p,e` and `es`
  both occur.
- Fixing it moves every sr/nr reference, so it waits for a planned
  regeneration.
- The 9 two-line program stubs (`oncvpsp*.F90`) have no copyright
  header. [Coding standard violations](coding_standard_violations.md)
  doesn't list them, and its note on the years in `linout_m` belongs to
  `psmlout_m`.

## Interfaces, data layout and idiom (I)

**I-1. Kinds and literals (Kinds and constants).**
- Literal styles:
  - `1.0d0` is the standard, but `0.d0`-style literals occur 193 times
    in 22 files and `_dp` literals 47 times in 10 files.
  - 3 literals have no kind.
  - `tfapot` writes 1/3 as `.33333333d0`.
- Kinds:
  - `input_text_m` and `input_toml_m` take `dp` from `iso_fortran_env`.
  - `qp => real128` is defined in 7 files.
- Conversions: `dble()` (159) and `real(x,dp)` (20) are both used.
- fortitude S291 (bare decimal, 194 hits, fixes itself) and C021
  (no real suffix, 3 hits) catch most of this.

**I-2. Declaration order.**
- 14 declarations in 11 files use an array bound before declaring it,
  e.g. `rr(mmax)` under `!Input variables` before `integer :: mmax`.
  The files include `der2exc_m.f90:58`, `lschfb_m.f90:66`,
  `lschvkbb_m.f90:39` and `wellstate_m.f90:101`.
- `-std=f2008` rejects them. This is §5 of Fortran best practices,
  which [Coding standard violations](coding_standard_violations.md)
  counts as closed ("Open: None").
- One implicit real-to-integer conversion is left:
  `main_setup_m.F90:205`.
- The `-std` pass stopped at these errors, so expect a few more.

**I-3. The kappa extent and the projector dimension.**
- The kappa extent is written `nkap` 76 times, `mxkap` 14 times and as
  a literal `2` 4 times; `ea` alone uses all three.
- The number of kappas for an l is computed four ways: `kap_count`, a
  duplicate `nkap_of` in ps_file_m, an inline `merge` (3 sites), and
  `do ikap=1,2`.
- In 4 loops in upfout, `kap_count` is passed a variable named `l1`
  that holds l.
- `mxprj` is a constants_m parameter, yet 18 routines take it as an
  argument. Declarations then mix the two, e.g.
  `evkb(mxprj,mxl+1,nkap)` at `diagout_m.f90:163`.

**I-4. `use` and `implicit none`.**
- 143 `use` statements sit at module level in 53 modules; in 20 of
  those modules nothing at module level needs them. The Files rule
  asks for them in the procedure that makes the call.
- `implicit none` is repeated in 101 module procedures and absent from
  291. The template in [Coding Standards](coding_standards.md) shows it
  in every routine.
- Proposed: module level only. fortitude S201 fixes it.

**I-5. Derived types.**
- `input_t` (51 components), `budget_t` (52) and `ae_entry_t` have no
  default values, while 20 other types are fully defaulted.
- The input's defaults are set in each reader, and both readers declare
  `inp` `intent(out)`.
- The reference configuration is stored twice in `input_t`, as
  `na/la/fa` and as `cnf(1)`, with a copy step at
  `oncvpsp_run_m.f90:237`.
- `budget_t` dummies that are only read are `intent(inout)`.

**I-6. Smaller items.**
- 33 assumed-shape dummies, which the Declarations rule doesn't allow.
  Proposed: allow them inside the library and keep explicit shapes at
  the entry points.
- Integer 0/1 flags (`vout`'s `mode`, `run_diag`'s `iso`) beside 58
  logical dummies.
- `deallocate` is explicit in 78 procedures, partial in 6 and absent
  in 10.
- 25 functions modify their arguments, and 9 take `ifail`.
- [Library plan](library_plan.md) names `oncvpsp_run` a C-callable
  entry point. It takes a logical, strings and derived types, which the
  entry-point rule forbids, and Phase 5 wraps it behind a handle anyway.
- `renorm_r` puts the optional `rnorm` after `ifail`.

## Scripts, tests, build and user-facing names (S)

**S-1. Tests registered without the helper's guards.**
- `well_rmax_r_14_Si` and `well_norm_r_confined_14_Si`
  (`tests/CMakeLists.txt:608`, `:634`) call `add_test` directly with
  `$<TARGET_FILE:oncvpspr>`. The `oncvpsp_test` helper would skip them
  when `ONCVPSP_WITH_FR` is off. As it stands, a configure without the
  fr executable should fail; read from the code, not reproduced.
- The 9 `well_tail`, `well_rmax` and `well_norm` tests don't test
  `oncvpsp_compares`. Without pyoncvpsp they fail on an ImportError
  instead of being left out.
- The cause is that tests run by Python don't go through the helper.

**S-2. Lint and CI.**
- `ruff.toml` has no explicit `select`. The pinned ruff 0.14 is clean,
  but ruff 0.17 finds 103 problems in `scripts/` with the same file, so
  the next autoupdate breaks the hook.
- The shellcheck exclude list names 5 files that no longer exist.
- The mkdocs navigation leaves out `numerics_audit.md` and
  `smoothness_plan.md`.
- `.github/workflows/cmake_ctest.yaml` runs ctest without a preset, so
  the 31 scans run there too. No CI job configures with FR or Libxc
  off, which would have caught S-1.

**S-3. Test names and labels.**
- Test names follow six patterns.
- Labels mix test kind (`toml`, `split`, `mesh`) with treatment (`sr`,
  `r`, `nr`). `ctest -L r` picks the 10 fr reference tests plus 6 TOML
  tests.
- `run` is missing from the label list in `developer_setup.md`.
- `occupation_ramp` is labelled `mesh`.
- Proposed: one registration helper that generates the name and a
  `kind;treatment` label pair.

**S-4. Relativity names.**
- Each treatment is spelled about four ways across the executables,
  targets, source files, program names, test names, labels, reference
  files and docs.
- `oncvpsp_nr.F90` and `oncvpsp_r.F90` don't match `oncvpspnr.x` and
  `oncvpspr.x`.
- The nr program's errors say `oncvpsp:`.
- `unit_run_*` uses 0, 1 and 2.
- In DIAGNOSTICS, `_r` means a radius (`ld_r`, `ghost_r`).
- 11 error prefixes name no routine, among them `run_plot_r`.

**S-5. Script code.**
- About 200 lines of Python sit in shell heredocs that no linter
  checks, 145 of them in `confine_test.sh`.
- The shell boilerplate is written several ways.
- `fortran_float` is written three times.
- Two TOML libraries are used.
- Five scripts patch `sys.path`.
- Tolerances live in five places and use three definitions of
  "relative".
- `kernel_diff.py` deliberately duplicates `compare_output.py`'s kernel
  comparison so that it needs only the standard library. Proposed: a
  test that checks the two agree.

**S-6. Fortran unit tests.**
- The same `check()` helper is copied into 8 files.
- All 20 take `dp` from `iso_fortran_env`, and 3 define `pi` locally,
  against Kinds and constants.

**S-7. Docs.**
- Six live docs name things that no longer exist:
  - `core_correction.md:57` `modcore.f90`;
  - `program_overview.md:44` `SR_SHOOTING_VPRIME.md`;
  - `known_issues.md:62` `vrel`, `run_config`;
  - `numerics_audit.md:189`;
  - `libxc_use.md:32` `compare.sh`;
  - `smoothness_plan.md:1822` `ctest -LE smooth`.
- About 90 other missing names are in historical records, on purpose.

## Proposed standards and their enforcement

Rules to add to [Coding Standards](coding_standards.md), and what would
check each one. "Checker" means one of the scripts proposed below.

| standard | check |
|---|---|
| Literals `1.0d0`; kinds `dp`/`qp` only from constants_m; `real(x,dp)` or `dble` (pick one) | fortitude S291, C021 (enable); checker for `=> real64/real128` outside constants_m |
| Integer bounds declared before the arrays they bound; standard F2008 | `-std=f2008 -Wconversion` in the Debug flags |
| `implicit none` at module level only; `use` in the procedure that needs it | fortitude S201 (enable); checker for module-level `use` not needed there |
| Compensated sums compiled without reassociation | per-file `-fno-associative-math`; unit test of a known sum |
| The mesh step and the index of a radius only from logmesh_m | checker: `log(rr(` outside logmesh_m |
| Tolerances, floors, iteration limits and repeated factors are named parameters | checker: `1.0d-NN` literal in a comparison outside a parameter; review |
| One Brent, one Gauss–Legendre, one j_l per precision, one SCF mixer | review; the duplicates removed |
| A glossary of argument names: canonical `@param` text, unit, type and rank | checker against `doc/glossary.md` |
| No one-letter names, type components included | checker |
| Public procedures carry their module's prefix; module `x_m` in `x_m.f90` | checker; the prefix table in the standards |
| Section labels match intent; `iout` is an input | checker (with a `--fix` for `iout`) |
| Every module procedure has a `!!` description and one `@param` per dummy, in order; `@return` dropped | checker |
| `!>` or `!<` on every module declaration and type component, one name each | checker |
| Error-message prefix is the routine, its host, or a public procedure of its module; `/`, `i0`/`f12.6`/`es12.4`, no string continued inside; 4.0.1 texts allowlisted | checker |
| `ierr` and `ifail` checked after every call; no `stop` in module procedures except rule-5 checks | checker |
| Copyright and GPL block in every file | checker |
| Kappa extent `nkap`; the kappa count by `kap_count`; `mxprj` from constants_m only | checker for literal `2` or `mxkap` in a kappa dimension; review |
| Arguments of new routines in the order mesh, potentials, state, sizes, `iout`, `ifail` | review only |
| Every test registered through one helper that applies the guards and generates name and labels | the helper; CI job with FR and Libxc off |
| ruff with an explicit `select`; shellcheck on every script | config |
| Docs name only files and routines that exist (historical records excepted) | `docrefs.py` as a ctest |
| Fortran unit tests use one assertion module and constants_m | review; the shared module |

**The checkers.**
- **Shape:** one script, `scripts/check_conventions.py`, built on a
  small statement parser of about 150 lines. The parser joins
  continuation lines, tracks procedures and their hosts, and skips
  interfaces and types.
- **Rules:** sections, headers, messages, errors, module documentation,
  license, names against the glossary, numerics literals, kinds and
  `use` placement, each a function of about 20–80 lines.
- **Runs:** in pre-commit on the changed files, and as a ctest test
  over `src/`. Each run takes well under a second.
- **Baseline:** a generated baseline file (`scripts/conventions_baseline.txt`)
  lets the checker go in before the fixes. It then fails only on new
  violations, and the baseline shrinks as code is fixed.
- **What it can't do:** neither fortitude nor gfortran checks
  documentation, names or messages. gfortran 11 has no
  `-Wimplicit-pure`.

**fortitude rules not to enable.** These conflict with the standards:
- MOD002, double-precision literals (2011 hits);
- S251, useless return, which conflicts with `return` before `end`;
- S262, mixed scalar and array declarations, which conflicts with the
  template;
- S221, function without a result variable (106 hits);
- S271, unsorted `use` statements, which is style only.

C031 (51 hits) and MOD031 (9 hits) are optional.

**Order of work.**
1. The bugs, which change no output: E-1, S-1, N-1 and M-3's sign.
2. M-1, with a full suite run.
3. The configuration: S-2, the fortitude rules with `--fix`, and
   `-std=f2008` with I-2.
4. The checker with its baseline. Then the mechanical fixes it lists
   (E-2, I-4, M-2, M-5, I-3, N-3), file by file.
5. The output-moving items (M-6's re-seeding, M-7, E-6, then M-8) in
   one batch, at the next reference regeneration.
