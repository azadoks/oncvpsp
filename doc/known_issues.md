# Known issues

Suspected bugs found while auditing dummy-argument `intent`s across `src/`,
then triaged against the code, the callers, the input checks in
`check_data.f90`, and the test inputs in `tests/data/`. The P1 and P2
issues have been fixed, and so have the P3 issues apart from dead code
that is kept on purpose. Line numbers refer to the code before the fixes.

The issues were found on `master`, and the fixes were then merged into
`oncvpsp-opt`, which had fixed some of them independently. Where the two
fixes differed, the entry says which was kept.

Later issues come from the audits in
[Fortran Best Practices](fortran_best_practices.md) and
[Coding Standard Violations](coding_standard_violations.md). Their line
numbers refer to commit `a672619`.

Priorities:

- **P1**: changes results in common use.
- **P2**: changes results with uncommon options, crashes on some compilers,
  or hides errors.
- **P3**: latent, unreachable in practice, dead code, or cosmetic.

"Refs" says whether fixing the issue would change the reference outputs in
`tests/`.

## Summary

| Priority | Issue | Refs | Status |
|---|---|---|---|
| P1 | `der2exc.f90:79` stencil sign | yes (all `icmod>=1` tests) | fixed |
| P2 | `m_psmlout.f90:434,969` deallocates a pointer to dummy memory | no | fixed |
| P2 | `run_config(_r).f90:152` stale `ii` | negligible (79_Au) | fixed |
| P2 | `modcore3.f90:414` stale `yy` in the `ircross` search | yes (34_Se) | fixed |
| P2 | `psatom.f90:216` non-convergence check can't fire | no | fixed |
| P2 | `psatom_r.f90:184` `lschvkbb` error masked | no | fixed |
| P2 | `run_phsft(_r).f90` `irphs` unset for large `rxpsh` | no | fixed |
| P2 | `input_toml_m.f90` required keys read without `stat` | no | fixed |
| P2 | `excwig.f90` `vxc`/`exc` unset for `rho<1e-15` | no | fixed |
| P2 | `vkbphsft.f90` `dnpi` can be undefined | no | fixed |
| P2 | `oncvpspnr.x` is scalar-relativistic | new tests (08_O_nr, 14_Si_nr) | fixed |
| P2 | psp8 file without the model core for `icmod=6`, `fcfact<=0` | yes (14_Si_bca) | fixed |
| P3 | build, dead code, silent errors, bounds, argument mismatches (see [P3 (fixed)](#p3-fixed)) | no | fixed |
| P3 | 5 latent bugs: `qroots`, `ldirac*` series starts, `icmod=4` grid, `_tmp_input`, exit status of errors | no | fixed |
| P3 | `functional_get_tau` (imported dead code) | no | kept |
| P3 | Adams 1/24 coefficient in `adams.inc` | last digits only; refs pass unchanged | fixed |
| P3 | input echoes that don't run again to the same output | echo lines only (5 refs) | fixed |
| P3 | 5 cosmetic items | no | fixed (unused arguments kept) |
| P3 | `rxpsh` inside a core radius not rejected | no | fixed |
| P3 | `oncvpspnr.x` banner says "scalar-relativistic version" | no (the nr refs already say non-relativistic) | fixed |

## P1 (fixed)

- **`der2exc.f90:79`**: the 4th stencil point is `rhotot - 2*hh*rho`, a
  copy of the 1st, instead of `rhotot + 2*hh*rho`. The ±2h terms cancel, so
  the result is 4/3 times a 2nd-order central difference. Every `d2exc`
  element is therefore about 33% too large. This goes back to upstream
  version 4.0.1.
  - Used by `modcore`, `modcore2` and `modcore3`, i.e. any `icmod>=1`.
  - For `icmod` 1–3 the error is confined to printed diagnostics:
    `d2excae`/`d2excps` and the "rms 2nd derivative error".
  - For `icmod=4` it scales the coarse-scan and Nelder–Mead objective.
    Because the scaling is nearly uniform, the optimised model core shifts
    only slightly.
  - Fix: `rhot(:)=rhotot(:)+2.d0*hh*rho(:,kk)`.
  - Refs: yes. The tests affected are 14_Si, 14_Si_UPF, 34_Se, 38_Sr_lxc,
    40_Zr, 60_Nd_GHOST, 79_Au_lxc, 83_Bi and 74_W (`icmod=4`). Only the
    printed `d2exc` tables and rms errors change; they are scaled by 3/4.
    At test tolerance the optimised model cores and psps are unchanged.

## P2 (fixed)

- **`m_psmlout.f90:434,969`**: `deallocate(chcore)` is called after
  `chcore => rhomod(:,1)`, and `chcore` is never `allocate`d. This is
  invalid Fortran.
  - gfortran happens to free the base of `rhomod` silently, and nothing uses
    `rhomod` afterwards.
  - Other compilers (ifort, nvfortran) can abort.
  - Reached with PSML output and `icmod>=1`. No test covers this
    combination.
  - Fix: delete both `deallocate`s.
- **`run_config.f90:152`, `run_config_r.f90:152`**: the "maximally ionized
  intermediate state" loop runs over `kk`, but it writes `fat(ii,2)` with a
  stale `ii`, so the loop does nothing.
  - `fat(:,2)` stays equal to the target occupations. The intended
    `min(reference, target)` gives a different value only when a target
    increases a reference-level occupation. In the tests, that only happens
    in 79_Au_lxc config 3.
  - Only the adiabatic path and the starting guess change, not the final
    self-consistent state, so this is a robustness issue.
  - Fix: `fat(kk,2)=dmin1(fat(kk,1),fat(kk,3))`.
- **`modcore3.f90`**: the `rcross` loops call `gg1cc(gg,yy)` without
  updating `yy=rr(ii)/x(2)`.
  - In the search loops (lines 170, 215, 278, 310, 346, 377), `rcross`
    only feeds an unused `tt`. That is dead computation.
  - Only the final loop (line 414), which sets `ircross`, matters.
    - For `icmod=4` the stale `yy` makes `gg=0`, which gives the same
      `ircross` as the correct formula.
    - For `icmod=3`, `yy=0`. In 34_Se this puts `ircross` at r≈0.99
      instead of r≈2.57. That moves the switch from analytic to numerical
      derivatives and the smoothing width `drint`. So the higher model-core
      derivatives (`rhomod(:,3:5)`) change slightly, but `rhomod(:,1)` does
      not.
  - Also, `ircross` is left undefined if no crossing is found. That needs an
    absurd `fcfact`.
  - Fix: set `yy` before each `gg1cc` call. Also initialise `ircross=0` and
    stop if it is never set.
  - Refs: 34_Se only. It is uncertain whether upstream intended this
    behaviour.
- **`psatom.f90:216`**: the non-convergence check is `if(it==400 ...)`
  inside `do it=1,100`, so it can never fire.
  - An unconverged test configuration returns `ierr=0`, and its excitation
    energies are printed as if valid. The pseudopotential itself is
    unaffected.
  - The tests converge in 3–68 iterations.
  - `psatom_r` uses `it==100` correctly.
  - Fix: `it==100`.
- **`psatom_r.f90:184`**: on an `lschvkbb` error, `exit` leaves only the
  `ikap` loop. The next `ii` iteration resets `ierr=0`, so the error is
  masked unless it occurred on the last state, and `rho` is missing that
  component.
  - Fix: add `if(ierr/=0) exit` after `end do !ikap`.
- **`run_phsft.f90`, `run_phsft_r.f90`**: if `rxpsh > rr(mmax)` (about 45
  bohr), `irphs` is never assigned.
  - On the first `l`, `irphs` holds garbage. On later `l`s it holds the
    previous value.
  - It is then used to index `rr(irphs)` and passed to
    `fphsft`/`vkbphsft`, so the read can go out of bounds.
  - The `irphs<=irc` check only prints an error, and the message prints
    `rr(irphs)` under the label "rxpsh".
  - Reached only when the user sets `rxpsh`. The default is -1 in both
    readers, and no test sets it.
  - Fix: set `irphs=0` before the search, and stop if it is still 0 or
    `<=irc`.
  - Fixed: it stops with `run_phsft: ERROR rxpsh beyond radial mesh`
    (test `error_rxpsh_beyond_mesh`). `oncvpsp-opt` had instead clamped
    `irphs` to `mmax-5`; the error was kept, because a log derivative at a
    radius other than the one asked for is easy to miss.
- **`input_toml_m.f90`**: required-key `get_value` calls have no `stat` or
  default.
  - toml-f's `val` is `intent(out)`, so a missing key silently leaves the
    value undefined.
  - `check_data` catches only some of these cases.
  - A missing `psfile` or `atsym` leads to deallocating an unallocated
    `tmp_str`.
  - Fix: pass `stat=` and stop on failure, or use a small `require_*`
    helper.
  - No test uses TOML input.
  - Fixed: the reader of `oncvpsp-opt` (`e083200`), which checks every key
    and has its own tests (label `toml`), replaced the `master` fix.
- **`excwig.f90`**: `vxc` and `exc` are not set where `rho<1e-15`, and
  that covers every density tail. `vout`'s `exca`/`vxcd` are not zeroed, so
  `eexc` and `difxc` pick up `rho*garbage`. If that garbage is NaN, `eexc`
  becomes NaN. `exchdl` and `excpzca` zero these values.
  - Affects every `iexc=1` (Wigner) run. No test uses `iexc=1`.
  - Fix: add an `else` branch that zeroes both.
- **`vkbphsft.f90`**: `dnpi` is set only if some scan energy lies within
  `depsh/2` of `epa(1,l)`. Otherwise, garbage is added to the whole printed
  pseudopotential log-derivative curve. This happens, for example, for a deep
  semicore projector below the default `epsh1=-2`. No test hits it.
  - Fix: set `dnpi=0` at entry, optionally with a warning.
  - Fixed: when `ep` is outside the scan, the curves are aligned at the
    nearer end of the scan (from `oncvpsp-opt`, which replaced the
    `dnpi=0` fix). Test `14_Si_lder` has its s projector below `epsh1`.

- **`oncvpspnr.x` is not non-relativistic.** CMake builds it from
  `oncvpsp.F90` with `RELATIVISTIC=0` (`src/CMakeLists.txt:121-124`), but no
  source reads that macro. `oncvpsp.F90:114` hard-codes `srel=.true.`, so
  the binary is identical to `oncvpsp.x`. The binary is installed, and
  anyone running it silently gets scalar-relativistic pseudopotentials. No
  test covers it.
  - The old non-relativistic driver, `oncvpsp_nr.F90`, sets
    `srel=.false.`, but CMake doesn't build it and it is stale (see
    [P3 (fixed)](#p3-fixed)).
  - Fixed: `oncvpsp_nr.F90` is now a wrapper that defines
    `ONCVPSP_NONREL` and includes `oncvpsp.F90`, which then sets
    `srel=.false.` (from `oncvpsp-opt`, which replaced the `RELATIVISTIC=0`
    fix on `master`). Both builds use it.
- **`linout.f90`, `linout_r.f90`: no model core in the psp8 file for
  `icmod=6` with `fcfact<=0`.** The psp8 writers wrote `fcfact` as
  `fchrg` and wrote the model core only when `fcfact>0`. Since
  `f098dc9` on `oncvpsp-opt`, `fcfact<=0` with `icmod=6` takes the cutoff
  from the pseudo-waves, so those psp8 files had `fchrg=0` and no core,
  although the pseudopotential was made with one. The UPF and PSML
  writers test `icmod>=1` and were right.
  - Found while checking the input echoes (2026-10-19).
  - Fixed: `fchrg` is `fcfact` when `icmod>=1`, or 1 when `fcfact<=0`, and
    the core is written when `fchrg>0`. Inputs with `icmod=0` or
    `fcfact>0` write the same psp8 file as before. The psp8 block of the 14_Si_bca ref, which had
    recorded the missing core, now has it. Tests `psp8_core_sr_14_Si` and
    `psp8_core_r_74_W`.

## P3 (kept, or fixed later)

### Dead code

- **`functionals.F90`**: in `functional_get_tau`, `n(1)`/`t(1)` are stale
  for `nspin=2`. The routine is never called. Kept, because the file is
  imported and is easier to keep in sync with upstream unchanged.

### Cosmetic

Fixed in the mechanical coding-standard pass (see
[Coding Standard Violations](coding_standard_violations.md#status)),
except the unused arguments:

- `run_diag_sr_so_r.f90:142`: the `npr=0` for `l==lloc` was redundant,
  because `oncvpsp_r.F90:539` already zeroes `nproj(lloc+1)`. Removed.
  Unlike `run_diag_r`, this routine doesn't stop on an `ldiracfb` error.
  `run_diag_r` runs first with identical calls, though.
- `vrel.f90`: the second of two consecutive `cycle` tests (`2.5*rc`)
  could never trigger. Removed.
- `lschpb.f90`: declared a local `real(dp) :: exp`. The name still
  resolved to the intrinsic. Removed.
- `sbf_basis_con.f90`: the unused local `tder(6)` is removed.
- Unused names:
  - locals `zz` in `lschkb`/`lschps`: removed;
  - arguments `nn` in `lschkb`/`lschps`, `icmod` and `rcfact` in
    `modcore`/`modcore2`, and `drdum`/`dqdum` in `eresid`: kept, because
    removing them would change the routines' calls for programs that use
    the library. Their comments now say they are not used.
- **Fixed:** `adams.inc` (formerly `aeo.f90`) wrote the Adams 1/24
  coefficient 4 times as `4.16666666667d-2`, a relative error of about
  8e-12. It is now `1.0d0/24.0d0`, fixed together with the introduction of
  `constants_m`.
  - This changes the last printed digits of most outputs. Quantities that
    are small differences move more: the largest relative changes are
    about 4e-3, in the 57_La excitation error (2.37e-6) and some 80_Hg_r
    convergence-profile entries.
  - Every reference comparison still passes, so the refs in `tests/refs`
    were kept as they are (they are the 4.0.1 outputs).
- **Fixed:** **`oncvpsp.F90`: `oncvpspnr.x` printed "scalar-relativistic version
  4.0.1"** in its banner. The banner chose its text with
  `#if (defined RELATIVISTIC) && (RELATIVISTIC == 0)`, but no build
  defines `RELATIVISTIC` (the non-relativistic build defines
  `ONCVPSP_NONREL`). It now uses `ONCVPSP_TREATMENT`, as the build line
  does. The 08_O_nr and 14_Si_nr refs, from upstream's non-relativistic
  executable, already say "non-relativistic version".
  - Found in the mechanical coding-standard pass (2026-10-19): the
    `#if` inside the continued `write` statement was the only code
    fortitude could not parse.

## P3 (fixed)

Fixed in batches 1–4 of the plan below. None of these changes the
reference outputs. All 28 reference tests pass with both the CMake and the
Make binaries. New error messages follow the format in
[Coding Standards](coding_standards.md#error-messages).

**Latent bugs that needed a decision**

- **`qroots.f90:92`**: `found_root=.true.` was set unconditionally, and
  `dlgd_last` took the value at the bisected root rather than at the grid
  point, which made a fake sign change on the next step. `dlgd_last` now
  takes the grid value, a failed bisection is skipped with
  `qroots: WARNING bisection failed, skipping`, and the final check counts
  the roots. No test prints the warning, and no output changes.
- **`ldiracfb`/`ldiracfs` series starts**: the sign of `ee` in the
  starting derivatives, and `fup(ii)=sc*gup(ii)` at `ldiracfb.f90:217`
  (fixed on `oncvpsp-opt`, `0b1c655`; FR outputs move by at most 4e-6).
- **`icmod=4` grid**: `modcore3` now takes the coarse-search grid from the
  input (text or TOML) in place of its hard-coded values, and scales the
  initial Nelder-Mead simplex to the grid steps. The default grid gives
  bit-identical output. `check_data` rejects a grid with a non-positive
  minimum or step, or a maximum below the minimum. The input echo
  (`write_icmod`) writes the grid when it is not the default, so a run
  from the echoed input uses the same grid. Tests `40_Zr_grid`,
  `toml_sr_40_Zr_grid` and `error_icmod4_bad_grid`.
- **`m_psmlout` `_tmp_input`**: gone; the PSML file is written into the
  output like the psp8 and upf files (`oncvpsp-opt`, `8b2a016`).
- **Text input from a pipe**: the text reader backspaced, so it failed
  with "Illegal seek" on a pipe. It now reads each line once
  (`oncvpsp-opt`, `17f5b9b`).
- **Exit status of errors**: the error `stop`s from 4.0.1 (66 of them) are
  now `error stop 1`, or `error stop '...'` where they had a stop code, so
  failures exit with status 1. The messages are unchanged, except that
  gfortran prints `ERROR STOP` instead of `STOP` before a stop code. The
  normal end of each driver still exits with 0.
- **Input echo**: an input echoed in the output or in a psp8, UPF or PSML
  file should run again to the same output. Checking every test input in
  all three formats found these problems:
  - No echo wrote `rxpsh`, so a rerun used the default radius. The echoes
    now write it when the input sets it (`input_echo_m`'s `epsh_lines`),
    which changes the two echoes in the 14_Si_lder ref.
  - With `icmod=4`, `modcore3` sets `fcfact=1`, which the psp8 header
    needs as `fchrg`. The psp8 and UPF echoes printed that value in place
    of the input's. They now echo the input value (`fcfact0`, like
    `rc0`). This changes one line of the psp8 echo in the 38_Sr_lxc,
    40_Zr, 40_Zr_grid and 74_W_r refs.
  - The PSML echo left out `rcfact`, the `icmod=4` grid, `rxpsh` and the
    log mesh. It wrote `iexc` as `i6`, which overflowed for Libxc codes
    (`******`), and the grid-adjusted `rc` to 2 decimals in place of the
    input `rc0`. It now uses the formats of the psp8 echo and shares the
    model core and log derivative sections with it (`input_echo_m`).

  The `echo` tests (`scripts/echo_rerun_test.sh`) check the output echo
  and the file echo for inputs that cover each of these. `08_O_wigner`
  can't be checked with UPF, because UPF output doesn't support `iexc=1`.

- **`rxpsh` inside a core radius**: `run_phsft` and `run_phsft_r`
  (lines 93 and 114) printed `ERROR: run_phsft: rxpsh=... <= rc(l)=...`
  for `rxpsh<=rc` and went on, evaluating the log derivatives inside rc;
  `run_phsft_r` called itself `run_phsft`, and the local-potential message
  printed `lloc+1` as the l. Found while documenting `rxpsh` (2026-10-19).
  - `check_data` now stops with
    `check_data: ERROR rxpsh not beyond every core radius` when a set
    `rxpsh` is not larger than every input rc (`rc(5)` too for `lloc=4`).
    Tests `error_rxpsh_inside` and `error_rxpsh_at_rc`.
  - `run_phsft(_r)` now give
    `run_phsft: WARNING rxpsh inside the core radius, l=..., rxpsh=..., rc=...`
    (or `inside the local potential core radius`) and go on, for callers
    of the library that skip `check_data`. The test is now `irphs<irc`: the
    mesh rc is the first point beyond the input rc, and the log derivative
    radius the first point at or beyond `rxpsh`, so with `rxpsh>rc` the
    programs never get the warning, while a radius between the input rc
    and the mesh rc (same mesh point) is allowed. Test `unit_run_phsft`.

**Build and dead code**

- **`src/CMakeLists.txt`**: removed the duplicate `fpovlp.f90` and
  `relatom.f90` entries.
- Deleted the uncalled `lschpse.f90` (which used `rr(ii)` in place of
  `vv(ii)`), `sbf_basis.f90` (shape and bounds errors) and `fpovlp_r.f90`
  (never compiled), and the commented-out `fpovlp_r` call in
  `oncvpsp_r.F90`.
- **`run_vkb_r.f90:78-83`**: removed the dead `l1==1` branch and the
  `l1=0` before it.
- **Legacy `src/Makefile`**: it didn't build. Checked on Linux with
  gfortran:
  - The committed `make.inc` held one developer's macOS paths, and its
    `FCCPPFLAGS = -ansi` makes GNU `cpp` reject the `//` in the URLs in
    the `.F90` licence headers. It is now `make.inc.example`, a blank
    template that explains each setting. `make.inc` itself is untracked and in
    `.gitignore`, and the Makefiles stop with a message if it is missing
    or unfilled.
  - `input_text_m.o` was not in any object list, so the drivers failed
    with "Cannot open module file 'input_text_m.mod'". It is now in
    `OBJS_ALL`, with explicit module dependencies so `make -j` works.
  - `oncvpspnr.x` was built from the stale `oncvpsp_nr.F90`, whose
    `run_phsft` call omitted `rxpsh`, shifting every later argument. It
    stopped with `ERROR: run_phsft: rxpsh= 0.00013 > rr(mmax)=****`.
    `oncvpsp_nr.F90` is now the `ONCVPSP_NONREL` wrapper of `oncvpsp.F90`,
    and the older `oncpsp_nr.f90` is deleted.
- **CMake libxc include directories**: a fresh configure against a system
  libxc (`/usr/include`, found via pkg-config) failed with "Cannot open
  module file 'xc_f03_lib_m.mod'". pkg-config drops system include
  directories, and `ONCVPSP_DEPENDENCY_INCLUDE_DIRS` was never applied to
  any target. The rewritten CMake build of `oncvpsp-opt` (`abf71da`)
  handles this, and replaced the fix on `master`.
- **`tests/CMakeLists.txt`**: the `nr` loop set its label on the `sr` test
  of the same name, and every loop used the property `LABEL` instead of
  `LABELS`, so `ctest -L` didn't work. Both fixed, and since replaced by
  the script-based tests of `oncvpsp-opt`.

**Errors that were silent**

- **Missing input file**: `oncvpsp.F90` and `oncvpsp_r.F90` checked the
  `open` status only after `read_input_text` had run, so `oncvpsp.x -i
  missing.dat` aborted with "Cannot read from file opened for WRITE". They
  now check right after `open` and stop with
  `oncvpsp: ERROR cannot open input file, file=...` and exit status 1. The
  old check after the read could no longer fire and was removed. Tests
  `error_{sr,nr,r}_missing_input` (label `errors`, using
  `scripts/error_test.sh`) check the message and the exit status.
- **`upfout.f90:327`, `upfout_r.f90:350`**: the unsupported-`iexc` message
  now uses `(t8,a,i0,a)`, so it prints instead of failing on the integer.
  The text is unchanged.
- **`run_plot_r.f90:178`**: `ldiracfs` errors now stop the run. `ldiracfs`
  only fails on bad `kap` or `zz`, and prints its own message first.
- **`run_ghosts.f90`**: `nbas` now starts at `nbmax`, so ghost testing
  isn't skipped when all `nbmax` barrier states are found.
- The messages added with the P2 fixes in `a672619` now follow the format
  and exit with status 1: `rxpsh` beyond the mesh in `run_phsft(_r)` and
  no model-core crossover in `modcore3`.

**Search results and index bounds**

- **`lschpsbar.f90:94`**: `imin=2` before the search.
- **`run_vkb(_r).f90`, `sr_so_r.f90`**: `ipk=1` before each peak search.
- **`upfout(_r).f90`**: `jj=max(jj,1)` after the norm search.
- **`psatom(_r).f90`**: `ierr` is set to 0 at entry and is now
  `intent(out)`.
- **`wellstate(_r).f90`**: `iumax=1` before the search.
- **`renorm_r.f90`**, **`modcore2.f90`**, **`eresid.f90`**: stop with a
  message if the wave function is zero (`nin`, `mmaxlin`) or the core
  charge never exceeds the valence pseudocharge (`rmatch`/`rhocmatch`).
- **`sbf_rc_der.f90`**: the guard is now `llin>4`.
- **`modcore.f90:142`**: stop if `ircc<23`.
- **`run_ghosts.f90`**: stop if `mbar+4>mmax`, which `lschpsbar` would
  write past.
- **`dpnint.f90`**: the `tt(1)==0` extrapolation now runs once after the
  loop (bit-identical result), and stops if `mm<npoly+2`.

**Argument mismatches**

- **`wellstate_r`**: dropped the unused `srel` argument that
  `oncvpsp_r.F90` didn't pass.
- **`debl`**: declared `debl(6)` in `linout_r`, `upfout_r` and
  `check_data`, matching the drivers.
- **`m_psmlout.f90`**: `f` in `write_psml_item` is no longer optional, and
  the `lmax = lmax` no-ops are removed so `lmax` is `intent(in)`.
- **`upfout(_r).f90`**: `nrl` is now a local copy of the new `intent(in)`
  argument `nrlin`, and the real→integer conversion is explicit.
- **`scripts/check_interfaces.sh`**: compiles each driver (all nine,
  including the ae/ps stages) with all external procedures as one file,
  so gfortran checks every call against its definition. It runs in CI
  after the build and reported exactly the three mismatches above before
  they were fixed.

## Next steps

Nothing is left open. The dead code in `functionals.F90` is kept on
purpose.

## Dismissed

Each of these was checked and turned out not to be a bug:

- `sbf_basis_con.f90` bounds: `check_data` enforces `3<=ncon<=5`, and
  `run_optimize` sets `ncon=ncon_in+iprj-1`.
- `fpovlp(_r).f90` `rr(101)`: `irc` is always far above 101 on the log
  mesh. `fpovlp_r` is never called.
- `optimize.f90` `nnull==0`: `check_data` guarantees `nnull>=2`.
- `linout_r.f90` `zero(20)`: `npr_sr` and `npr_so` are ≤10.
- `vrel.f90` edge reads: the rc window keeps `ii` in the interior.
- `run_ghosts.f90` `hmev` past `nbas`: `hmev` has size `nbmax` and is
  zeroed. `dsyev` with N=0 is a legal quick return.
- `qroots.f90` `nroot>=4`: `check_data` guarantees it.
- `input_toml_m.f90` icmod range: `check_data` rejects values outside 0–4
  for both readers.
- `linout(_r).f90` `ixc_abinit`: `check_data` rejects `iexc>4`.
- `sratom`/`relatom` `ierr=100` unchecked: the drivers already stop with
  "ERROR all-electron reference atom not converged" when `it>=100`, and
  every other failure inside `sratom`/`relatom` stops there.
- `wellstate_r.f90:154` commented-out `ldiracfb` check: the well search
  retries with a new well radius when `ldiracfb` fails, and the outer
  convergence test catches a bad `et`. The sibling `wellstate` does stop,
  but restoring the check could stop runs that now recover, so it is left
  as is.
- Leaks in `vrel`, `run_diag_sr_so_r`, `run_config(_r)` and `modcore`:
  these are local allocatables, which Fortran deallocates automatically on
  return.
