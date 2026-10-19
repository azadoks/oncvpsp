# Known issues

Suspected bugs found while auditing dummy-argument `intent`s across `src/`,
then triaged against the code, the callers, the input checks in
`check_data.f90`, and the test inputs in `tests/data/`. The P1 and P2 issues
have been fixed. The P3 issues are open. Line
numbers refer to the code before the fixes.

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
| P3 | everything else below | no | open |

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

- **`oncvpspnr.x` is not non-relativistic.** CMake builds it from
  `oncvpsp.F90` with `RELATIVISTIC=0` (`src/CMakeLists.txt:121-124`), but no
  source reads that macro. `oncvpsp.F90:114` hard-codes `srel=.true.`, so
  the binary is identical to `oncvpsp.x`. The binary is installed, and
  anyone running it silently gets scalar-relativistic pseudopotentials. No
  test covers it.
  - The old non-relativistic driver, `oncvpsp_nr.F90`, sets
    `srel=.false.`, but CMake doesn't build it and it is stale (see
    [Build](#build)).
  - Fixed: `oncvpsp.F90` now sets `srel=.false.` when built with
    `RELATIVISTIC=0`, and the banner, `--help` and `--version` text say
    "non-relativistic". The check is `defined RELATIVISTIC &&
    RELATIVISTIC == 0`, so a build without the macro stays
    scalar-relativistic. New tests `08_O_nr` and `14_Si_nr` run
    `oncvpspnr.x` against `tests/refs/*_nr.out`.

## P3 (open)

### Latent or unreachable in practice

- **`ldiracfb.f90:217`**: `fup(ii)=sc*gup(ii)` should be `sc*fup(ii)`.
  Each iteration recomputes everything, and no caller reads `up(:,2)`, so
  the typo has no effect.
- **`ldiracfb.f90:142-146,186-189`, `ldiracfs.f90:107,110`**: the series
  starts near r=0 use `+cci*(ee+vv)` and `(2c²-ee-vv)`. The
  predictor/corrector uses `-cci*(ee-vv)` and `(2c²+ee-vv)`, i.e. the sign
  of `ee` is flipped. In `ldiracfb` the inward start has the same
  problem.
  - The effect is negligible: near the origin |V|≫E, and the inward error
    decays.
  - It may cause last-digit changes in the `_r` refs; this hasn't been
    checked.
- **`qroots.f90:92`**: `found_root=.true.` is set unconditionally, so the
  "failed to find root" check is dead.
  - Do **not** just delete this line. `dlgd_last=dlgd` (line 99) stores
    the value at the bisected root, not at the grid point. That creates a
    fake sign change on every following step, and bisection then "fails"
    several times per `l` in every test. Deleting line 92 would make every
    test stop.
  - Fix: save the grid value before bisecting, and treat a failed bisection
    as a skip.
- **`upfout.f90:103`, `upfout_r.f90:111`**: the backward `jj` search could
  end at `jj=0` and read `rr(0)`. That needs a normalised orbital whose norm
  never exceeds 1e-6, which can't happen. Fix: `jj=max(jj,1)`.
- **`sbf_rc_der.f90`**: the guard is `llin>5`, but `llin=5` would write
  `sbfad(1,11)`. `check_data` caps `lmax<=3`. Fix: change the guard to
  `llin>4`.
- **`modcore.f90:142`**: `ii=ircc-25+3*jj` needs `ircc>=23`. `ircc` sits
  near rc, at an index above 1000.
- **`lschpsbar.f90`**: the barrier start writes `uu(mbar:mbar+4)` without
  checking it against `mmax`. That only matters for rc > ~15 bohr.
- **`lschpsbar.f90:94`**: `imin` is undefined if no point up to `mbar` lies
  below `emax`. It is then clamped to ≥2. Fix: set `imin=2` at entry.
- **`run_ghosts.f90`**: `nbas` is set only on the `lschpsbar` error path.
  If all `nbmax` solves succeed, `nbas` stays 0 and ghost detection is
  silently skipped. `nbmax` is sized so the error path normally fires. Fix:
  set `nbas=nbmax` first, and bound the `jj` loops by `nbas`.
- **`oncvpsp_r.F90:484`**: calls `wellstate_r` with 13 arguments, but the
  routine takes 14; the last is `srel`. This is harmless because `srel` is
  unused, and the compiler can't catch it because the routine has no
  explicit interface. Fix: drop `srel` from the signature.
- **`input_toml_m.f90`**: `fcfact_min/max/step` and `rcfact_min/max/step`
  are declared `intent(out)` but never set, which gfortran warns about. No
  code reads them: `modcore3` hard-codes the `icmod=4` grid at lines
  161–165. That also means the text reader's explicit-grid form of
  `icmod=4` is silently ignored. Fix: set them to defaults, or wire the
  grid into `modcore3`.
- **`m_psmlout.f90`**: `write_psml_item` declares `f` optional but uses it
  unconditionally. Every caller passes it. Fix: make `f` non-optional.
- **`dpnint.f90`**: the `tt(1)==0` extrapolation block sits inside the `jj`
  loop. That wastes work, but the result is bit-identical. With `mm<9` it
  reads `ss` beyond `mm`.

- **`sratom.f90:198`, `relatom.f90:209`**: these return `ierr=100` when
  the all-electron SCF doesn't converge, but `oncvpsp.F90:283` and
  `oncvpsp_r.F90:278` never check it. A warning is printed and the run then
  continues with the unconverged atom. Some other solver `ierr`s are also
  ignored: `ldiracfs` in `run_plot_r.f90:178` and `ldiracfb` in
  `wellstate_r.f90:154`, where the check has been commented out. In
  `wellstate_r` the outer convergence loop still catches a bad `et`. Fix:
  stop on `ierr/=0` after the AE solve.
- **`upfout.f90:327`, `upfout_r.f90:350`**: the unsupported-`iexc` error
  writes a string, an integer and a string with the format `(t8,a)`. On
  format reversion the integer `iexc` meets an `a` descriptor, so gfortran
  should abort with a runtime I/O error instead of printing the message.
  This hasn't been tried. It is reachable with an unsupported libxc
  functional and `psfile=upf`, and the run would stop either way. Fix:
  `(t8,a,i0,a)`.
- **Missing input file**: `oncvpsp.F90:171` opens with `iostat=ios` but
  only checks `ios` after `read_input_text` has run. Tested: `oncvpsp.x -i
  missing.dat` aborts with "Fortran runtime error: Cannot read from file
  opened for WRITE" and exit code 2. The TOML reader reports the error
  clearly. The same applies to `oncvpsp_r.F90`. Fix: check `ios` right
  after `open`.
- **Possibly used before being set.** Each case needs degenerate data:
  - `ipk` in `run_vkb.f90:298`/`run_vkb_r.f90:331`: needs an identically
    zero projector.
  - `nin` in `renorm_r.f90:61`: needs an identically zero wave function.
  - `iumax` in `wellstate.f90:78`/`wellstate_r.f90:81`: needs `|uu|` to be
    monotonic over the whole grid.
  - `rmatch`/`rhocmatch` in `modcore2.f90:89`: needs the core charge to
    never exceed the pseudo valence charge, e.g. `nc=0`.
  - `mmaxlin` in `eresid.f90:95`: needs `|uu|<=1e-10` everywhere.

  Fix: initialise each variable before its search.
- **`m_psmlout.f90:138,621`**: the PSML writer echoes the input through a
  fixed-name scratch file, `_tmp_input`, in the working directory. It
  doesn't check for I/O errors and never deletes the file. Two PSML runs in
  the same directory can race. Fix: use `status='scratch'` or an internal
  buffer.

### Dead code

- **`lschpse.f90:78`**: uses `rr(ii)+0.5*sls/rr(ii)**2` where every sibling
  routine uses `vv(ii)`. The routine is compiled but never called.
- **`run_vkb_r.f90:78-83`**: `l1=0` is set right before
  `else if (l1==1)`, so that branch is dead. For `lloc=0` the `else` branch
  computes the same `vloc`.
- **`sbf_basis.f90`**: has shape and bounds errors, but nothing calls it.
  Delete it.
- **`functionals.F90`**: in `functional_get_tau`, `n(1)`/`t(1)` are stale
  for `nspin=2`. The routine is never called.

### Build

- **`src/CMakeLists.txt`**: `fpovlp.f90` and `relatom.f90` are each listed
  twice (lines 56/71 and 60/75). CMake ignores the duplicates.
  `fpovlp_r.f90` is never compiled, and its only call site is commented
  out.
- **Legacy `src/Makefile`** (fixed): it didn't build. Checked on Linux
  with gfortran:
  - The committed `make.inc` held one developer's macOS paths, and its
    `FCCPPFLAGS = -ansi` makes GNU `cpp` reject the `//` in the URLs in
    the `.F90` licence headers. It is now `make.inc.example`, a blank
    template that explains each setting and recommends
    `-traditional-cpp -P`. `make.inc` itself is untracked and in
    `.gitignore`, and the Makefiles stop with a message if it is missing
    or unfilled.
  - `input_text_m.o` was not in any object list, so the drivers failed
    with "Cannot open module file 'input_text_m.mod'". It is now in
    `OBJS_ALL`, with explicit module dependencies so `make -j` works.
  - `oncvpspnr.x` was built from the stale `oncvpsp_nr.F90`, whose
    `run_phsft` call omitted `rxpsh`, shifting every later argument. It
    stopped with `ERROR: run_phsft: rxpsh= 0.00013 > rr(mmax)=****`. It is
    now built from `oncvpsp.F90` with `-DRELATIVISTIC=0`, as in CMake, and
    `oncvpsp_nr.F90` is deleted.
  - All 28 `tests/refs` comparisons pass with the Make binaries (libxc
    build), and the build without libxc (`exc_libxc_stub.o`) also links.
- **CMake libxc include directories** (fixed): a fresh configure against
  a system libxc (`/usr/include`, found via pkg-config) failed with
  "Cannot open module file 'xc_f03_lib_m.mod'". pkg-config drops system
  include directories, and `ONCVPSP_DEPENDENCY_INCLUDE_DIRS` was never
  applied to any target. `CMakeLists.txt` now falls back to pkg-config's
  `INCLUDEDIR`, and `src/CMakeLists.txt` passes each directory to
  `oncvpsp` as a Fortran `-I` option. They are passed as options rather
  than as include directories, which CMake would drop for `/usr/include`.
- **`src/oncpsp_nr.f90`**: deleted. It was not built by CMake or Make; it
  was a stale copy of `oncvpsp_nr.F90` from before the input module
  existed.

### Cosmetic

- `run_diag_sr_so_r.f90:142`: the `npr=0` for `l==lloc` is redundant,
  because `oncvpsp_r.F90:539` already zeroes `nproj(lloc+1)`. Unlike
  `run_diag_r`, this routine doesn't stop on an `ldiracfb` error.
  `run_diag_r` runs first with identical calls, though.
- `m_psmlout.f90:239,723`: the no-op `lmax = lmax` forces
  `intent(inout)`.
- `upfout.f90:112`, `upfout_r.f90:121`: `nrl` is enlarged in place, with an
  implicit real→integer conversion. PSML doesn't use `nrl`, so `both` mode
  is unaffected. Use a local copy instead.
- `vrel.f90`: the second of two consecutive `cycle` tests (`2.5*rc`) can
  never trigger.
- `lschpb.f90`: declares a local `real(dp) :: exp`. The name still resolves
  to the intrinsic.
- `aeo.f90`: `aii` lacks `implicit none`, although all its names are
  declared.
- `psatom.f90`, `psatom_r.f90`: `ierr` is not reset at entry. Callers
  always pass in 0.
- `sr_so_r.f90`: `ipk` is unset only if a whole column is NaN. Fix: set
  `ipk=1` before each search.
- `sbf_basis_con.f90`: the local `tder(6)` is unused.
- Unused names:
  - arguments `nn` in `lschkb`/`lschps`;
  - locals `zz` in `lschkb`/`lschps`;
  - arguments `icmod` and `rcfact` in `modcore`/`modcore2`;
  - local `srel` in `run_config_r`;
  - `drdum`/`dqdum` in `eresid`.

- `oncvpsp_r.F90:753,763` pass `debl(6)` to `linout_r`/`upfout_r`, which
  declare `debl(6,2)`. This works by sequence association, and only
  `debl(l1,1)` is read, which lies within the actual array.
  `check_data.f90:37` declares `debl(5)` and reads at most `lmax+1<=4`.
  Fix: declare `debl(6)` in both routines.
- `aeo.f90`: the Adams 1/24 coefficient is written 4 times as
  `4.16666666667d-2`, a relative error of about 8e-12. Writing it as
  `1.0d0/24.0d0` would change the last digits of the refs.

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
- Leaks in `vrel`, `run_diag_sr_so_r`, `run_config(_r)` and `modcore`:
  these are local allocatables, which Fortran deallocates automatically on
  return.
