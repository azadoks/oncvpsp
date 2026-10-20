# Fortran best practices

Language-level deviations from modern Fortran practice (F2008/F2018) in the
82 files that `src/CMakeLists.txt` compiles (including the optional
Libxc/TOML files and the two main programs). The vendored
`src/xmlf90-wxml/` is excluded. Some have since been fixed: see the
[remediation order](#recommended-remediation-order) and the
[Coding Standard Violations](coding_standard_violations.md). Line numbers
refer to commit `7e9480a`, before any Known Issues fixes. File names are
relative to `src/`. The `foo`/`foo_r` pairs and the two drivers it cites
have since been merged into single routines taking `irel` or `nkap` (see
[Merging sr and fr](fr_merge_plan.md)), so its `_r` files and
`oncvpsp_r.F90` lines are as they were then.

Related pages:

- [Known Issues](known_issues.md) lists suspected bugs. Its items are not
  repeated here, only cross-referenced. Those since fixed are no longer
  listed there; their fixes are in the commit messages.
- [Coding Standard Violations](coding_standard_violations.md) covers layout,
  naming, comments and consistency (`.eq.`, `character*N`, `enddo`,
  commented-out code, header blocks). Those items are only summarised here.

Items marked **(new)** are possible bugs that are not in Known Issues.

## How the counts were made

Counts come from scripts that strip comments and strings, plus three
gfortran 11.4 passes over all compiled sources, concatenated into one
translation unit per executable (`oncvpsp.F90` + library and
`oncvpsp_r.F90` + library). gfortran then checks every call against the
callee's definition, as an explicit interface would:

1. `-fsyntax-only -Wall -Wextra -Wconversion -Wimplicit-interface
   -Wcompare-reals -Wintrinsic-shadow -Wcharacter-truncation` (and
   `-Wconversion-extra` separately);
2. `-O2 -Wmaybe-uninitialized`;
3. `-std=f2018`.

The modules and `exc_off_pbe.f` were compiled separately with the same flags.
The diagnostics in both builds were merged, and each `file:line` was counted
once.

## Summary

Ordered by risk: correctness, then portability, then maintainability, then
style.

| # | Category | Scope | Risk |
|---|---|---|---|
| 1 | External procedures, no explicit interfaces | 87 of 118 procedures; 323 implicit-interface call sites; 3 argument mismatches found | correctness |
| 2 | Build and preprocessor | `RELATIVISTIC` macro unused; `oncvpspnr.x` is scalar-relativistic | correctness |
| 3 | Error handling | 98 `stop`s (54 in library routines, most exit with status 0); 14 `ierr` results ignored; unchecked I/O | correctness |
| 4 | Variables set on only some paths | 19 `-Wmaybe-uninitialized` sites in 14 files | correctness |
| 5 | Non-standard declaration order | 20 files use an array bound before declaring it | portability |
| 6 | Implicit typing | 8 external procedures, 4 module procedures | correctness / portability |
| 7 | Hard-coded limits and magic numbers | 219 literal-shape dummy arrays, 110 fixed-size locals, `c=137.036` in 7 places | correctness / maintainability |
| 8 | Floating-point equality | 31 comparisons in 20 files | correctness (mostly benign) |
| 9 | Implicit real→integer conversion | 13 assignments in 11 files | correctness (mostly intended) |
| 10 | Explicit-shape and assumed-size dummies | 616 array dummies, 4 assumed-shape | maintainability |
| 11 | Kind handling | 4 `dp` schemes; 189 lines with `d`-specific intrinsics | maintainability / portability |
| 12 | Allocation, pointers, `save` | 175 `allocate`s, 0 with `stat=`; 5 aliasing pointers | maintainability |
| 13 | Modules, `use`, global namespace | 6 modules; 13 `use` without `only` | maintainability |
| 14 | Code duplication | 13 `foo`/`foo_r` pairs at 64–87% identical lines | maintainability |
| 15 | Dead code | 8 uncalled procedures, 4 uncompiled files | maintainability |
| 16 | Unused names and shadowing | 264 unused locals, 39 unused/unset dummies, 16 locals named like intrinsics | style |
| 17 | Fixed form and legacy constructs | 1 fixed-form file; no `GOTO`, `COMMON`, `EQUIVALENCE` or `ENTRY` | style |
| 18 | Hard-coded I/O units | 834 `write(6,…)` | style |

Missing `intent`s: none remain. All 1458 dummy arguments, in all 119
procedure definitions (counting `exc_libxc` and its stub separately), have an
`intent`, including those in `exc_off_pbe.f` (declared
with standalone `intent(in) ::` statements).

## 1. External procedures without explicit interfaces

**Why it matters.** An external procedure (one not in a module and not
`contains`ed) has an implicit interface at every call site. The compiler
cannot check the number, type, kind or rank of the arguments, or `intent` or
`optional`. A mismatch compiles and links silently, then corrupts memory or
gives wrong results at run time. Moving procedures into modules makes these
mismatches compile errors.

**Scope.**

- **Program units.** 87 external procedures plus 2 main programs. There are 7
  modules (`m_psmlout`, `m_libxc_list`, `parser`, `xc_messages`,
  `functionals_m`, `input_text_m`, `input_toml_m`), holding 29 module
  procedures and 2 internal ones. `exc_libxc` and `exc_libxc_stub` are
  counted once.
- **Lines.** The module files are about 3,560 of 20,650 compiled lines
  (17%). Everything numerical is external: the radial solvers, the
  optimizer, the XC drivers and the `run_*` drivers.
- **Call sites.** 323 implicit-interface call sites, to 82 distinct
  procedures. The most called are `dpnint` (30), `aeo`/`aio` (18 each),
  `vout` (16), `gg1cc` (16), `vpinteg` and `der2exc` (15 each), and the LAPACK
  routines `dsyev` (13), `dgesv` (2) and `dgesvd` (1).
- **Function types.** The external functions `aeo`, `aio`, `aei`, `aii`,
  `gg1cc`, `gp1cc`, `gpp1cc`, `sbf8` and `tfapot` are typed by hand in each
  caller (e.g. `lschfb.f90:57`
  `real(dp) :: aei,aeo,aii,aio,als !functions in aeo.f90`).
- **Global names.** All 87 names are in the global linker namespace. Generic
  names such as `derivs` (`exc_libxc.f90:118`), `vout`, `optimize`, `pspot`
  and `qroots` could collide with a host code if ONCVPSP is ever linked as a
  library.

**Mismatches found** by compiling the whole program as one unit:

- **`oncvpsp_r.F90:484`**: `wellstate_r` is called with 13 of 14 arguments.
  This is already in [Known Issues](known_issues.md).
- **(new) `oncvpsp_r.F90:753,763`**: `debl` is declared `debl(6)` in the
  caller (`oncvpsp_r.F90:72`). It is passed to `linout_r` and `upfout_r`,
  which declare it `debl(6,2)` (`linout_r.f90:80`, `upfout_r.f90:80`).
  gfortran rejects this with "Actual argument contains too few elements for
  dummy argument 'debl' (6/12)" when it can see both sides. The callees only
  read `debl(l1,1)`, so it works through sequence association, but the
  declared shape is wrong. `check_data.f90:37` declares the same array as
  `debl(5)`.
- **(new) `oncvpsp_nr.F90:620`**: `run_phsft` is called without the `rxpsh`
  argument that was added in #34 (`run_phsft.f90:19`), so every later
  argument is shifted by one. CMake does not compile this file, but the legacy
  `src/Makefile` does (`OBJS_NR`). See [section 2](#2-build-and-preprocessor).

Other shape inconsistencies are legal but fragile. `nproj` is `nproj(6)` in
the callers (`oncvpsp.F90:56`) and `nproj(5)` in `psatom.f90:50`,
`psatom_r.f90:49`, `run_config.f90:56` and `run_config_r.f90:56`.

About 60 call statements pass an array element as the start of a column,
e.g. `call dpnint(rr,vpuns(1,l1),mmax,rl,vpl(1,l1),nrl)`
(`linout_r.f90:120`). This relies on sequence association, which only works
for explicit-shape or assumed-size dummies. It has to become a section
(`vpuns(:,l1)`) when the dummies become assumed-shape.

## 2. Build and preprocessor

**Why it matters.** Preprocessor macros and build lists decide what is
actually compiled. A macro that the source never tests, or a stale file that
one build system still uses, silently produces a different program from the
one intended.

- **`RELATIVISTIC` is never tested (new).** `src/CMakeLists.txt:109,124,139`
  defines `RELATIVISTIC=1/0/2` for the three executables. No source file
  references it.
  - `oncvpspnr.x` is built from `oncvpsp.F90`, which hard-codes
    `srel=.true.` (`oncvpsp.F90:114`). So the "non-relativistic" executable
    is scalar-relativistic.
  - The real non-relativistic driver, `oncvpsp_nr.F90`, differs only in
    `srel=.false.` and the banner text. It is not compiled, and it has the
    stale `run_phsft` call [above](#1-external-procedures-without-explicit-interfaces).
  - Coding Standard Violations notes the same thing. The fix is a
    `#if RELATIVISTIC == 0` around `srel`, and then deleting
    `oncvpsp_nr.F90`.
- **CMake source list.** `fpovlp.f90` and `relatom.f90` are listed twice
  (`src/CMakeLists.txt:56,71` and `:60,75`), and `fpovlp_r.f90` is never
  compiled. The fully-relativistic block probably meant `fpovlp_r.f90`, but
  `oncvpsp_r.F90:522` has the `fpovlp_r` call commented out in favour of
  `fpovlp`.
- **Legacy `src/Makefile`.**
  - It still builds `oncvpspnr.x` from `oncvpsp_nr.o`.
  - It does not list `input_text_m.o` or `input_toml_m.o`, so it can no
    longer link the current drivers.
  - Either remove it or regenerate it from the CMake list.
- **Preprocessor use elsewhere.**
  - `WITH_TOML` (3 blocks each in `oncvpsp.F90`/`oncvpsp_r.F90`) is used
    cleanly.
  - `functionals.F90` has 80 directive lines, nearly all
    `#if XC_MAJOR_VERSION<=4 … #else … #endif`, to support Libxc 3–7. This is
    the usual pattern, but 22 of the 26 `#if` blocks are inside procedures; a thin
    version-specific wrapper module would localise them.
  - `m_libxc_list.F90` has a `.F90` extension but no directives.
- **Compiler flags.** Release uses `-O3 -ffast-math`
  (`CMakeLists.txt:37`). `-ffast-math` allows reassociation and assumes no
  NaN/Inf. That interacts badly with the exact-equality convergence tests in
  [section 8](#8-floating-point-equality). Debug has `-fcheck=all`, but no
  build type turns on `-Wimplicit-interface` or `-std=`.

## 3. Error handling

**Why it matters.** A `stop` inside a library routine ends the process with
no chance for the caller to clean up, retry or report context. A plain `stop`
exits with status 0, so scripts and the CTest harness see a failure as
success. Error codes that are returned and then ignored hide failures.

### `stop` in library routines

- **Count.** 98 `stop` statements in 36 files: 66 bare, 25 with an integer
  code, 7 with a string. There are 2 `error stop`s (`oncvpsp.F90:150`,
  `oncvpsp_r.F90:151`).
- **Library routines.** 54 of these are outside the main programs and the
  input readers, in 34 files. For example: `dpnint.f90:48,60,64`,
  `vout.f90:123,140,172`, `run_vkb.f90:145,199,240`,
  `sr_so_r.f90:123,179,187`, `wellstate.f90:152,202`, `modcore.f90:93,204`,
  `relatom.f90:142`, `sratom.f90:129`, `qroots.f90:95,103` and
  `functionals.F90:60,396,640,646`.
- **Exit status.** All 54 are bare (`stop` or `stop 'text'`), so they exit
  with status 0. Only `input_toml_m.f90` and the two main programs use
  `stop 1`. Checked with gfortran: a bare `stop` returns 0.
- **`m_psmlout.f90:1138,1161`**: `add_zero_r` and `get_unit` stop with a
  message but no context.

### Ignored error codes

53 calls pass `ierr`. In 14 of them `ierr` is not tested in the following
statements:

- **(new) `oncvpsp.F90:283`, `oncvpsp_r.F90:278`**: `sratom`/`relatom`
  return `ierr=100` when the all-electron SCF does not converge
  (`sratom.f90:197-199`). The drivers never test it. Only a "WARNING failed
  to converge" line is printed, and the run continues with an unconverged
  reference atom.
- `fphsft.f90:69`, `fphsft_r.f90:67`, `run_diag.f90:98`,
  `run_diag_r.f90:109`, `run_diag_sr_so_r.f90:133`, `run_plot.f90:153,161`,
  `run_plot_r.f90:170,178`, `vkboutwf.f90:60,83`, `wellstate_r.f90:154`:
  results of the scattering-state solvers (`lschfs`, `ldiracfs`, `lschps`,
  `lschkb`) and of `lschvkbb`/`ldiracfb` are used without a check.
  - `lschfs`, `lschps` and `lschkb` always return `ierr=0`.
  - `ldiracfs` can return 2, 3 or 60, and `lschvkbb`/`ldiracfb` can fail, so
    `run_plot.f90:153`, `run_plot_r.f90:170,178`, `fphsft_r.f90:67`,
    `run_diag_r.f90:109`, `run_diag_sr_so_r.f90:133` and
    `wellstate_r.f90:154` can use the output of a failed solve. At
    `wellstate_r.f90:156-161` the check is commented out.
- Already in Known Issues: `psatom.f90:216` (check can't fire),
  `psatom_r.f90:184` (masked).

### Unchecked I/O

- **`oncvpsp.F90:171,185`, `oncvpsp_r.F90:171,184`**: `open(…, iostat=ios)`.
  `ios` is not tested until after `read_input_text`/`read_input_toml` has
  already read from the unit (`oncvpsp.F90:203`). A missing input file
  therefore shows up as a misleading "Read error, input data file line 1"
  from `read_error`, or as a toml-f load error.
- **`input_text_m.f90:477`**: `cmtskp` reads with no `iostat`, so an
  end-of-file while skipping comments crashes with a runtime error rather than
  a message. The dead copy at `oncvpsp_r.F90:798` is the same.
- **`m_psmlout.f90:138,206,621,690`**:
  - The input is echoed to a scratch file `_tmp_input` in the working
    directory and then read back.
  - `open` and `close` have no `iostat`, and the file is never deleted.
  - Two runs in the same directory overwrite each other's file.
  - `get_unit` (`m_psmlout.f90:1147`) searches units 10–99 by hand.
  - Fix: use `open(newunit=…)` and `status='scratch'`, or build the string in
    memory.
- **(new) `upfout.f90:327-328`, `upfout_r.f90:350-351`**: the "unsupported
  iexc" message writes `'…ERROR iexc = ',iexc,' is …'` with format `(t8,a)`.
  Format reversion applies the `a` descriptor to the integer `iexc`, so
  gfortran prints garbage bytes where the value should be (checked). Use
  `(t8,a,i0,a)`. `check_data` normally rejects such an `iexc` first.
- **Positive.** `input_text_m.f90` checks `iostat` on 36 of its 37 reads.

## 4. Variables set on only some paths

**Why it matters.** A variable set only inside an `if` or a search loop is
undefined when the condition never holds. The result is then compiler- and
optimisation-dependent garbage, and `-fcheck=all` does not catch it.

`-O2 -Wmaybe-uninitialized` reports 19 sites in 14 files. The ones not
already in Known Issues:

- `run_vkb.f90:298-301`, `run_vkb_r.f90:334`: `ipk` (the peak index) is set
  only when some `|vkb|>0`. This is the same pattern as `sr_so_r.f90` in
  Known Issues.
- `renorm_r.f90:61,70`: `nin` is set only if some `uu(ii,1)≠0`.
- `wellstate.f90:78,84`, `wellstate_r.f90`: `iumax` is set only if `|uu|`
  ever decreases going inward.
- `modcore2.f90:89,192`: `rmatch` and `rhocmatch` are set only if
  `rhoc>rhotps` somewhere. When they are not, `ircc=0` is caught earlier,
  but the compiler cannot see that.
- `modcore2.f90:251`, `modcore3.f90:483`: `iint` is set inside the
  smoothing-mesh search.
- `eresid.f90:95,101`: `mmaxlin` is set inside a search loop and then used
  as an allocation size.
- `lschfb.f90:164`, `lschfs.f90:124`, `fpovlp.f90:67`: `gamma` is set by
  `if(ll==0)` / `if(ll>0)`. This is a false positive, because `ll>=0`, but
  `if/else` would make it provable.

Already in Known Issues: `lschpsbar.f90:94` (`imin`), `modcore3.f90:75`
(`ircross`), `vkbphsft.f90:158` (`dnpi`), `qroots`/`eresid` (`found_root`).

The usual fix is to initialise the variable to a sentinel before the search
and test the sentinel afterwards.

## 5. Non-standard declaration order

**Why it matters.** The standard requires a variable that appears in a
specification expression, such as an array bound, to be declared earlier in
the same scope or to be implicitly typed. With `implicit none`, a bound used
before it is declared is a GNU extension. gfortran accepts it by default, but
`-std=f2008` rejects it, and other compilers may too.

- **Scope.** 22 declarations in 20 files: `der2exc`, `exc_libxc`, `excwig`,
  `fpovlp`, `ldiracfb`, `ldiracfs`, `lschfb`, `lschpb`, `lschpsbar`,
  `lschvkbb`, `lschvkbbe`, `renorm_r`, `upfout`, `upfout_r`, `vkboutwf`,
  `vpinteg`, `vrel`, `wellstate`, `wellstate_r` and `wf_rc_der`.
- **Examples.**
  - `lschfb.f90:41` declares `rr(mmax),vv(mmax)` and `:44` declares
    `mmax`.
  - `upfout.f90:71` declares `uupsa(mmax,nv)`, and `nv` is declared at
    line 80.
- **Cause.** The convention of grouping declarations by Input/Output/Local
  section sometimes puts an array before its size argument. This predates the
  `intent` annotation. Declaring size arguments first in each section fixes
  it.

## 6. Implicit typing

**Why it matters.** Without `implicit none`, a misspelt name becomes a new,
implicitly typed variable. Names starting with I–N become integers, so a
misspelt real such as `ekin` → `ekn` silently truncates.

- **`exc_off_pbe.f`**: 7 routines (`exchpbe:201`, `corpbe:274`,
  `gcor2:417`, `exchpw91:436`, `corlsd:481`, `gcor:516`, `corpw91:536`) use
  `IMPLICIT REAL*8 (A-H,O-Z)` (e.g. line 232).
  - The integers `lgga` and `lpot` are integer only because of their first
    letter.
  - `excggc:24` in the same file has `implicit none`.
  - This is imported Perdew/Burke code. Converting it is low-risk but has to
    be checked against the PBE reference outputs.
- **`aeo.f90:62`**: `aii` has no `implicit none`. All its names happen to be
  declared. This is in Known Issues.
- **`m_psmlout.f90:23`**: the module has no `implicit none` at module scope.
  - `psmlout` and `psmlout_r` (and their internal `do_configuration`) have
    their own.
  - The module procedures `my_add_attribute:1118`, `add_zero_r:1128`,
    `get_unit:1147` and `write_psml_item:1164` do not, so they rely on
    implicit typing.
  - One `implicit none` after `module m_psmlout` covers all of them.
- **`m_libxc_list.F90:1`** also has no module-level `implicit none`. It
  contains only types and parameters, so this is harmless.
- **Positive.** All other 79 external procedures and the 6 other modules have
  `implicit none`.

## 7. Hard-coded limits and magic numbers

**Why it matters.** Fixed dimensions that are repeated by hand in every
routine must all change together. A mismatch, such as `debl(5)`/`(6)`/`(6,2)`
or `nproj(5)`/`(6)`, compiles silently (see
[section 1](#1-external-procedures-without-explicit-interfaces)). Physical
constants that are copied around drift apart.

### Array limits

- **Dummy arrays.** 219 of 616 dummy arrays, in 32 files, have literal
  shapes:
  - `30` is the maximum number of states (`fa(30)`, `na(30)`, `ea(30,2)`).
  - `6` is the maximum number of l channels + 1 (`rc(6)`, `ep(6)`, `nproj(6)`,
    `debl(6)`).
  - `5` is the maximum number of test configurations (`facnf(30,5)`,
    `nvcnf(5)`).
  - `4` is the maximum l + 1 (`vkb(mmax,mxprj,4)`).
  - These are spread over `check_data.f90:37-41`, `input_text_m.f90:36-109`,
    `linout*.f90`, `upfout*.f90`, `psatom*.f90` and `run_config*.f90`.
- **Fixed-size locals.** 110 local arrays in 28 files. Notable ones:
  - `const_basis.f90:51`: `usvd(10,10)`, `ss(10)`, `work(100)` and
    `amat(10,nbas)`. `dgesvd` is called with `lwork=64` (`:67-68`), although
    `work` has 100 elements.
  - `run_ghosts.f90:75-77`: `lpgh(20)`, `epgh(20)` and friends, guarded by
    `npgh<20` (`:265`).
  - `modcore2.f90:70`, `modcore3.f90:72`: `rint(20)`, `fint(20)`. Only 8
    elements are used.
  - `linout_r.f90:89` `zero(20)`, `gnu_script*.f90:44` `gln(100)`,
    `sbf_rc_der.f90:41` `sbfad(5,10)`.
- **Program parameters.** `mxprj=5` (`oncvpsp.F90:243`,
  `oncvpsp_r.F90:248`) is assigned at run time, not declared as a
  `parameter`.
- **Fix.** Put these limits in one constants module (`integer, parameter ::
  nstate_max=30, nl_max=6, ncnf_max=5, nproj_max=5`) and use them for every
  declaration.

### Numerical constants

- **Speed of light.** `137.036d0` appears in 7 places: `ldiracfb.f90:64`,
  `ldiracfs.f90:65`, `renorm_r.f90:54`, `lschfb.f90:85`, `lschfs.f90:74`,
  `fpovlp.f90:57` and `vrel.f90:82`. Some write it as `cc`, some as
  `fss=(1/137.036)**2`. CODATA gives 137.035999…. One named constant would
  keep these consistent.
- **Truncated rational coefficient.** The Adams–Bashforth/Moulton
  coefficient 1/24 is written `4.16666666667d-2` (`aeo.f90:31,44,57,69`).
  That is accurate to about 8e-12 relative, not to double precision, and every
  outward/inward integration uses it. Write `1.0_dp/24.0_dp`. Also
  `tfapot.f90:45`: `**.33333333d0` for 1/3.
- **`pi`.** Defined locally 25 times in four ways: 19 times as a 40-digit
  literal, plus `4*datan(1)` (`exc_off_pbe.f:55`), `real(8)` and
  `parameter(pi=…)`.
- **Radial mesh.** `amesh=1.006d0`, `rr1=0.0005d0/zz` and the 45 bohr outer
  radius are literals in `oncvpsp.F90:239-248`.
  - The mesh step is recomputed from `r(101)` in `exc_off_pbe.f:51`
    (`al = 0.01d0*dlog(r(101)/r(1))`), and the same expression appears in
    `lschfb.f90:67` and its siblings.
  - These assume the log mesh and an `mmax>101`.
- **Version strings.** The PSML `creator` attribute says `ONCVPSP-2.1.2+psml`
  (`m_psmlout.f90:200,443,684,979`), while the program is 4.0.1. The UPF
  header says `06/20/2107` (`upfout.f90:176`, `upfout_r.f90:205`).

## 8. Floating-point equality

**Why it matters.** `x == 0.0d0` and `x /= y` on computed reals depend on
rounding, optimisation level and flags such as `-ffast-math`. They are safe
only for values that were assigned exactly, such as an input sentinel.

`-Wcompare-reals` reports 31 comparisons in 20 files.

- **Input sentinels (benign).** Examples:
  - `oncvpsp.F90:299` (`rc(l1)==0.0d0` means "use the default");
  - `psatom.f90:150`, `psatom_r.f90:151`, `sratom.f90:118`,
    `run_config.f90:180,277`, `run_config_r.f90:181` (zero occupation);
  - `linout_r.f90:194` (`soscale==1.d0`);
  - `run_plot.f90:202`, `run_plot_r.f90:230,237` (`cvgplt` slot unused).
- **Convergence by exact equality.** `sratom.f90:133`, `relatom.f90:146`,
  `psatom.f90:170` and `psatom_r.f90:188` set `convg=.false.` if
  `ea(ii)/=et`. The comment at `sratom.f90:132` says this relies on
  `lschfb` returning a bit-identical eigenvalue once its own `eps` is met. It
  works, but it is fragile under `-ffast-math` or a change to the solver. Use
  `abs(ea(ii)-et) > tol`. (`sratom` and `relatom` now also require at least
  two iterations and a screening-potential residual below 1e-8 Ha: the
  eigenvalues alone stood still at iteration 1 for H on a mesh starting at
  1e-8 bohr, 0.29 Ha from self-consistency. `psatom` is unchanged.)
- **Computed values.** These are the ones most worth rewriting:
  - `oncvpsp.F90:585`, `oncvpsp_r.F90:715`: `rho(ii)==0.0d0` finds the
    density tail.
  - `vrel.f90:128`: `vr(ii)==0.0d0` patches singular points.
  - `modcore.f90:234`: `a0max==0.0d0`.
  - `lschfb.f90:106`: `ee==0.0d0`.
  - `lschpsbar.f90:84`: `emin==0.0d0`.
  - `optimize.f90:220`: `qout(kk)==0.0d0`.
  - `vout.f90:159`: `rhoc(1) .ne. 0.0d0`.
  - `dpnint.f90:101,119,126`.
  - `functionals.F90:412,679`.
  - `run_ghosts.f90:141`: `emax/=ebcut(l1)`.

## 9. Implicit real→integer conversion

**Why it matters.** Assigning a real expression to an integer truncates
toward zero with no visible `int()`. A reader cannot tell whether truncation
or rounding was intended, and a negative argument rounds the wrong way for
`floor` semantics.

`-Wconversion` reports 13 sites in 11 files. All look like intended mesh-index
arithmetic:

- `oncvpsp.F90:248`, `oncvpsp_r.F90:239`: `mmax=dlog(45.0d0/rr1)/al`.
- `eresid.f90:90,95,142`: `irclin=rc/dr + 1.5d0`, i.e. rounding by adding
  0.5.
- `ldiracfb.f90:178`, `lschfb.f90:192`, `lschpb.f90:141`,
  `lschpsbar.f90:169`, `lschvkbb.f90:152`: `nin=mch+2.3d0/al`.
- `run_optimize.f90:161`: `nqout=2 + 0.75d0*qmax/dqout`.
- `upfout.f90:112`, `upfout_r.f90:121`: `nrl` (in Known Issues).

Write `int(…)`, `nint(…)` or `floor(…)` explicitly.

Integer division: `-Winteger-division` flags only `dpnint.f90:89-93`
(`npoly/2`, intended).

## 10. Explicit-shape and assumed-size dummy arrays

**Why it matters.**

- An explicit-shape dummy (`a(n)`) needs its size passed separately. The
  compiler trusts the caller to pass a consistent size.
- An assumed-size dummy (`a(*)`) has no size at all, so it cannot be
  bounds-checked, even with `-fcheck=bounds`, and cannot be used in
  whole-array expressions.
- Assumed-shape dummies (`a(:)`) carry their bounds, can be checked, and
  remove the size arguments. They need an explicit interface
  ([section 1](#1-external-procedures-without-explicit-interfaces)).

**Scope.** Of 616 dummy arrays:

- 278 have variable bounds, e.g. `rr(mmax)`;
- 219 have literal bounds ([section 7](#7-hard-coded-limits-and-magic-numbers));
- 119 mix the two, e.g. `vkb(mmax,mxprj,4)`;
- 14 are assumed-size;
- 4 are assumed-shape, all in `m_psmlout.f90` (`add_zero_r`,
  `write_psml_item`).

**Assumed-size (all 14):**

- `aeo.f90:28,41,54,66` (`yy(*)` in `aeo`/`aio`/`aei`/`aii`);
- `dpnint.f90:32,36` (`xx`, `yy`, `tt`, `ss`);
- `exc_off_pbe.f:32-33` (`rho`, `r`, `vxc`, `exc`);
- `vkbphsft.f90:50` (`vkb(mmax,*)`, `evkb(*)`).

**Redundant size arguments.** 62 procedures take `mmax` only to declare
`rr(mmax)` and friends. `mxprj`, `nv` and `nc` are passed the same way. With
assumed shape, `size(rr)` replaces them. Where a routine deliberately works on
a prefix, as in `call vpinteg(…,ircut,…)` (`der2exc.f90:91`), the caller
would pass `rr(1:ircut)`.

## 11. Kind handling

**Why it matters.** A single, centrally defined working-precision kind makes
the code's precision explicit and changeable. Non-standard or obsolescent
spellings (`real*8`, `real(8)`, `dble`, `dsqrt`) tie the code to compilers
that use kind 8 for double precision, and hide the kind at the call site.

- **`dp` defined four ways:**
  - `integer, parameter :: dp=kind(1.0d0)`, re-declared locally 81 times
    (in every external procedure);
  - `selected_real_kind(10,100)` in `m_psmlout.f90:28`;
  - `iso_fortran_env` `real64` in `input_text_m.f90:2` and
    `input_toml_m.f90:2`;
  - none at all in `functionals.F90` (30 `real(8)` lines, e.g. 29, 107,
    369-392) and `exc_off_pbe.f` (`real*8` at 32-35, plus `IMPLICIT REAL*8`).
  - These all give the same kind on every current compiler, but nothing
    enforces it.
  - `double precision` is used in `m_psmlout.f90:1131-1135`.
- **Specific-name intrinsics.** These are obsolescent and some are
  non-standard:
  - `dlog` (59), `dsqrt` (47), `dabs` (38), `dmin1` (22), `dexp` (20) and
    `dmax1` (17), plus `datan` (1): 189 lines in 47 files. Worst:
    `exc_off_pbe.f` (20), `lschfb.f90` (11), `lschpsbar.f90` (9) and
    `run_plot*.f90` (9 each).
  - `dfloat` (non-standard GNU extension): 14 lines in 13 files, e.g.
    `lschpb.f90:177` `sn=ro**(2*ll+3)/dfloat(2*ll+3)`.
  - `idint`: 10 lines in 3 files (`run_ghosts.f90:108,117`,
    `run_plot.f90:79-86`, `run_plot_r.f90`).
  - `dble()`: 32 lines in 12 files, e.g. `linout.f90:97`,
    `oncvpsp_r.F90:407`.
  - Replace these with the generic `log`, `sqrt`, `abs`, `min`, `max`, `exp`,
    `real(…,dp)` and `int`.
- **Literal constants.** This is in good shape. 1557 lines write literals with
  a `d` exponent (`1.0d0`), and only 15 default-real literals exist, in 7
  files.
  - Most are exact (`0.0`, `0.5`, `2.`), e.g. `m_psmlout.f90:923`
    `jval = ll(i) + (3-2*jj)*0.5`, `exc_off_pbe.f:585,601` `2.*C3*RS`, and
    `run_phsft.f90:76,91` `rxpsh < 0.0`.
  - The one inexact case is `functionals.F90:419` `abs(a-b) < 1.0e-8` (a
    tolerance, harmless).
  - `linout.f90:164` and `linout_r.f90:187` write a default-real `0.0` into
    the psp8 header.
  - The `1.0d0` → `1.0_dp` change is a style decision; see Coding Standard
    Violations.
- **Character kinds.** `functionals.F90:332,343` assign a 1024-character
  Libxc string into `character(len=1000)` (`-Wcharacter-truncation`).

## 12. Allocation, pointers and `save`

**Why it matters.** Allocatables are freed automatically and cannot alias;
pointers can dangle, leak and alias. `allocate` without `stat=` aborts with a
compiler-specific message. Implicit `save`, where an initialised local keeps
its value between calls, is a common source of surprises.

- **Allocation.** 175 `allocate` statements in 59 files, 0 with `stat=`.
  The arrays are small (`mmax`≈2000), so this is low risk. A wrapper or
  `stat=`/`errmsg=` on the large ones in the main programs would still give a
  clear message.
- **Pointers.** `m_psmlout.f90` has 9 aliasing pointers into `target`
  dummies: `r`, `chval`, `chcore`, `vps` and `vlocal` in `psmlout`
  (`:106-107`), and the same names except `vps` (which is allocatable there)
  in `psmlout_r` (`:589-590`).
  - They exist only to give columns short names. `associate` or plain
    sections would do the same without the `target` attribute.
  - The `deallocate(chcore)` at `:434,969` is in Known Issues.
  - The `toml_table`/`toml_array` pointers in `input_toml_m.f90` are required
    by the toml-f API.
- **Explicit deallocation.** The 138 `deallocate`s of local allocatables are
  redundant: automatic deallocation already happens on return. They are
  harmless.
- **Implicit `save`:**
  - `gg1cc.f90:48-51`: `c21…c34` are initialised in their declarations, which
    implies `save`. They are never modified, so they should be `parameter`.
  - `m_psmlout.f90:116,599`: the local `lsymb` is initialised, so it is
    implicitly saved. It also shadows the module-level `lsymb` (`:29`).
- **Explicit `save`.** `exc_libxc.f90:37` (`integer, save :: functprinted`)
  is hidden global state in an external procedure. It is not thread-safe, and
  it makes the printed output depend on call history.
- **Automatic arrays.** 43 local arrays are sized by arguments (e.g.
  `psatom.f90:68` `uua(mmax,nv)`, `run_vkb.f90:65-68`). They live on the
  stack. They are fine at current sizes, but `allocatable` is safer with
  large `mmax`.

## 13. Modules, `use` and the global namespace

**Why it matters.**

- A module without `private` exports every helper and every imported name.
- `use` without `only` imports everything, which hides where names come from
  and invites clashes.

**Default accessibility:**

- `m_psmlout`, `functionals_m`, `input_text_m` and `input_toml_m` have
  `private` by default and an explicit `public` list. This is good.
- `parser` and `xc_messages` (`functionals.F90:21,41`) have no `private`, so
  all their procedures are public.
- `m_libxc_list.F90` marks some items `public` but has no default
  `private`, so the helper parameters `EXCH`, `CORR` and `EXCH_CORR` are
  public too.

**`use` without `only`** (13 of 24 `use` statements):

- `functionals.F90:89-100` (`parser`, `xc_messages`, `iso_c_binding`,
  `xc_f90_lib_m`/`xc_f03_lib_m`, `xc_f03_funcs_m`);
- `exc_libxc.f90:24` (`functionals_m`);
- `m_psmlout.f90:64,65,547,548,1171` (`xmlf90_wxml`, `m_libxc_list`);
- `input_toml_m.f90:3` (`tomlf`).

`iso_fortran_env` is always imported with `only`.

**Module coverage.** Only I/O, XC wrappers and input parsing are in modules
(see [section 1](#1-external-procedures-without-explicit-interfaces)). The
project standard ([Coding Standards](coding_standards.md), rule 3) forbids
"data in modules", but modules that hold only procedures and constants are
compatible with it.

## 14. Code duplication

Since resolved for the `foo`/`foo_r` pairs and the drivers, by the
merge described in [Merging sr and fr](fr_merge_plan.md).

**Why it matters.** A bug fixed in one twin stays in the other. Many of the
Known Issues exist in one sibling but not the other, e.g. `psatom` `it==400`
versus `psatom_r` `it==100`.

| Pair | Lines | Identical lines (ratio) |
|---|---|---|
| `wellstate` / `wellstate_r` | 216 / 227 | 193 (0.87) |
| `fphsft` / `fphsft_r` | 85 / 83 | 71 (0.85) |
| `psatom` / `psatom_r` | 250 / 268 | 217 (0.84) |
| `run_vkb` / `run_vkb_r` | 321 / 354 | 262 (0.78) |
| `run_config` / `run_config_r` | 294 / 342 | 240 (0.75) |
| `run_phsft` / `run_phsft_r` | 129 / 143 | 102 (0.75) |
| `linout` / `linout_r` | 275 / 329 | 222 (0.74) |
| `upfout` / `upfout_r` | 547 / 675 | 445 (0.73) |
| `run_diag` / `run_diag_r` | 155 / 164 | 113 (0.71) |
| `run_plot` / `run_plot_r` | 211 / 248 | 155 (0.68) |
| `fpovlp` / `fpovlp_r` (not compiled) | 81 / 88 | 55 (0.65) |
| `oncvpsp` / `oncvpsp_r` | 664 / 829 | 478 (0.64) |
| `gnu_script` / `gnu_script_r` | 263 / 317 | 185 (0.64) |
| `read_input_toml` / `_r` (in module) | 382 / 382 | ratio 0.97 |
| `read_input_text` / `_r` (in module) | 233 / 224 | ratio 0.91 |
| `psmlout` / `psmlout_r` (in module) | 458 / 580 | ratio 0.81 |

Other near-copies:

- `sratom`/`relatom` (0.80), `lschpb`/`lschpse` (0.79), `lschvkbb`/`lschvkbbe`
  (0.72), `sbf_basis`/`sbf_basis_con` (0.69), `lschfb`/`ldiracfb` (0.60) and
  `modcore2`/`modcore3` (0.50).
- The near-origin series start `sn=ro**(2*ll+3)/dfloat(2*ll+3)` is repeated in
  7 solvers.
- `oncvpsp.F90` and `oncvpsp_nr.F90` differ only in `srel` and banner text.

**Approach.** The fully-relativistic variants mostly add a `kap`/`ikap` loop
and a trailing dimension of 2. Options are:

- pass `srel`/`nkap` and loop over `1:nkap`, with `nkap=1` for the
  scalar-relativistic build;
- or move the shared bodies into module procedures parameterised by the
  number of spin-orbit components.

The two input readers are the easiest target (91–97% identical).

## 15. Dead code

**Why it matters.** Dead code is still compiled, still has to be kept
building, and can be called by mistake. Several Known Issues are in dead
routines.

- **Never called but compiled** (8 procedures):
  - `lschpb.f90:19` and `lschpse.f90:19` (whole files). `lschpse` carries the
    Known Issues `vv`/`rr` bug.
  - `sbf_basis.f90:19` (whole file).
  - `exc_off_pbe.f` `exchpw91:436`, `corlsd:481` and `corpw91:536` (and
    `gcor:516`, which only `corlsd` calls): about 190 lines of PW91 code.
  - `functionals.F90:76` `messages_experimental` and `:715`
    `functional_get_tau`.
  - `oncvpsp_r.F90:786,807`: `cmtskp`/`read_error`, superseded by
    `input_text_m`.
- **Not compiled.** `oncvpsp_nr.F90`, `oncpsp_nr.f90`, `oncvpsp_cpp.f90`,
  `m_libxc_list_cpp.f90` and `fpovlp_r.f90`.
- **Commented-out code.** About 106 lines; see Coding Standard Violations.

## 16. Unused names and intrinsic shadowing

**Why it matters.** Unused dummies and locals hide real problems among the
warnings. A local named after an intrinsic makes the intrinsic unusable in
that scope.

- **Unused locals.** 264 (`-Wunused-variable`) in 50 files. Worst:
  `oncvpsp.F90` (38), `oncvpsp_r.F90` (29), `vrel.f90` (22),
  `modcore.f90` (20), `eresid.f90` (13), `ldiracfs.f90` (11).
- **Unused or unset dummies.** 39 (`-Wunused-dummy-argument`). Beyond those
  in Known Issues:
  - `psatom.f90:22`, `psatom_r.f90:22` (`rcmax`);
  - `modcore2.f90:27`, `modcore3.f90:27` (`irps`);
  - `renorm_r.f90:19` (`ll`);
  - `functionals.F90:165` (`ndim`, `nel`);
  - `m_psmlout.f90:35,37,518,520` (`epstot`, `nrl`, `rc0`).
  - `exc_libxc_stub.f90:20` (all arguments, intentional) declares `vxc`/`exc`
    `intent(out)` and never sets them, before it `stop`s.
- **Unused parameters.** 10 (`-Wunused-parameter`), e.g. `pi` in
  `m_psmlout.f90:69,552`, `eps` in `oncvpsp.F90:77`, and `eps`/`blend` in
  `modcore3.f90:65-66`.
- **Locals named like intrinsics.** 16:
  - `nint` in 9 solvers (e.g. `lschfb.f90:62`, `ldiracfb.f90:62`,
    `lschvkbb.f90:65`);
  - `gamma` in 4 (`lschfb.f90:59`, `lschfs.f90:58`, `fpovlp.f90:50`,
    `vrel.f90:57`);
  - `sum` (`dpnint.f90:39`) and `exp` (`lschpb.f90:52`, in Known Issues).
  - None of these scopes calls the shadowed intrinsic today. A later call to
    the intrinsic in the same scope would not compile.
  - Dummy names `type` and `set` (`m_psmlout.f90:1164`) are legal, but they
    are keywords.

## 17. Fixed-form source and legacy constructs

**Why it matters.** Fixed form is obsolescent, has column-sensitive syntax
and silently truncates at column 72. Constructs such as `GOTO`, `COMMON` and
`EQUIVALENCE` defeat modern analysis.

- **Fixed form.** Only `exc_off_pbe.f` (629 lines) is fixed form. It has
  `c`/`C` comment lines, column-6 continuation (`&` and `1`), upper-case
  keywords, bare `END` (7) and tab-indented code at lines 245 and 247
  (non-standard). `-std=f2018` also rejects its `REAL*8`.
- **Absent from all compiled sources** (verified by script): `GOTO` (none,
  computed or assigned), arithmetic `IF`, labeled `DO`/shared `DO`
  termination, `CONTINUE`, `COMMON` (a commented-out `common/GAS/` at
  `exc_off_pbe.f:39` only), `EQUIVALENCE`, `ENTRY`, `DATA`, statement
  functions and labeled `FORMAT`.
- **Still present (obsolescent or old style):**
  - `character*N`: 31 lines, reported 47 times by `-std=f2018`;
  - `.eq.`-style operators: 89 lines;
  - `return` immediately before `end` (84 of 119 procedures).
  - These are style items; see Coding Standard Violations.

## 18. Hard-coded I/O units

**Why it matters.** Unit numbers 5 and 6 are conventional, not guaranteed.
`iso_fortran_env` provides `input_unit`, `output_unit` and `error_unit`.
Error messages on standard output get mixed into the main output file, which
the tests compare.

- **Unit 6.** 834 `write(6,…)` lines in 49 files. Errors go to unit 6 as
  well, except in `input_toml_m.f90` and parts of `oncvpsp*.F90`, which use
  `stderr`.
- **Unit 5.** `read(5,*)` appears only in the dead `oncvpsp_r.F90:798`.
- **`get_unit`.** `m_psmlout.f90:1147` hand-searches units 10–99; use
  `open(newunit=…)`, as the drivers already do.

## Recommended remediation order

1. **Fix the build.** Make `oncvpspnr.x` really non-relativistic (test
   `RELATIVISTIC`, or set `srel` from it), and remove `oncvpsp_nr.F90` and
   `oncpsp_nr.f90`. Fix the duplicated or missing CMake entries. Drop or
   regenerate `src/Makefile`.
2. **Fix the argument mismatches.** `debl` shape in `linout_r`/`upfout_r` and
   `check_data`, `nproj(5)`, and `wellstate_r` (Known Issues). Add a CI job
   that compiles the concatenated sources with gfortran (as described
   [above](#how-the-counts-were-made)) and fails on argument mismatches. That
   gives most of the benefit of explicit interfaces before any refactor.
3. **Error handling.** Check `ierr` after `sratom`/`relatom` and the
   unchecked solver calls. Test `open` `iostat` before reading. Replace
   library `stop`s with `error stop` plus a non-zero code (a one-line change
   per site), so failures are visible to scripts and CTest. Fix the
   `upfout` format.
4. **Initialise search results** ([section 4](#4-variables-set-on-only-some-paths))
   and add `implicit none` to `aii` and `m_psmlout`.
5. **Fix declaration order** ([section 5](#5-non-standard-declaration-order))
   so the code builds with `-std=f2008` (apart from `exc_off_pbe.f`). Add
   `-std=f2008 -Wimplicit-interface` to the Debug flags.
6. **Introduce a constants module.** It holds `dp` (from `iso_fortran_env`),
   `pi`, `c_light` and the array limits (30/6/5/4/`mxprj`). Replace the 81
   local `dp` definitions, the 7 copies of 137.036 and the literal shapes.
   Fix the 1/24 coefficient (this changes the reference outputs in the last
   digits).
   - **Done:** `src/constants_m.f90` holds `dp`, `pi`, `c_light` and
     `mxprj`. Every routine except the imported `functionals.F90` and
     `m_libxc_list.F90` takes `dp` from it. `pi` and 137.036 are no longer
     defined locally, except for the literals in the imported
     `exc_off_pbe.f`. These changes leave the output bit-identical. The
     1/24 coefficient is fixed.
   - **Done:** the array limits. `constants_m` also holds `mxstate=30`
     (states), `mxcnf=5` (reference and test configurations), `mxl=3`
     (largest l; `mxl+1` projector channels, `mxl+2` semi-local potentials
     with the local one last), `mxlin=5` (largest lmax the readers accept;
     `mxlin+1` slots of the per-l input arrays), `mxkap=2`, `nrhomod=5`
     (model core and 4 derivatives), `nlev=3` (levels of a configuration)
     and `nerrlv=7` (convergence-profile error levels), and every
     declaration with those meanings uses them. The extents are unchanged
     except where they disagreed: `nproj(5)` in `psatom` and `run_config` is
     `mxlin+1` like its actual, and `linout`'s `npr_sr(5)`, `npr_so(5)` are
     `mxl+1`, the channels they count.
   - **Left:** local scratch arrays (`usvd(10,10)`, `lpgh(20)`, ...),
     small fixed tuples (`fcgrid(3)`, `log_mesh(3)`, `eterm(4)`), the pairs
     `uu(mmax,2)`, `fat(:,2)` and `lschvkbs`'s two projectors, and
     `tfapot.f90`'s `**.33333333d0`.
7. **Move procedures into modules**, starting with the leaf routines (`aeo`
   family, `dpnint`, `vpinteg`, the `lsch*`/`ldirac*` solvers). Then switch
   their dummies to assumed shape and drop the redundant size arguments.
   Replace the sequence-association calls (`a(1,l1)` → `a(:,l1)`) at the
   same time.
   - **Superseded** by [Library and wrappers](library_plan.md), phase 4:
     the dummies stay explicit-shape, so that the routines remain callable
     from C, and modules come after error codes and the driver split.
8. **Merge the `foo`/`foo_r` twins**, beginning with the input readers and
   `psmlout`/`psmlout_r`.
9. **Mechanical cleanup.**
   - **Done** (with fortitude; see
     [Coding Standard Violations](coding_standard_violations.md)):
     generic intrinsics in place of `dlog`, `dsqrt`, `dabs`, `dexp`,
     `dmin1`, `dmax1`, `idint` and `dfloat`; `only` on every `use`;
     `intent` on every argument; `implicit none` in `m_psmlout`; the
     unused locals. All except in the imported files.
   - Add explicit `int()`/`nint()`.
   - Replace exact-equality convergence tests with tolerances.
   - Remove the dead routines.
   - Make `gg1cc`'s coefficients `parameter`s.
   - Convert `exc_off_pbe.f` to free form with `implicit none`.
