## ONCVPSP CODING STANDARDS

The coding standards I've tried to follow and will want for any contributed additions are pretty basic.
They should also be obvious from examininga few of the source routines.

1) Fortran 95 free-format.

2) One subroutine per file.

3) All data is passed as subroutine scalar or array variables.  There are to
be no common blocks, data in modules, or derived data types.

4) All variables are to be explicitly typed (eg. "implicit none").

5) No one-letter variable names.

6) Input, output, and local variables should be grouped.

7) All local arrays should be allocated and deallocated (unless they're tiny).

8) Each subroutine should have a brief statement of what it does immediately
before or after the subroutine statement.

9) All subroutine arguments should be defined in comments following the
subroutine statement.

10) There should be comments throughout to indicate what is going on.
Any new math should be described in a text or pdf contribution to the doc directory.

11) I'd like to stick to one input file read on standard input and one
output file to standard output, with shell scripts to extract any new
sections.

12) If you add new input data, the process of reading and checking it should
be done in such a way that older files still run.  This could also be
accomplished by building a separate executable for the new feature.

There are probably other conventions I've been following but you get the idea.

## Conventions

Hamann's list ends with "There are probably other conventions I've been
following". These are those conventions, read off the majority of his
files, and they apply to all the code except the imported files (see
[Linting](#linting)) and the exceptions kept in
[Coding Standard Violations](coding_standard_violations.md). A routine
looks like this:

```fortran
!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! ... (the GPL notice)
!
 subroutine shift_pot(rr,vv,rcut,vshift,mmax,nshift,iout)

!! adds a constant to a potential inside a radius
!!
!! @param rr  log radial mesh
!! @param vv  potential (Ha), shifted on output
!! @param rcut  radius inside which vv is shifted (bohr)
!! @param vshift  the shift (Ha)
!! @param mmax  size of the log radial mesh
!! @param nshift  number of mesh points shifted
!! @param iout  unit of the text output

 use constants_m, only: dp

!Input variables
 integer, intent(in) :: mmax,iout
 real(dp), intent(in) :: rr(mmax),rcut,vshift

!Input/output variables
 real(dp), intent(inout) :: vv(mmax)

!Output variables
 integer, intent(out) :: nshift

!Local variables
 integer :: ii

! shift the points inside rcut
 nshift=0
 do ii=1,mmax
   if(rr(ii)<rcut) then
     vv(ii)=vv(ii)+vshift
     nshift=nshift+1
   end if
 end do
 if(nshift==0) then
   write(iout,'(/2a,f12.6)') 'shift_pot: WARNING no mesh point inside rcut', &
&        ', rcut=',rcut
 end if
 return
 end subroutine shift_pot
```

**Files.**

- Every procedure is in a module; there are no external procedures, and
  the only program units outside modules are the main programs. This
  replaces rules 2 and 3 on this branch: the compiler checks every call
  against its definition, and the names stay out of the global
  namespace.
- One module per file, `<name>_m` in `<name>_m.f90`. A module around a
  single routine is named after it (`lschfb_m` holds `lschfb`; `modcore`
  is in `modcore1_m`, as `modcore_m` holds the helpers the model cores
  share), one that gathers several after their subject (`boundsolve_m`).
  It starts `implicit none`, `private` and a `public` list, and callers
  `use` it with an `only` list in the procedure that makes the call.
- `.f90`, or `.F90` when the file uses the preprocessor.
- Helpers used by one routine only are internal procedures after
  `contains`, private procedures of its module, or functions in an
  include file (`adams.inc`).
- The external libraries are called through interfaces: LAPACK and BLAS
  through `lapack_m`, which declares the routines ONCVPSP calls as the
  reference LAPACK does, libxc through its own Fortran module.
- No global state, which is what is left of rule 3: a module holds
  parameters (`constants_m`), types and procedures, but no variables, so
  everything that changes is passed as arguments, and runs can share a
  process in threads. `tests/unit/test_run_threads.f90` checks it: runs
  at once in OpenMP threads give the results of each alone, bit for bit.
  The procedures are compiled recursive (`-frecursive`), so no local
  array is static either.
- Derived types are allowed for values that belong together and are
  passed together (`piecewise_t` in `pieces_m`, the integrator state in
  `outwf_m`, a configuration of the atom, `config_t` in `input_m`, the
  results of a run, `results_t` in `results_m`). A type is named `_t`,
  documents each component with `!>`, gives it a default value where one
  makes sense, and a procedure that needs only part of a larger type
  takes that part (`fill_t` in `modcore_budget_m`). Older code keeps
  separate arrays with a common index, like the reference configuration's
  `na`, `la` and `fa`. Types that a library defines (toml-f's
  `toml_table`, xmlf90's `xmlf_t`) are used where its calls need them.
- The non-, scalar- and fully-relativistic codes share their routines.
  A routine whose work depends on the equation takes `irel` (0 non-, 1
  scalar-, 2 fully-relativistic), tested at run time, rather than a
  logical or a preprocessor macro. Arrays indexed by l have a last index
  for kappa of size `nkap` (1, or 2 for `irel=2`): 1 for j=l+1/2, 2 for
  j=l-1/2. Loops over it run to `mkap`, which is 1 when `nkap=1` or l=0.
  `boundstate`, `scatterstate` and `aeatom` choose the Schroedinger or Dirac
  solver. `_r` is left only on routines that exist for the
  fully-relativistic code alone (`sr_so_r`, `renorm_r`). See
  [Merging sr and fr](fr_merge_plan.md).
- The entry points a C wrapper will call (see [Library
  plan](library_plan.md)) take integer and `real(dp)` scalars and
  explicit-shape arrays, with their sizes as arguments: no logical,
  optional, assumed-shape, allocatable, pointer or derived-type dummies.
  Routines inside the library may take derived types.

**Header.**

- The `!` + copyright + GPL block heads every file; a new file carries the
  GPL notice under its own authors' copyright.
- A description, as `!!` comments, right after the `subroutine` or
  `function` statement (rule 8).
- Then, after an empty `!!` line, one `!! @param name  description` line
  for each dummy argument (rule 9), name and description separated by
  two spaces, with units (bohr, Ha) where they apply; continuation lines
  are indented under the description. A function result is
  `!! @return description`. Arguments described together are named
  `a,b` or `a..b` (a through b in argument order).
- The description and the arguments are one block of `!!` lines, with no
  blank or plain `!` line inside it, and nothing between it and the
  statement but blank lines. These are the doc comments that fortls,
  the language server of the Modern Fortran extension, shows on hover
  (see [Developer setup](developer_setup.md)); a plain `!` comment is
  not shown, and a blank line ends the block (fortls keeps only the last
  block). Text that would not read well run together, such as a table,
  goes in plain `!` comments after the block.
- Module constants (and any module variable) have a `!>` comment above
  the declaration, or a `!<` comment after it on the same line, with one
  name per documented declaration. A module may have a `!>` description
  above its `module` statement, as `constants_m` does.

**Declarations.**

- `use` statements with an `only` list. `implicit none` is the
  module's, after its own `use` statements, and is not repeated in its
  procedures (rule 4 holds through it); a main program has its own.
- The four section labels, as full-line comments in this order, each
  followed by the declarations it covers (rule 6): `!Input variables`
  (`intent(in)`), `!Input/output variables` (`intent(inout)`), `!Output
  variables` (`intent(out)`, and a function result) and `!Local
  variables` (no intent, including parameters). A label with nothing
  under it is left out; `!Output variables - printing only` marks a
  routine whose only output is printed.
- Every dummy argument has an `intent`, and every declaration has `::`.
- Several names may share a declaration (`real(dp) :: de,emax,emin`).
- Reals are `real(dp)` with `dp` from `constants_m` (see
  [Kinds and constants](#kinds-and-constants)); literals are written
  `1.0d0`, with a digit on each side of the point and the kind always
  given.
- Within each section, the integers that bound arrays are declared before
  the arrays (`integer, intent(in) :: mmax` before `rr(mmax)`): standard
  Fortran 2008, which the Debug build checks (`-std=f2008`).
- Explicit-shape (`rr(mmax)`) or allocatable arrays; no assumed-size
  (`*`) dummies, except in `adams.inc` and the LAPACK interfaces of
  `lapack_m`.

**Names.**

- No one-letter names (rule 5). Scalars double the letter: `ii`, `jj`,
  `ll` (angular momentum), `rr` (radius), `zz` (atomic number); `l1` is
  `ll+1`, the angular-momentum index.

**Layout.** `scripts/reindent.py` sets the indentation; the rest is by
hand.

- The `subroutine`/`function`/`program`/`module` statement and the
  statements directly inside it in column 2.
- Each block body 2 columns right of its opener. `else`, `else if` and
  `case` are in the opener's column, and `case` in the column of
  `select case`.
- Full-line comments in column 1. A comment line describes what follows
  it; prose comments start with `! ` and lower case, and doc comments
  (see Header) with `!! `, `!> ` or `!< `.
- Continuation lines start with `&` in column 1, and the line before ends
  with `&`. The text after the `&` is aligned for reading.
- At most 80 columns.
- Lower-case keywords; `end do`, `end if`, `else if`, `end select` and
  `end subroutine name`; `==`, `/=`, `<` etc., not `.eq.`; `if(`,
  `write(`, `allocate(` without a space before the parenthesis.
- One statement per line, apart from short initialisations
  (`uu(:)=0.0d0; up(:)=0.0d0`).
- `return` before the `end` of every routine, and before `contains` in a
  routine with internal procedures.

**Output.**

- `write(iout,'(…)')` with the format inline; no labeled `FORMAT`
  statements, and no list-directed (`*`) output in new code.
- A routine that writes, or calls one that does, takes the unit of the
  text output as an argument, `iout`, last but for `ifail`; internal
  procedures use their host's. The programs, `oncvpsp_main` and the unit
  tests pass 6 (`integer, parameter :: iout=6`). Only the command line
  (`cli_m`, `parse_command_line`) writes to unit 6 itself, and the checks
  that stay `error stop` (rule 5 of [Error messages](#error-messages)) to
  `error_unit`.
- Errors and warnings as in [Error messages](#error-messages).

**Comments and dead code.**

- No commented-out code. A real alternative that a user may want (such
  as the other built-in meshes in the main programs) is described in
  prose.

## Kinds and constants

Take the real kind and shared constants from `constants_m`
(`src/constants_m.f90`), with an `only` list:

```fortran
 use constants_m, only: dp, pi
```

Don't define `dp`, `pi` or the speed of light locally. Write rational
coefficients as quotients (`1.0d0/24.0d0`), not as truncated decimals.

## Input echo

The input is read into an `input_t` (`input_m`), which the readers
(`read_input_text`, `read_input_toml`), `check_data` and the writers take
whole. A new input value is a component of `input_t`, documented there,
read by both readers. The input is echoed in the output and in the psp8,
UPF and PSML files, and each echo must run again to the same output (the
`echo` tests check this). A new input value goes into all four echoes:
the main programs, `linout` (`_r`), `upfout` (`_r`) and `echo_input` in
`psmlout_m`, or into a section
of `input_echo_m`, which they share. An optional value is echoed only when
the input sets it, so that older inputs echo as before (rule 12 above).
Echo the value as read, not one the program has changed since (see `rc0`
and `fcfact0`).

## Error messages

These rules were added after the list above. Messages that were in the
4.0.1 release keep their text, because users' scripts and the reference
outputs may depend on it, but failures now exit with status 1. New
messages follow these rules.

A routine reports a failure it cannot recover from like this, and returns
(see [Library and wrappers](library_plan.md), phase 1):

```fortran
 if(irphs==0) then
   write(iout,'(/2a,f12.6,a,f12.6)') 'run_phsft: ERROR rxpsh beyond radial mesh', &
&        ', rxpsh=',rxpsh,', rr(mmax)=',rr(mmax)
   ifail=1
   return
 end if
```

and each caller passes it on:

```fortran
 call run_phsft(...,iout,ifail)
 if(ifail/=0) return
```

1) Write the message to `iout`, the output file (unit 6 from the
programs), so it appears next to the output that led up to it. This holds
for every message; command-line errors go to unit 6.

2) Start the message with a blank line (`/` at the start of the format).

3) The text is `<routine>: ERROR <description>`. `<routine>` is the name
of the subroutine or program that makes the check, exactly as in its
`subroutine` or `program` statement, with no line number. `ERROR` is in
capitals. The description is lower case and has no final period.

4) After the description, give the values that explain the error as
`name=value` pairs separated by `, `. Use the names from the source or the
input file, e.g. `rxpsh=` or `rr(mmax)=`. Use `i0` for integers, `f12.6`
for reals of known scale (radii in bohr, energies in Ha) and `es12.4` for
anything else. (`f0.N` would drop the leading zero, and `es0.N` needs
Fortran 2018.) Use `a` with `trim` for strings. To split a long message,
end one string and start another (`'...', &` / `&  '...'`) rather than
continuing inside a string, which would copy the indentation into the
message.

5) Fail with `ifail=1` and `return`, `ifail` being the routine's last
argument (0 on success). Only the programs, the command line and
`oncvpsp_main` stop, with `error stop 1`, not `stop`, so that the exit
status is non-zero and scripts and CTest see the failure. A check of a
condition the code has just established, or that `check_data` has already
validated, is a bug if it fails: it stays `error stop 1`, with its message
on `error_unit` (`use, intrinsic :: iso_fortran_env, only: error_unit`).

6) If a called routine reports a failure through `ierr`, the caller that
decides to stop prints the message. It names the called routine and
includes `ierr`, e.g. `oncvpsp: ERROR sratom failed, ierr=100`.

7) A routine that already has an `ierr` argument should set `ierr` and
return when its callers can recover from the failure; `ifail` is for the
failures they cannot. The two are separate arguments.

8) A warning has the same form with `WARNING` in place of `ERROR`, and the
program continues.

## Linting

The [pre-commit](https://pre-commit.com) hooks in `.pre-commit-config.yaml`
check every commit, and CI runs them on every file. Install them once in a
clone with

```bash
pre-commit install
```

and run them on everything with `pre-commit run --all-files`. They trim
trailing whitespace and end files with one newline, re-indent the Fortran
with `scripts/reindent.py` (all three fix the files in place), check the
YAML and TOML files and the scripts' shebangs, and run

- [fortitude](https://fortitude.readthedocs.io) on the Fortran, with the
  settings in `fortitude.toml`;
- [shellcheck](https://www.shellcheck.net) on the shell scripts written
  since 4.0.1;
- [ruff](https://docs.astral.sh/ruff) on the Python, with `ruff.toml`.

`.editorconfig` tells editors the same whitespace rules (LF line endings,
no trailing whitespace, one final newline, spaces not tabs, 2-column
Fortran indentation, 80 columns) so that files come out of the editor
already conforming. Most editors honour it natively or with the
[EditorConfig](https://editorconfig.org) plugin. [Developer
Setup](developer_setup.md) has the installation of the tools and the VS Code
setup that runs `reindent.py` as the formatter and fortitude in the editor.

fortitude uses its default rules and its preview rules, without the three
that conflict with the rules above: double-quoted strings (S241), sorted
`use` statements (S271), and `implicit none (type, external)` (C003),
which is Fortran 2018. It adds a literal's kind (C021) and `implicit none`
only on the module (S201); with the preview default S291, a literal has a
digit on each side of the point (`0.0d0`, not `0.d0`). The code keeps to Fortran 2008
in that, and gets the check C003 asks for from gfortran instead:
`-Werror=implicit-interface`, in every build type, makes a call of a
procedure without an explicit interface a compile error. fortitude's
rules against external procedures (C091, C092) are on, so a routine
outside a module or an `external` declaration fails the check;
`lapack_m.f90` alone may have assumed-size arrays (C071). Among other
things fortitude wants `intent` on every argument, `only` on every `use`,
`==` rather than `.eq.`, `end do` rather than `enddo`, `character(len=N)`,
`::` in every declaration, generic intrinsics (`log`, not `dlog`), no
assumed-size (`*`) arrays, and lines of at most 80 columns (set in
`fortitude.toml`). `fortitude check --fix` fixes some of these itself.

gfortran adds two checks in the Debug build: `-std=f2008`, so the code is
standard Fortran 2008 (a bound, for one, is declared before the array it
bounds), and `-Wconversion`, a warning for a conversion that may change a
value (write `int()` or `real(,dp)`). The imported files are compiled
without them.

The imported files (`exc_off_pbe_m.f`, `functionals.F90`,
`gg1cc`/`gp1cc`/`gpp1cc` and `m_libxc_list.F90`) are excluded from the
hooks, so that they can be kept in sync with their sources, and so are the
reference outputs in `tests/refs`, the patches in `cmake/patches` and
`COPYING`.
