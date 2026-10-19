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

## Kinds and constants

Take the real kind and shared constants from `constants_m`
(`src/constants_m.f90`), with an `only` list, before `implicit none`:

```fortran
 use constants_m, only: dp, pi
 implicit none
```

Don't define `dp`, `pi` or the speed of light locally. Write rational
coefficients as quotients (`1.0d0/24.0d0`), not as truncated decimals.

## Input echo

The input is echoed in the output and in the psp8, UPF and PSML files, and
each echo must run again to the same output (the `echo` tests check this).
A new input value goes into all four echoes: the main programs, `linout`
(`_r`), `upfout` (`_r`) and `echo_input` in `psmlout_m`, or into a section
of `input_echo_m`, which they share. An optional value is echoed only when
the input sets it, so that older inputs echo as before (rule 12 above).
Echo the value as read, not one the program has changed since (see `rc0`
and `fcfact0`).

## Error messages

These rules were added after the list above. Messages that were in the
4.0.1 release keep their text, because users' scripts and the reference
outputs may depend on it, but they now stop with `error stop` so that
failures exit with status 1. New messages follow these rules.

A fatal error is reported like this:

```fortran
 if(irphs==0) then
   write(6,'(/2a,f12.6,a,f12.6)') 'run_phsft: ERROR rxpsh beyond radial mesh', &
&        ', rxpsh=',rxpsh,', rr(mmax)=',rr(mmax)
   error stop 1
 end if
```

1) Write the message to unit 6, the output file, so it appears next to the
output that led up to it. This holds for every message, including
command-line errors.

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

5) Stop with `error stop 1`, not `stop`, so that the exit status is
non-zero and scripts and CTest see the failure.

6) If a called routine reports a failure through `ierr`, the caller that
decides to stop prints the message. It names the called routine and
includes `ierr`, e.g. `oncvpsp: ERROR sratom failed, ierr=100`.

7) A routine that already has an `ierr` argument should set `ierr` and
return, not stop, when its callers can recover from the failure.

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
trailing whitespace and end files with one newline (fixing the files in
place), check the YAML and TOML files and the scripts' shebangs, and run

- [fortitude](https://fortitude.readthedocs.io) on the Fortran, with the
  settings in `fortitude.toml`;
- [shellcheck](https://www.shellcheck.net) on the shell scripts written
  since 4.0.1;
- [ruff](https://docs.astral.sh/ruff) on the Python, with `ruff.toml`.

fortitude uses its default rules, without the four that conflict with the
rules above: double-quoted strings (S241), and three that ask for modules
or `implicit none (external)` in place of external procedures (C003, C091,
C092). Among other things it wants `intent` on every argument, `only` on
every `use`, `==` rather than `.eq.`, `end do` rather than `enddo`,
`character(len=N)`, `::` in every declaration, generic intrinsics (`log`,
not `dlog`), no assumed-size (`*`) arrays, and lines of at most 100
columns. `fortitude check --fix` fixes some of these itself.

The imported files (`exc_off_pbe.f`, `functionals.F90`,
`gg1cc`/`gp1cc`/`gpp1cc` and `m_libxc_list.F90`) are excluded from the
hooks, so that they can be kept in sync with their sources, and so are the
reference outputs in `tests/refs`, the patches in `cmake/patches` and
`COPYING`.
