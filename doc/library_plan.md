## A library and Python/Julia wrappers

**Status: to do. Nothing here has started.**

The goal is to call ONCVPSP from Python and Julia: pass in an input, get
back the pseudopotential and the test results as arrays, and never take the
interpreter down. Today the code is a program. `oncvpsp_main(irel,istage)`
parses the command line, reads the input file, computes, writes to standard
output and the output files, and ends with `error stop` when anything fails.
This is the plan for turning it into a library with a C interface, and
wrappers on top of that.

Moving routines into modules is part of the plan but not its first step.
Wrappers never call the internal routines, module or external: compilers
name their symbols differently (`vrel_`, `__m_MOD_vrel`, `m_mp_vrel_`), so
the only portable entry points are `bind(C)` routines. What blocks wrapping
is the program-style behaviour (`stop`, unit 6, files), not the lack of
modules.

### Rules

- **Output is unchanged.** The rules of [Merging sr and fr](fr_merge_plan.md)
  apply to every step: `ctest -j8`, `scripts/check_interfaces.sh`,
  line-for-line comparison of the outputs with those
  of the commit before the step, and one commit per step. The executables
  keep their command lines, messages and exit codes.
- **No new `#if`.**
- **Interfaces that C could call.** Dummies stay integer and `real(dp)`
  scalars and explicit-shape arrays with their sizes passed in. This is
  what `bind(C)`, ctypes, cffi, f2py and Julia's `ccall` all handle without
  copies, so the internal routines are not switched to assumed shape (as
  item 7 of the [remediation order](fortran_best_practices.md#recommended-remediation-order)
  once suggested).
- **No data in modules** ([Coding Standards](coding_standards.md)), so
  that two calculations in one process cannot see each other.

### Phase 1: errors returned, not stopped

There are 56 `stop`/`error stop` lines in 34 files. Those in
`oncvpsp_main.F90` (the command line and the input file) stay. Every other
one becomes an error code:

- [ ] Each routine that can fail gets an `ierr` argument (0 for success) and
      returns at the failure, after writing its message as now.
- [ ] Each caller tests `ierr` and passes it up, as `boundstate` and
      `lschvkbb` already are at `oncvpsp_main.F90:630,835`.
- [ ] `oncvpsp_main` turns a non-zero `ierr` into `error stop 1`, so the
      executables behave as before.
- [ ] A list of the codes, by routine, in this page or in the headers.

Order: the leaf routines first (`dpnint`, `qroots`, `vout`, `sbf_*`), then
upwards. The imported `functionals.F90` keeps its `stop`s; `exc_libxc`
checks the functional before calling it.

### Phase 2: output to a chosen unit

598 `write(6,…)` lines in 45 files. A library has to be able to send that
text elsewhere, or nowhere.

- [ ] Decide how the unit reaches the routines. Options: an `iout` argument
      on every routine that writes (long argument lists, but explicit), or
      one `integer, parameter` in `constants_m` that is opened on a file or
      `/dev/null` (simpler, but a global unit).
- [ ] Error messages go to the same unit as now, so the reference outputs
      do not change.

### Phase 3: a driver split in three

`oncvpsp_main.F90` (1090 lines) is split into:

- [ ] `oncvpsp_run(...)`: the calculation, from the input values (as
      arguments) to the results (as arrays), with `ierr`. No command line,
      no input file, no output files.
- [ ] the readers, which fill the input arguments (`read_input_text`,
      `read_input_toml`, already separate);
- [ ] the writers (`linout`, `upfout`, `psp8`, `psmlout`), called by the
      driver after `oncvpsp_run`;
- [ ] `oncvpsp_main`, left as the command line, the readers, `oncvpsp_run`
      and the writers.

Open question: the results whose size is known only during the run (the
number of projectors, the mesh). Either the caller passes maximum sizes
(`mxprj`, `mmax`), or a first call returns the sizes.

### Phase 4: routines into modules

Optional for the wrappers, but it makes every call checked by the compiler
and takes the routine names out of the global namespace, which matters once
the library shares a process with other Fortran code.

- [ ] Decide whether to do it at all: it touches nearly every file, which
      makes merging upstream fixes from `master` much harder.
- [ ] If so, one group per commit, leaves first: the mesh and quadrature
      utilities (`aeo` family, `dpnint`, `vpinteg`, `der2exc`), the XC
      routines, the radial solvers (`lsch*`, `ldirac*`, `boundstate`), the
      pseudization and optimizer, the `run_*` drivers.
- [ ] Each group updates `src/CMakeLists.txt`.
- [ ] Once all are done, `scripts/check_interfaces.sh` is no longer needed.

### Phase 5: the C interface

- [ ] A module `oncvpsp_c` of `bind(C, name="oncv_…")` routines over
      `oncvpsp_run`, using `iso_c_binding` kinds (`c_int`, `c_double`,
      checked equal to `dp`), with character input as `c_char` arrays and
      every result in caller-allocated arrays.
- [ ] A C header for it, and a shared library target (`liboncvpsp`) in
      CMake.
- [ ] A C test in `tests/` that runs one of the reference inputs through
      the interface and compares with the reference output.
- [ ] Check that two runs in one process give the same results (the
      `save` in `exc_libxc.f90` changes what is printed, not the numbers).

### Phase 6: wrappers

- [ ] Python: ctypes or cffi over `liboncvpsp`, built with
      scikit-build-core, returning NumPy arrays; tests against the
      references.
- [ ] Julia: `ccall` over the same library, later packaged as a JLL with
      BinaryBuilder.
- [ ] Both take the input as a dictionary or keyword arguments, in the
      layout of the TOML input.
- [ ] A user page on the wrappers, in the `nav` of `mkdocs.yml`.
