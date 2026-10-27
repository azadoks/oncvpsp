## A library and Python/Julia wrappers

**Status: in progress. Phases 2 and 4 are done, and phase 1 but for the
stops that move with phase 3.**

The goal is to call ONCVPSP from Python and Julia: pass in an input, get
back the pseudopotential and the test results as arrays, and never take the
interpreter down. Today the code is a program. `oncvpsp_main(irel,istage)`
parses the command line, reads the input file, computes, writes to standard
output and the output files, and ends with `error stop` when anything fails.
This is the plan for turning it into a library with a C interface, and
wrappers on top of that.

Moving routines into modules was part of the plan, though not needed for
wrapping, and is done (phase 4). Wrappers never call the internal routines, module or external: compilers
name their symbols differently (`vrel_`, `__m_MOD_vrel`, `m_mp_vrel_`), so
the only portable entry points are `bind(C)` routines. What blocks wrapping
is the program-style behaviour (`stop`, unit 6, files), not the lack of
modules.

### Rules

- **Output is unchanged.** The rules of [Merging sr and fr](fr_merge_plan.md)
  apply to every step: `ctest -j8`,
  line-for-line comparison of the outputs with those
  of the commit before the step, and one commit per step. The executables
  keep their command lines, messages and exit codes.
- **No new `#if`.**
- **Interfaces that C could call.** The dummies of the entry points (the
  `bind(C)` routines and `oncvpsp_run`) are integer and `real(dp)`
  scalars and explicit-shape arrays with their sizes passed in. This is
  what `bind(C)`, ctypes, cffi, f2py and Julia's `ccall` all handle without
  copies, so the internal routines are not switched to assumed shape (as
  item 7 of the [remediation order](fortran_best_practices.md#recommended-remediation-order)
  once suggested).
- **No data in modules** ([Coding Standards](coding_standards.md)), so
  that two calculations in one process cannot see each other.

### Phase 1: errors returned, not stopped

A routine that can fail returns, instead of stopping, with a last argument
`ifail`: 0 on success, 1 after a failure whose message it has written
(unit 6, as before). Every caller tests it right after the call and returns
with it set, so the failure reaches `oncvpsp_main`, which ends with
`error stop 1` as before: the output and the exit status do not change.

- `ifail` is a new argument, not the existing `ierr`. The solvers'
  `ierr` (`lschfb`, `ldiracfb`, `lschvkbb`, `scatterstate`, ...) carries
  codes that some callers accept and go on from (a state that did not
  converge, 100; a ghost below the node count, 200); folding the failures
  into it would let those callers carry on after one. `ifail` has one
  meaning: the calculation cannot go on.
- A procedure with optional arguments gets `ifail` last, passed by
  keyword (`ifail=ifail`), so that it cannot land in an optional slot.
- A reader that fails after opening its file closes it first.
- Checks of conditions the code has just established, or that
  `check_data` has already validated, stay `error stop 1`: they are bugs,
  not failures, and returning from them would put `ifail` on every call of
  the XC routines and of the uniform-grid quadratures. They are `vout`'s
  `mode`, `exc_dispatch`'s `iexc`, `sbf_rc_der`'s and `sbfq_rc_der`'s
  `llin`, `vploc_fit`'s `lpopt`, the segments of `pw_alloc`,
  `pw_product_integral` and `outwf_setup`'s sources, `add_zero_r`'s sizes,
  the `ug_*` node counts, grid positions, brackets and `ug_fmin`'s
  iterations, `diagout`'s bracket,
  `eresid`'s q grid, `gauss_legendre`'s Newton iterations, and
  `dirac_outward`'s series at the origin, which converges at any energy
  of the solvers' brackets.
- A build without libxc refuses a libxc functional in `check_data`
  (`have_libxc`), so `exc_libxc`'s stub, which stops, is not reached.
- The imported `functionals.F90` keeps its `stop`s.

Steps, one commit each, every one with the outputs of the reference
inputs and of a corpus of failing inputs (stdout, and stderr up to the
backtrace, and the exit status) unchanged:

- [x] The readers: `read_input_text`, `read_input_toml`, `check_data`,
      the all-electron and pseudopotential data files (`ae_file_m`,
      `ps_file_m`, `load_pair`).
- [x] The radial solvers and their callers: `match_point`, `bs_solve`
      (its `match` callback), `lschfb`, `ldiracfb`, `boundstate`,
      `outwf_m`, `lschvkbb`, `renorm_r`, the well states and free tails,
      and up through the atoms and configurations. A routine of a later
      step stops on the `ifail` of one of these until its own step.
- [x] The pseudization: `run_optimize`, `optimize`, `qroots`, `pspot`,
      `const_basis`, `sbf_basis_con`, `eresid`, `ode_rc_der`, `run_vkb`,
      `sr_so_r`.
- [x] The model cores: `modcore`, `modcore2`, `modcore3`,
      `poly_core_fit`, `teter_scale`, `modcore_budget` (its steps,
      `no_model`, and `zeroin` and its callers) and
      `modcore_budget_table`. Nothing the searches call back can fail:
      `set_model`'s three checks stay error stops (a radius beyond the
      table, which covers the searched interval; `dhseqr` on the tail's
      companion matrix; `dgglse` on constraints of full rank by
      construction).
- [x] The tests and the writers: `run_diag`, `run_plot`, `run_ghosts`,
      `run_phsft`, `phsft_m`, `outmesh_m`, `linout` and `upfout`. The
      modules of `oncvpsp_confine` alone (`confine_ws_m`, `ws_basis_m`,
      `confine_wall_m`) are the program's, and stop on an `ifail` of the
      library: their phase searches call back into `phsft_m`.
- [ ] `oncvpsp_main`'s own stops on the results (the command line and the
      input file stay): they move with phase 3.

### Phase 2: output to a chosen unit

A library has to be able to send the text output elsewhere, or nowhere.

- [x] Every procedure that writes, or calls one that does, takes the unit
      as an argument, `iout`, before `ifail` (by keyword where there are
      optional arguments), and writes to it: 140 procedures. Internal
      procedures use their host's `iout`. A module variable would have
      been shorter, but it is the global state that the rules exclude.
- [x] `oncvpsp_main`, the programs and the unit tests pass 6
      (`integer, parameter :: iout=6`), so the output does not change.
      `oncvpsp_confine`'s two callbacks of the root searches (`phase`,
      `u_rs`), whose interfaces have no room for it, do the same.
- [x] Unit 6 stays for the command line (`cli_m`, `parse_command_line`,
      `write_help`) and for the messages of the checks that stay
      `error stop` (phase 1): they are bugs, not output.
- [ ] A caller that wants no output opens `/dev/null` (or `NUL`) and
      passes its unit: with phase 3.

### Phase 3: a driver split in three

`oncvpsp_main_m.F90` (1090 lines) is split into:

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

Done on smooooth: every procedure is a module procedure, so every call is
checked by the compiler and the routine names are out of the global
namespace (see [Coding standards](coding_standards.md)).

- [x] One group per commit, leaves first: the XC, mesh and Bessel-basis
      routines; the radial solvers; the atom, pseudization and model-core
      routines; the `run_*` drivers and the output writers. Each `foo.f90`
      became `foo_m.f90`, and the output is byte-identical.
- [x] LAPACK and BLAS called through the interfaces of `lapack_m`.
- [x] `scripts/check_interfaces.sh` removed; fortitude's rules against
      external procedures (C091, C092) are on.

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
      `save` in `exc_libxc_m.f90` changes what is printed, not the numbers).

### Phase 6: wrappers

- [ ] Python: ctypes or cffi over `liboncvpsp`, built with
      scikit-build-core, returning NumPy arrays; tests against the
      references.
- [ ] Julia: `ccall` over the same library, later packaged as a JLL with
      BinaryBuilder.
- [ ] Both take the input as a dictionary or keyword arguments, in the
      layout of the TOML input.
- [ ] A user page on the wrappers, in the `nav` of `mkdocs.yml`.
