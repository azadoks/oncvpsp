## A library and Python/Julia wrappers

**Status: in progress. Phases 1 to 4 are done; phases 5 (the C
interface) and 6 (the wrappers) are next.**

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
  `bind(C)` routines; `oncvpsp_run`, which takes derived types, is
  wrapped behind a handle in phase 5) are integer and `real(dp)`
  scalars and explicit-shape arrays with their sizes passed in. This is
  what `bind(C)`, ctypes, cffi, f2py and Julia's `ccall` all handle without
  copies, so the internal routines are not switched to assumed shape (as
  item 7 of the [remediation order](fortran_best_practices.md#recommended-remediation-order)
  once suggested).
- **No data in modules** ([Coding Standards](coding_standards.md)), so
  that two calculations in one process cannot see each other.

### Phase 1: errors returned, not stopped

A routine that can fail returns, instead of stopping, with a last argument
`ifail`: 0 on success, or the kind of the failure, whose message it has
written (unit 6, as before, or `error_unit` for a check the code should
never fail). The kinds are the `fail_*` constants of `failure_m` (input,
convergence, not_found, construction, numerical, io, internal), one
exception class each in the wrappers. Every caller tests `ifail` right
after the call and returns with it unchanged, so the failure reaches
`oncvpsp_main`, which ends with `error stop 1` as before: the output and
the exit status do not change.

- `ifail` is a new argument, not the existing `ierr`. The solvers'
  `ierr` (`lschfb`, `ldiracfb`, `lschvkbb`, `scatterstate`, ...) carries
  codes that some callers accept and go on from (a state that did not
  converge, 100; a ghost below the node count, 200); folding the failures
  into it would let those callers carry on after one. `ifail` has one
  meaning: the calculation cannot go on.
- A procedure with optional arguments gets `ifail` last, passed by
  keyword (`ifail=ifail`), so that it cannot land in an optional slot.
- A reader that fails after opening its file closes it first.
- No library routine stops, not even on a check of a condition the code
  has just established, or that `check_data` has already validated: a bug
  must not take a wrapper's interpreter down either. Those checks return
  `fail_internal`; they are `exc_dispatch`'s `iexc`, `sbf_rc_der`'s and
  `sbfq_rc_der`'s `llin`, `vploc_fit`'s `lpopt`, the segments of
  `pw_alloc`, `pw_product_integral` and `outwf_setup`'s sources,
  `add_zero_r`'s sizes, the `ug_*` node counts, grid positions and
  brackets, `diagout`'s bracket, `eresid`'s q grid, `ldiracfs`'s and the
  bound-state solvers' arguments, and `set_model`'s radius. The
  iteration limits (`root_brent`, `root_fmin`, `lm_root`,
  `quad_gauss_legendre`, `bs_box_energy`, the origin series of
  `sr_outward` and `dirac_outward`) return `fail_convergence`, the LAPACK
  failures (`set_model`'s `dhseqr` and `dgglse`, `ws_basis_m`'s `dsyev`)
  `fail_numerical`. So the XC routines, the uniform-grid quadratures and
  the radial integrators take `ifail` too, and the callbacks of
  `bs_solve` and `bs_box_energy` are subroutines with `ifail`.
- A build without libxc refuses a libxc functional in `check_data`
  (`have_libxc`), so `exc_libxc`'s stub, which returns `fail_input`, is
  not reached; with libxc, `check_data` refuses a functional libxc cannot
  evaluate (`exc_libxc_check`), so `exc_libxc_init`'s `fail_input` is not
  reached either.
- The statuses that pass between modules, of the radial solvers
  (`ierr_box_state`, `ierr_iterations`, `ierr_not_converged`,
  `ierr_ghost`), of `psatom` (`ierr_ps_scf`) and of Libxc
  (`ierr_lxc_*`), are named in `failure_m`, with the values of 4.0.1
  that some messages print.

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
      `modcore_budget_table`. What the searches call back reports through
      `budget_t%ifail`: `set_model`'s three checks (a radius beyond the
      table, which covers the searched interval; `dhseqr` on the tail's
      companion matrix; `dgglse` on constraints of full rank by
      construction), its XC and `der2exc`, after which the search's
      functions return 0 so that the searches end at once.
- [x] The tests and the writers: `run_diag`, `run_plot`, `run_ghosts`,
      `run_phsft`, `phsft_m`, `outmesh_m`, `linout` and `upfout`. The
      modules of `oncvpsp_confine` alone (`confine_ws_m`, `ws_basis_m`,
      `confine_wall_m`) return `ifail` as the library does, and the
      program stops on it: their phase searches call back into `phsft_m`.
- [x] `oncvpsp_main`'s own stops on the results (the command line and the
      input file stay): they moved with phase 3 into `oncvpsp_run`.
- [x] The kinds of failure (`failure_m`), every `ifail=1` made the kind of
      its failure, and the checks that stayed `error stop` returned
      (above), with `write_reference_results`' check of the reference
      atom's iterations.

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
      `write_help`). The messages of the checks that stay `error stop`
      (phase 1) go to stderr (`error_unit`), beside the stop's own: they
      are bugs, not output, and a library's unit 6 is its host's.
- [x] The libxc functionals' description goes to `iout`, once per run,
      before the all-electron atom (`exc_libxc_info`, from `oncvpsp_run`,
      and from `oncvpsp_confine` after it loads the atom), not to unit 6
      once per process from the first `exc_libxc` call.
- [x] A caller that wants no output opens `/dev/null` (or `NUL`) and
      passes its unit, as `tests/unit/test_run.f90` does.

### Phase 3: a driver split in three

`oncvpsp_main_m.F90` (895 lines) is split into:

- [x] `oncvpsp_run(irel,istage,inp,ae_filename,ps_filename,skip_ld,psp,iout,ifail)`
      (`oncvpsp_run_m`): the calculation, from the input as read (an
      `input_t`, which it completes as the writers echo it) to the
      results (a `results_t`, `psp`, of `results_m`), with the tests' text on
      `iout` and `ifail`. It ends with `ifail=1` where the driver
      stopped (the input checks before `check_data`, the reference atom).
      The staging files of `_ae.x`/`_ps.x` and `oncvpsp_ld.x` stay in it,
      by name, blank for none: a library caller passes `istage=0` and
      blank names.
- [x] the readers, which fill the `input_t` (`read_input_text`,
      `read_input_toml`, already separate);
- [x] the writers (`linout`, `upfout`, `psmlout`), called from `psp`
      by `write_psp_files(inp,psp,iout,ifail)` (`psp_files_m`), which
      writes the files that `psfile` asks for;
- [x] `oncvpsp_main`, left as the command line, the readers,
      `oncvpsp_run` and `write_psp_files`.
- [x] `results_t` holds the pseudopotential (`psp%pot`, a `pseudo_t`
      of `pseudo_m`: the atom, the mesh, the projectors, the local and
      semi-local potentials, the valence and model-core charges), what
      the writers take besides (the all-electron charges and
      eigenvalues, the PSML potentials, the orbitals), the input values
      as read for the echoes (`psp%echo`, an `echo_t`) and the tests
      (`psp%tests`).
- [x] The test results in `psp` as arrays, not only as text: a
      `tests_t` (`test_results_m`) of one type per test, `logder_t` (the
      log derivatives), `ghosts_t`, `configs_t` (the test
      configurations), `convergence_t` (the crossings of the error
      levels, the curves, their tails and the optimization's
      diagnostics), `weighted_cut_t` (the model core budget's weighted
      cutoff) and `vloc_cut_t`, plus the B matrices' eigenvalues
      (`psp%tests%bbeig`). `diagout` writes the DIAGNOSTICS block from
      it. `tests/unit/test_run.f90` checks
      them against that block, which prints them exactly (nr, sr, fr,
      and the weighted cutoff). The tables of `run_diag` (the
      eigenvalues and log derivatives at rc of each projector's energy)
      and the plot data of `run_plot` stay text only.

The results whose size is known only during the run (the number of
projectors, the mesh) are allocatable components of the types, so the
Fortran caller passes no sizes. The C interface (phase 5) cannot take a
derived type with allocatable components: it holds them behind a handle,
with a struct of the sizes and calls for the arrays.

The speed of light was a module variable (`relativity_m`'s `c_light`,
set from the input by `set_relativity`), so two runs at once in threads
with different values used one. It is an argument, `c_light`, of the
solvers and the relativistic exchange and of the routines above them,
beside `irel` or `iexc`: `oncvpsp_run` takes it from the input
(`c_light_of`), `oncvpsp_ld.x` and `oncvpsp_confine.x` from the files
(`load_pair`), and the types that carry `irel` or `iexc` to their
procedures carry it too (`phsft_t`, `budget_t`).

The weights of the log mesh's stencils and cell rule (`logmesh_m`) were
computed on first use, a first computation that raced in threads; they
are parameters, rational numbers that `scripts/logmesh_weights.py`
rounds once to double precision in `src/logmesh_weights.inc`, and the
test `unit_logmesh_weights` checks that the file is what it writes.

The root finds of `oncvpsp_confine`'s sphere scan (`confine_ws_m`) took
their functions' data (the solver, the radius, the target phase, the
poles) from module variables; the functions are internal procedures of
the routines that search, which pass the data down.

So the library keeps no state between or across runs, and runs can share
a process in threads. `tests/unit/test_run_threads.f90` (`unit_run_threads`)
runs five cases, nr, sr and fr with two values of `c_light` (and, with
Libxc, its PBE in icmod 4 and 6), each alone and then all of them twice
over in OpenMP threads at once, and compares the results and the text
output bit for bit; before the speed of light was an argument, 85 of its
checks failed. The test reads the inputs before the threads start and
gives each run a scratch file: gfortran does not open one file, even
`/dev/null`, on two units at once, which failed the Debug build at
random. A threaded BLAS (MKL, OpenBLAS) uses all its threads outside the
parallel region and one inside it, and its sums then differ in their last
bits (07_N_opt's budget, through LAPACK), so the test holds it to one
thread. Only the test is compiled with OpenMP. The procedures
are compiled recursive in every build type (`-frecursive`, and
`-recursive`, `-Mrecursive`; `cmake/OncvpspCompilerOptions.cmake`):
gfortran makes a local array over 64 kB static otherwise (there is none
now), and its run-time check of recursion in the Debug build stops two
threads in one procedure.

The cache of all-electron atoms that the test configurations reuse is an
`ae_cache_t` (`config_ae_m`) that `oncvpsp_run` owns and passes down, so
a second run starts with an empty one; `tests/unit/test_run.f90` checks
that run bit for bit.

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

### Phase 4b: the stages of the run

`oncvpsp_run` is split into stages that a caller can run on their own
and repeat, each from the results of the ones before, which it does not
change; `oncvpsp_run` calls them in order, with the same text output. The
order follows what the optimizer would reuse first.

- [x] `complete_input(inp,cnf,iout,ifail)` (`oncvpsp_run_m`): the input as
      the stages take it, the configurations, `check_data`, and the
      default core radii.
- [x] `solve_ae` (`ae_state_m`): the all-electron atom, solved or read
      (istage 2), in an `ae_state_t` that holds what the all-electron
      data file holds, and the cache of the test configurations' atoms.
- [x] `optimize_channel` (`channel_m`): one l's channel (both kappas)
      from the `ae_state_t` and the channel's inputs alone, in a
      `channel_t`: the projectors' states, the optimized pseudo-waves,
      the semi-local potentials and the convergence curves.
      `tests/unit/test_stages.f90` solves the atom once, runs every
      channel twice from it, and checks them against each other and the
      whole run bit for bit (nr, sr, fr, and a local channel).
- [x] `run_vkb` (`run_vkb_m`) split into `local_potential`, in a
      `local_t` (its inside pieces, on the mesh and as pieces), and
      `channel_vkb`, one channel's projectors from its `channel_t`
      arrays and the `local_t`: with `lloc=4` a channel changed alone
      changes only its own projectors. `run_vkb` calls them;
      `test_stages` checks them against the whole run.
- [ ] The coupled block, from the valence densities and the model core to
      the tests, re-runnable on the earlier stages with another model
      core (icmod, fcfact, rcfact).

### Phase 5: the C interface

A `bind(C)` derived type has C's layout, but its components can only be
scalars of interoperable kinds, arrays of fixed size, `c_ptr`,
`c_funptr` and other `bind(C)` types: no allocatable or pointer
components, no `character(len>1)`, no `extends` or type-bound
procedures. `results_t`, `input_t` and the test types are made of
allocatable arrays, so they cannot cross as structs. On the other side,
Julia's `ccall` (and Python's ctypes and cffi) map an `isbits` struct of
`Cint`, `Cdouble`, `NTuple` and `Ptr` fields to a C struct one to one,
but check nothing: a struct that differs from the Fortran type after a
field is added reads garbage, with no error. The C descriptors of
Fortran 2018 (`CFI_cdesc_t`) pass allocatable and assumed-shape arrays
as arguments, not as components, and their layout is the compiler's
(`ISO_Fortran_binding.h` differs between gfortran, ifx and nvfortran),
so they stay out of the interface.

So the interface is a hybrid:

- the results stay Fortran's, behind an opaque handle: a `type(c_ptr)`
  to a run object allocated by the library, which holds the `results_t`
  and the completed `input_t`, and is freed by the caller;
- one small `bind(C)` struct, `oncv_dims`, carries what a caller needs
  to allocate: the sizes (mesh points, projectors per l, lmax, test
  configurations, log-derivative energies, …) and an interface version,
  so a wrapper refuses a layout it does not know rather than misread it;
- getters copy each array into memory the caller allocates, with its
  size checked against the dimensions, and return an `ierr`. No pointer
  into the run object is handed out, so nothing outlives it.

- [ ] A module `oncvpsp_c` of `bind(C, name="oncv_…")` routines over
      `oncvpsp_run`, using `iso_c_binding` kinds (`c_int`, `c_double`,
      checked equal to `dp`; `c_bool` for logicals).
- [ ] The input as the TOML text, a `c_char` array with its length,
      read by `read_input_toml`: the wrappers take a dictionary in the
      TOML layout (phase 6), and the parser already checks it. Setters
      on an input handle only if a caller needs them.
- [ ] `oncv_run(text, n, iout_path, h, ierr)`, `oncv_free(h)`,
      `oncv_dims(h, dims)` and one `oncv_get_…(h, buf, n, ierr)` per
      array (the potentials, projectors and energies, the charges, the
      tests of `tests_t`), and `oncv_write_files` over
      `write_psp_files`. The text output goes to a file the caller
      names, `/dev/null` for none; the error messages to stderr.
- [ ] Stage entry points beside `oncv_run`, as handles on the same
      pattern (dims and getters), over the stage splits of `oncvpsp_run`:
      `oncv_ae(text)`, an all-electron handle that `oncv_run` takes in
      place of `ae_filename`; `oncv_optimize_channel(ae, ..., l)`, with
      the pseudo-waves, the convergence curves (ecurve, wcurve) and the
      node check; the per-channel projectors; and the coupled block (model
      core to the tests) again on a run handle's earlier stages with
      another model core.
- [ ] Every number the optimizer's metrics read (oncvpsp-opt's
      `diagnostics.py`, `cutoffs.py`, `log_derivative.py`, the ghosts,
      configurations, `weighted_cut_t`, `vloc_cut_t`) from `tests_t` and
      `psp%pot`, with no text; check that the model-core checks
      (`core_truncation.py`, `model_core.py`) need only the charges and
      potentials. `run_diag`'s tables may stay text only.
- [ ] The text output optionally into a buffer the caller reads, beside
      a file or `/dev/null`, so the text parser can check the arrays
      during the switch-over without temporary files.
- [ ] Later, `oncv_confine(ae, run, radii, walls, channels)` returning the
      fields of `confine_summary_t`, in place of a third process.
- [ ] The failure's context, beside the integer `ifail`, which stays the
      code the bindings read: a small `bind(C)`-compatible record
      {kind, msg_id, stage, l, ikap, iproj, ierr} and a fixed-length message,
      filled at the failure site and kept in the run object, with a
      getter on the handle; ierr is the solver's status (`ierr_ghost`,
      `ierr_not_converged`, ...), so that failures of one message stay
      apart without its text. The bindings raise from it exceptions that
      map onto pyoncvpsp's Stage, Scope and Attribution. Not a
      polymorphic `error_t`: the required `ifail` keeps the compiler
      checking every call.
- [ ] A rule in `scripts/check_conventions.py` that flags any statement
      other than `return` right after `if(ifail/=0)`, to catch a caller
      that does not pass a failure on.
- [ ] A C header for it, and a shared library target (`liboncvpsp`) in
      CMake.
- [ ] A C test in `tests/` that runs one of the reference inputs through
      the interface and compares with the reference output.
- [ ] Check that two runs in one process give the same results.
      `tests/unit/test_run.f90` does so from Fortran, bit for bit, and
      `tests/unit/test_run_threads.f90` for runs in threads; the C test
      repeats it through the interface, two handles at once.

### Phase 6: wrappers

- [ ] Python: ctypes or cffi over `liboncvpsp`, built with
      scikit-build-core, returning NumPy arrays; tests against the
      references.
- [ ] Julia: `ccall` over the same library, `oncv_dims` as an `isbits`
      struct, the arrays allocated in Julia and filled by the getters;
      later packaged as a JLL with BinaryBuilder.
- [ ] Both check the interface version of `oncv_dims`, and free the
      handle with a finalizer on the result object.
- [ ] Both take the input as a dictionary or keyword arguments, in the
      layout of the TOML input.
- [ ] A user page on the wrappers, in the `nav` of `mkdocs.yml`.
