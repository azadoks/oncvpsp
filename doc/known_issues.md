# Known issues

Suspected bugs in `src/`, found by auditing the code (dummy-argument
`intent`s, [Fortran Best Practices](fortran_best_practices.md) and
[Coding Standard Violations](coding_standard_violations.md)) and triaged
against the callers, the input checks in `check_data_m.f90` and the test
inputs in `tests/data/`. The issues that have been fixed are no longer
listed; their fixes and the reasons for them are in the commit messages
(up to `264a1a4`).

Priorities for new entries:

- **P1**: changes results in common use.
- **P2**: changes results with uncommon options, crashes on some compilers,
  or hides errors.
- **P3**: latent, unreachable in practice, dead code, or cosmetic.

Give each entry its file and line, the inputs that reach it, and whether
fixing it would change the reference outputs in `tests/`.

## Open

None.

## Kept on purpose

- **Unused arguments**: `icmod` and `rcfact` in `modcore`/`modcore2`.
  Removing them would change the routines' calls for programs that use the
  library. Their comments say they are not used.
- **The sr/fr divergences**: where the scalar- and fully-relativistic
  codes differed for no evident reason before they were merged, both
  behaviours are kept behind `irel`; they are listed under
  [Divergences](fr_merge_plan.md#divergences).

## Dismissed

Each of these was checked and turned out not to be a bug:

- `sbf_basis_con_m.f90` bounds: `check_data` enforces `3<=ncon<=5`, and
  `run_optimize` sets `ncon=ncon_in+iprj-1`.
- `fpovlp_m.f90` `rr(101)`: `irc` is always far above 101 on the log mesh.
- `optimize_m.f90` `nnull==0`: `check_data` guarantees `nnull>=2`.
- `linout_m.f90` `zero(20)`: `npr_sr` and `npr_so` are ≤10.
- `vrel.f90` edge reads (vrel is now `ae_delta_m`): the rc window keeps
  `ii` in the interior.
- `run_ghosts_m.f90` `hmev` past `nbas`: `hmev` has size `nbmax` and is
  zeroed. `dsyev` with N=0 is a legal quick return.
- `qroots_m.f90` `nroot>=4`: `check_data` guarantees it.
- `input_toml_m.f90` icmod range: `check_data` rejects values outside 0–4,
  6 and 7 for both readers.
- `linout_m.f90` `ixc_abinit`: `check_data` rejects `iexc>4`.
- `sratom`/`relatom` `ierr=100` unchecked: the drivers already stop with
  "ERROR all-electron reference atom not converged" when `it>=100`, and
  every other failure inside `sratom`/`relatom` stops there.
- `wellstate_m.f90` with `irel=2` doesn't check the `ldiracfb` error in
  its well search: the search retries with a new well radius when
  `ldiracfb` fails, and the outer convergence test catches a bad `et`.
  With `irel<2` it does stop on an `lschfb` error, but checking for
  `irel=2` could stop runs that now recover, so it is left as is.
- Leaks in `vrel` (now `ae_delta_m`), `run_diag`, `run_config` (now
  `config_ae_m` and `run_config_m`) and `modcore`:
  these are local allocatables, which Fortran deallocates automatically on
  return.
