## Merging the scalar- and fully-relativistic code

The non-relativistic (nr) and scalar-relativistic (sr) programs already
share one code path: `oncvpsp_nr.F90` includes `oncvpsp.F90`, and `srel`
selects the Schroedinger equation with or without the scalar-relativistic
terms. The fully-relativistic (fr) program is a copy: a `foo_r` for most
routines `foo`, and an `_r` entry point in five modules, about 6000 lines
in all. This is the plan for merging the two, and the record of how far it
has got. Tick a box when its commit is in.

### Rules

- **Output is unchanged.** Every step passes `ctest -j8` (the reference
  comparisons allow no numerical difference), `scripts/check_interfaces.sh`,
  the legacy `make` check and the comparison of the split programs
  (`_ae.x` + `_ps.x`) with the single ones. The reference comparisons
  ignore text, so the outputs are also compared line for line with those
  of the code before the merge (53bd8ac), apart from the build line and
  dates, as are the help, version and argument-error messages of every
  executable. One commit per step.
- **Little preprocessor.** No new `#if` blocks. Which equation is solved
  is an integer argument, `irel` (0 nr, 1 sr, 2 fr), tested at run time.
- **Interfaces that C could call.** A merged or new routine takes integer
  and `real(dp)` scalars and explicit-shape arrays with their sizes passed
  in: no logical, optional, assumed-shape, allocatable or pointer dummy
  arguments. `irel` replaces the logical `srel` in each routine as it is
  merged. The character arguments already there (`atsym`, `psfile`) stay
  for now.
- **The fr layout throughout.** Arrays indexed by l get a last index for
  kappa of size `nkap` (1 nr/sr, 2 fr): `ikap=1` is j=l+1/2, `ikap=2`
  j=l-1/2, as the fr code has them now.
- **Differences kept.** Where sr and fr differ for no evident reason
  (listed under [Divergences](#divergences)), the merged code keeps both
  behind `irel`, so the references don't change. Unifying them is a
  separate decision.

### Phase 0: tests

The fr path had 3 reference tests to the sr path's 28.

- [x] Add fr references for the inputs covering what the 3 did not:
      no core (03_Li), lloc<4 (07_N_dl), iexc=1 (08_O_wigner), icmod 1
      (14_Si), 2 (79_Au_lxc), 3 (34_Se), 4 with a grid (40_Zr_grid) and 6
      (14_Si_bc), rxpsh (14_Si_lder), lpopt 1 (20_Ca), a single projector
      (38_Sr_lxc), three projectors (73_Ta), ghosts (60_Nd_GHOST); nr
      references for 03_Li, 07_N_dl and 14_Si_lder; TOML input tests for
      the new fr inputs

### Phase 1: one array layout

- [x] Give the per-l and per-level arrays of `oncvpsp.F90` the kappa
      index (`ea`, `ep`, `epa`, `vp`, `vkb`, `evkb`, `pswf`, `uua`,
      `uupsa`, `cvgplt`, `ecurve`, `wcurve`, `vr`). With `nkap=1` an sr
      routine's explicit-shape dummy takes the whole array, by sequence
      association, until the routine is merged

### Phase 2: choosing the equation in one place

- [x] Bound state (`boundstate`): one routine calling `lschfb` or
      `ldiracfb`, on `uu(mmax,2)`, `up(mmax,2)` (callers that need the
      large component normalized still call `renorm_r`)
- [x] Scattering state (`scatstate`): the same for `lschfs` and `ldiracfs`
- [x] All-electron atom (`aeatom`): the same for `sratom` and `relatom`

### Phase 3: merging the pairs

- [x] `read_input_text` / `read_input_text_r`
- [x] `read_input_toml` / `read_input_toml_r`
- [x] `run_config_ae` / `run_config_ae_r`, and their warnings
- [x] `write_ae_file`, `read_ae_file` / their `_r`
- [x] `fphsft` / `fphsft_r`
- [x] `wellstate` / `wellstate_r`
- [x] `run_config_ps` / `run_config_ps_r`
- [x] `psatom` / `psatom_r`
- [x] `run_phsft` / `run_phsft_r`
- [x] `run_diag` / `run_diag_r`
- [x] `run_plot` / `run_plot_r`
- [x] `run_vkb` / `run_vkb_r`
- [x] The projector loop of the driver (from `npa` to `run_optimize`)
      and the pseudo-wave-function loop after `run_vkb`: the same text in
      both drivers, for phase 4 to keep one copy

### Phase 4: one driver

- [x] Move the body of `oncvpsp.F90` into a subroutine taking `irel` and
      the stage (whole, all-electron only, pseudopotential only), with each
      executable a short program calling it; drop `oncvpsp_r.F90` and the
      `#if`s of the driver: `oncvpsp_main(irel,istage)` in
      `oncvpsp_main.F90`, an external procedure so that
      `check_interfaces.sh` checks its calls; the nine program files call
      it, and only the `WITH_TOML` `#if`s are left
- [x] The fr-only steps (`sr_so_r`, `run_diag_sr_so_r`, the second ghost
      test, the averaged unscreening) under `irel==2`: `oncvpsp.F90` and
      `oncvpsp_r.F90` differ only in their headers and the setting of
      `irel`

### Phase 5: output

The psp8, UPF and PSML files differ in form for fr. Where the difference
is a few sections, the two writers merge as the rest did; otherwise they
stay separate and what they share moves to common routines.

- [x] `gnu_script` / `gnu_script_r`
- [x] The input echo of `linout`, `linout_r`, `upfout`, `upfout_r`:
      `write_input_echo` in `input_echo_m`
- [x] `linout` / `linout_r`, merged: they differ only in the projector
      section and a few header lines
- [ ] The shared parts of `upfout` / `upfout_r`
- [ ] The shared parts of `psmlout` / `psmlout_r`

### Phase 6: documentation

- [ ] `coding_standards.md`: replace the `foo`/`foo_r` pairing with
      `irel`/`nkap` and the interface rules above
- [ ] `program_overview.md`, `relativistic.md` (its list of fr
      routines), `coding_standard_violations.md`, `known_issues.md` (the
      divergences, if kept); `fortran_best_practices.md` cites the `_r`
      files as they were

### Divergences

Differences between `oncvpsp.F90` and `oncvpsp_r.F90` with no evident
physical reason, kept behind `irel` until decided:

- sr makes the number of points of the linear output mesh even ("PWSCF
  wants an even number"); fr does not
- the first point of the default log mesh, `rr1`, is computed differently
- `lschvkbb` searches 0.75 to 1.25 of the energy in sr, 0.90 to 1.10 in fr
- when `wellstate` finds the well state with more nodes than asked for,
  fr keeps its n in `npa` and sr the n asked for; sr stops if `lschfb`
  fails in the well, fr goes on
- the error-message prefixes and some numeric formats (`e` and `d`);
  `run_vkb` prints kappa with the B-matrix Hermiticity errors in fr only
- the plot data of `run_plot`: fr writes every point of the potentials,
  starts the wave functions a step further out, labels j=l+1/2 with -l,
  and averages the convergence profiles over kappa
- the test-configuration results table of `run_config_ps`: sr leaves out
  empty valence levels, fr lists every level once for each kappa

### Fixed on the way

Behaviour changed by a merge, in paths that no test reaches:

- `run_config_ae_r` searched for an empty level's j=l+1/2 bound state with
  kappa=l+1, which `ldiracfb` rejects, so the search always failed; the
  merged `run_config_ae` uses -(l+1). It also records the error code of a
  failed sr search in `ierrcnf`, as fr did.
- `run_vkb` declared `pswf(mmax,mxprj,5)` for an actual of 4 columns of
  l; the merged routine declares 4. Its fatal overlap-matrix message is
  fr's wording ("sovl matrix eigenvalue ERROR") in sr too.
