# Known issues

Suspected bugs found while auditing dummy-argument `intent`s across `src/`.
None of these have been fixed yet. Line numbers refer to the commit that
added this file.

## Likely wrong results

- **`der2exc.f90:57`** — the 4th point of the finite-difference stencil uses
  `rhotot - 2*hh*rho` (same as the 1st point) instead of `rhotot + 2*hh*rho`,
  so the ±2h terms cancel and the derivative estimate is wrong.
- **`ldiracfb.f90:217`** — the rescaling loop sets `fup(ii)=sc*gup(ii)`;
  almost certainly meant `sc*fup(ii)`. The small-component derivative is
  overwritten with the large-component one.
- **`run_config.f90:152`, `run_config_r.f90:152`** — the "maximally ionized
  intermediate state" loop runs over `kk` but writes
  `fat(ii,2)=dmin1(fat(ii,1),fat(ii,3))` with a stale `ii`, so only one element
  is ever set.
- **`modcore3.f90:170,180,215,...`** — the `rcross` search loops call
  `gg1cc(gg,yy)` with a `yy` left over from an earlier pass instead of
  `rr(ii)/xt(2)`, so `rcross`/`ircross` are computed from the wrong value.
- **`lschpse.f90:78`** — the `emin` search uses `rr(ii)+0.5*sls/rr(ii)**2`;
  `vv(ii)` was probably intended.
- **`ldiracfb.f90`, `ldiracfs.f90`** — the series start near r=0 uses
  `cci*(ee+vv)` and `(2c^2 - ee - vv)` while the predictor/corrector uses
  `-cci*(ee-vv)` and `(2c^2 + ee - vv)`. The sign on `ee` is inconsistent
  (small effect, near the origin only).
- **`run_vkb_r.f90:78-83`** — `l1=0` is set immediately before
  `else if (l1==1)`, so that branch is dead. Probably meant `lloc==0`.
- **`run_diag_sr_so_r.f90:142`** — `npr=0` for `l==lloc` is set after `nvkbt`
  has already been computed from `npr`, so the local channel still uses its
  projectors. Unlike `run_diag`/`run_diag_r` it also doesn't stop on an
  `ldiracfb` error.

## Silent failures / missing error handling

- **`psatom.f90:216`** — the non-convergence check is `if(it==400 ...)` inside
  `do it=1,100`, so it can never fire; `ierr=101` is never set and an
  unconverged run exits silently. `psatom_r` uses `it==100` correctly.
- **`psatom_r.f90:155`** — on an `lschvkbb` error, `exit` only leaves the
  `ikap` loop; the next `ii` iteration resets `ierr=0`, masking the error.
- **`psatom.f90`, `psatom_r.f90`** — if every occupancy `fa` is zero, `ierr`
  is tested without having been set in the routine.
- **`qroots.f90:83,92`** — `found_root=.true.` is set unconditionally after
  the bisection loop, so the "failed to find root" check is dead code. The
  routine also assumes `nroot>=4`.
- **`run_phsft.f90`, `run_phsft_r.f90`** — if `rxpsh >= rr(mmax)`, `irphs` is
  never set; the radius check prints an error but execution continues.
- **`input_toml_m.f90`** — `read_input_toml`/`read_input_toml_r` declare
  `fcfact_min/max/step` and `rcfact_min/max/step` `intent(out)` but never set
  them (gfortran warns). The text reader sets defaults
  (1.5/6.0/0.5, 1.0/1.9/0.1), so TOML input with `icmod=4` passes undefined
  values. Also, required-key `get_value` calls have no `stat`/default, so a
  missing key leaves the output undefined without an error, and the `icmod`
  check accepts 0–4 only, while the text reader accepts `icmod>=4`.

## Invalid memory operations

- **`m_psmlout.f90:434,969`** — `deallocate(chcore)` after
  `chcore => rhomod(:,1)` deallocates through a pointer to dummy-argument
  memory; this is invalid and likely to fail at runtime when `icmod>=1`.
- **`upfout.f90:103`, `upfout_r.f90:111`** — the backward `do jj=...,1,-1`
  search can complete without `exit` (e.g. an all-zero orbital), leaving
  `jj=0` and reading `rr(0)`.
- **`sbf_rc_der.f90`** — guards `llin>5`, but `llin=5` still writes
  `sbfad(1,11)`, outside `sbfad(5,10)`.
- **`sbf_basis.f90`** — `tder(6)` is assigned from `orbasis_der(:,jj)` of
  length `nconmx` (shape mismatch unless `nconmx==6`), and the `nconmx` loop
  reads `sbfder(5)` out of bounds for `nconmx>5`. (The routine is not called
  anywhere.)
- **`sbf_basis_con.f90`** — likewise requires `ncon_in<=5` and
  `ncon>=ncon_in+iprj-1`.
- **`fpovlp.f90:53`, `fpovlp_r.f90:51`** — use `rr(101)` while `rr` is
  dimensioned `rr(nn)`; out of bounds if `nn<101`.
- **`modcore.f90`** — the derivative loop starts at `ii=ircc-25+3*jj`, out of
  bounds for `ircc` below ~20.
- **`lschpsbar.f90`** — the barrier start writes `uu(mbar:mbar+4)` without
  checking against `mmax`.
- **`run_ghosts.f90`** — `hmev(jj)` can be read past `nbas` when ghost states
  push `jj` up, and `dsyev` can be called with `nbas=0`.
- **`optimize.f90`** — `nnull==0` indexes `eresiddot_ev(1)` out of bounds.
- **`linout_r.f90`** — local `zero(20)` is indexed up to `npr_sr`/`npr_so`
  (up to `2*nproj`), with no bounds check.
- **`vrel.f90`** — `vr(ii-1)`/`vr(ii+1)` could go out of bounds at the array
  edges (prevented in practice by the rc window).

## Undefined values

- **`lschpsbar.f90:94`** — `imin` is only set inside `if(emin==0)` (and only
  if some point lies below `emax`); otherwise `imin=max(imin,2)` reads an
  undefined value.
- **`modcore3.f90`** — `ircross` is undefined if no crossing is found;
  `icmod>=5` would reach the final section with `xx(:,1)` unset.
- **`excwig.f90`** — `vxc`/`exc` are not set where `rho<1e-15`, so `vout`'s
  `vxc` and local `exca`/`vxcd` hold stale or garbage values there when
  `iexc==1`.
- **`vkbphsft.f90`** — `dnpi` is undefined if no scan energy falls within
  `depsh/2` of `ep`, yet it is still added to all of `pshp`.
- **`sr_so_r.f90`** — `ipk` can be used uninitialized if a projector column
  is all zeros.
- **`linout.f90`, `linout_r.f90`** — `ixc_abinit` is never set for `iexc>4`
  but is written out.

## Interface and consistency issues

- **`oncvpsp_r.F90:484`** — calls `wellstate_r` with 13 arguments; the
  routine takes 14 (the last is `srel`). Harmless only because `srel` is
  unused.
- **`upfout.f90:112`, `upfout_r.f90:121`** — `nrl` is enlarged in place
  (`nrl = 1 + uurcut/drl`, implicit real→integer truncation). Because
  `psmlout` is called afterwards when `psfile=='both'`, PSML output gets a
  different `nrl` in "both" mode than in PSML-only mode.
- **`m_psmlout.f90:239,723`** — no-op `lmax = lmax` forces `lmax` to be
  `intent(inout)`; removing it would allow `intent(in)`.
- **`m_psmlout.f90`** — `write_psml_item` declares `f` optional but uses it
  unconditionally (every caller passes it).
- **`functionals.F90`** — `functional_get_tau` uses `n(1)`/`t(1)` instead of
  `n(is)`, and passes literal `1` instead of `xc_one` for libxc ≤ 4 (routine
  is unused).

## Resource leaks and cleanup

- **`vrel.f90:70`** — returns early when `.not. srel` without deallocating
  `cf`, `dv`, `fr`, `frp`.
- **`run_diag_sr_so_r.f90`** — never deallocates `ur`, `urp`, `evkbt`, `vkbt`.
- **`run_config.f90`, `run_config_r.f90`** — `psatom` error-return paths
  don't deallocate `uu`/`up`.
- **`modcore.f90`** — `vxct` is allocated but never deallocated.

## Minor / cosmetic

- `dpnint.f90` — the `tt(1)==0` extrapolation block sits inside the `jj` loop
  and reads `ss(2:9)` before they're computed on early iterations (final
  result is correct, just wasted work).
- `vrel.f90` — two consecutive `cycle` tests duplicate each other (the
  `2.5*rc` bound is ineffective).
- `lschpb.f90` — declares a local `real(dp) :: exp`, shadowing the intrinsic
  name (still resolves to the intrinsic).
- `aeo.f90` — `aii` lacks `implicit none`, unlike its siblings.
- Unused arguments/locals: `nn` in `lschkb`/`lschps`; `zz` local in
  `lschkb`/`lschps`; `icmod`, `rcfact` in `modcore`/`modcore2`; `srel` in
  `run_config_r`; `drdum`/`dqdum` in `eresid`.
- `src/oncpsp_nr.f90` is not compiled and appears to be a stale duplicate of
  `oncvpsp_nr.F90`.
