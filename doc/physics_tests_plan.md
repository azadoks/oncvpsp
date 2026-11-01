## Unit tests of the physics

The unit tests in `tests/unit` mostly check the code against itself: a
stage run alone is bit for bit the stage of a whole run, a solver's answer
does not depend on its starting guess, an analytic derivative equals its
finite difference. The reference outputs check the numbers against those of
an earlier build. Neither says whether the numbers are right. This plan adds
tests that do: each one states a physical invariant, a relation the exact
answer must satisfy, and checks that it holds across a variety of cases
(atoms light and heavy, nr, sr and fr, LDA and GGA, configurations neutral
and ionized). A test that fails says which law was broken, in which case,
by how much.

### Rules

- **One invariant per check, stated in the test's header** as an equation,
  with the cases it runs over and the tolerance and where that comes from:
  the mesh's discretization error (from the Adams order and `amesh`, as
  `test_adams` measures it), the solvers' convergence tolerances, or the
  round-off of a finite difference. A tolerance is never set from the
  observed error alone.
- **The cases are a table in the test**, not one atom: at least a light and
  a heavy atom, each `irel` the invariant applies to, and a GGA beside the
  LDA where the functional enters.
- **Analytic inputs where possible.** The solver and functional tests use
  model potentials and densities built in the test; the all-electron tests
  call `aeatom` with the configuration in the test; only the
  pseudopotential tests read an input file.
- **Smaller units where the code hides the quantity.** If an invariant needs
  a quantity a routine computes but does not return (an orbital, a term of
  the potential), the routine is split so that it can be called alone, in
  a module, with no change to any output (the reference tests are the
  check). Each such split is its own commit.
- The existing tests that already check one of these invariants (hydrogen
  in `test_solver_path`, the Adams orders in `test_adams`, the functional
  derivatives in `test_fxc`) are not repeated; a new test extends their
  cases where they are narrow.

### Step 1: the all-electron atom (`test_ae_physics`)

`aeatom` on a mesh made by `make_log_mesh`, for He, O, Si, Cu and Au
(neutral) and O+, Si+ and Cu(0.5+), nr, sr and fr, the VWN LDA (iexc 6)
and PBE (iexc 4).

1. **Charge.** `int r**2 rho dr = sum f` over all states, and
   `int r**2 rhoc dr` = the core's `sum f`.
2. **Coulomb limits of the potential.** `r vfull(r) -> -Z` as `r -> 0`, and
   `r vfull(r) -> -(Z - N)` where the density has died off (the cations
   test the Hartree term's total charge; for the neutral atoms the limit
   is 0).
3. **Virial theorem (nr).** For a Kohn-Sham ground state of any density
   functional, the scaling `n -> g**3 n(g r)` makes `E` stationary at
   `g = 1`, giving
   `2 T_s + E_ne + E_H = int n r dv_xc/dr d3r`
   (the right side is `-E_x` for exchange alone). The terms on the left are
   `eterm`; the right needs `v_xc` of the density (`exc_dispatch`).
   It checks the kinetic, nuclear, Hartree and exchange-correlation
   energies against each other.
4. **Janak's theorem.** `dE/df_i = eps_i`: the central difference of `etot`
   in the occupation of a valence state (`f +- 0.01`) equals its
   eigenvalue, for the highest state and a deeper valence state of each
   atom, nr and sr. For fr, the share of `f` between the two j levels is
   fixed (`kap_j_share`), so the derivative is the shares' average of the
   two eigenvalues.
5. **Hellmann-Feynman.** With an external potential `lambda W`
   (`aeatom`'s `vext`, a smooth well), `dE/dlambda = int n W d3r`.
6. **The non-relativistic limit.** At `c` 10 and 20 times the physical
   value, sr and fr give the nr eigenvalues and total energy, and the
   difference falls by 4 when `c` doubles (it is `O(1/c**2)`).
7. **Spin-orbit order.** For `l > 0`, fr's `j = l - 1/2` level lies below
   its `j = l + 1/2` level, and the splitting falls by 4 when `c` doubles.
8. **Shell order.** Within each `l`, the eigenvalues rise with `n`.

### Step 2: the bound-state solvers on model potentials (`test_bound_physics`)

`lschfb` and `ldiracfb` on potentials whose spectra are known, which
`test_solver_path` uses only for hydrogen's 1s, 2p and 3d. Each quantity
is computed on the program's mesh scaled by `1/Z` and on the mesh of half
its step; the difference from the exact value falls by `2**4`, so the
Richardson value checks the physics well below the mesh's error.

1. **Hydrogen-like spectra.** `E = -Z**2/(2 n**2)` (nr) for `Z` 1, 6, 26,
   79 and `n` to 5, every `l < n`; the Dirac energies
   `E = c**2 ((1 + (Z a/(n - |k| + sqrt(k**2 - Z**2 a**2)))**2)**(-1/2) - 1)`
   for every `kappa` (point nucleus).
2. **sr is Dirac without spin-orbit.** For `l > 0` the sr energy is the
   `(2j+1)`-weighted average of the two Dirac levels to `O((Z/c)**4)`; for
   s states the sr equation is Dirac's for the large component of
   `kappa = -1`, so the energies are equal.
3. **Scaling.** `E(Z) = Z**2 E(1)` (nr) on meshes scaled by `1/Z`.
4. **Expectation values.** For hydrogen-like states,
   `<r> = (3 n**2 - l (l+1))/(2 Z)` and `<1/r> = Z/n**2` from the returned
   function.
5. **Nodes and orthogonality.** A state has `n - l - 1` nodes, and two
   states of one `l` (one `kappa`) are orthogonal.
6. **A non-Coulomb potential.** The isotropic harmonic oscillator,
   `E = omega (2 n_r + l + 3/2)`.
7. **Hellmann-Feynman for one state.** `dE/dlambda = <u|W|u>` for
   `V + lambda W`, nr and fr.

### Step 3: the exchange-correlation functionals (`test_xc_physics`)

`exc_dispatch` on uniform densities and on model atoms, and PBE's exchange
and correlation alone (`pbe_x`, `pbe_c`, split out of `excggc` for this).
The LDAs' exchange is Slater's, so their correlation is `eps_xc - eps_x`.

1. **High-density limit.** `eps_c -> A ln r_s + B` with Gell-Mann and
   Brueckner's `A = (1 - ln 2)/pi**2` (PZ, VWN, PBE's PW92).
2. **The potential is the derivative of the energy.**
   `v_xc = d(n eps_xc)/dn` over `r_s` from 0.1 to 20 for the LDAs, and
   `dE_xc/dt = int v_xc drho` on a model atom for every built-in
   functional, PBE included.
3. **Correlation.** `eps_c <= 0`, rising with `r_s`, and `r_s eps_c`
   falling (Levy and Perdew's scaling inequality for a uniform density);
   PZ's two fits join at `r_s = 1` only to their coefficients' rounding.
4. **Lieb-Oxford.** `E_xc >= -1.679 int n**(4/3)` on model atoms, and PBE's
   `1 <= F_x(s) <= 1.804`, rising.
5. **Uniform scaling of exchange.** `eps_x(g**3 n, g**4 n') = g eps_x`
   (PBE).
6. **PBE's uniform limit and linear response.** At zero gradient PBE is
   Slater plus PW92; its second-order gradient terms cancel
   (`mu = beta pi**2/3`), leaving the next order.
7. **Relativistic exchange** (iexc 7) is weaker than Slater's, by
   `O((k_F/c)**2)`.

### Step 4: the pseudopotential (`test_ps_physics`)

A whole run (`oncvpsp_run`) and its stages (`optimize_channel`), over the
inputs of `test_stages` and more (08_O nr, 14_Si sr and fr, 07_N_dl,
74_W, 60_Nd_GHOST, 03_Li fr).

1. **Norm conservation.** `int_0^rc phi**2 = int_0^rc u**2` for every
   projector's pseudo-wave and its all-electron state.
2. **Matching.** `phi = u` beyond `rc`.
3. **Reference eigenvalues.** The pseudo-atom's reference states have the
   all-electron eigenvalues (the first configuration of the configuration
   tests).
4. **Log derivatives at the projector energies.** The pseudo and
   all-electron log derivatives agree at each projector's energy.
5. **Valence charge.** `int r**2 rho_ps dr = zion`, and `rho_ps` is the
   all-electron valence density beyond `rcmax`.
6. **The ionic tail.** `r vloc(r) -> -zion` where the core has died off.
7. **Projectors are local.** `vkb(r) = 0` beyond each channel's `rc`.
8. **The model core is the core beyond `r_core`**, and is positive.
9. **Missing kinetic energy falls with the cutoff.** Each convergence
   curve is non-increasing in `q` (its integrand is `q**4 |phi(q)|**2 >= 0`).
10. **fr's two projector sets.** The scalar-relativistic and spin-orbit
    projectors recombine to the two `j` sets.

### Step 5: transferability and consistency

1. **Excitation energies.** The pseudo-atom's total energy differences
   between the test configurations follow the all-electron ones to the
   size of the reference eigenvalue errors.
2. **Mesh convergence of physical results.** The all-electron total energy
   and eigenvalues converge as `amesh -> 1` at the Adams order (extending
   `test_adams` from single states to the self-consistent atom).

### Status

- [x] Step 1 (2026-10-31, `tests/unit/test_ae_physics.f90`). All eight
  invariants hold for nr and fr: Janak to 8e-9 Ha, Hellmann-Feynman to
  4e-10 relative, the virial to 2e-9 of `T`, the charges to 5e-14, the
  ionic tail to 6e-12. What the checks found:
  - **PZ (iexc 3) is not the derivative of its energy at `rs = 1`**, where
    its potential steps (`excpzca` adjusts its constants for the
    continuity of `exc` and `vxc` there, but `vxc` still steps by a
    little): Janak's theorem misses by 2e-6 Ha for He and the virial by
    2e-7 of `T` for Si. The test uses VWN for its LDA instead.
  - **The sr total energy is not the variational energy of the sr
    eigenvalues.** The density is that of the large components, each
    normalized alone (`sratom`: `rho = sum f (u/r)**2`), and the
    scalar-relativistic equation is not derivable from a Hamiltonian, so
    `dE/df` misses the eigenvalue and `dE/dlambda` misses `<W>` at
    `O((Z/c)**2)` relative: 2e-5 Ha for O's 2s and 2p, 3e-3 Ha for Au's
    5s. nr and fr satisfy both theorems to round-off. The test bounds sr
    by `(Z/c)**2` relative. Whether to change this (a density with the
    small components, as fr has) is open: it would change every sr
    output.
  - At `c` 100 times the physical value the relativistic shifts of the
    valence levels (~1e-8 Ha) are at the eigenvalues' noise; the limits
    are taken at 10 and 20 times, with the next order `(Z/c)**2` in the
    tolerance.
  - The Coulomb limit at the nucleus holds to 2e-8 of `Z` only: `vfull`
    is the last input potential, which at the first mesh points is that
    far from the output (most for sr and fr, whose mixing filters the
    residual near the nucleus).
  The test takes 38 s (Au fr dominates).
- [x] Step 2 (2026-10-31, `tests/unit/test_bound_physics.f90`). All seven
  invariants hold, 860 checks in 0.7 s. No routine needed splitting: the
  solvers return the energies and functions the invariants are about.
  - On the program's mesh (amesh 1.006) the energies are off by up to
    1.5e-8 relative (hydrogen 5s) and the oscillator's by 5e-8 (`n_r` 2,
    `l` 3): the whole difference is the 4th-order mesh error (it falls by
    15.9 when the step halves). The Richardson values give the exact
    energies to 5e-12, the moments to 8e-11.
  - The sr s states are the Dirac `kappa = -1` states exactly (to 4e-12);
    for `l > 0` the sr energy differs from the Dirac j-average by
    `C (Z/c)**4 Z**2/n**3` with `C` the same for every `Z` (0.06 for p,
    0.003 for d, 0.0006 for f), the order the derivation gives.
  - Scaling with `Z` holds to round-off (2e-15), orthogonality to 6e-10,
    Hellmann-Feynman to 3e-10 (nr) and 3e-9 (fr, Au 2p1/2).
- [x] Step 3 (2026-10-31, `tests/unit/test_xc_physics.f90`; the split
  `pbe_x`/`pbe_c` of `excggc`, its output unchanged to the bit). All hold,
  in 0.02 s:
  - The LDAs' potentials are their energies' derivatives to 4e-13; on a
    model atom the functional derivative holds to 8e-11 for the LDAs and
    5e-10 for PBE, whose divergence form leaves the mesh's error (8e-9 at
    amesh 1.012, 17 times less at 1.006).
  - PZ's Gell-Mann-Brueckner slope is off by 2.6e-4 (its `A` is rounded
    to 0.0311), and its fits step at `r_s = 1` by 3.2e-5 Ha in `eps_xc`
    and 2.8e-5 Ha in `v_xc`: the step of step 1's finding, measured.
  - PBE's linear-response residual is `C s**2` with the `C` the next order
    of `H` and `F_x` gives (15.9 at `r_s` 0.1), to 0.2% at `s` 1e-2 (at
    1e-3 the round-off, up to 0.8%).
  - The model atoms reach at most 0.93 of the Lieb-Oxford bound (PBE);
    Wigner (1) and Hedin-Lundqvist (2) do not have the Gell-Mann-Brueckner
    limit and are not checked for it.
- [ ] Step 4
- [ ] Step 5
