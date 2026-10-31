## Smooth dependence on the inputs

The generator's results should be smooth functions of its continuous
inputs. Today they are not, in ways that are hard to see: the core radii
are rounded up to the log mesh, `qcut` to a q-grid, model-core crossovers
to mesh points; solvers stop at loose tolerances; searches return grid
points. A finite-difference derivative of anything with respect to `rc`
is 0 or huge, and an optimizer driving the generator sees a staircase.
This is the plan to remove every such discontinuity that is not
intrinsically discrete, with numerics chosen for accuracy and for
smoothness together, and the record of how it goes. Branch `smooooth`.

### Contents

- [Definitions and scope](#definitions-and-scope)
- [Where it stands](#where-it-stands)
- [Goals](#goals)
- [Numerical treatment](#numerical-treatment)
- [Inventory and decisions](#inventory-and-decisions)
- [Tests](#tests)
- [Phases](#phases)
- [Decisions](#decisions)
- [Bugs found on the way](#bugs-found-on-the-way)

### Definitions and scope

**Continuous inputs.** The core radii `rc(l)` and `rc(5)`, `ep`, `qcut`,
`debl`, `dvloc0`, `fcfact`, `rcfact`, `slack`, the `icmod=4` search grid,
`rxpsh`, `epsh1`, `epsh2`, `depsh`. **Discrete inputs**: everything
integer (`lmax`, `ncon`, `nbas`, `nproj`, `lpopt`, `icmod`, `iexc`, the
configurations).

**Discretization controls** (`amesh`, `rr1`, `rmax`, `drl`, `rlmax`,
`upf_mesh`) are not in scope for smoothness: changing them changes the
discrete problem, and what is required is convergence as they are
refined, not smoothness in them. Results must not get less accurate.

**Smooth** means: for fixed discrete inputs and discretization, every
written quantity (the nonlocal kernel per l and j block, V_loc, the model
core, the valence density, the pseudo-wave functions) and every printed
diagnostic that is used as an objective (eigenvalues and energies of the
test configurations, log derivatives, ghost energies, residual kinetic
energies, the ecut estimates) is C-infinity in each continuous input up
to a stated noise floor, except at the **truly discrete events** listed
under [Inventory and decisions](#inventory-and-decisions) (a node
appears, a ghost appears, a bound state ceases to exist, a different
local minimum becomes global). Never the projectors themselves: they are
not unique (see `scripts/kernel_diff.py`).

**Not smooth by construction, and fine.** A function built in pieces that
are matched at a radius `rc` is only C^k there (k set by the matching).
Its value at a *fixed* radius r is then only C^k in `rc` as `rc` passes
r. That is a property of the pseudopotential, not of the numerics; the
tests leave those radii out (`smooth_scan.py`).

### Where it stands

`scripts/smooth_scan.py` steps one input through 41 values (since
2026-10-22 21 to 41, by the degree of the fit) over two
log-mesh cells (radii) or a small span (others), fits each feature with a
degree-6 polynomial and reports the largest residual relative to the
feature's scale, and the ratio of the largest step to the median step
("jump": about 1 for a smooth function, infinite for a staircase). The
features are read from the UPF file written on the log mesh itself
(`upf_mesh=600`), so no interpolation sits between the generator and the
measurement. Baseline at 4352a49 (degree 6, 41 points):

| program | input | key | kernel | V_loc | core | jump |
|---|---|---|---|---|---|---|
| sr | 14_Si | `rc[0]` | 1.1e-2 | 1.0e-4 | – | inf |
| sr | 14_Si | `rc[1]` | 3.9e-3 | 1.3e-5 | – | inf |
| sr | 14_Si | `rc[2]` | 8.5e-3 | 6.3e-9 | – | inf |
| sr | 14_Si | `rcloc` (lpopt 5) | 9.7e-2 | 3.6e-3 | – | inf |
| sr | 20_Ca | `rcloc` (lpopt 1) | 9.7e-3 | 1.7e-3 | – | inf |
| sr | 14_Si | `qcut[0]` | 5.5e-3 | 4.3e-5 | – | inf |
| nr | 08_O | `rc[1]` | 6.0e-3 | 1.5e-4 | – | inf |
| fr | 52_Te | `rc[1]` | 4.2e-3 | 2.1e-5 | – | inf |
| sr | 14_Si | `fcfact` (icmod 1) | 2e-16 | 1.7e-3 | 1.8e-2 | inf |
| sr | 79_Au_lxc | `fcfact` (icmod 2) | 2e-16 | 6.4e-4 | 1.7e-2 | inf |
| sr | 34_Se | `fcfact`, `rcfact` (icmod 3) | 2e-16 | < 1e-13 | < 1e-12 | ~1 |
| sr | 40_Zr_grid | `fcfact_min` (icmod 4) | 2e-16 | 3.0e-3 | 1.1e-1 | 7e2 |
| sr | 14_Si_bc | `fcfact` = ecut (icmod 6) | 2e-16 | 7.0e-5 | 2.4e-4 | 14 |
| sr | 14_Si_hm | `slack` (icmod 7) | 2e-16 | 2.3e-3 | 8.2e-3 | inf |
| sr | 14_Si | `ep[2]` (well state) | 1.3e-8 | 2e-16 | – | 1.0 |
| sr | 14_Si | `dvloc0` | 2.2e-9 | 1.2e-13 | – | 1.2 |
| sr | 14_Si | `debl[0]` | 1.5e-9 | 5e-14 | – | 1.0 |
| fr | 52_Te | `dvloc0` | < 1e-12 | < 1e-12 | – | ~1 |

("–": the input has no model core, or it does not move. The pseudo-atom
test eigenvalues, printed to 8 decimals, show the same pattern at
1e-9 to 1e-5.) icmod 3 is smooth in its own parameters in the values; its
discontinuities are in the psp8 derivative columns (M-5) and in rc (M-1).

Three regimes show:
snapped inputs are piecewise constant (residuals 1e-5 to 1e-1, infinite
jump); never-snapped inputs sit on a noise floor of 1e-9 to 1e-8 from
solver tolerances (eigen-solvers stopping at 1e-8, bisections at 1e-10);
the icmod 4 Nelder–Mead search jumps outright.

The table is kept as the baseline. Every row in it, and every row added
since, has been brought to G1 or G2 or to the noise floor recorded with it
under [Phases](#phases); what remains is in
[Order of the remaining work](#order-of-the-remaining-work-2026-10-21-after-m-12).
Of that list only D-5 (deferred) and steps 4 and 5 of the test-suite
item are left (status 2026-10-31). The confined atom (item 7 of the
diagnostics asked for by the oncvpsp-opt session, `oncvpsp_confine.x`)
is done. Since api 9 the optimizer moves from icmod 6 to icmod 9 with
r_s among its search parameters (M-9).

### Goals

Each goal is checked by a test in [Tests](#tests); the numbers are
relative to the feature's scale, as `smooth_scan.py` reports them.

- **G1 No staircases.** Every scan in the matrix (T2) has jump ratio
  below 3 and residual at most **1e-10** for the kernel and V_loc and
  **1e-12** for the model core; the test eigenvalues at most **1e-11**
  once printed to full precision (G6).
- **G2 Lower noise floor.** The never-snapped inputs (`dvloc0`, `debl`,
  `ep` of a well state) reach residual at most **1e-12** for the kernel:
  solvers converge to round-off and their results do not depend on the
  iteration path.
- **G3 No loss of accuracy.** (a) With every radius and `qcut` placed
  exactly where the old code rounded them to, the new code reproduces
  the old output to the accuracy of the old numerics: kernels and V_loc
  within **1e-6**, every difference above 1e-8 explained by a named
  change of method. (b) Under mesh refinement (`amesh` 1.012, 1.006,
  1.003) the kernels converge at least as fast as today's, measured on
  14_Si, 52_Te and 08_O before and after.
- **G4 Kinks are where they should be, of the order they should be.**
  Every breakpoint (a radius where a written function is built from two
  pieces) is listed with its continuity order, and no stencil,
  quadrature cell or interpolation crosses one (T1, T6).
- **G5 Discrete events are explicit.** Each truly discrete event is
  listed here, and where an objective might care about its approach,
  the output prints a continuous margin (for example the ghost energy
  margin `E_basis - E_ref`).
- **G6 Full precision out.** Everything in the scope above is written
  with 17 significant digits: UPF (done, phase 0), psp8 (`e21.13`, 14
  digits, unchanged by format; its values are interpolated), and a
  machine-readable block of the printed diagnostics in the output.
- **G7 Unchanged otherwise.** Runtime within 1.25x of today's on the
  reference inputs; every existing test passes, with references changed
  only in documented commits.

### Numerical treatment

The log mesh is r_i = rr1·exp(al·(i−1)): uniform in x = ln r with step
al = ln(amesh). Every operation below works in x, where the mesh is
uniform and every function the generator handles is smooth away from
the origin and from its breakpoints.

#### Interpolation: Floater–Hormann in ln r

For a value at an arbitrary radius the generator used `dpnint`, an
8-point Lagrange interpolant (no longer used; the value at the origin
of the output files is FH in r, see the follow-ups) whose stencil jumps as the
point crosses a mesh node: its value is continuous there but its
derivative is not.
Splines are smooth but only C^(k−1) with a global solve, and their high
derivatives are poor. The choice here is the **Floater–Hormann
barycentric rational interpolant** with blending degree d:

    r(x) = sum_i w_i f_i/(x − x_i)  /  sum_i w_i/(x − x_i),
    w_i = (−1)^(i−d) sum_{k in J_i} C(d, i−k),  J_i = {k : 0<=k<=n−d, k<=i<=k+d}

It has no real poles. It reproduces polynomials of degree d in x, with
error O(al^(d+1)). It is C-infinity in x, and stable in barycentric form.
At a node it returns the node value exactly. It blends all local degree-d
interpolants p_k with weights λ_k that fall off as |x − x_k|^−(d+1), but
the far interpolants grow as |x − x_k|^d, so its cardinal functions fall
off only as 1/|x − x_k|. Anything rough anywhere in the array therefore
leaks everywhere at that rate: a kink, or values a solver left unset past
its last point (that is how this was found: `lschfb`'s `up` beyond `nin`
moved u′ at rc by 1e-4). Cutting the sum off near x (the first version of
this plan proposed about 60 cells) would instead make jumps of that size.

The interpolant used is **windowed**: each λ_k is multiplied by a
C-infinity bump ω of the distance from x to the window's centre, 1 within
10 steps and 0 beyond 20. In barycentric form the node weights become
w_i(x) = (−1)^i Σ_k ω(|x − c_k|) C(d, i − k), smooth in x. The result is
still C-infinity, interpolatory, and exact for polynomials of degree d;
ω does not increase away from x, so Floater and Hormann's alternating-sum
argument still excludes poles; and a node more than about 27 steps away
has no influence at all (tested bit for bit). Locality is not a
substitute for clean data: a consumer checks what it reads near the
radius (`ode_rc_der` stops if u′ from the solver's du/dn and from
stencils on u disagree beyond the solver's own error, which catches
garbage within reach), and the solvers' arrays are fixed at the source. **d = 6** (order 7;
measured: the error of r²e^−r falls 191-fold from al = 0.1 to 0.05, and
is 8e-16 at al = ln 1.006; 5 µs per value). At al = ln 1.006 this is round-off for every
smooth function the generator handles at the radii of interest. The unit
tests check reproduction, order, node exactness, and continuity of the
slope on rough data (which a switching Lagrange stencil fails).

Interpolation is only ever applied to values. Derivatives are never
taken of an interpolant; they are fields made on the mesh first (next
section) and then interpolated.

#### Derivatives at a radius

- **Generic smooth functions on the mesh** (V_AE, the core charge):
  fixed centred stencils in x of 8th order for d^j/dx^j, j ≤ 4, at
  every mesh point, converted to d/dr exactly with the Stirling numbers
  of the first kind, r^k f^(k) = sum_j s(k,j) D_x^j f. The fields are
  then FH-interpolated to the radius. The stencils are fixed, so the
  result is smooth in the radius. Round-off grows as eps/al^k: about
  1e-13, 1e-11, 1e-9, 1e-7, 1e-5 relative for k = 0..4 at al = 0.006
  (the tolerances of `test_logmesh`).
- **All-electron wave functions at rc** (the matching constraints): from
  the radial equation, not from stencils on u. `sr_flux_m` integrates
  the first-order system y = (u, q), y' = A(r) y, with A depending on
  V only (nr and sr; the Dirac system for fr has the same form). Then

      y^(n+1) = sum_{j=0..n} C(n,j) A^(j) y^(n−j)

  gives every derivative from y at the radius and A, A', A'' and A''',
  that is, from V and its first three derivatives. y itself comes by FH
  interpolation of (u, q) from the mesh. The derivatives of V come from
  the fields above. For a smooth V this replaces eps/al^4 ≈ 1e-5
  round-off in the 4th derivative of u by eps/al^3 ≈ 1e-9 in V'''.
  For a GGA potential, whose own noise is above round-off, the two are
  comparable. Measured (phase 1, `ode_rc_der`, with u′ from the solver's
  du/dn = al·r·du/dr, not from a stencil):
  - hydrogen 1s solved by `boundstate` (nr), against the analytic state,
    at amesh 1.006 / 1.003 / 1.0015: every derivative carries the
    solver's own error in u (1.0e-10 / 6.5e-12 / 4.0e-13), the 4th up to
    1.1e-10 / 8.7e-11 / 1.8e-9 (V‴ round-off); the stencils of
    `wf_rc_der` reach 7.9e-7 / 2.6e-6 / 2.1e-5 in the 4th, growing as
    the mesh refines;
  - 14_Si (sr, PBE), bound 3s and 3p at a fixed rc: the ODE route agrees
    across the three meshes to 1e-8 relative in the 3rd and 4th
    derivatives; the stencils differ from it by up to 2e-6 in the 4th.
  The ODE route is kept, with `wf_rc_der` as the cross-check (T1).

#### Integrals

- **To any limit:** the cumulative integral C_i = ∫_0^{r_i} f dr is
  computed once on the mesh. Each cell is integrated exactly over the
  degree-7 Lagrange interpolant of f·r on a fixed centred 8-point stencil
  in x, so the order is 8; at the ends the stencil is one-sided. Below
  rr(1), f is taken as the power law r^m with the integrand's known
  exponent. Then ∫_a^b f dr = C(b) − C(a), with C FH-interpolated. This is
  C-infinity in both limits, exact at the nodes, and consistent: one
  object, not a sum of full cells plus a partial-cell correction whose
  rule changes as the limit crosses a node.
- **In closed form wherever the integrand is analytic.** Inside rc the
  pseudo-wave function is a sum of spherical Bessel functions. Their
  overlaps on [0, rc] are Lommel's integrals:

      ∫_0^a r² j_l(αr) j_l(βr) dr = a²[β j_l(αa) j_l'(βa) − α j_l'(αa) j_l(βa)]/(α² − β²)
      ∫_0^a r² j_l(αr)² dr = (a³/2)[j_l(αa)² − j_{l−1}(αa) j_{l+1}(αa)]

  Where these cancel, a series takes over: when both αa, βa ≤ 1.5 their
  leading powers cancel exactly, and the integral's own power series
  (a fixed number of terms) is used; when (β − α)a ≤ min(1, αa/4) the
  numerator is expanded in (β − α) about αa, with the Taylor
  coefficients of j_l from the Bessel equation, so there is no
  cancellation. Both switches are smooth to round-off (measured ≤ 8e-15
  relative). The equal-argument form follows from the second. The same formula gives the Fourier–Bessel
  transform of the inside piece at any q, so the residual kinetic energy
  needs no r-mesh inside rc at all. The kinetic energy inside rc is
  q_j²/2 times the overlap.
- **Products of piecewise functions** (projector overlaps, B_ij) are
  split at every breakpoint. Each piece is integrated as its own smooth
  function (with its smooth extension across the breakpoint, see Kinks),
  never as one sampled array across a kink.
- **Fixed node counts.** Every quadrature whose extent depends on an
  input either keeps its node count and scales its nodes (the Gauss
  panels of modcore_budget, the partial cell at qcut), or keeps its
  spacing and gains nodes only where the integrand is negligible: the
  r-grid of eresid at its end, where u is at machine precision (Q-2), and
  its q-grid at q = 0, where the integrands vanish as q⁴ (C-10). This
  removes integer truncations like `nq = qmax/dq + 1`.

#### Kinks and breakpoints

A function built from pieces (pseudo inside a radius, all-electron
outside) is represented by its pieces, each a smooth function with a
smooth extension past the breakpoint:

- A Bessel sum is analytic everywhere.
- The vploc polynomial is a polynomial.
- V_AE and u_AE are smooth through rc.
- The Teter function and polynomials of the model core are analytic.

Every operation is done per piece on that piece's extension: values
(choose the piece by side), derivatives (one-sided by construction),
integrals (split), and interpolation onto an output mesh (each output
point from its own piece). No stencil, quadrature cell or interpolant ever
spans a breakpoint, so the error does not depend on where the
breakpoint sits relative to the mesh.

Where a product form avoids a division, it is used. The projector seed
is χ = (ε − T − V_loc)φ = (ε − V_loc)φ + ½φ'' − l(l+1)φ/(2r²). It needs
no V_l = ε − Tφ/φ, which is singular where an extension of φ has a node.

The breakpoints and their continuity orders:

| breakpoint | function | order |
|---|---|---|
| rc_l | pseudo-wave function | C^(ncon−1) |
| rc_l | V_l | C^(ncon−3) (nr); sr: see S-1 |
| rc_l | χ, β | C^(ncon−3) |
| rc(5) | V_loc, lpopt 1, 2 | C² |
| rc(5) | V_loc, lpopt 3, 4, 5 | C³ |
| 0.95·rc_l | the vr window of run_vkb | C¹ (kept, D-4); quadratures split there |
| r* | model core, icmod 1 | C⁴ |
| r_s | model core, icmod 6/7 | C² (C⁴ with `match_order = 4`, D-5, not built) |
| 3/b … 4/b | Teter taper, icmod 2–4 | C-infinity (`cinf_taper`, M-5; was a C¹ cut at 3/b) |
| r_w | wellstate well onset | C-infinity (R-12; was C²) |
| ra, rb | ionic tail switch | C-infinity (`cinf_taper`) |
| r(r_s = 1) | V_xc of PZ81 (iexc 3) | a jump of about 1e-5 Ha (the published coefficients, F-1, D-6) |

#### Roots and optima

Each one is bracketed on the mesh or a grid (a discrete step, but one that
only picks the bracket) and then polished to round-off on the
FH interpolant by Brent's method, or Newton's where the derivative is at
hand. By the implicit function theorem the root or minimizer is then as
smooth as the function, except where the bracket choice itself changes:
a new crossing or a new basin, which are the discrete events. This
applies to:

- the model-core crossover r*;
- the icmod 6/7 radius r_s;
- the wellstate width;
- the qroots wave vectors;
- the optimize secular equation;
- the icmod 3/4 minimum;
- the ecut crossings.

#### Iterations

Every iterative solve ends with its last correction applied and then one
more Newton or secant step. Its result is then at the round-off of the
fixed point, whatever the path or the iteration count. Tolerances are
relative and near round-off: 1e-13 for eigenvalues; the SCF loops test
both eigenvalues and the potential residual explicitly. An iteration
that ends on its iteration limit is an error, not a result (C-6; the
qroots skip is gone, see Q-5; the other limits of audit C-15 since
2026-10-27: `lm_root`, `root_brent` (then `ug_brent`) where the caller does not take
`converged`, the box state's bracket and bisection, the
scalar-relativistic origin series, the phase shifts' halving depth and
icmod 4's Nelder-Mead, which warned).

As done in phase 2:

- The eigen-solvers (lschfb, ldiracfb, lschvkbb, lschvkbbe, lschpsbar)
  stop after two passes in a row whose perturbative correction is below
  1e-12·max(|E|, 0.2). The second pass is at the energy with the first's
  correction applied, and its energy and function are returned. Their
  correction is Newton's step, so the second one is far below
  round-off.
- A Newton polish (qroots, optimize) stops when a step was within the
  round-off of the residual, or within a few ulps of the unknown: at a
  root of slope g′ the residual cannot get below ulp(x)·g′. The unknown
  is then evaluated once more at the moved point.
- A fixed point (ionic_screening) is iterated to a tolerance 20 times
  above its measured floor (1e-13 Ha against 1e-17 to 5e-15). With a
  contraction of 7–10 per step, the result is then within 2e-14 of the
  fixed point.
- A bracketed scalar equation (wellstate) is solved by Illinois regula
  falsi. It stops at a few ulps of the residual, or when the bracket has
  closed to a few ulps of the unknown. A closed bracket with a residual
  above 1e-10 means a jump in the function, and the program stops.

#### Policy branches

A branch on a continuous quantity that picks between methods (use the
bound state or a scattering state, drop a term below 1e-3) is either
replaced by a continuous formulation or made to depend on discrete facts
and inputs only. Each one is listed below with its decision.

### Inventory and decisions

From four audits of the code at 4352a49 (construction of the
pseudo-wave functions; local potential, projectors and output; model
core; solvers and diagnostics). Category: **a** a radius rounded to the
mesh, **b** an input rounded, **c** an integral ending at a mesh index,
**d** a kink crossed by a stencil, quadrature or interpolant, **e** a
stencil chosen by an input, **f** a tolerance or iteration count, **g** a
threshold or policy branch, **h** a discrete search, **i** a grid whose
size depends on an input. **W** marks those that change the written
pseudopotential (the others change printed diagnostics only). Phase
numbers refer to [Phases](#phases).

#### Core radii (phase 3)

| id | where | cat | W | decision |
|---|---|---|---|---|
| R-1 | `oncvpsp_main` 462–465: rc truncated to 6 decimals | b | W | done in phase 3: removed |
| R-2 | `oncvpsp_main` 575–588: rc → first mesh point beyond | a | W | done in phase 3: removed; `irc` is the last mesh point ≤ rc, for the routines that index the mesh (`pspot`'s pieces by side of rc, `wellstate`, `vrel`, `run_vkb`'s cut, the diagnostics); `run_optimize`, `fpovlp` and `run_vkb`'s windows take rc |
| R-3 | `rcmax` = max of snapped rc | a | W | done in phase 3: max of the given rc |
| R-4 | `wf_rc_der`: 7-point stencils applied 4 times at `irc` | a, e | W | done in phase 3: ODE-recursion derivatives at rc (`ode_rc_der`) |
| R-5 | `fpovlp`, `sbf_basis_con`: overlaps to `irc` with the (23,28,9)/24 end rule | c | W | done in phase 3: `lm_integral` for all-electron overlaps; `sbf_overlap` (Lommel) for Bessel ones, in quad precision |
| R-6 | `pspot`: inside/outside by mesh index; V_l joins V_AE at irc+1 | a, d | W | done in phase 3: `pieces_m` (a function of up to 4 segments, each piece held on the whole mesh as its smooth extension, integrals per segment by `lm_integral`; unit test `pieces`); `pspot` returns the inside piece and its seed (ε − T)φ̃ on the nodes to past rc; `irc` is the last node ≤ rc, so the node values are by side of rc |
| R-7 | `vploc`: v′, v″, v‴ at irc(5) by three nested 5-point stencils, polynomial for ii ≤ irc(5); round-off grows as al⁻³ (Fe: V_loc near 0 moves 3.6e-6, then 1.8e-5, as the mesh is halved twice) | a, e | W | done in phase 3 (`vploc_m`, unit test `vploc`): V, V′, V″, V‴ at rc(5) from `lm_derivs` (smooth in rc(5)); the polynomial evaluated analytically at any radius, its extension past rc(5) a piece of the projectors' integrals. The bump-weighted least-squares fit was tried: of degree 8 on the nodes within 0.25 bohr its V‴ was good to 2e-6 for V = −2/r + e^−r at 1.7 bohr (the truncation of a fit of fixed width for a potential with the derivatives of −Z/r), `lm_derivs` to 1e-9; its al⁻³ round-off (2e-11 on the default mesh) is a fixed pattern on the nodes, smooth in rc(5) through the interpolant |
| R-8 | `run_vkb`: support cut at irc; window `20(r/rr(irc) − 1)`, ff = (1−x²)², C¹; ff·vr·φ̃ with vr = Δ/2u | a, d | W | done in phase 3: the projectors are `piecewise_t` with segments split at 0.95 rc, rc(5) (and 0.95 rc(5) for lloc<4) when inside, and rc; on each segment χ = seed − V_loc^piece φ̃ − f Δ/2 (D-4, `ae_delta`), every term its smooth extension; the B matrix, normalizations and overlaps by `pw_integral`/`pw_product_integral`; the orthonormal combinations act on the pieces; `vkbpw` returned beside the node values. On the T3 inputs the fr kernels move by ≤1.1e-4 (the Δ form), sr/nr by ≤1e-6 |
| R-9 | `vpinteg` to `irc`; points 1–3 skipped | c | W | done in phase 3: `run_vkb` (R-8), `sr_so_r` and `run_ghosts` on the pieces (`pw_product_integral`, `pw_integral`), `vkboutwf` per segment on its grid (R-10); `der2exc` on `lm_integral` in phase 5, `vpinteg` removed in phase 7 |
| R-10 | `vkboutwf`, `lschvkbb`: rc found as the last nonzero vkb. Since e809a29 the overlaps end at irc and the Adams steps past irc use the source-free u″ at their back points (no O(al²) from the kink), both assuming the support ends on a node | a, e | W | done in phase 3 (`outwf_m`, unit test `outwf`): no integration on the nodes is smooth, since a step whose stencil spans a kink has an error that depends on where in the cell the kink falls (a restart on the nodes at the first node past the kink still jumps by one local truncation error, 1e-10, when the kink crosses a node); so each segment between kinks is integrated on its own grid, uniform in ln r with the mesh step, starting exactly at the segment's lower bound and running 32 points past its upper one: the first on the nodes from the power series at the first node, each later one from (u, u′) interpolated from the grid before; every segment primed by 3 Runge-Kutta steps in 4 sub-steps (so the Adams scheme starts from data exact to its own error), V and the projectors interpolated onto the grid from their pieces, the overlaps per segment by the cumulative 8th-order rule to the bound, the solution interpolated back to the nodes. Every operation then has the same form wherever the kinks sit. The method stays 4th order (the log derivative at r = 3 of the unit test's model falls 15.6, 15.4, 14.7 times on three halvings; truncation 2e-7 at the default step). The potential's pieces (`run_vkb`) and the projectors' are threaded to every solver; `lschvkbb`'s window stays beyond the last kink of either. The weights of the stencils and the cell rule are computed in quad once (they were recomputed per integral, 60 % of the runtime); the first segment's overlaps, most of the work, are dot products with the weights of `lm_integral` to its bound (`lm_integral_weights`, computed once per setup), in a compensated sum of four lanes or, with `ONCVPSP_BLAS_DOT`, one `dgemm`; the interpolants' weights are shared by u and u′ and by V and the sources. Against e3a587f the runs take 74_W_r 8.7 → 10.2 s (9.5 with BLAS), 14_Si_r 2.8 → 3.0 s, 32_Ge 2.5 → 3.1 s, 08_O_nr 1.0 → 1.2 s. More accurate: in the nr references the reference configuration's pseudo-atom eigenvalues meet the all-electron ones to 1.1e-9 Ha (03_Li, was 2.6e-8), 8.3e-10 (14_Si, was 2.0e-9), and 32_Ge's 3d to 4.1e-7 (was 1.1e-6); none got worse |
| R-11 | `vrel` window 0.5–1.5 rr(irc); 1e-10 threshold on u; vr = 0 points replaced by the average of their neighbours | a, g | W | done in phase 3: `ae_delta` (unit test `ae_delta`), Δ = u″ − [l(l+1)/r² + 2(V−ε)]u at the nodes about the window from the first-order system (u, u′ = up/(al r), V′ by `lm_derivs`), no division by u; `vrel` removed. The local potential of lloc<4 keeps Δ/(2u) for its own taper (S-1), with an error if u has a node within reach of rc |
| R-12 | `wellstate`: well onset at rr(irc+5), x³ onset (C²), inside the ±12-point stencils of `wf_rc_der` (uord 4, 5 of well states O(al²), O(al) off); thresholds 0.75, 1.5 rr(irc) override `ep`; negative `ep` reset to 0.25 | a, d, g | W | done in phase 3: onset at rc·amesh⁵ (continuous), then at rc·`well_onset`, an input of default 1.006⁵ that does not depend on the mesh (T4 at the end of phase 3); the well cwell x³/(1+x³) e^(−0.01/x), the well of before times a C-infinity factor within 5 % of 1 from x = 0.2 on, so that the states found are those of before (with e^(−1/x³), a different shape, 14_Si's kernels moved 7 % and 11 inputs found another n at the wall's lower limit; with e^(−0.1/x), still 11; with e^(−0.01/x) only 40_Zr's l=1, whose wall sat at the limit). The `rc[0]` scan 8.9e-6 → 4.4e-12. The overrides stay by default; `ep_as_given` (D-1, revised) |
| R-13 | `ionic_window`: ra = rcmax (snapped); rcore = first mesh point where `|vion+zion/r|` < 1e-10; rb = max(rcore, 2ra) | a, g, h | W | done in phase 3: ra the largest given rc; rcore the root of the interpolant of (V_ion + zion/r)² − tol² between the mesh points where it crosses (`lm_root`; the square, not the log, as V_ion + zion/r may cross zero there); rb the C-infinity maximum of rcore and 2ra within 0.05 bohr |
| R-14 | `run_diag` `mchf = max(irc)+5`; `run_phsft` rxpsh → irphs (X-2, done); `run_ghosts` mbar (X-3, done); `run_plot` n1..n4 | a | – | phase 6 |
| R-15 | `upfout` `nrlproj`, `cutoff_radius_index`; `linout` rc+2·drl format switch | a, g | W | β exactly 0 beyond rc (phase 4); indices then harmless |
| R-16 | `outwf_m`: two breakpoints that cross (rcloc and 0.95 rc(l) or rc(l) of a channel with rcloc inside its rc). The thin segment between them is integrated on its own grid, restarted at both ends from the grid before, whose stencil reaches `og_ext` points past the segment into its piece's extension; that piece is a different one on either side of the crossing (in 20_Ca, V_ae without the taper below 0.95 rc(s) = 1.805, the polynomial V_loc with it above), so the restart error goes to the same value from both sides with different slopes: a kink in the l = 0 eigenvalues where rcloc = 0.95 rc(s) (3s slope 20 % off, 1.7e-11 Ha per 1.35e-4 bohr), carried into the pseudo valence density, V_loc beyond rcloc (1.8e-11 relative, 1.9 to 3.9 bohr) and the test configurations (1.0e-10 at 161 points, 5.8e-11 at 21; 1.2e-12 with rc(s) = 1.91, which moves 0.95 rc out of the window). The B matrix (`pw_integral`, per segment on the nodes) is smooth through it (2026-10-22) | d | W | resolved (2026-10-27, the user's decision): within tolerance, below G1; kept as is |

#### Wave vectors, residual kinetic energy, optimisation (phase 3)

| id | where | cat | W | decision |
|---|---|---|---|---|
| Q-1 | `eresid` 268–274: qcut → nearest node of a q-grid with dq = 2π/(150 rc) | a | W | done in phase 3: E_r at exactly qcut, the value at the next node of the q grid plus the partial cell by 8-point Gauss-Legendre with the transforms at its nodes |
| Q-2 | `eresid`: nq = qmax/dq + 1, top 3 intervals dropped; irclin, mmaxlin integer, tail cut at \|u\| < 1e-10 | i, c, g | W | done in phase 3: the inside piece by `sbf_overlap` at each q (exact); the outside piece on a linear mesh from rc of spacing 2π/(75 qmax), continuous in rc and qmax, to the last nonzero point of u (the node count changes where u is at machine precision; a spacing set by the end point moved every node, and the 1e-12 interpolation error of u onto them jumped together), u by `lm_interp`, integrated by the 8th-order rule (`ugrid_m`); the q grid q_i = qmax − (nq − i) dq, dq = 2π/(150 rc) exactly, a node entering at q = 0 (C-10; it was of spacing qmax/(nq−1), nq an integer), accumulated from qmax by the 8th-order cumulative rule (Kahan), from the small values up |
| Q-3 | `run_optimize` 158–164: qmax = dq·int(2·qroot(nbas)/dq), max with 20, branch on qcut | i, g | W | done in phase 3: qmax = smax(smax(2 qroot(nbas), 20), 2 qcut), smax the C-infinity maximum (a+b+√((a−b)²+1))/2, within 0.5 bohr⁻¹ of the larger. The truncation at qmax stays, as a smooth bias. Considered and not taken: Parseval, E_r = T − (1/π)∫₀^q q'⁴F² (PseudoPotentialGenerators.jl `plane_wave.jl`), which needs no qmax and F below q only, but leaves E_r as the difference of two O(1 Ha) numbers, a round-off floor of ε·T/E_r ≈ 1e-11 at E_r = 1e-5 Ha (PPG's test notes a 1e-9 Ha floor); the direct sum from qmax down has the floor of E_r itself. An analytic tail beyond qmax (the transform's asymptotic form from the C^(ncon−1) matching) would remove the bias; phase 7 if wanted |
| Q-4 | `run_optimize` cvgplt/ecurve: linear interpolation assuming spacing dqout while eresid moved the q's | e, h | W (icmod 6/7) | done in phase 3: the curve of the optimized wave on the q grid (its coefficients returned by `optimize`), the crossings by Brent on the FH interpolant of ln E in the grid's index (`ug_cross_down`), in `run_optimize` and in `modcore_budget_table` (whose table is the whole curve since M-13). `dqout` is no longer used |
| Q-5 | `qroots`: scan dq = 0.02, pole filter \|f·f_last\| < 1, bisection to \|f\| < 1e-12 | h, f | W | done in phase 2: the bisection's root Newton-polished on g(q) to round-off; the scan stays (D-3). 2026-10-23: the scan and the bisection on the pole-free g(q) = ulgd j_l − q j_l′ itself, whose sign changes are exactly the roots, without the pole filter: for a large ulgd each root sits beside a pole of f, and the filter took every root for one (a Pb l=3 well state, ulgd = 65, from the optimizer's batches: the basis from roots near q = 100, or none). The bisection runs to a few ulps and cannot fail, so the skip (and its warning, whose format crashed the run) is gone. The reference outputs are unchanged |
| Q-6 | `optimize`: bisection to a 1e-10 bracket | f | W | done in phase 2: Newton polish of the norm constraint after the bisection |
| Q-7 | `optimize`: hard-case sign of x₁ | h | W | intrinsic (two global minima); keep, warn |
| Q-8 | `const_basis`: absolute 1e-8 singular-value threshold, values below left unscaled | g | W | done in phase 3: the rows of the constraint matrix (derivatives of orders 0–4 at rc, overlaps) are scaled to unit norm with their constraints, which changes nothing but the decomposition's accuracy (they differed by (q rc)⁴ and more, so the particular solution carried round-off at ε times that, 1e-11); an error if the smallest singular value is below 1e-10 of the largest |
| Q-9 | `wellstate`: rwell by ×/÷0.8 then bisection to \|E − ep\| < 1e-8; ep := E | f, h | W | done in phase 2: Illinois regula falsi on E(rwell) − ep once bracketed, to round-off; a bracket closed on a residual above 1e-10 is an error |
| Q-10 | `sbf8`, `eresid`'s `rsbf_lin`: the closed forms of j_1…j_3 cancel at small x (relative loss 3ε/x², 45ε/x⁴, 1575ε/x⁶; from `eresid`'s sin and cos, rotated from one q to the next and so good only to an absolute δ, δ(2l−1)!!(2l+1)!!/x^(2l+1), j_0 included); `sbf8` with nm = 3, 4 leaves j_1 unset; its recursion starts at nm+1 below x = 1/15 (top orders off by 1e-9, a jump where the start moves); argument reduced by mod(x, 2π) in double | – | W | done in phase 3: `sbf_m` (unit test `sbf`): `sbf_trig`, j_0…j_3 from x, sin x, cos x, the power series below x = 0.5, 1.5, 2.5, 3.5 (where no term of the series exceeds its first, the sum is ≥0.48 of it, and the closed form has lost under 2ε), 2.6ε against quad precision relative to the envelope √(j_l² + j_(l+1)²), and 2.7e-14 from sin and cos off by 1e-14; `sbf_trig_rising`, r·j_l(qr) on rising radii for `eresid`, each form in its own loop (no runtime cost). `sbf8` takes its l ≤ 3 from `sbf_trig` (every order set, no `mod`) and starts the recursion at nm + 1.36x + 15 at every x: 0.62ε·max(10, x), the recursion's own round-off over its x oscillatory orders. The kernels of the references move by ≤3e-12 (74_W_r); nothing printed moves beyond its tolerance |

#### Solvers (phase 2)

| id | where | cat | W | decision |
|---|---|---|---|---|
| S-1 | sr: V_l inside rc reproduces V+vr, outside V: V_l jumps by vr(rc) at rc | d | W | done on oncvpsp-opt (d1393c5): for lloc < 4 run_vkb applies the ff·vr window to the semi-local local potential, so V_loc(rc) = V_AE(rc) and every χ(rc) = 0; phase 3 moves it to the Δ form of D-4 with the rest of R-8 |
| S-2 | `lschvkbb`, `lschpsbar`: eps 1e-8; result is the guess before the last correction | f | W (rho, unscreening) | done: two passes in a row below 1e-12·max(\|E\|, 0.2), the second at the corrected energy (none if the first correction is at round-off); T5 to 3e-15 |
| S-3 | `lschvkbbe`, `ldiracfb`, `ldiracfs`: eps 1e-10, same | f | W | done, as S-2 (`ldiracfs` integrates outward only: its eps was unused, removed) |
| S-4 | `lschfb`: eps 1e-12, same | f | W | done, as S-2 |
| S-5 | `psatom` 204: SCF converged when the eigenvalue guess comes back unchanged | f | – | done: each eigenvalue moves by less than 1e-12 Ha and the screening residual (`ae_residual`) is below 1e-12, at the target occupations only (the ramp is now exact at its end); Anderson from the second iteration at them. `relatom` takes `sratom`'s criterion (\|ΔE\| < tolres), which unchanged eigenvalues needed |
| S-6 | `ionic_screening`: fixed point to 1e-10, 50 iterations | f | W | done: plain iteration to 1e-13 Ha; it contracts 7–10 a step to a floor of 1e-17…5e-15 (the references), so no Anderson was needed. Moves `coulomb_tail.r_coulomb`, the 1e-10 radius, in every reference |
| S-7 | `lschvkbb` nin = mch + 2.3/al, truncation of the density tail | a, b | W (~0) | done: `inward_start` in `lschvkbb` and in `lschpsbar` (inside the barrier) |
| S-8 | matching point `mch` at the turning point index; `lschfb`'s mch and `inward_start`'s nin are discrete in E (about 1e-12 in the all-electron eigenvalues) | h | W | done in phase 3: it was kept on the grounds that the all-electron atom does not depend on the scanned inputs, but well states do (their well width carried a sawtooth of 1.4e-13 Ha in E(rwell), in step with mch, into the ep scans). `match_m`: matched at the continuous turning point (an FH root), the energy correction the bump-weighted average over ±4 points of the pointwise one, the two solutions joined by a C-infinity step, in `lschfb`, `ldiracfb`, `lschvkbb` (kept beyond the projectors by a smooth maximum); `lschpsbar`'s barrier matching is phase 6 (R-14). The solvers also converge to round-off now (S-2 to S-4): the second pass below eps left its unapplied correction, up to 1e-14 relative, as noise. Unit test `matching`: E(Z) of hydrogen smooth to 4e-15 (was 3e-14 to 1e-13). nin stays discrete: the function there is at machine precision |
| S-9 | `sr_flux_m`: series/RK4 switch at rs, RK4 step count by ceiling | g, i | W (~1e-13) | done: rs = a/(8 max(1, \|b\|)), constant in E for all bound and scattering states (\|b\| ≤ 1), so the series/RK4 branch and the step count are too |
| S-10 | `exc_libxc`: ρ′, ρ″ of ρ_v + ρ̃_c by mesh stencils (`derivs`) across the kinks of the model core (r_s C², r* C⁴, the Teter cutoff C¹) and of ρ_v at rc | d, e | W (GGA) | dropped (the user's decision, 2026-10-27): RadialAtomicSolvers.jl showed that the stencil derivatives across the joints do not matter (planned: ρ′, ρ″ from the pieces) |

#### Exchange-correlation

| id | where | cat | W | decision |
|---|---|---|---|---|
| F-1 | `excpzca`: PZ81's two branches meet at r_s = 1 with the published, rounded C = 0.0020 and D = −0.0116, so ε_xc and V_xc jump there (V_xc by about 1e-5 Ha); the radius where r_s = 1 moves with the density. The solvers integrate across the jump: the all-electron eigenvalue of a state matched near it depends on the matching point by up to 1e-8 Ha (3p of Ti, matched 1 to 3 points inside its turning point). The output meshes' interpolants (`outmesh_m`) reach across it, so V_loc and the kernels of the psp8 file and of the UPF file on the linear mesh carry it over the interpolant's reach (14_Si, fcfact scanned: V_loc's 2nd difference spikes at the node it passes, 2.5e-5 off the fit there; on the linear mesh 0.79 to 1.03 bohr, up to 1.2e-5 Ha at 0.9 and 3e-8 relative at the sample at 1.0; 34_Se rcfact 9e-8; 52_Te with every rc moved, 1e-8, the spin-orbit kernel 5e-9). The smooth rows check those files' V_loc only under inputs that leave r(r_s = 1) in place (2026-10-22). The pseudo-atom's test configurations carry it too: in 08_O's 2s1 2p4 the pseudo valence density at the origin is 0.212, below n(r_s = 1) = 0.2387, so r_s = 1 falls just outside it where the density is flat, and moves fast with rc: a sawtooth of 1.2e-10 Ha in the 2s eigenvalue each time it passes a node (8 times as rc moves by 0.02; residual 4.1e-11 with PZ81 in the code or Libxc, 1.3e-13 with PW92, 4.4e-12 with PBE) | d | W | kept (D-6): noted, not changed |
| F-2 | `excpzca`, `excvwn`, `excwig`, `exchdl`, `exc_relx`: V_xc and ε_xc set to 0 below ρ = 1e-20 (1e-15 for Wigner), a jump of up to 1e-7 Ha at the radius where ρ crosses it; `exc_off_pbe` below 1e-18, about 1e-6 Ha; libxc its own thresholds | g | W (tails) | kept (D-6): noted, not changed. Measured: computing the LDA functionals at every ρ > 0 moves the kernels by at most 3.3e-9 and the AE terms by 1e-10 (`exchdl` would need ln(1+1/x) and its ε_c by a series at large x, where the closed form cancels); PBE and libxc overflow as ρ → 0 and would need a smooth damping of f(n, ∇n) with its potential |

#### Model core (phase 5)

| id | where | cat | W | decision |
|---|---|---|---|---|
| M-1 | `modcore` ircc, `modcore2` ircc, `modcore3` ircc: crossover rhoc = fcfact·rhotps at a mesh point | a | W | done in phase 5: `lm_root` of ln rhoc − ln(fcfact·rhotps); match at r* |
| M-2 | `modcore`: derivatives at ircc by stencils | a, e | W | done in phase 5: `lm_derivs` at r* |
| M-3 | `modcore`: a0 by bisection to 1e-7, monotonicity checked only at mesh points | f, h | W | done in phase 5: closed form: p' is affine in a0; critical points by companion matrix |
| M-4 | `modcore2`: x by bisection to 1e-7 | f | W | done in phase 5: Newton on the monotone x·g'/g |
| M-5 | `gg1cc`, `gpp1cc`: Teter function cut to 0 at x = 3 (C¹); columns 4, 5 by stencils across it; `dpnint` patches near the origin (zeroing rhomod(1,3:5)) | d, e | W | done in phase 5: analytic derivatives (h = √g a sum of five shifted sincs, series near their centres); delete the patches; revised: not "no cutoff" (the x⁻¹⁰ tail is far fatter than a core charge: V_ion then differs from −zion/r to 1e-10 out to 22.7 bohr in 79_Au_lxc instead of 6, and the outermost crossing of that tolerance hops between the lobes of sin²), but g tapered to 0 between x = 3 and 4 by the C∞ switch of `cinf_taper` |
| M-6 | `modcore3`: grid argmin then Nelder–Mead to an f-spread of 1e-4; returns a vertex | h, f | W | done in phase 5: keep grid and NM to choose the basin; then Newton (trust region) on ∇f = 0 to round-off; revised: to the noise of f instead, ~5e-13 in ln a (f carries ~1e-17 of der2exc's round-off on 1.2e-4, and the valley's soft curvature is 5e-4), by central differences along the Hessian's eigenvectors with steps ∝ 1/√λ (axis differences: ~1e-11) |
| M-7 | `modcore3` ircross, dead rcross loops | a | W | done in phase 5: removed (analytic derivatives everywhere) |
| M-8 | `modcore_budget` ifrom: first q of a 0.1 grid above the budget; the only way slack enters | a, b | W | done in phase 5: tested cutoffs start at q_b exactly; max over q by grid plus Brent; revised: the error sphere tapered by `cinf_taper` over 6…10 bohr (its hard edge rippled err(q) with period π/10 and the argmax hopped between ripples, kinking r_s); cutoffs to q_b + 1 first, extended to the peak + 1 when a later cutoff errs more at that r_s, so each r_s is a root of the max over all of them |
| M-9 | `modcore_budget` r_s by integer bisection; irtop, ilo from snapped rc | h, a | W | done in phase 5: bracket on the mesh, Brent on g(r_s) = max_q err − errmax; r_s in [r_top/50, r_top], r_top = min(rc) or the last radius with core charge (the real rc, not its index); revised (api 8): g has several roots on some inputs (13 of 36 Cs samples of oncvpsp-opt), and the bisection's was neither the smallest nor the largest, so the window of the minimum tail starts at q_b (one root on each sample) and r_s is the largest root, by steps of 3 % down from r_top, then Brent. Where an interior peak of the error decides r_s (more often with the window from q_b), the largest error hops between the crests of its ripple (~0.3 bohr⁻¹), and r_s is C⁰ with kinks of ~0.5 % in slope (14_Si_bc, slack 0.89 to 1; its scan moved to 0.85). And the excess is not monotone in r_s: |ΔE_xc| at q_b changes sign as r_s falls, and the lobe of the other sign meets errmax over a range of slack, so the largest root jumps when the lobe drops below it (14_Si_bc: r_s 0.426 at slack 0.88, 0.387 at 0.89; a p-norm over q moves the jump, it does not remove it). Kept, and icmod 6/7 deprecated (api 9) for icmod 8/9, which take r_s as an input: the model is smooth in r_s and slack (core 1e-15; the error's 16-mean and signed errors 1e-10, PBE; scans 14_Si_bc9_rs, 14_Si_bc9_slack, 14_Si_bc8_rs) |
| M-10 | `modcore_budget` icmod 7 candidates snapped with nint, chosen by sorting | b, h | W | done in phase 5: continuous 1-D minimisation of the mismatch on [r_s6, r_top]; revised: a ladder r_top/1.025ⁱ down to r_s6 picks the basin, then Brent on the central-difference slope of the mismatch in ln r_s (step 1e-2); the mismatch itself (`der2exc`) still reads the model core on the mesh, C² in r_s (S-10) |
| M-11 | `modcore_budget`: panel count ceiling(qhi−qlo), nxg = 20 + ceiling(qhi·r_s) | i | W | done in phase 5: fixed counts: 80 window panels × 8 nodes, 400 fill nodes, 8 per mesh interval beyond r_s, GL5 per 0.1 of q (doubling any changes nothing printed); the fill in a Legendre basis in 2t − 1 (monomials had alternating 1e5 coefficients and 1e-11 noise), monotone exactly (p′ at the ends and the roots of p″) |
| M-12 | `modcore_budget` outer integral from mesh index `is`; forward transform a plain sum across the C² joint | c, d | W | done in phase 5: split at r_s, `lm_integral`; transforms split there too (fill by Gauss–Legendre on [0, r_s], beyond it the interpolant of ρ_c by Gauss–Legendre per mesh interval, the suffix sums tabulated once). The error's integral below r_s then sampled the fill, a polynomial of degree 16 in r/r_s, on the fixed log mesh, so its quadrature error rippled as r_s slid over the mesh (4e-12 on an error of 1e-4, ~1e-9 in r_s); fixed in phase 5: that integral on a log mesh of the same step ending at r_s (the functionals need a log mesh; Gauss nodes will not do for GGA), b_q there from the rotated sines, ρ_v interpolated; +1.3 s on 14_Si_bc |
| M-13 | `modcore_budget_table`: linear interpolation in q with an assumed spacing, last crossing | e, h | W | done: the root by Q-4; the table now holds the whole q-grid curve (`ncurve` 4000, an error if a curve is longer) instead of every nstep-th node, whose interpolant changed when nstep = 1 + (nq−2)/999 did (nq is 829…1051 in 07_N_opt) |
| M-14 | `oncvpsp_main` irct (PSML matching radius) | a | W (metadata) | done: the joint, r* or r_s (with no joint, as before: the first point beyond which the model is the core charge); `number-of-continuous-derivatives` 2 for icmod 6/7, and the annotation by icmod |

#### Output (phase 4)

| id | where | cat | W | decision |
|---|---|---|---|---|
| O-1 | `linout`, `upfout`: `dpnint` (switching 8-point Lagrange) of kinked arrays onto the linear mesh; ringing beyond rc of about 1e-4 relative | d, e | W | done in phase 4: each output point from its piece by FH on that piece's extension; β ≡ 0 beyond rc |
| O-2 | `oncvpsp_main` 882–891 + `linout`: fr psp8 drops sr/so terms below 1e-3 of the first | g | W | done in phase 4: write all; see D-2 |
| O-3 | UPF numbers e20.10, PP_R f10.4 | b | W | 17 digits (done, phase 0) |
| O-4 | `upfout` nrl extended to a tail-norm radius (1e-6) | g, i | W | kept (user decision): only how many points are written |
| O-5 | `oncvpsp_main` 997: vpuns = −zion/r where rho == 0 | g | W (~1e-10) | done in phase 4: switch over (r_max/2, r_max) with `cinf_taper` |
| O-6 | printed diagnostics at f12.6–f14.8 and d10.2 | b | – | machine-readable block at 17 digits: done (`diagout`, the DIAGNOSTICS block: the Ecut crossings at every level and the log derivatives with their radius), parsed by pyoncvpsp (`diagnostics`) |

#### Diagnostics (phase 6)

| id | where | cat | decision |
|---|---|---|---|
| X-1 | `run_diag`, `run_diag_sr_so_r`: matching at mchf = max(irc)+5; renorm to the last nonzero point | a, c | done (2026-10-27): r_m = max(rc_l, rc_loc)·`diag_match`, 1.006⁵ independent of the mesh (amesh⁵ would put a first-order mesh error in the diagnostic, phase 3's T4); u and al r u′ by FH at r_m (`lm_interp`), a scattering state normalized to r_m (`log_norm_to`, `renorm_r`'s `rnorm`), `lschvkbbe` matching at r_m; rcore printed as the input rc. 07_N_opt rc(s) 1.20–1.22: jump ≤ 1.2, residuals ≤ 1e-12 except the l = 0 slope test and well-state energy shift, a mesh-period ripple of 3e-9 and 8e-10 from the interpolated derivative arrays (the table prints 1e-7) |
| X-2 | `run_phsft`: rxpsh → irphs; npsh = int(range/depsh); 2π unwrapping on any decrease; `vkbphsft` heuristic thresholds and nearest-grid alignment | a, i, g, h | done (`phsft_m`): (u, al r u′) interpolated to rxpsh exactly, the default rc·amesh² (rc of the local potential beyond lmax), the old node whenever rc is one; the phase followed in its doubled angle (sign-invariant) by continuation, the energy step halved until no step turns it by more than π/2, so no thresholds and no π or 2π corrections; aligned at ep itself. npsh kept: it depends on epsh1, epsh2 and depsh only, which the optimizer fixes |
| X-3 | `run_ghosts`: barrier at mesh index mbar; basis = states below ebcut (count jumps) | a, i | done: the wall at 3 rc exactly, a node of the log mesh scaled by s ∈ [1, amesh) (V_loc and the projectors interpolated there from their pieces, `pw_value`, `pw_resample`); each barrier state's overlaps with the projectors weighted by `cinf_taper` of its energy, 1 to ebcut and 0 from 1.5 ebcut, so states enter at weight 0; the margins (e_basis against e_ref, the least ⟨r⟩/rc) in the DIAGNOSTICS block |
| X-4 | `run_plot`: index ranges and thinning | a, h | done (2026-10-27): mostly resolved already, as the plotted values are mesh-node values of functions phase 3 made smooth, and the thinning from rr(n1) by drl kept a fixed set of nodes whenever drl < rr(end)/200. What was left: each plot's spacing min(drl, rr(end)/200) with its own end, and the ends from irc. Now one spacing min(drl, rr(n2)/200) for all plots (n2 from drl and nrl only) and the ends at rc + 1 and 1.1 rc of the input rc; the line count at an end can still change by one (plots only) |
| X-5 | `run_config`/`psatom`: see S-5; occupation ramp may exit early | f, g | done: the explicit criterion of S-5 (phase 2), at the target occupations only; the ramp ends on them exactly (audit C-1) |
| X-6 | `coulomb_tail`: α_z and its cut changes by the rectangle rule on the nodes; the 1e-10 radius the node past the last one above it, which sits at the noise floor of the unscreening (R-10 moved it 4.71 → 4.58 bohr in 14_Si nr as max \|V_loc − V_ion\| beyond rcmax fell 2.0e-9 → 5.8e-10) | c, a, g | done (2026-10-27): `lm_integral`, cut at the radii exactly; the radius by `lm_root` on (V_loc + zion/r)² − 1e-20 past the last node above it; the noise floor (max \|V_loc − V_ion\| beyond rcmax) is the next line, as before. 07_N_opt rcloc 1.18–1.20: α_z smooth to 1e-10 of 5 (degree 10), the cuts at their 1e-12 noise, the radius independent of rcloc; in the references only the radius moves (5.1574 → 5.1474 in 14_Si) |

#### Truly discrete (kept)

- Node counts in every solver, and so the principal quantum number of a
  reference. With projectors a node count is not discrete in this sense:
  the projectors' terms start as r**(l+1) like the homogeneous solution,
  and where they cancel at small r a node enters from r = 0 without an
  eigenvalue passing, then moves out with the energy. The count from
  max(0.1 rsrc, 0.05) bohr (Hamann's cutoff, vkboutwf) jumped where such
  a node crossed it: the well-state diagnostic of a Pb input failed when
  the free tail moved one from 0.23 to 0.38 bohr. `lschvkbb` and
  `lschvkbbe` count eigenvalues below E on (0, R] instead (`outwf_index`,
  by the inertia of the Schur complement of the projectors), which
  changes only where an eigenvalue passes, and `run_diag` takes its
  target from the pseudo-atom at epa rather than from the all-electron
  node count (unit test `outwf_index`: the index against finite
  differences). No solver counts nodes with the cutoff any more.
- A bound state exists or not; `boundstate`'s box-state fallback
  (flagged, `ierr=10`).
- The bound/scattering branch of the diagnostics at epa = 0.
- A ghost exists or not (the margin is printed, G5).
- The π multiple by which run_phsft aligns the pseudo phase to the
  all-electron one (changes only where they differ by π/2).
- B or S singular in `run_vkb` (a true pole of the kernel).
- The hard case of `optimize` (Q-7).
- The model-core fill type (minimum tail, Martins, none) and Martins
  admissibility. Their boundary is made exact (roots of p').
- Which basin or local minimum a search lands in (M-6, M-10).
- `pspot`'s error for a node inside rc.
- Integer and mesh inputs.
- Sign and order of individual projectors; the kernel does not see them.

### Tests

Written before the code they test. A test of a phase's code is written
and committed (red) at the start of the phase, and the phase is done
when it is green. All of these carried the ctest label `smooth` until the
plan was complete, so `ctest -LE smooth` kept the suite green on the
branch in between; phase 7 dropped it: the scans are labelled `scan`
and run by their own test preset, `ctest --preset scan`.

- **T1 Unit tests of the primitives** (`tests/unit/test_logmesh.f90`,
  written red against stubs of `src/logmesh_m.f90` and
  `src/sbf_overlap_m.f90`; green since phase 1). They check:
  - FH reproduction of degree-6 polynomials in ln r;
  - node exactness;
  - order 7 in al;
  - error below 1e-13 on the default mesh;
  - continuity of the slope across nodes for random data;
  - derivatives at off-mesh radii to 1e-13…1e-5 by order;
  - cumulative integrals and integrals between off-mesh radii to 1e-13;
  - continuity of d/db ∫^b for random data;
  - roots to round-off, and no bracket giving `ierr=1`;
  - Lommel overlaps for l = 0..4, including α = β and α = β(1+1e-9), to
    1e-13 against a 200-point Gauss rule.
- **T2 Scan matrix** (`scripts/smooth_scan.py`, ctest
  `smooth_<program>_<input>_<key>`, in `tests/CMakeLists.txt` with each
  row's tolerance and the reason for any above its goal). Every row of
  the baseline table, and with them:
  - each rc of 14_Si;
  - `rcloc` with lpopt 1 and 5;
  - `qcut`;
  - `ep` of a bound and a well state;
  - `debl`, `dvloc0`;
  - `fcfact` for icmod 1, 2, 3;
  - `rcfact` for icmod 3;
  - the icmod 4 grid;
  - `slack` for icmod 6 and 7;
  - the model core of the psp8 file with its four derivatives, for
    icmod 1, 3 and 6 (added in phase 5);
  - rc for nr 08_O and fr 52_Te;
  - `rxpsh` and the test eigenvalues, added in phase 6 with the
    full-precision diagnostics they need (`rxpsh` done; the test
    eigenvalues are still read from their 8-decimal print);
  - the optimizer's setup, 07_N_opt (`sr_07_N_opt_*`), with the phases,
    Ecut crossings and ghost margins of the DIAGNOSTICS block.

  Each scan straddles mesh nodes. Thresholds are those of G1 and G2; a
  row's threshold is tightened to its goal when its phase is done.
  The harness has its own self-test (a smooth function passes, a jump of
  1e-7 and a kink of slope 1e-6 are caught).
- **T3 On-node equivalence** (G3a), at the start and end of phase 3: the
  reference inputs with rc, rc(5) and qcut set to the values the old
  code rounded them to, run by the code before and after the phase, are
  compared with `kernel_diff.py`. The differences are recorded in the
  phase's commit, with the method change that explains each.
- **T4 Mesh convergence** (G3b), at the start and end of every phase
  that changes numerics: kernels at amesh 1.012, 1.006, 1.003 against a
  1.0015 run; the convergence rates are recorded and must not get worse.
- **T5 Solver path independence** (phase 2, written at its start): each
  eigen-solver from two start guesses on either side of the eigenvalue
  gives the same eigenvalue to 1e-13 relative and the same normalized
  function to 1e-12; hydrogen-like potentials, where the answer is known.
- **T6 Piecewise operations** (phase 3, at its start):
  - all-electron derivatives at rc against hydrogen's analytic functions
    (nr; to 1e-9 for orders ≤ 3, 1e-7 for order 4), and against the
    stencil route for sr and fr;
  - an integral of a C¹-kinked analytic function, split at the kink,
    exact to 1e-13 wherever the kink sits relative to the mesh;
  - the Wronskian identity ∫_0^rc u_i u_j = [u_i u_j' − u_i' u_j]/(2(E_j − E_i))
    for two hydrogen states, as a check of `lm_integral` against
    boundary values only.
- **T7 Existing suite**: every reference test passes. References are
  regenerated (with `scripts/ref_patch.py`) only at the end of a phase,
  in their own commit, which records the largest kernel and V_loc
  change per input.

### Phases

Each phase is: tests first (red), then the code in small commits,
`ctest -LE smooth` green throughout, the phase's T2 rows green at the
end, references patched in a final commit with the changes recorded.
Solver hygiene comes before continuous radii so that the scans of the
later phases can see 1e-10.

#### Order of the remaining work (2026-10-21, after M-12)

Reordered by what the optimizer (oncvpsp-opt, the GA/BO search and the
planned Gauss–Newton on the log-derivative residuals) actually varies and
reads. It varies rc, qcut and debl (or ep) per l, and rcloc. It always
runs Libxc PBE (iexc −101130), icmod 6 (budget at the EOS ecutrho;
icmod 9 with r_s searched since api 9, 2026-10-31), lpopt 5,
scalar-relativistic, with UPF on the log mesh; rxpsh is fixed at a mesh
point just past the rc bound. It reads the Ecut convergence
tables (log-interpolated at 1e-2…1e-5), the arctan log-derivative
phases at rxpsh, the ghost test's eigenvalues (a hard reject),
rho_core from the UPF, and the projector signs, all from stdout and the
UPF. It uses no PSML or psp8. The gradient plan takes forward
differences, h ≈ 1e-2 on [0, 1] now, limited by the f12.6 phases, and
≈ 1e-4 with O-6. It uses a Jacobian only, no Hessians.

The T2 rows so far scan kernels, V_loc and the core, and none in that
setup: the icmod 6/7 rows are LDA, and the Libxc rows use other cores.
So the order is:

1. T2 rows in the optimizer's setup: its N input (PBE, icmod 6,
   lpopt 5, UPF log mesh) scanned in rc, qcut, debl and rcloc, with
   smooth_scan.py reading the phases at rxpsh, the Ecut crossings, the
   ghost eigenvalues and UPF rho_core as well. This measures what is
   kinked there, and is the acceptance test for the items below.
   Done (sr_07_N_opt_* rows). Kernels ≤ 5.7e-12, V_loc ≤ 4.3e-12, core
   ≤ 1.6e-11 (the noise of icmod 6's search). The LD phases (1.2–1.5e-6)
   and Ecut (1.2–1.3e-4 relative) sit at their print resolution. The
   bound-ghost test's basis energy is a staircase in rc (3.4e-4 Ha steps
   for l = 1, every ~19 points, about one mesh step: its basis is cut at
   a mesh point), red until X-3. rcloc moves no Ecut. The default has
   rc[0] = rcloc, and rc < rcloc is refused (test_data), so a step down
   in rc[0] or up in rcloc is infeasible there. Running the input also
   showed the AE SCF regression fixed in b3a7c7d: a fixed 1e-11
   residual is below the floor of this mesh.
2. O-6: a machine-readable block at 17 digits with the phases at rxpsh,
   the convergence tables and their crossings at 1e-2…1e-5 per
   projector (Brent on the interpolant, as in Q-4), parsed by pyoncvpsp.
   Done: `diagout` writes the block after the gnuplot script, from
   `cvgplt` (all seven levels, 1e-2 down by half decades) and the phases
   `run_phsft` now returns; pyoncvpsp's `diagnostics` reads it, and
   smooth_scan.py takes ld and ecut from it. In the N rows ld falls from
   1.2–1.5e-6 to ≤ 1.2e-10 and ecut from 1.2–1.3e-4 to ≤ 3.7e-10, so their
   tolerances are 3e-10 and 1e-9. The references gained the block only.
   A discrete event that no item removes (G5): where a convergence curve
   has a shoulder, the last crossing of a level jumps when the shoulder
   passes it (fully-relativistic 14_Si, l = 1: 3.85 Ha at 3.2e-4, then
   7.29 Ha at 1e-4). M-13 finds the same last crossing more exactly; it
   does not remove the jump. Not seen in the N rows. A continuous margin
   (the shoulder's height above the level) would be the G5 remedy.
3. M-13: the budget table's crossings, smooth in the requested cutoff
   (runs on every trial). Done: the root-finding was Q-4's already; what
   was left is the table of every nstep-th node of the curve, which
   moved the crossing by up to 3e-8 Ha (52_Te; 2e-9 in 07_N_opt, where
   nq of l = 1 p2 is 1051) and jumped when nstep did. The table holds the
   whole curve now: the budget's Ecut is the same root on the same curve
   as the diagnostics block's. Kernels unchanged; V_loc ≤ 5e-12, core
   ≤ 7e-10 (14_Si_bca, 14_Si_hm, 14_Si_hm_r: the inputs with
   fcfact ≤ 0).
4. X-3: the ghost test continuous, with a printed margin. Done: the
   ghost group of the N rows falls from 6.5e-4 (rc[1], the staircase) to
   ≤ 9.5e-13, tolerance 3e-12, and a scan of 2e-6 in rc[1] centred where
   the wall lands on a node shows no step (3.2e-13). The block gains
   `ghost l ikap n e_ref e_basis` per row of the bound test and
   `ghost_r l ikap e r` (the basis state of least ⟨r⟩/rc more than 0.05 Ha
   above the reference, GHOST(+) below 1), parsed by pyoncvpsp. Only the
   ghost tables change: the basis eigenvalues move by up to 4e-3 Ha
   (60_Nd_GHOST l = 0, the states between ebcut and 1.5 ebcut; most by
   ≤ 1e-4), and 52_Te_r's l = 0 now reports GHOST(+) at ⟨r⟩/rc 0.995
   (the wall and the ratio at rc rather than at the mesh point below it).
5. X-2: the phase shifts, with rxpsh honoured exactly. Done: with rxpsh
   at the old node of every input, the phases are the old ones to 4e-15
   modulo π (the first now in (−π/2, π/2], an all-electron curve can move
   by π); the exceptions are two old π artifacts, where the all-electron
   curve sat π above the pseudopotential's as its sign flipped (79_Au_lxc
   l = 3 below −3 Ha, 34_Se l = 2 at its last energy), now aligned. At the
   default radius (now rc·amesh², up to a mesh step from the old node) the
   phases move by up to 0.27. A new N row scans rxpsh across a node (span
   0.004, degree 10): ld 1.8e-11 (the old code: a staircase, jump ∞). Only
   the log-derivative tables change.
6. References, then phase 7. References done: patched for phases 5
   and 6 together; the kernels do not move, V_loc moves only where the
   model core does, the log derivatives everywhere (X-2).
7. Low, for the optimizer: D-5 (deferred 2026-10-27, no current
   motivation); X-1, X-4 and X-6 (done 2026-10-27); X-5 (done by
   S-5 and audit C-1), the M-12 ripple (done: r_s smooth to 5e-12), M-14
   (PSML) and the psp8 model-core pieces (both done: the psp8 core smooth
   to 1e-10). The UPF on the log mesh writes ρ_c's values, which are
   smooth already.
8. S-10: dropped (the user's decision, 2026-10-27; RadialAtomicSolvers.jl
   showed that it does not matter). The analysis as planned: GGA ρ′, ρ″
   from the pieces, last (a small effect). Item 1 shows
   no kink above 1e-11 in the optimizer's setup: a 5-point stencil across
   a C² joint errs by a truncated power of the joint's offset in the
   cell, which is C² in the joint's radius, so finite-difference
   Jacobians do not see it. Scope: internal PBE (iexc 4), Libxc PBE and
   PBEsol. Sketch from reading the code:
   - ρ_v′, ρ_v″ from the solver's (u, al r u′) and u″ from the equation,
     u″ = [l(l+1)/r² + 2(V − e)]u + 2 Σ_j evkb_j ⟨β_j|u⟩ β_j (main's
     unscreening, `ionic_valence`, `psatom`); ρ̃_c′, ρ̃_c″ from the model
     core's columns;
   - the potential pointwise, without the stencil divergence of the
     flux: internal PBE through `exchpbe`/`corpbe` with lpot = 1 (the
     PW86(24) form, its U, V, W from the exact ρ′, ρ″); Libxc with fxc,
     v = v_ρ − 2[v_σ(ρ″ + 2ρ′/r) + ρ′(v_ρσ ρ′ + 2 v_σσ ρ′ρ″)];
   - `vout` with optional derivative arguments (into a module); the
     stencil route kept where the density is all-electron. In
     `ionic_tail` both V_xc of ρ_v^AE + ρ_c and of ρ_v^AE + ρ̃_c take the
     same route, so that they cancel beyond the core radius; in
     `der2exc` (icmod 6/7) the trial core needs its derivative columns.

#### Phase 0: infrastructure

- [x] UPF numbers to 17 digits (`upfout`): `e24.16` for every array and
      header value, the linear mesh too; references with UPF files
      patched (14_Si_UPF, 32_Ge, 52_Te_r, 74_W_r, 80_Hg_r)
- [x] `scripts/smooth_scan.py` with its self-test; the baseline above
- [x] T1 written against a stub (`logmesh_m`, `sbf_overlap`), red
- [x] T2 registered in ctest (label `smooth`): 17 of 20 red, as the
      baseline says; the three green rows guard what is already smooth
- [x] Test comparisons through pyoncvpsp and upf-tools instead of
      fldiff (on oncvpsp-opt, b4368a6); `smooth_scan.py` reads outputs
      through them too (the baseline unchanged)

#### Phase 1: primitives

- [x] `logmesh_m`: `lm_interp`, `lm_derivs`, `lm_cumint`,
      `lm_integral`, `lm_root` (FH d = 6, summed over the whole mesh)
- [x] `sbf_overlap` (Lommel, with the small-argument and near-degenerate
      series)
- [x] ODE-recursion derivatives of (u, q) for irel 0, 1, 2
      (`ode_rc_der`); measured against the stencil route on hydrogen and
      14_Si (PBE), see Derivatives at a radius
- [x] T1 green (`unit_logmesh`, `unit_ode_rc_der`)

#### Phase 2: solvers to round-off

- [x] T5 written (`tests/unit/test_solver_path.f90`), red: from guesses
      on either side, lschfb differs by up to 6e-13, ldiracfb 2p1/2 by
      1e-10, lschvkbb by 7e-9, lschvkbbe by 3e-12, lschpsbar by 2.5e-9
- [x] S-2…S-4 in lschvkbb, lschpsbar, lschvkbbe, lschfb, ldiracfb (see
      Iterations); T5 green, to 3e-15
- [x] S-5 psatom criterion (and relatom's); S-6 ionic_screening; S-7;
      S-9
- [x] Q-6 optimize, Q-5 qroots polish, Q-9 wellstate width (Illinois)
- [x] T2 rows `dvloc0`, `debl`, `ep` (well) at G2, reached in phase 3
      (R-5 in quad: `debl[0]` 5.6e-14, `ep[2]` 2.4e-13). At the end of
      phase 2: `dvloc0` was at
      2.4e-14 (14_Si) and 1.6e-14 (52_Te_r). Its 2.2e-9 at degree 6 was the
      fit's own truncation (the kernel bends in dvloc0), so the G2 rows fit
      with degree 10. `debl[0]` (5e-10) and `ep[2]` (8e-10) stay red until
      phase 3: their noise is in the l block whose energy moves and is
      already in the AE derivatives at rc, u‴/u at 1.5e-10 and u⁗/u at
      1.8e-9 against u′/u at 1e-14. That is the round-off of the
      `wf_rc_der` stencils (eps/al³, eps/al⁴), which R-4 replaces by
      `ode_rc_der`.
- [x] References regenerated on oncvpsp-opt's exact ones (eecbf77), each
      keeping its header; the commit lists the changes per input. At the
      code before and after, on the 51 reference inputs:
      - kernels move by at most 6.6e-8 (60_Nd_GHOST l=3: 1.9e-7);
      - `coulomb_tail.r_coulomb` moves in every reference (S-6);
      - 60_Nd_GHOST's truncation changes move by 2e-5;
      - nothing else moves beyond the tolerances.
      Runtime is 1.05–1.19 times the old (80_Hg_r 1.19, 74_W_r 1.15).

#### Phase 3: continuous radii, qcut, projectors

- [x] T6 written (`unit_ode_rc_der`, `unit_piecewise`; green, since it
      tests the phase 1 operations on pieces: a C¹ kink split is exact to
      8e-16 wherever it sits, one sampled array is off by 5e-12 to 1e-9).
- [x] T3 baseline recorded. `scripts/t3_inputs.py` writes the 51 reference
      inputs with rc, rc(5) and qcut where the code before this phase rounds
      them. Run by 7507bb3 on those inputs, the kernels equal the references
      exactly, and only the echoed inputs differ.
- [x] T4 baseline recorded (`scripts/mesh_convergence.py`, the psp8 kernels
      at amesh 1.012, 1.006, 1.003 against 1.0015). It is erratic, because
      each mesh snaps rc to a different point:

      | input | block | 1.012 | 1.006 | 1.003 |
      |---|---|---|---|---|
      | 14_Si | l=0 | 7.5e-3 | 7.5e-4 | 5.8e-3 |
      | 14_Si | l=1 | 2.3e-2 | 9.9e-3 | 8.8e-3 |
      | 14_Si | l=2 | 1.7e-2 | 2.7e-3 | 1.1e-2 |
      | 08_O | l=0 | 7.7e-2 | 2.4e-2 | 9.1e-3 |
      | 08_O | l=1 | 1.1e-2 | 1.6e-2 | 1.1e-2 |
      | 52_Te_r | l=0 | 1.7e-2 | 4.7e-3 | 4.0e-3 |
      | 52_Te_r | l=1 | 2.5e-2 | 1.2e-2 | 3.2e-3 |
      | 52_Te_r | l=2 | 1.5e-2 | 4.5e-3 | 1.1e-3 |
      | 52_Te_r | l=3 | 1.7e-2 | 1.1e-2 | 5.1e-3 |
      | 52_Te_r | so l=1 | 9.0e-4 | 3.9e-4 | 2.5e-4 |
      | 52_Te_r | so l=2 | 1.2e-4 | 3.3e-5 | 6.9e-6 |
      | 52_Te_r | so l=3 | 1.8e-5 | 1.1e-5 | 5.2e-6 |
- [x] R-1, R-2, R-3: rc as given; `irc` the last point ≤ rc. On the T3
      inputs (rc at the old rounding) the kernels move by ≤1e-9 against
      170bc99, so the method is unchanged; on the references rc moves by
      up to one mesh step (0.6 %), and with it the kernels (table in the
      commit). The `rc[0]` scan goes from a staircase (1e-2, jumps) to
      8.9e-6 with jumps of 2; `vloc` 2e-10; `rcloc` stays a staircase
      (`vploc`, R-7). What remains indexes the mesh: `vpinteg` to `irc`
      (R-9), `run_vkb`'s cut and the window term (R-8), `pspot`'s pieces
      (R-6), `vploc` (R-7), `vrel` (R-11), `wellstate` (R-12), the model
      core (phase 5)
- [x] R-4 matching data at rc from the ODE (`ode_rc_der` in
      `run_optimize`, in vfull, which every well potential equals out to
      rr(irc+5)). Against the stencils: 1e-9 to 1e-7 relative on bound
      states, up to 1e-4 on the 4th derivative of well states (the well's
      C² onset inside the stencils, R-12). Kernels move by at most 5.7e-6
      (74_W_r); the consistency test's "vpsp" at rc for 74_W_r's third
      projectors moves 1e-4, the stencils' error in (u/r)″ (2.4e-6) over a
      small u/r at rc (-0.013). G2: `debl[0]` 5e-10 → 4.9e-11, `ep[2]`
      8e-10 → 5.0e-10; the rest is S-8 (see there)
- [x] AE SCF (`sratom`, `relatom`) to a screening residual of 1e-11 Ha
      and eigenvalue changes of 1e-11 relative, 3 times the residual's
      round-off floor (3-4e-12 for Nd, Bi), as S-5 did for `psatom`. At
      1e-8 the energy terms depended on the SCF path to 3e-7 Ha, so any
      1e-13 change in a solver moved them by that much. Total AE
      iterations 6431 → 8516 on the 51 inputs, at most 40 in one SCF.
      The AE totals move by ≤4e-10 Ha, the terms by ≤3.1e-7 (73_Ta_r),
      the eigenvalues by ≤1e-7, the kernels by ≤5.4e-8
- [x] R-5 overlaps: `sbf_basis_con` on the closed-form overlaps of
      `sbf_overlap` (the overlap matrix, and the overlaps with the previous
      pseudo-waves from their basis coefficients, so it no longer needs
      the mesh); `fpovlp` on `lm_integral`. The basis is nearly dependent
      (the extra wave vectors of `qroots`: the overlap matrix of 14_Si l=2
      has eigenvalues 4.8e-9 … 0.81), and the mesh integrals used before
      were accurate to 1e-8, so its last orthonormal function was mostly
      their error. With Q-8: `ep[2]` 5.0e-10 → 6.5e-12. What was left was
      ε·κ of the matching solve (κ ≈ 350 after scaling) amplified by
      1/√λ_min ≈ 1.4e4 into the near-dependent direction of the basis,
      harmless in φ̃ (the combination nearly cancels) but not in φ̃″. So
      the linear algebra of the basis inside rc is in quad precision
      (`sbf_q_m`: Bessel values, derivatives at rc and Lommel overlaps;
      `qlinalg_m`: Jacobi eigen-solver; `sbf_basis_con`, `const_basis`,
      and the inside transforms of `eresid` from the quad coefficients;
      unit test `quadbasis`): `ep[2]` 1.8e-13, `debl[0]` 6.3e-14, `qcut[0]`
      5.6e-13, all at G2. The kernels move by ≤1.8e-10 (what double
      precision had lost); runtime unchanged (14_Si_r 2.7 s, 74_W_r 8.7 s:
      the many overlaps of `eresid` share the Bessel values at rc, and
      Lommel's form is kept in quad until it has cost 13 digits)
- [x] R-6 pieces (`pieces_m`)
- [x] Q-1…Q-4 residual kinetic energy at exact qcut, fixed quadratures
      (`ugrid_m`, unit test `ugrid`, `eresid` rewritten). T2 `qcut0` at G1
      passes (kernel 1e-11), `debl0` at G2 passes (5.8e-12). References
      regenerated: the kernels move by up to 1.2e-2 (19_K_st), typically
      2e-3, V_loc by up to 2e-2 (74_W_r, lloc 4 variants), which is the
      exact qcut: on the T3 inputs (qcut already at the old rounding) the
      kernels move by ≤4e-6, 60_Nd_GHOST 2e-5, 07_N 7e-5, 08_O 2e-5
      (qmax by the smooth maximum where 2 qroot ≈ 20, and the exact
      inside transforms). Runtime is unchanged (14_Si_r 2.7 → 2.5 s, 74_W_r
      9.5 → 9.1 s, after `fh_value` evaluates each window once). `ep[2]`
      stays at 5.0e-10, unchanged by this, so it is not in E_resid
- [x] R-7 vploc at continuous rc(5); R-8, R-9 run_vkb on the pieces;
      R-11 vrel → `ae_delta`; S-1 (the local potential's taper in the Δ
      form, Δ/(2u) on the window)
- [x] R-10 vkboutwf, lschvkbb, sr_so_r, run_ghosts on the pieces: the
      pseudo-atom's outward integration per segment on a grid from the
      segment's bound (`outwf_m`). The `vloc` rows: 14_Si `rc[0]`
      6.8e-10 → 5.6e-14, 08_O `rc[1]` 4.5e-8 → 5.8e-14. The kernels of
      the references move by ≤2.8e-10, V_loc by ≤6.7e-6 Ha (32_Ge, the 3d
      density inside rc), the test eigenvalues by ≤2.6e-6 Ha (73_Ta_r)
- [x] R-12 wellstate onset; D-1 behind `ep_as_given`; R-13 ionic window.
      Scans: `rc[0]` kernel 4.4e-12, `rc[2]` 8.2e-11, `rcloc` 9.3e-12 (its
      row at degree 14: the kernel changes by 70 % across the window, so
      the degree-6 fit's truncation was 1e-5), 08_O `rc[1]` 9.2e-12, all at
      G1; `vloc` 6.8e-10 (14_Si), 4.5e-8 (08_O), from the pseudo-atom
      (R-10)
- [x] Q-10 (`sbf_m`, unit test `sbf`): the low-order Bessel functions to
      a few ε at every x. The kernels of the references move by ≤3e-12,
      nothing printed beyond its tolerance; runtime unchanged
- [x] T2 rows of phase 3 at G1, kernel / vloc: 14_Si `rc[0]` 4.4e-12 /
      5.6e-14, `rc[1]` 8.0e-12 / 5.9e-13 (fitted at degree 14, as `rcloc`:
      degree 6 leaves 1.3e-9 of truncation, 10 leaves 3.3e-11), `rc[2]`
      8.2e-11 / 6.1e-14, `rcloc` 9.3e-12 / 1.0e-14, `qcut[0]` 8.0e-13 /
      1.4e-14; 20_Ca `rcloc` 4.7e-13 / 6.2e-12; 08_O `rc[1]` 9.2e-12 /
      5.8e-14; 52_Te `rc[1]` 2.7e-11 / 3.5e-14. At G2: 14_Si `dvloc0`
      7.8e-14, `debl[0]` 5.6e-14, `ep[2]` 2.4e-13, 52_Te `dvloc0` 1.2e-14
- [x] T3 at the end of the phase: the T3 inputs run by each commit of
      the phase and diffed against the commit before. Largest kernel move
      per commit:

      | commit | items | kernel |
      |---|---|---|
      | 79b1cd8 | R-4 | 5.7e-6 (74_W_r) |
      | a9e131c | AE SCF | 5.4e-8 |
      | 8eae4dd | S-8 | 6.1e-7 (60_Nd_GHOST_r) |
      | 6c310b9 | Q-1…Q-4 | 7.2e-5 (07_N) |
      | bc00d9b | R-5, Q-8 | 5.9e-7 |
      | 170bc99 | R-5 in quad | 2.1e-10 |
      | 54c5267 | R-1…R-3 | 0 (rc already on its node) |
      | e3a587f | R-6…R-13 | 1.9e-2 (40_Zr), 1.6e-3 (14_Si) |
      | dc24b82 | R-10 | 1.4e-10 |
      | 0971ea2 | Q-10 | 3.3e-12 |

      Above G3a's 1e-6, each with its change of method:
      - R-12: the well's factor exp(−0.01/x). With it removed (HEAD,
        onset rr(irc+5) as before on these inputs), 40_Zr's 1.9e-2 (its
        l=1 well state, n 5 → 6) and 14_Si's 1.6e-3 vanish against
        170bc99.
      - R-11: Δ from the first-order system instead of `vrel`, which took
        Δ/u by stencils with thresholds and neighbour averages. What
        remains without the well's factor is 2.1e-7 at most in the nr
        inputs, where Δ ≡ 0, 1.2e-5 in sr (60_Nd_GHOST's l=3) and 1.1e-4
        in fr (73_Ta_r), where Δ carries the spin-orbit term.
      - Q-1…Q-4: E_r at the exact qcut, its inside part exact. The
        kernels move by up to 7.2e-5 (07_N), 2.2e-5 (60_Nd_GHOST), 1.8e-5
        (08_O_wigner_r). V_loc moves only through the model core of
        icmod 6 and 7, by 2.9e-2 Ha (14_Si_bca) and 8.2e-3 (14_Si_hm):
        qmax went from 20 to 20.06 by the smooth maximum (Q-3), and with
        it the E_r from which the ecut budget sets the core. 14_Si_bc,
        whose fcfact is given, moves by 1e-7.
      - R-4: 5.7e-6, the stencils' error in the matching data.
- [x] T4 at the end of the phase. The kernels now converge monotonically
      in every block, where at the start each mesh snapped rc to a
      different point. At 1.003 every block is 2.5 to 320 times closer
      to the 1.0015 run than before. Only 14_Si l=1 at the coarse meshes
      is 6–9 % worse (2.51e-2 at 1.012, against 2.31e-2):

      | input | block | 1.012 | 1.006 | 1.003 | order (1.006→1.003) |
      |---|---|---|---|---|---|
      | 14_Si | l=0 | 1.8e-4 | 6.5e-5 | 1.9e-5 | 1.75 |
      | 14_Si | l=1 | 2.5e-2 | 1.1e-2 | 3.5e-3 | 1.61 |
      | 14_Si | l=2 | 3.9e-3 | 1.5e-3 | 4.8e-4 | 1.65 |
      | 08_O | l=0 | 8.3e-3 | 4.2e-4 | 1.0e-4 | 2.06 |
      | 08_O | l=1 | 6.2e-4 | 9.7e-5 | 3.4e-5 | 1.53 |
      | 52_Te_r | l=0 | 4.0e-4 | 1.8e-4 | 6.1e-5 | 1.57 |
      | 52_Te_r | l=1 | 2.2e-3 | 8.5e-4 | 2.7e-4 | 1.68 |
      | 52_Te_r | l=2 | 1.5e-4 | 4.9e-5 | 1.4e-5 | 1.86 |
      | 52_Te_r | l=3 | 4.7e-4 | 1.5e-4 | 4.1e-5 | 1.89 |
      | 52_Te_r | so l=1 | 1.5e-4 | 6.2e-5 | 2.0e-5 | 1.62 |
      | 52_Te_r | so l=2 | 2.1e-6 | 7.5e-7 | 2.2e-7 | 1.75 |
      | 52_Te_r | so l=3 | 1.6e-6 | 7.1e-7 | 2.3e-7 | 1.61 |

      Measured against the 1.0015 run, an order p shows as
      log2(2^p + 1) here, so p is only about 0.9 to 1.7. The solvers are
      4th order and the quadratures 8th, so something first order
      remains.

      It is not `dpnint`; it is the well's onset. This was measured on
      nested meshes, amesh = 1.012^(2^-k) for k = 0..4 (rr1 does not
      depend on amesh, so each mesh holds every point of the coarser
      ones). The UPF file was written on the log mesh points
      (`upf_mesh`), with no interpolation, and compared at the points
      of the coarsest mesh. The order is log2(e_k/e_(k+1)), where e_k
      is the change from mesh k to mesh k+1:
      - The UPF and psp8 kernels change by the same amounts to three
        digits, in every block of 14_Si, 08_O and 52_Te_r. Both
        converge at order 1.0 to 1.3, tending to 1 (14_Si l=1: e =
        1.46e-2, 7.08e-3, 3.44e-3, 1.69e-3).
      - With one projector per l, 14_Si's l=0 and l=1 kernels change
        by at most 2.5e-6. l=2 still converges at order 1 (its ep=0.1
        state comes from `wellstate`). The first-order term is in the
        kernels of well states: every second projector, and every
        positive ep.
      - The well state's half-point radius converges at order 1 too
        (14_Si l=0: 7.1190, 7.1335, 7.1405, 7.1440, 7.1458), because
        its onset rc·amesh⁵ (R-12) is rc(1 + 5 al) to first order.
        With the onset fixed at rc·1.012⁵ on every mesh, the 14_Si
        kernels converge at order 4.0–4.1 and then 3.3–3.5, to 1e-7.
        The 08_O kernels converge at order 2 to 4 below its coarsest
        mesh, which is off by 9e-3 with either onset. Probably the well
        reaches the kernels through the tail beyond rc that `eresid`
        transforms, but this was not checked.
      - Below about 1e-6 the changes are erratic (1e-8 to 1.4e-6, the
        largest at the finest mesh), in the bound-state kernels too: a
        floor, not an order.
      - With the onset fixed, `dpnint` shows: 08_O l=1 changes by
        2.4e-6 in psp8 and 9.6e-8 in UPF. O-1 still has a target,
        three orders down.
      - V_loc changes by about 1e-9 in 08_O and 52_Te_r. In 14_Si it
        jumps (e = 2.5e-4, 5e-10, 6.3e-5) because the icmod 1 model
        core does (`nlcc` likewise): the crossover still falls on a
        mesh point (phase 5).

      The fix is to make the onset rc·c with c independent of the mesh.
      With c = 1.006⁵, a run on the default mesh changes only by
      round-off.
      X-1's matching radius max(rc)·amesh⁵ has the same form.
- [x] The well's onset an input: `well_onset` in [pseudopotentials]
      (TOML; the `.dat` input on its lmax line since the follow-ups),
      the onset as a multiple of rc, at least 1, default 1.006⁵. 1 is
      the natural limit: the well is flat to every order at its onset,
      so it can start at rc.
      Set, it is echoed in the psp8 and UPF files as a comment, like
      `ep_as_given`. The kernels still depend on the onset, on the
      well's height and width, and on its node count (phase 8 removes
      that).
      Tests, written first and red:
      - `well_onset_{sr_14_Si,nr_08_O,r_52_Te}`: the well states on amesh
        1.006 and 1.012 agree to 5e-4 bohr in half-point radius (they
        differed by 1e-2 to 9e-2), for the default and for 1; 1 differs
        from the default; the echo appears only when set.
      - `toml_input_well_onset_range` (below 1 is an error) and
        `toml_input_well_onset_default`.
      - The scan `smooth_sr_14_Si_well_onset` over [1.015, 1.025]:
        kernel 3.1e-11, V_loc 8e-15. The kernel depends strongly on
        the onset there (0.47 per unit ln c), so the window is half the
        usual span; at 0.02, degree 6 leaves 4.9e-10 of fit truncation
        (9.8e-11 at degree 10), with no jump.
      Outputs: every reference is on the default mesh, where the onset
      moves by round-off only (rc·1.006⁵ in place of rc·amesh⁵ with
      amesh from the mesh), so no reference changes. Against 86589e6,
      the kernels of the 51 references:

      | input | kernel |
      |---|---|
      | 74_W_r | 1.5e-11 |
      | 52_Te_r | 4.7e-12 |
      | 73_Ta | 5.6e-13 |
      | 34_Se | 5.1e-13 |
      | 08_O_wigner_r | 4.4e-13 |
      | others | <= 3.2e-13 |

      V_loc moves by at most 7e-9 relative (07_N_dl_r, whose V_loc is
      the semilocal l=2 potential of a well state, inverted from its
      pseudo wave function), and otherwise by at most 1.5e-10
      (79_Au_lxc_r). The fourth derivative of 14_Si_hm's model core at
      r=0 moves by 1e-6 (a finite difference over 0.01 of a change of
      1e-13). The eigenvalues and energies do not move beyond their
      tolerances.
- [x] References patched (R-10: `scripts/ref_patch.py` now takes the lines
      of the quantities with tolerances near their printed precision, the
      pseudo-atom energies and test eigenvalues, whenever they differ)

#### Phase 4: output

- [x] Tests written, red (`ctest -L smooth -R 'psp8|upf|projector'`):
      - `smooth_scan.py --output psp8 | upf` takes the features from the
        psp8 file or the UPF file on the linear mesh. For a radius it now
        also leaves out 0.95 times the window (the taper's start moves
        with rc: 52_Te's sample at 3.0 bohr sat in it, on the linear mesh
        only), and samples the 8 points of the file's mesh on each side of
        each window, which an interpolant across the radius reaches. The
        log-mesh rows stay green with them. T2 rows, kernel / V_loc:
        14_Si `rc[0]` psp8 5.9e-11 / 5.8e-14 (green: its points next to
        rc lie in the taper's window), UPF 2.5e-10 with an infinite jump
        (`nrlproj`, R-15, cuts the projectors at 4 + int(rc/drl) points);
        08_O `rc[1]` psp8 1.4e-8; 52_Te `rc[1]` psp8 7.0e-9 (l=1 at 3.09
        bohr), UPF 7.7e-9.
      - `projector_output_{sr_14_Si,r_14_Si,nr_08_O,r_52_Te}`
        (`scripts/projector_output_test.py`): every projector of the psp8
        and UPF files 0 beyond its rc (today 2 points of psp8 and 4–10 of
        UPF carry `dpnint`'s ringing, up to 1.9e-7 in 14_Si and 3.8e-5 in
        52_Te, 9.4e-4 in a spin-orbit projector of 14_Si_r); the fr psp8
        file with nproj scalar-relativistic terms for l = 0 and 2 nproj of
        each kind for l > 0 (14_Si_r writes 2 of 4 for l = 1 and 2, and 3
        of 4 spin-orbit terms for l = 2).
- [x] O-1 output meshes from the pieces; β ≡ 0 beyond rc. `outmesh_m`
      (new): at each output radius the windowed FH interpolant
      (`lm_interp`) of a smooth function, of the piece of the radius's
      segment (`vkbpw`, `vsrpw`, `vsopw`; 0 beyond the last bound), or of
      V_loc's piece of `vlocpw` plus the smooth remainder
      vpuns − node values (the screening and the ionic tail, as in
      `outwf_setup`). `linout` and `upfout` take the pieces in place of
      `vkb`, `vsr`, `vso`; densities, model core and PP_CHI by FH, 0
      beyond a function's last nonzero node (a model core of compact
      support: FH left 1e-157 there, which `e21.13` cannot write).
      The origin keeps `dpnint`'s limit (replaced by FH in r 2026-10-27,
      see the follow-ups). On the log
      mesh (`upf_mesh` > 0) the values are the node values, unchanged.
      R-15: `linout`'s format switch and `nrlproj` now only shorten zeros.
      T2 rows, kernel: 14_Si `rc[0]` psp8 / UPF 4.4e-12 / 4.4e-12, 08_O
      `rc[1]` psp8 9.2e-12, 52_Te `rc[1]` psp8 / UPF 2.7e-11 / 2.7e-11,
      as on the log mesh. `projector_output_*`: 0 beyond rc everywhere;
      the fr ones red on the term counts only (O-2). Outputs against
      611d1da: kernels 2e-8 (14_Si, 22_Ti, 29_Cu) to 3e-5 (60_Nd) relative,
      the ringing beyond rc removed; V_loc ≤ 4e-6, densities ≤ 1e-7,
      PP_CHI ≤ 1e-9 absolute; the model core's derivative columns, which
      jump at its crossover (phase 5), up to 4e-2 / 1e-1 of their largest
      value in 14_Si_bc / 14_Si_bca (the 4th derivative), ρ_mod itself
      ≤ 2e-7
- [x] O-2 (D-2): every term of the fr decomposition written (no 1e-3
      cut), so the diagnostics of `run_diag_sr_so_r` use them all too.
      `linout` writes 2 nproj scalar and spin-orbit terms per l > 0,
      nproj scalar for l = 0. `projector_output_r_*` green. Only the fr
      psp8 files move (2 or 3 more scalar and 1 or 2 more spin-orbit
      terms per l); against 9796df1 their kernels change by the dropped
      terms, 5.7e-7 (08_O_wigner_r) to 6.7e-4 (34_Se_r) relative; every
      other output identical
- [x] O-4: dropped, the extension kept. The UPF mesh from rlmax only
      would cut the pseudo-wave functions (PP_CHI) there: of their norm,
      14_Si_UPF's 3p would lose 5.8e-2, 32_Ge's 4p 1.7e-1, 52_Te_r's 5p
      1.0e-2 (the extension to where 1e-6 is left is what keeps them). The
      extension changes only how many points are written, not a value at
      any radius
- [x] O-5: the unscreened potentials switched to −zion/r by
      `cinf_taper` over the outer half of the log mesh, (r_max/2,
      r_max), not set to it where ρ underflows to 0 (43.7–43.9 bohr of
      44.5, by mesh index; V_loc + zion/r is ≤ 6e-17 there). The written
      files are bit-identical; α_z's tail beyond 10 bohr moves by ≤ 1e-10
      (19_K), the excitation errors of the pseudo-atom tests by ≤ 1.7e-11
      (29_Cu)
- [x] T4 at the end of the phase. O-1 was the order-1 error that phase 3
      left: from 1.012 to 1.006 the blocks converge at order 3 to 4.1
      (52_Te_r l=2, already at 4.6e-8, at 1.4; 13 and 14 in 08_O, whose
      1.012 run is off by 9e-3 as before), and at 1.006 every block is
      within 1.1e-6 of the 1.0015 run (it was 7e-7 to 1.1e-2). From 1.006
      to 1.003 five blocks keep order 3.4–4.5; the others sit at a floor
      of 2e-8 to 8e-8 and wander on it (order −0.7 to 1.7). At 1.003
      every block is 150 to 70 000 times closer to the 1.0015 run than
      at the end of phase 3:

      | input | block | 1.012 | 1.006 | 1.003 | order (1.006→1.003) |
      |---|---|---|---|---|---|
      | 14_Si | l=0 | 3.5e-7 | 3.7e-8 | 2.0e-8 | 0.91 |
      | 14_Si | l=1 | 1.9e-5 | 1.1e-6 | 5.0e-8 | 4.48 |
      | 14_Si | l=2 | 9.5e-6 | 5.9e-7 | 2.7e-8 | 4.44 |
      | 08_O | l=0 | 9.2e-3 | 4.7e-7 | 6.8e-7 | −0.54 |
      | 08_O | l=1 | 5.2e-4 | 7.5e-8 | 2.3e-8 | 1.71 |
      | 52_Te_r | l=0 | 9.9e-7 | 8.7e-8 | 5.0e-8 | 0.81 |
      | 52_Te_r | l=1 | 5.1e-7 | 6.6e-8 | 7.7e-8 | −0.22 |
      | 52_Te_r | l=2 | 4.6e-8 | 1.7e-8 | 2.7e-8 | −0.65 |
      | 52_Te_r | l=3 | 2.8e-7 | 2.5e-8 | 2.2e-8 | 0.18 |
      | 52_Te_r | so l=1 | 2.9e-7 | 1.9e-8 | 1.7e-9 | 3.43 |
      | 52_Te_r | so l=2 | 3.3e-8 | 2.1e-9 | 9.6e-11 | 4.47 |
      | 52_Te_r | so l=3 | 8.5e-8 | 5.0e-9 | 2.6e-10 | 4.25 |

      The floor's source is not yet known.
- [x] References patched. `scripts/ref_patch.py --files` takes every
      changed line of the psp8 and UPF files: the kernels are compared to
      1e-6, and lines each kept within 1e-5 left five of them (32_Ge,
      60_Nd_GHOST, 52_Te_r, 74_W_r, 38_Sr_lxc_r) up to 2e-6 off

#### Phase 5: model core

- [x] Unit tests: crossover root, closed-form a0, Teter derivatives
      against series and finite differences, written first
      (`unit_modcore`, label smooth, against the stubs of `modcore_m`):
      - Teter function g and its derivatives to the 4th with no cutoff,
        against the closed form and `gg1cc`, `gp1cc`, `gpp1cc` away from
        x = 0, ½, 1, against their series constants there, each
        derivative against an 8th-order difference of the one below, and
        smooth on 41-point windows of width 0.04 from 0 to 3.5;
      - the crossover as the outermost root of ln ρ_c − ln(fcfact ρ̃_v),
        for one and three crossings, and smooth in fcfact;
      - the icmod 1 polynomial: the match at r*, p′/r² ≤ 0 on [0, r*]
        and > 0 somewhere with a0 lowered by 1e-8, a0 against a
        bisection on a fine grid, smooth in r*. The test is on p′/r²,
        not p′: where the polynomial touches 0 at the origin (p′ ≈ δ r²
        there), a grid of p′ cannot see it;
      - the icmod 2 scale x at a given x g′/g, and ierr out of range.
- [x] T2 rows added, red: `slack` for icmod 6 (14_Si_bc, about 0.9:
      V_loc 1.7e-3, core 5.9e-3, infinite jump); the model core of the
      psp8 file with its four derivatives (`smooth_scan.py`: from the
      psp8 file the `core` group has each derivative too, scaled by its
      own largest value) for fcfact of icmod 1 (14_Si, 3.1e-2 in the 4th
      derivative), rcfact of icmod 3 (34_Se, 2.1e-5, no staircase: the
      stencils and `dpnint` patches near the origin), and the cutoff of
      icmod 6 (14_Si_bc, 7.2e-2 in the 3rd). T4 at the start of the
      phase is the one at the end of phase 4.
- [x] M-1–M-5, M-7 (icmod 1–4). `modcore`, `modcore2`, `modcore3` match
      at the root r* of ln ρ_c − ln(fcfact ρ̃_v) (`core_crossover`), with
      the derivatives there interpolated (`lm_derivs`); icmod 1's a0 in
      closed form; icmod 2's scale by Newton; the Teter function and its
      four derivatives analytic everywhere, tapered to 0 over x = 3…4
      (M-5 revised, above; the switch's e^(−1/t) counts as 0 below 1e-30,
      as a model core of 1e-260 was written without its E); no ircross,
      blend, stencils or `dpnint` patches. Beyond r*, icmod 1's columns
      are `lm_derivs` of ρ_c at the mesh points. `der2exc` integrates to
      a radius (`lm_integral`): max(rcmax, r*) for icmod 1, the mesh's end
      otherwise. `gg1cc`, `gp1cc`, `gpp1cc` moved to tests/unit as the
      references of `unit_modcore`. T2: fcfact of 14_Si (icmod 1) and
      79_Au_lxc (icmod 2), rcfact of 34_Se and its psp8 core, green at
      1e-13 or better; still red: 40_Zr_grid (M-6), the icmod 6/7 rows
      (M-8…M-12), and the psp8 core of 14_Si's fcfact, 1.5e-4 (was 3.1e-2)
      in the 4th derivative, which is only C⁰ in r* (a C⁴ match) and is
      interpolated across the joint by the writers. The model core moves
      by the snap of r* to the mesh (up to a step below the root; 34_Se:
      rmatch 1.2470 → 1.2538, ρ_c there −3.9 %): kernels ≤ 7.6e-7, V_loc
      ≤ 2.8e-3, model core 1e-3 to 6e-2 relative (largest in 34_Se)
- [x] M-6 (icmod 4): after the grid and Nelder–Mead, Newton on the
      central-difference gradient in (ln a, ln s), each difference along
      an eigenvector of the Hessian (step 5e-3) with step 1e-3·√(λ₁/λₖ),
      in a trust radius of 0.1; 8–14 steps to a step of 3e-12. The
      stationary point no longer depends on the grid: 40_Zr and
      40_Zr_grid (different grids) land on the same (2.3788, 1.5525). The
      rms error falls in every icmod 4 input (40_Zr_grid_r 1.2097e-4 →
      1.2057e-4); against 3b812b9: V_loc ≤ 1.1e-3, model core ≤ 7.9e-2
      relative (40_Zr_grid_r), kernels unchanged. T2: fcfact_min of
      40_Zr_grid green (V_loc 2.3e-14, core 9.0e-13; was 1.8e-3, 7.6e-2)
      at a core tol of 3e-12, the noise above; smooth_scan.py no longer
      counts a jump ratio in steps under the tol (pure noise has a ratio
      ~3, and this row's output does not depend on its input at all)
- [x] M-8–M-12 (icmod 6/7, `modcore_budget`): cutoffs from q_b exactly,
      max over q by grid and Brent; r_s by Brent on the excess over
      [r_top/50, r_top] (rcmin passed as a radius); icmod 7 by Brent on
      the mismatch's slope in the basin of a ladder; fixed node counts;
      transforms and the error's integral split at r_s; the fill in a
      Legendre basis; the error sphere tapered by `cinf_taper` over
      6…10 bohr (`ionic_switch` renamed, in its own file, as it tapers
      more than the ionic tail). smooth_scan.py drops the radii the
      reported r_s sweeps, as for radius keys. T2: 14_Si_bc fcfact
      (degree 14) core 9.6e-11, V_loc 3.4e-11; 14_Si_hm slack (degree 14)
      4.3e-11, 1.2e-11 (xc round-off); 14_Si_bc slack 2.7e-9, 1.5e-9 (the
      quadrature ripple, M-12), each at a tol over its floor (was ~1e-6,
      kinked). Against 8974b5c, kernels unchanged: 14_Si_bc r_s 0.370954
      → 0.367395, Ecut 23.81 → 23.44 Ha, V_loc 2.7e-3, core 7.9e-3;
      14_Si_bca 0.670688 → 0.655583, 14.05 → 14.13, 6.6e-3, 2.4e-2;
      14_Si_hm 0.852000 → 0.836009, 23.46 → 23.76, 1.0e-2, 3.4e-2 (core:
      the value column relative; its 3rd–4th derivative columns move by
      O(1) at the C² joint, which the writers' pieces address: done below)
- [x] M-13: the budget table holds the whole convergence curve
- [x] The model core's pieces and M-14: `core_pieces` (modcore_m) keeps
      a model with a joint (icmod 1 at r*, 6 and 7 at r_s) as two pieces
      per column, the fill and the core charge's `lm_derivs`, each
      continued lm_fh_reach() + 2 nodes past the joint; the psp8 and UPF
      writers take each output point from its piece (`om_piece_values`),
      as for the projectors (O-1). The fill's columns in 14_Si_bc's psp8
      file are now the exact polynomial's derivatives to 1e-9 (the old
      interpolant across the jump: 4.5e-2, a 4th derivative of −7.2e7
      for 1.6e5). icmod 1 prints r* (`model radius r*`), and
      smooth_scan.py drops the radii it sweeps, as for r_s. T2, psp8
      core: 14_Si fcfact 5.4e-13 (was 1.5e-4), 34_Se rcfact 8.5e-13,
      14_Si_bc fcfact 1.2e-10 at degree 24 (was 3.0e-2, kinked; the
      floor is r_s's round-off in the steep 3rd derivative at 0.36 bohr,
      the same at span 0.02). PSML: the matching radius is the joint
      (M-14). References: only the psp8 core lines within a few points of
      the joint, and the r* line
- [x] S-10 GGA gradients from the pieces: dropped (2026-10-27),
      RadialAtomicSolvers.jl showed the effect does not matter
- [ ] `match_order` (D-5): input key, minimum tail and Martins fill to
      C⁴, T2 row `slack` with `match_order = 4`. Deferred (the user's
      decision, 2026-10-27): low priority, with no current motivation;
      S-10, the GGA across the C² joint, was dropped, and C⁴ has shown
      no practical benefit
- [x] T2 rows fcfact, rcfact, icmod 4 grid, slack: icmod 1–3 at G1;
      icmod 4, 6 and 7 at the floors measured above (core 3e-12 to 1e-10,
      V_loc 3e-12 to 1e-10), each tolerance and its reason in
      `tests/CMakeLists.txt`
- [x] References patched, with phase 6's (after X-2): only the inputs
      with a model core change beyond the log derivatives; V_loc by
      1.3e-3 to 2.7e-2 Ha, rho_mod by 3.2e-2 to 6.0 (the derivative
      columns more, at the C² joint), the kernels not at all
- [x] Revisit the icmod 6/7 quadrature ripple (M-12, known limit): fixed.
      The xc error's integral below r_s runs on its own log mesh, ri(k) =
      r_s exp((k - kin) al), from at least rr(1) to the points its
      interpolant reads beyond r_s; ρ_v is interpolated there, b_q taken
      from the rotated sines (the spectrum's points are rr(1:nx), then
      ri). Dense scan of 14_Si_bc's slack (81 points, degree 12): r_s
      residual 3.0e-9 → 5.7e-12 (noise; the ripple's period was ~8
      points), UPF core near the joint 4.7e-10 → 9.7e-13. T2: 14_Si_bc
      slack at degree 10, core 2.7e-9 → 2.9e-12, V_loc 1.5e-9 → 8.6e-13
      (degree 6 had hidden behind the ripple; the old code leaves 1.6e-9
      at degree 10); 14_Si_bc fcfact had no ripple above its fit, now at
      degree 18, core 9.6e-11 → 5.4e-12, V_loc 3.4e-11 → 1.9e-12;
      14_Si_hm slack 3.6e-11 → 2.1e-11. Tolerances 1e-11/3e-12 and
      2e-11/1e-11. r_s moves by ~1e-9, the needed cutoffs not at all;
      14_Si_bc 4.7 → 6.0 s, 14_Si_hm 3.3 → 3.8 s

#### Phase 6: diagnostics

- [x] X-1: `run_diag` matching at r_m = max(rc_l, rc_loc)·1.006⁵
      exactly (2026-10-27)
- [x] X-6: `coulomb_tail` by `lm_integral` and `lm_root` (2026-10-27)
- [x] X-4: `run_plot` on one fixed thinning, its ends at the input rc
      (2026-10-27)
- [x] X-5: S-5 and audit C-1
- [x] X-2: the phases at rxpsh exactly, continuous in the energy (`phsft_m`)
- [x] X-3: the ghost test on a wall at 3 rc exactly, a smooth basis
      cutoff, its margins in the DIAGNOSTICS block (≤ 9.5e-13 in the N
      rows)
- [x] Machine-readable diagnostics block (O-6), parsed by pyoncvpsp: the
      Ecut crossings and the log derivatives to 17 digits; the N rows'
      ld and ecut at ≤ 3.7e-10
- [x] T2 row rxpsh (sr_07_N_opt_rxpsh, ld 1.8e-11)
- [x] T2 rows of the test eigenvalues at G1: smooth_scan.py's tests
      group reads the DIAGNOSTICS block's test configurations (17 digits;
      diagnostic 4 below), each one's pseudo-atom eigenvalues and its
      total-energy difference from the reference, and is in the 19 rows
      whose input has test configurations, at 1e-10 (2026-10-22). Measured
      ≤ 3.9e-12, but 08_O rc 3.5e-11 and 20_Ca rcloc 9.0e-11, two
      discrete defects below G1, explained since (2026-10-22): in 08_O a
      sawtooth in the pseudo 2s of 2s1 2p4, PZ81's jump crossing nodes
      (F-1; gone with PW92, contrary to the first look), and in
      20_Ca a kink where rcloc crosses 0.95 rc(s) = 1.805, two breakpoints
      of the outward integration (R-16). Left out of the rows that
      move r(r_s = 1) through the pseudo-atom's density (14_Si and 34_Se
      fcfact, 34_Se rcfact, 52_Te rc[*]): there PZ81's jump (F-1) steps
      the energy differences, residuals 2.8e-7 to 7.0e-7 (14_Si fcfact
      with PW92 or PBE: 7e-15 and 1.1e-14)
- [x] References patched (with phase 5's): the log-derivative tables
      of every input

#### Phase 7: wrap-up

- [x] The numerics above moved into the permanent docs
      (`program_overview.md`, "Numerics: smooth dependence on the
      inputs", 2026-10-27); this page kept as the record
- [x] Dead code removed (`lschpb`, `lschps`, `lschkb`, `vpinteg`: no
      caller left)
- [x] The bugs below fixed where the phases did not already
- [x] `smooth` label dropped; T2 part of the default suite
      (2026-10-27): the scans carry the label `scan`; `unit_modcore` is a
      unit test and the projector-output tests are labelled as the other
      TOML tests. The dev loop is `ctest -LE "scan|nist"`. 2026-10-28:
      the scans (about 4700 s·proc, three times the rest) have their own
      test preset, `scan`, which the default one leaves out

#### Phase 8: the kernels free of the well (follow-up)

Revised 2026-10-23 (the user's decision): the free tail is for ep >= 0
only. For ep < 0 its windowed tail grows as exp(6 kappa L), and the
round-off of E_resid with its peak squared, from kappa L ~ 1.5 (an I d
well state of oncvpsp-opt at -1.61 Ha, kappa L 3.6: E_resid 2.5e3 Ha);
a well state below 0 with `well_tail = "free"` is an input error, and
"confined" or "well" are the modes for it. See
doc/changes_from_4.0.1.md, well states.

A follow-up, planned for after phases 4 to 7 and done (2026-10-22) with
parts of phases 5 to 7 still open. It changes the method,
not only the numerics, so it needs its own review (D-7).

**Why.** A well state is the atom's state at a positive energy ep (every
second projector, and every first projector with a positive ep),
confined by the well of `wellstate` so that it can be normalized. On
[0, rc] the well changes nothing: it is flat to every order at its onset
(≥ rc), so the state there is the atom's at ep, up to a factor, and the
kernels do not see the factor. Beyond rc the well shapes the tail, and
`eresid` transforms that tail (Q-2) for the residual kinetic energy that
the optimization minimizes. So the kernels depend on the well's onset,
height and width, and on its node count. Measured on 14_Si at amesh
1.006 (onset rc·c, scratch build):
- d‖K‖/d ln c is 0.45–0.48 for l=1 with today's tight box (n=4,
  half-point radius 4.5 bohr), and 0.01–0.03 one node up (n=5, 7.9 bohr).
- Where the well state takes another node, the kernels jump: 7 % for
  l=1 at c ≈ 1.08, and 3e-3 between the larger boxes. The well's node
  count is a hidden discrete variable, which changes with rc, ep and
  debl (40_Zr's l=1 in T3).
- With the tail removed from `eresid`, the optimization fails (a node in
  the first pseudo-wave function): the tail is what makes the junction
  at rc count.

`well_onset` fixes the onset at a multiple of rc that does not depend on
the mesh. The arbitrariness remains, and so do the jumps.

**Method.** In `eresid`, the tail of a well state is the atom's free
solution at ep, with no well, under a window that takes it to zero far
from rc:

    t(r) = u_free(r) w(r),  w(r) = erfc((r − R0)/L)/2,  r ≥ rc

- u_free solves the radial equation in vfull at ep, outward from rc,
  starting from the well state's u and u' at rc. Since the well state
  solves the same equation on [rc, onset], u_free equals it there.
  It is integrated on `eresid`'s own linear tail mesh (spacing
  2π/(75 qmax)), with vfull by `lm_interp` (and −zion/r beyond the
  log mesh's rmax, 45 bohr by default). On that mesh the integration
  is accurate and does not depend on amesh; on the log mesh, a 4 bohr
  wavelength at 60 bohr would get about 11 points.
- erfc is analytic, so the window's spectrum falls as a Gaussian. Its
  share of E_r above qcut is about exp(−((qcut − k)L)²/4), where k is
  the state's wave number √(2ep). With (qcut − k)L ≥ 11 that is
  below 1e-13.
- R0 = rc + 6L. There w differs from 1 by erfc(6)/2 ≈ 1e-17, and its
  derivatives by as little, so the tail joins the inner function as
  the well state does. The tail ends at R0 + 6L.
- L must be a smooth function of qcut and ep, with no branch:
  - For ep < 0 (a negative ep with `ep_as_given`, or a second projector
    below 0), u_free grows as e^(κr), κ = √(−2ep). Over the window
    (R0 + 3L − rc = 9L) the Gaussian beats the growth once
    (Δq L)²/4 − 9κL ≥ 30, with Δq = qcut − k; for ep ≥ 0 the same bound
    with κ = 0 is the leakage bound above.
  - So L is the positive root of (Δq L)²/4 − 9κL = 30, with k and κ
    from a smooth split of ep: k = √(√(ep² + δ²) + ep) and
    κ = √(√(ep² + δ²) − ep), with δ small. A max(0, ·) would put a kink
    at ep = 0.
  - qcut ≤ k (ep > qcut²/2, 12.5 Ha at qcut 5) is an input error: the
    state then oscillates faster than the cutoff.
- The well states are normalized to 1 inside rc instead of over all
  space. The kernels do not see the change; the printed residual
  energies and "Fraction of norm inside rc" of the well states do.
- What stays:
  - The well, for the energy overrides (the 0.75·rc and 1.5·rc tests,
    D-1), which use the bound state, and for the printed half-point
    radius.
  - In the fr code, the well for the j = l−1/2 energy. It is the
    eigenvalue of the j = l+1/2 well (D-7c), so `well_onset` still
    matters there, smoothly.
- Cost: the tail runs 12L from rc, about 44 bohr at Δq = 3, against 15
  to 30 bohr today, so `eresid` takes about 2.5 times as long for these
  states. To measure against G7.

**Tests, written first:**
- Unit (`eresid_tail`): the transform of a windowed sin(kr), l = 0, whose
  exact transform is δ(q − k):
  - (1/π)∫_qcut^qmax q⁴|F|² below 1e-13 of its total, for (qcut − k)L
    at the bound. This is the window's leakage.
  - The same for a growing e^(κr), with L from the ep < 0 bound.
  - L, R0 and the end point smooth in (qcut, ep) across ep = 0.
- Independence of the well: with the new tail, the kernels of 14_Si
  and 08_O change by less than 1e-10 as `well_onset` goes from 1 to 2.5.
  That range crosses the node switches, where today's kernels jump by
  7 %. 52_Te_r's do too if D-7c takes the input ep for both j.
- Independence of the window: R0 and L doubled, the kernels change by
  less than 1e-10.
- Smoothness across a former node switch: a smooth scan, at G1, of the
  ep or debl input under which today's well state changes n. It is
  found first and recorded with its jump today.
- Mesh convergence: on nested meshes (the T4 measurement of phase 3),
  the kernels of the well states converge at order ≥ 3, like those of
  the bound states.
- The rest of the T2 rows stay green.

**Steps:**
1. The tests above (red).
2. `free_tail`, with `eresid` taking the tail from it for well states
   behind the D-7 input (`well_tail`). Bound states keep their own tail.
3. The new normalization of the well states, and the printed
   diagnostics that follow from it.
4. Validation on every reference, new tail against old:
   - the kernels' change per block (the size of today's arbitrariness);
   - the log derivatives at the analysis radius;
   - the ghost test (60_Nd_GHOST, and no new ghosts);
   - the test-configuration errors;
   - the convergence profiles.
   None may get worse beyond its noise.
5. The default decided (D-7a); then the references are patched with the
   change table.

**Checklist:**
- [x] D-7 decided (2026-10-22, below): `well_tail` with "well" the
      default; the windowed growing tail for ep < 0; the j = l+1/2 well's
      eigenvalue for j = l−1/2; the well states normalized inside rc, so
      their convergence tables are on that scale
- [x] Unit test `free_tail` (the plan's `eresid_tail`): the tail of
      sin(kr) to 1.5e-12; the window's share of E_resid above qcut,
      L against 2L, 2e-27 of the tail's own there, and 5e-21 for
      sinh(κr) (L against 1.1L: at 2L the tail reaches e^34 and the
      transforms' round-off would dominate); L without a kink at ep = 0;
      `tail_kinetic` against ∫t(−t''/2) to 1e-13. The independence
      tests (`well_tail_sr_14_Si`, `well_tail_nr_08_O`,
      scripts/invariance_test.py, `well_onset` 1, 1.03, 1.09, 2.5,
      tol 1e-12) are red with the well's tail (kernels 6.5e-2 and
      4.6e-3, ecut 0.28 and 0.29 relative). T2 rows `sr_07_N_opt_rc1_free`,
      `sr_14_Si_rc1_free`, `sr_14_Si_ep2_free` (smooth_scan.py `--set`).
      The window's independence is the unit test's (no input scales L).
      The former node switch: the ep[2] row crosses ep = 0, where the
      well found no state at all below 0.005 Ha; with the free tail there
      is no n left to switch (nr, sr)
- [x] `free_tail_m` and the window; normalization inside rc. As built:
      - `well_tail = "well" | "free"` in [pseudopotentials] (TOML; the
        default "well" until D-7a), echoed in the psp8 and UPF files.
      - With "free", a well state is the atom's solution at exactly the
        energy asked for: `wellstate` only resolves that energy (its
        replacements, D-1) and, in the nr and sr codes, builds no well;
        `scatterstate` integrates outward past the reach of everything
        taken at or near rc (irc + pw_reach + 16 + the interpolants'
        reach); normalized to 1 inside rc; the windowed free tail
        beyond, on the log mesh (`free_tail_log`) and on eresid's mesh
        (`free_tail`, Runge-Kutta, 2 steps per step). Two earlier
        versions left 1e-9 in 14_Si's l = 1 kernel: the well state
        itself, whose eigenvalue the well's bracket can leave 1e-10 from
        ep, and Runge-Kutta values beyond rc next to the solver's inside
        (each node changes source as rc crosses it: 3.9e-10 in the rc[1]
        scan, 3.8e-12 now).
      - The window ends at R0 + x_e L, x_e the root of
        x² − κLx − 6κL = 39 (6.24 for κ = 0), where a growing tail is
        1e-17 of its value at rc; R0 + 6L left e^(12κL) 1e-17, order 1
        for κL ≥ 3.
      - The r-space kinetic energy of the consistency test is that of
        the same function: ∫₀^rc of the Bessel sum, and the tail's
        t t'/2 + ∫(t'² + l(l+1)t²/r²)/2 (`tail_kinetic`); ratio
        1.00000000.
      - `run_diag`, `run_diag_sr_so_r`, `run_plot` and `modcore_probes`
        took ep < 0 to mean a bound state: now a flag per projector
        (`wellprj`), as a well state can have ep < 0 (ep_as_given);
        `run_plot` still plots the fr j = l−1/2 partner of a well state
        below 0 as a bound state, as before. The references do not
        change.
      - Cost: 14_Si 2.0 s against 1.8 s, the 07_N_opt rows the same.
      - L had √(81κ² + 60Δq²) under its root, the bound at 60 rather
        than 30 (√2 too long for κ = 0; found by the numerics audit,
        doc/numerics_audit.md C-5): fixed, the validation numbers below
        the same at their printed precision.
      - The `.dat` input: `ep_as_given`, `well_onset` and `well_tail` as
        optional values of the lmax line (follow-ups, below).
- [x] Validation on the references recorded (every input of
      tests/data/toml, nr, sr and fr, "free" against "well"):
      - Kernels: by the size of today's arbitrariness, per block 1e-5 to
        7e-2 (14_Si l = 1 6.8e-2, 56_Ba l = 2 5.7e-2, 60_Nd_GHOST l = 1
        up to 9.1e-2 fr), 0.16 at most (74_W nr l = 3).
      - Test-configuration errors: the same to 3 digits everywhere.
      - Log derivatives at rxpsh (RMS of the phase error over the
        energies): the same within 1 % of itself, 74_W better
        (0.194 → 0.181 sr).
      - Ghosts: the same count everywhere; no new ones.
      - Convergence tables: those of the bound-state projectors do not
        move; those of the well states rise with their normalization
        inside rc (a well state normalized over all space had only its
        fraction inside rc there), at 1e-4 Ha by 1 to 20 Ha, 119 Ha in
        73_Ta.
      - Mesh convergence (T4, mesh_convergence.py `--set`): 14_Si's well
        blocks at orders 4.0–4.5 as before; 08_O the same as before.
- [x] References: none change (the default is "well")

#### Follow-ups (once everything above is done)

- [x] The `.dat` input brought to parity with the TOML input: `c_light`,
      `upf_mesh` and the log mesh were already optional values; now the
      lmax line is `lmax [ep_as_given [well_onset [well_tail]]]` (0 or 1;
      0 or less for the default onset; 0 "well", 1 "free"), read like
      the atom line's c_light (a line that ends early leaves the rest at
      their defaults, a bad value is an error). The echoes (stdout, psp8,
      UPF, PSML) write the three only when one is not its default, so a
      re-run reproduces the input (echo tests `*_14_Si_wt`); pyoncvpsp's
      `PseudoPotentialsInput` has the fields, in its dat and TOML forms
- [x] The value at the origin of the psp8 and linear-mesh UPF files:
      `outmesh_m` keeps `dpnint`'s limit, the polynomial of degree 7
      through the next 8 output points (0.01–0.08 bohr). Not yet trusted:
      check it against the functions' known behaviour at r = 0 (r β → 0,
      ρ'(0) = 0, V_loc finite) and against a fit on the log mesh. The
      model core's fill and the Teter cores are closed forms with
      analytic derivatives (`poly_fill`, `tail_derivatives`,
      `martins_derivatives`, `teter_core_derivs`): their columns could be
      evaluated at the output radii, r = 0 included, instead of
      interpolated from the nodes. Only the core charge beyond the joint
      is tabulated and needs `lm_derivs`;
      done 2026-10-27: `lm_origin` (logmesh_m), the FH interpolant in r
      (not ln r, where the origin is at infinity) of blending degree 6
      through 10 log-mesh nodes a factor of about 2 apart from rr(1),
      at r = 0, of the function or of the piece of its first segment;
      each window's weight is about 2^-7 of the one before, so the
      nearest dominates and there is no pole. `dpnint`'s limit and the
      override of V_loc by vpuns(1) (sr) are gone. Against a degree-7
      fit in r on the log-mesh UPF (nodes at ratio 2 from 1e-3 and 2e-3
      bohr; spread of the two fits, relative to the function's maximum,
      in parentheses), 9 inputs sr and fr: r β 3e-8 → 6e-10 (6e-10),
      PP_CHI 2e-11 → 3e-13, the model core 2e-9 → 3e-11 (3e-11), its odd
      derivatives at 0 now 1e-15 to 1e-12 (were up to 67, 14_Si_bc); the
      all-electron core density at 0 follows the e^(-2Zr) cusp (14_Si
      2.27e4 → 2.41e4). V_loc is not determined better than about 1e-8
      Ha absolute (4e-10 relative): below about 1e-4 bohr the log-mesh
      V_loc has its own structure (a slope at rr(1) in 29_Cu and 79_Au,
      GGA noise of 1e-10 per node in 07_N_opt), and the fit's value
      moves by that much with its first node; FH from rr(1) follows the
      mesh's values (07_N_opt 1.7e-8 → 5e-9 relative, 29_Cu 3e-12 →
      2e-10). In fr runs the all-electron densities diverge at r = 0
      (r^(2γ-2)); their r = 0 value is the interpolant's, not a limit.
      References: the r = 0 line of the psp8 files only (23 files).
      Then (same day): PSML's r = 0 value (`add_zero_r`, a line through
      nodes 2 and 3) is `lm_origin`'s too, so psp8, UPF and PSML agree
      there (14_Si to 13 digits); and the model core's columns are
      evaluated in closed form wherever the model is one (icmod 1's
      polynomial, the Teter function of icmod 2 to 5, the fill of 6 and
      7 inside r_s): a piece may carry a closed form (`pw_set_form`,
      pieces_m), which `pw_value` and `outmesh_m` evaluate between the
      nodes and at r = 0 in place of the interpolant; only the core
      charge beyond a joint is interpolated. The odd derivatives of
      icmod 1 at r = 0 are now exactly 0; in the Teter taper (79_Au_lxc,
      3.03 bohr) the interpolant had been 1e-8 off in the value and 4e-4
      in the 4th derivative. Then (same day) the rest of audit C-16: the
      fully relativistic PSML potentials (vp − vo per kappa) take the
      switch to −zion/r over the outer half of the mesh that vpuns has
      (O-5, `coulomb_switch`); UPF's `nrlproj` is clamped to the mesh
      size; and an odd linear-mesh size, at `rlmax` or at the orbitals'
      extent, is made even downward where upward would pass rr(mmax)
      (`outmesh_m` refused such a mesh: 14_Si at rlmax 44.695, drl 0.02,
      where rr(mmax) = 44.6995, now runs). Nothing printed moves: beyond
      the core vp − vo is V_ion, already −zion/r to 1e-12, below PSML's
      13 digits (14_Si and 80_Hg fr), and `nrlproj` stays below the
      linear mesh, which `upfout` widens to the orbitals' extent
- [x] Weighted wave cutoff for the icmod 6/7 budget (asked for by the
      oncvpsp-opt session; approved and done 2026-10-22, see
      doc/core_correction.md; T2 row sr_07_N_opt_rc1_wcut, Ecut_w smooth
      to 1.2e-11):
      `budget_cutoff = "max" | "weighted"` in `[model_core_charge]`, an
      optional 5th field after slack on the `.dat` icmod line (0 max, 1
      weighted), the default the current maximum. Ecut_w is the root of
      Σ w_i E^r_i(Ecut) = errmax Σ w_i, the weights the occupation of
      each occupied valence wave, 2(2l+1) for the first wave of an empty
      channel, 0 otherwise; Ecut_w and the weights in the DIAGNOSTICS
      block; budget slack·Ecut_w (times −fcfact if fcfact < 0), qlo from
      Ecut_w. pyoncvpsp's `ModelCoreChargeInput` gets the field and its
      dat and TOML writers
- [x] Diagnostics asked for by the oncvpsp-opt session (2026-10-22), as
      transferability proxies for the optimizer's guardrails; opt-in or
      additive DIAGNOSTICS output, approved 2026-10-22, in its order
      (all seven done: 1, 2, 3, 4 and 6 on 2026-10-22, 5 and 7 below;
      see doc/program_overview.md, step 10; 2 as the cutoffs at which a Cauchy-Schwarz bound on the
      truncation of the local energy of the pseudo valence density falls
      to 1e-2..1e-5 Ha per electron, `ecut_vloc`, smooth to 1.0e-12 to
      2.4e-12 in the 07_N_opt rows rc, qcut, debl and rcloc (1.6e-11 in
      rc1_wcut), which the ecut group now checks; 1 with the levels, plus `ekin_check`, the r-space less the
      Fourier kinetic energy to q_infinity, the tail the curve leaves out
      plus the r-space integral's error: of either sign in the references
      (195 of 421 negative), |t| ≤ 5e-10 Ha in 14_Si, up to 5.4e-7 in the
      f channels of 60_Nd and 73_Ta, so the error dominated; since the
      r-space integral is split at rc (the Parseval item), 107 of 421
      negative, ≤ 9.4e-11 in 14_Si, up to 1.4e-7 in 60_Nd f, mostly the
      energy above q∞):
      1. Longer, denser Ecut tables: levels down to 1e-6 and 1e-7, or
         the whole E^r(Ecut) curve every 1–2 Ha to E^r < 1e-7.
      2. A Fourier convergence table for V_loc (its residual, or that of
         its energy, against Ecut).
      3. Projector conditioning per channel: the condition number of B,
         the KB energies, the largest projector coefficient.
      4. Test configurations in DIAGNOSTICS, machine-readable: AE − PS
         total-energy difference against the reference configuration,
         and the eigenvalue errors.
      5. Log derivatives at a second radius (optional rxpsh2) or over a
         wider energy range. Done 2026-10-22 as a separate executable, at
         any radii (the user's request): `oncvpsp_ld.x` reads the
         all-electron file and the pseudopotential file of `--ps-file`
         (`ps_file_m`) and prints the DIAGNOSTICS `ld_r`/`ld` lines per
         l, j and radius, by default the run's to the last digit (sr, nr
         and fr, every input of tests/data); `--no-log-derivatives` takes
         them out of the run (7 to 15 % of it: 73_Ta's `_ps.x` 7.05 s to
         6.0 s, `oncvpsp_ld.x` 1.15 s for one radius, 1.6 s for 20).
         `phsft_m` takes a set of radii from one solution per energy,
         each radius continued on its own, so it gives the same phases
         as alone (to 7e-14). Found on the way: the pseudopotential's
         solution now runs past the projectors' support, without which a
         radius inside rc took their projections only to it (52_Te s and
         p at 2.5 bohr, rc 3.0 and 3.1: 3.5e-3; not reachable from the
         run, whose rxpsh check_data keeps beyond every rc); and the
         continuation misses a resonance narrower than the energy step,
         where the doubled angle turns by 2 pi in one step (73_Ta's
         all-electron 4f near -0.95 Ha at 2.6 bohr: the phases below it
         are off by pi unless a halving happens to land in it, the
         pseudopotential's too, through the alignment). Fixed 2026-10-28:
         each step is checked against the Pruefer angle (atan(r phi'/phi)
         less pi per node of u inside the radius), which a node crossing
         the radius carries across any resonance; a step whose change
         differs from it by a multiple of pi takes the Pruefer angle's
         for the all-electron potential (the same function up to a
         constant there), and is halved, to the deepest halving, for the
         pseudopotential, whose Pruefer angle also jumps by pi where a
         node appears inside the projectors (73_Ta 4f at 0.05 Ha steps:
         4.4315 below the level, 1.2899 before; tests/unit/test_phsft).
      6. Model-core extent: where the model core and the AE core fall
         below 1e-2 and 1e-4 of their value at rc; optionally its
         minimum and the r of its maximum (a hollow core).
      7. (Raised 2026-10-23, the user's request, larger; done 2026-10-23,
         the sphere and the compression below) AE
         against PS in a confined atom (harmonic or soft wall, a few
         strengths): eigenvalue and energy shifts, a stand-in for
         compression; and, in the same confinement, the cutoff that a
         channel with no occupied state needs at the energies of the
         valence bands, in place of the cutoffs at the valence energies
         (3eff290), removed 2026-10-23. Those continued the channel's
         scattering state beyond rc into a model interstitial (the local
         potential stepped to a flat level below the window), and the
         number measured that model as much as the channel: with the step
         8/qcut wide and centred 3 widths out, the tail grew in the
         classically forbidden stretch before it, by 50 to 4500 in
         amplitude in the references (2e9 in an oncvpsp-opt Li p row at
         qcut 1.25, E_resid 5e13 Ha); a C-infinity ramp from rc removes
         the growth (1.1 to 82) but its own spectrum sets the cutoff, the
         1e-5 cutoff falling with the ramp's width, 56, 31, 25, 20 Ha for
         0.5 to 4 bohr in Si d (15.6 before; Li p 94 to 22, Au f 119 to
         39), with no plateau. A confinement with a physical boundary
         condition (a Wigner-Seitz sphere: u' = 0 and u = 0 for the
         bonding and antibonding edges) needs no invented tail.
         The sphere is done (2026-10-23), as a post-processing
         executable, the user's request: `oncvpsp_confine.x` reads the
         files of `oncvpsp_ld.x` and gives, per l, j and radius rs, the
         Wigner-Seitz edges (r phi'/phi = 0, -l-1 and infinite at rs:
         bottom, centre and top; on phi = u/r, the rule's, not u), AE and
         PS at frozen potentials, and each PS edge's cutoff: its
         variational energy in the r j_l(q r) of the same condition, which
         are the free Hamiltonian's eigenstates in the sphere
         (`ws_basis_m`). The edges are roots of atan(r phi'/phi) - pi
         (nodes inside rs) (`phsft_node_phase`), monotonic beyond the core
         radii for the PS too (the projectors' terms cancel in the
         Wronskian of u and du/dE), so no scan step can miss one; the
         continuation of phsft_m missed Si 3s at 5.4 bohr with 0.02 Ha
         steps, where the phase turns by more than pi/2 within one (the
         same aliasing as the 73_Ta 4f of item 5). The PS phases are
         shifted by the multiple of pi they differ from the AE ones by at
         most energies, not at ep: Li's s reference is its 1s, whose band
         at 6 bohr is 1e-7 Ha wide, and the 2.6e-7 Ha between the AE and
         PS 1s put the two phases pi apart at ep. 14_Si at 1.25 rc: valence
         edges within 1e-5 Ha of the AE ones, 1e-5 Ha cutoffs 7 to 12 Ha
         for s, 8 to 28 Ha for d from the bottom to the top of the 3d band.
         The basis adds a q^-5 tail from V'(rs) at the bottom and centre,
         none at the top, where the tails of the PS's kinks are as slow.
         The compression is done too (2026-10-23): aeatom and psatom
         self-consistent in a wall 1/2 k x^2 exp(-0.5/x), x = r - r0
         beyond r0 = rcmax, C-infinity so outwf_m takes it into the local
         potential, ramped in with psatom's occupations and with the
         eigenvalues bracketed around their last values (the wall can
         raise them above 0: Li 2s at k = 0.5). The free atoms reproduce
         the run's AE_ref and PS_ref to the last digit. 14_Si at k = 0.01,
         0.03 and 0.1 Ha/bohr^2: dE_AE - dE_PS = 2.7e-6, 9.0e-6 and 3.3e-5
         Ha, 3p 1.1e-5 to 8.6e-5 Ha; 03_Li fr at 0.5: 1.2e-5 Ha, its 1s
         semicore 4.2e-3 Ha.
- [x] (Proposed 2026-10-22 by the oncvpsp-opt session; done 2026-10-22,
      the user's request) A
      regularized pseudo-wave optimization, so the wave is determined for
      every input. At large qcut·rc E_resid(qcut) is flat to round-off in
      the null space (6 random N rows, rc0 2.7–3.8, qcut0 8.5–13.9: e1 ~
      1e-18, |f1| ~ 1e-18–3e-16, |x1| ~ 0.7–0.86), and `optimize` takes
      the sign of x1 from f1, so a ± pair of waves, one with a node, is
      chosen by the quadrature's rounding (859c36c and d0a99ca chose
      differently). The objective gains ε times a well-conditioned term,
      for a unique wave smooth across E_resid(qcut) → 0 (not a sign
      tie-break): E_resid at a lower q (a fraction of qcut, or the wave's
      own crossing), or the squared distance from a reference wave (pswf0,
      or the AE wave in r < rc), whichever is better conditioned and
      smoother. ε an input (per channel if natural in the optimize call,
      else global) in the TOML and as an optional trailing `.dat` field,
      like budget_cutoff; ε = 0 the default and bit for bit today's
      result. In DIAGNOSTICS per wave: E_resid(qcut) at the optimum, the
      regularization term's value, e1 and |f1|. Also asked for, separately:
      a warning when e1 and |f1| are at round-off.
      Done: `wave_regularization` = w ([pseudopotentials], or the fifth
      value of the lmax line), global; `optimize` minimizes E_resid(qcut)
      + w E_resid(qcut/4) with the same secular solution, E_resid at
      qcut/4 from eresid's third q (w = 0: two q, today's path, outputs
      the same but the new lines, in 90 runs of the inputs, sr, fr and nr,
      bar the last digits of two fr ghost lines that move with the memory
      layout as the old binary's do). The distance to a reference wave was
      not taken: pswf0 is orthogonal to the null space, so it adds nothing,
      and the AE wave has nodes, which the flat null space can reproduce.
      qcut/2 and 0.35 qcut are still nearly flat in the 6 rows (regularized
      e1 5e-16 and 6e-12 at w = 1e-6), qcut/4 is not (2e-9). Warning
      `optimize: WARNING E_resid flat to round-off` when e1 and |f1| are
      both below 1e-12 Ha (6e-19–2e-17 and 2e-18–3e-16 in the 6 rows, at
      least 9e-9 in the 181 waves of the references). The block's `opt l
      ikap iproj e reg e1 f1`: E_resid(qcut), the term, e1 and |f1| of
      E_resid(qcut) alone. Measured:
      - The 6 rows: every w from 1e-12 to 1e-4 gives a nodeless wave
        (E_resid(qcut) 4e-14–2e-11 at 1e-10, up to 4e-8 at 1e-4).
      - 07_N_opt with qcut0 raised: the node error is not only the round-
        off tie. At qcut0 17–18 e1 and |f1| are ~1e-11, determined, and
        the minimum of E_resid(qcut) itself has the node; from 22 on it is
        flat (e1 2.5e-14). With w = 1e-10 qcut0 16–30 run, but a qcut
        scan at 17 meets the node at 16.98; w = 1e-8 runs everywhere.
      - Smoothness (span 0.04 in qcut0, the optimizer's groups): at large
        qcut the kernels have a floor of 1e-9 without the regularization
        too (qcut0 14, w = 0: 1.2e-9, steps even, no jump), the round-off
        of E_resid(qcut) over a weak curvature. With w = 1e-8 the scans at
        qcut0 16, 17 and 22 are at that floor (kernels 1.7e-9–2.3e-9, ld
        1.0e-9–3.6e-9, jump ratios ~1), where w = 0 has a jump of 158 in
        ld at 16 and stops at 22. In N row 22 (rc0 3.44) the rc0 scan
        falls from 8e-4 (w = 1e-10) to 2e-9 (1e-8) and 2e-10 (1e-6) in the
        kernels; ld and ecut stay steep there (ecut with shoulder jumps,
        G5).
      - Cost on ordinary inputs: the kernels move in proportion to w,
        2.8e-8 relative at 1e-10 and 2.8e-6 at 1e-8 (14_Si l = 0), 7.9e-9
        and 7.9e-7 (07_N_opt). So w ~ 1e-8 for the optimizer's high-qcut
        region, 0 elsewhere. Test: 07_N_flat (07_N_opt at qcut0 22, w =
        1e-8; with w = 0 it warns and stops on the node), with its TOML
        twin, an echo test, and errors for a negative value.
- [x] (Proposed 2026-10-22; decided 2026-10-22: the curve to infinity,
      not by Parseval, not in the optimization) E_resid by Parseval: the
      residual kinetic energy as T_r − ∫₀^q t(k) dk, the r-space kinetic
      energy less the Fourier one below q, which is the residual to true
      infinity, in place of ∫_q^q∞; in the curve and in the matrices of
      `eresid` that the optimization minimizes (the basis' r-space
      kinetic matrix less its part below qcut, in quad precision as the
      null basis). It would drop q∞ and its smooth maximum (Q-3).
      `ekin_check` bounds what it changes, and shows its obstacle: the
      r-space kinetic energy is good to 1e-10–5e-7 Ha in the references,
      no better than the 1e-7 level, so T_r would need a more accurate
      quadrature first.
      Measured 2026-10-22 (14_Si, 60_Nd_GHOST, 73_Ta, sr and fr; amesh
      1.012 to 1.00075, and copies with dq halved or q∞ doubled):
      `ekin_check` is three terms. (a) T_r's mesh error, the largest at
      the default mesh: O(h⁴) in 14_Si, about O(h³) in Nd and Ta, 3e-7
      Ha in the Nd f channels at amesh 1.006 (it has the wrong sign
      there). The jump at rc is in the third derivative of u·Tu (the
      matched derivatives), so the cell of rc matters less than first
      thought. (b) The q quadrature, dq = 2π/(150 rc): a mesh-independent
      floor, 8.8e-11 in Si l=0 #2 and 8.4e-9 in the third Ta s
      projector, gone with dq halved. (c) The energy above q∞, which is
      the term Parseval would add: 9.3e-8 and 1.4e-7 Ha in the Nd f #2
      waves, up to 1.2e-8 in the Ta f waves, below 1e-12 in Si (each gone with
      q∞ doubled). So the truncation at q∞ biases the 1e-7 Ecut level in
      Nd f, and Parseval would fix it only with T_r good to 1e-8, which
      needs amesh 1.003 or a better quadrature at the default mesh.
      Done 2026-10-22, the better quadrature (the user's request):
      `pspot` sums no longer over the mesh, switching functions at irc
      inside the cell of rc, but splits at rc. To rc, the Bessel sum's
      ψ Tψ = ψ(εψ − χ) by `lm_integral` on its smooth extension past rc
      (as the free tail's path already did); beyond, the AE wave's
      gradient form, ½∫(u′² + l(l+1)u²/r²) + ½u(rc)u′(rc), the
      boundary term from ψ: the AE mesh's u′(rc) differs from the
      matched one by the solver's own O(h⁴) error (2e-8 Ha in Nd f), and
      with it T_r converged only as h⁴. With ψ's, term (a) goes as about
      h⁶, at amesh 1.006 6.6e-12 Ha in 14_Si, 1.1e-10 in 60_Nd and
      5.9e-10 in 73_Ta (fr; it was 3e-7). What `ekin_check` shows now is
      (b) and (c), and T_r is good enough for Parseval at the default
      mesh.
      Decided 2026-10-22 (the user's): not in the optimization. The tail
      it adds, ≤ 1.4e-7 Ha, is far below the residuals the optimization
      is asked for (1e-5–1e-3 Ha); T_r − ∫₀^q has T_r's 1e-10 floor where
      the direct ∫_q^q∞ has E_r's own (e1 ~ 1e-18 in the flat rows of
      the regularization item, which it would turn into noise); the q
      quadrature (b) stays; and the matrices of `eresid` would need an
      r-space kinetic matrix of the basis. In the diagnostics, yes:
      `run_optimize` finds the DIAGNOSTICS levels below the tables'
      (3.2e-6 to 1e-7) on the curve plus `ekin_check`, by Parseval the
      curve to true infinity, good to (b), 1e-8 Ha. The tables, the plots,
      `ecurve` (the model core's budget) and the levels to 1e-5 keep the
      curve to q∞, whose tail is ≤ 1.4% of 1e-5. Below 1e-8 a level
      would need dq halved.
      Then, the same day (the user's: are both curves right, and why
      not the tail everywhere?). The curve to infinity is the quantity
      meant; the one to q∞ is low by the tail. But T_r − ∫₀^q carries
      what ∫_q^q∞ does not: the q quadrature below q, (b), up to 8.4e-9
      Ha, and T_r's mesh error. With (b) larger than the tail, in the
      107 negative channels, the correction added noise (73_Ta p #3's
      3.2e-7 level fell 6.3 Ha). dq = 2π/(300 rc) takes (b) to 2.1e-11
      (against 600; +15% CPU), but T_r's mesh error remains: it moves
      as rc crosses a cell of the mesh, and at 1e-5 Ha the curve's slope
      turns 1.8e-12 Ha of it into 4e-7 Ha of Ecut (the smooth rows of
      07_N_opt failed, ecut residual 1.1e-8 against 1e-9), and in
      74_W_r's third projectors it is −5e-9 Ha at the default mesh
      (−1.2e-11 at amesh 1.003). So the tail is not taken by Parseval.
      It is the optimized wave's own E_resid from q∞ to 3q∞, by `eresid`
      with no null functions, the same transforms as the curve: smooth
      (that scan's residual back to 1.4e-11) and independent of the
      r-space mesh (74_W_r d #3: 9.899e-9, 9.897e-9, 9.897e-9 Ha at
      amesh 1.006, 1.003, 1.0015). Above q∞ the integrand has no
      structure finer than 2π/rc: 16 dq and 25 points per period of 3q∞
      in the tail's linear mesh move it by 4e-13 Ha from dq and 75
      points; beyond 3q∞ there remains 5e-11 Ha (to 4q∞ and 6q∞), beyond
      2q∞ 8e-10, in 60_Nd f. Its cost, the 54 references 12 runs at a
      time: +3.0% of CPU (to 2q∞ on the curve's own grids it was +42%).
      With the tail in q space, dq = 2π/(300 rc) moves no Ecut level by
      more than 7e-5 Ha to 1e-5, 0.03 Ha at 1e-7, so dq stays. The curve
      with its tail is now in the tables, the plots, `ecurve` (the model
      core's budget) and every DIAGNOSTICS level (`ekin_tail`, ≤ 1.44e-7
      Ha); the optimization and its "Residual kinetic energy error" line
      keep E_resid to q∞. `ekin_check` is T_r less the Fourier energy to
      3q∞: (b) and T_r's error, −5.0e-9 to 8.4e-9 Ha in the references.
- [ ] Optimize the test suite (the user's request, 2026-10-22; approved
      2026-10-22, next after the Parseval item; steps 1 to 3 done, 4 and
      5 open as of 2026-10-31: the scan preset has grown with the icmod
      8/9 rows, so step 4 is the one that pays): more coverage for less
      wall time, by better choices of test cases and of what each
      computes. The user's aim for the smooth label: the fewest inputs
      that keep at least today's smoothness coverage (a map of what each
      row exercises first, then the smallest set that covers it all). Where the time
      goes now (`ctest -LE smooth -j 16`, the label now `scan`, 310 tests, 186 s wall, about
      2700 s·proc): toml 760, echo 556, mesh 490, nist 388, split 195,
      r 187, sr 182 s·proc; the smooth label about 100 s per row at 4
      jobs, before parsing the outputs in parallel (smooth_scan.py).
      1. Measure the coverage first (gfortran --coverage, gcov or gcovr
         per label), so that the choice is by what each test reaches,
         not by habit.
      2. Find the tests that rerun the same atom only to check another
         output or input path (echo, toml, split, mesh): share one run's
         outputs among the checks that read them, or use a light atom
         (H, Li, N) where the heavy one adds nothing to what is tested.
      3. Keep the heavy atoms (f shells, fully relativistic, GHOST)
         where they cover what nothing lighter does, and say in
         tests/CMakeLists.txt what each case is for.
      4. Fewer scan points or rows where the fit allows it, and the rows
         ordered so that ctest -j fills the cores (the longest first).
      5. Add what the coverage shows is missing (error paths, rarely
         taken branches of the solvers) at the time saved.
      The user widened the aim the same day: the fewest inputs for the
      broadest code coverage, first matching today's coverage, then
      filling its gaps. The nist label is left out of the coverage: it
      checks accuracy against reference all-electron data, and is rerun
      when the all-electron or solver numerics change.
      Done 2026-10-22 (step 1 and the first cut):
      - Coverage per test from a -O0 --coverage build (gcov, every
        test with its own GCOV_PREFIX): the non-smooth suite without
        nist runs 88 % of the lines of src/ and 62 % of the branches
        that gcov counts.
      - Each test is put in a class by what it checks (ref against
        tests/refs, toml against the text input, echo, split, log and
        UPF mesh, pipe, psml, psp8 core, well onset and tail,
        projector output, unit), and each class must keep its lines
        and branches. Line coverage alone would not do: it drops the
        14_Si reference test because the echo and mesh tests run the
        same lines, but those do not compare the numbers with a
        reference.
      - The fewest inputs keeping every class's coverage are 29 of 37,
        and on them 163 of 260 tests (greedy, by -O2 cost). Kept as
        well: every error, toml_input and ae_file case (each a
        documented message, some of it in toml-f, which is not
        instrumented; 0.1–2 s each) and 73_Ta's reference (the example
        of doc/users_guide.md). Removed: the inputs 07_N, 14_Si_st,
        17_Cl, 19_K, 19_K_st, 22_Ti and 56_Ba, and 28 references no
        test reads.
      - `ctest -LE smooth -j 16`, the same parser for both:

        | | tests | wall | s·proc | without nist |
        |---|---|---|---|---|
        | before | 264 | 185 s | 2729 | 2356 |
        | after | 179 | 112 s | 1743 | 1360 |

        Without nist (`ctest -LE "smooth|nist" -j 16`): 257 tests in
        150 s wall before, 172 in 88 s after. With nist the wall time is
        set by nist_RLDA_heavy (90 s).
- [x] Configurable Adams–Moulton order in the radial ODE solvers (the
      user's idea, listed 2026-10-22 at their request through another
      session; done 2026-10-28, steps 1 to 3, with step 4 checked and
      step 5 measured on 14_Si; the plan as written follows the note).
      **Done.** `adams_m` holds the Adams–Moulton coefficients of order
      k = 4 to 8 (exact rationals); `adams_order` is the optional fourth
      value of the log-mesh line, or `[log_mesh] adams_order`, default
      4, echoed only when it is not 4. A `type(adams_t)` is passed down
      every call chain from oncvpsp_run (no module state; its default
      initialization is order 4), and held by `outwf_t`, `phsft_t` and
      modcore_budget's state. Each step's implicit equation is solved
      exactly (`am_solve`, 2×2), so there is no predictor; every
      integration starts from k points (the series at the origin, the
      exponential at nin, nin at most mmax − k, Runge-Kutta priming of
      k − 1 steps in outwf_m), and vout's Hartree potential is the
      cumulative rule of order k. oncvpsp_ld and oncvpsp_confine use
      order 4. Unit test `adams`: order 4 is the default, every order is
      exact on polynomials of degree k, and the hydrogen-like 2p and 3d
      errors fall from amesh 1.024 to 1.012 by 15, 32, 44 and 122 at
      orders 4 to 7, order 8 reaching round-off at 1.012. At order 4 the
      exact solve moves the references within their tolerances, except
      the barrier states of the positive-energy ghost test above 50 Ha
      (`ghost_localization`: <r>/rc by up to 2.4e-4, the energy by 3e-5
      relative, in 15 references). At 58 to 222 Ha such a state turns
      by 0.2 to 0.3 radians per mesh step, where order 4's own error is
      larger: 40_Zr's l = 2 entry is another state (122 Ha) at orders 6
      and 8, which differ from each other by 1.6e-3 in <r>/rc. The
      07_N_opt rcloc scan's ld residual is 3.0e-10 (2.7e-10 at order 6,
      1.9e-10 at 8; jump ratio 1.03), its tolerance raised to 4e-10.
      **Measured (step 5, 14_Si all-electron eigenvalues, against amesh
      1.003 at order 8; the same in nr):**

      | amesh | order | 1s | 2p | 3s | 3p |
      |---|---|---|---|---|---|
      | 1.012 | 4 | 1.1e-7 | 5.9e-8 | 3.6e-9 | 2.9e-9 |
      | 1.012 | 6, 8 | 4.7e-8 | 4.3e-8 | 5.0e-9 | −6.8e-11 |
      | 1.006 | 4 | 1.8e-8 | 1.4e-8 | 1.4e-9 | 1.7e-10 |
      | 1.006 | 6, 8 | 1.4e-8 | 1.3e-8 | 1.5e-9 | −2.0e-11 |
      | 1.003 | 4 | 2.4e-10 | 6.4e-11 | −4.6e-12 | 1.2e-11 |

      Orders 6 and 8 agree to 1e-11 at every mesh, so the ODE's error is
      gone, but an error independent of the order remains, 4e-8 Ha at
      1.012 and 1.3e-8 at 1.006 (about al^1.7). It is not the start at
      the origin (rr1 = 1e-7 leaves it), nor relativistic. This is the
      "about 1e-7 where the differences stop falling smoothly" of step 7.
      **Found (2026-10-28): the functional, not the numerics.** 14_Si is
      iexc 3, Perdew–Zunger, whose two branches meet at rs = 1 with
      4-digit CC and DD (`excpzca`): ε_c jumps there by 3.2e-5 Ha and
      v_c by 2.8e-5 Ha. The step sits at a fixed radius, so the error
      depends on where the mesh falls, not on the step: moving rr1
      through one cell at amesh 1.006 (order 8, nr) gives a sawtooth in
      1s of ±4e-8 Ha (linear, with one jump per cell; rmax ±0.03 moves
      nothing). libxc's PZ (−001009) gives the same sawtooth to 1e-10;
      PZ_MOD (−001010, PZ's continuous coefficients) ±4e-9, from the
      kink left; PW92 (−001012) 3e-12. With PW92, 14_Si against amesh
      1.003 at order 8 (sr and nr alike, 1s the largest):

      | amesh | order 4 | order 6 | order 8 |
      |---|---|---|---|
      | 1.048 | 1.4e-5 | 2.3e-7 | 1.5e-8 |
      | 1.024 | 9.3e-7 | 3.2e-9 | 4.7e-11 |
      | 1.016 | 1.9e-7 | 2.7e-10 | 1e-11 |
      | 1.012 | 6.0e-8 | 5.8e-11 | 1.4e-11 |
      | 1.006 | 3.8e-9 | 1.6e-11 | 1.5e-11 |

      So order 8 at amesh 1.024 (a quarter of the default's points) is
      within 5e-11 Ha, and order 6 at 1.016 within 3e-10, for a smooth
      functional. For iexc 3 (and libxc PZ) the floor is PZ's own step,
      about 4e-8 at 1.006 and 1.3e-7 at 1.024, whatever the order; the
      coefficients are the xc implementation's, left as they are.
      The plan as written: AB4 predictor and two
      fixed-point AM4 corrector passes per mesh step (`adams.inc`, in
      sr_flux_m, lschfb, lschpsbar, outwf_m, vout, dirac_m, ldiracfb,
      lschvkbb), global error O(al⁴), while the interpolation,
      derivatives and quadrature are of order 7–8 (audit §6a).
      1. The implicit AM step solved exactly (a 2×2 linear solve; the
         equations are linear): the true implicit method, cheaper, no
         iteration-count artefact, still smooth (rational in the
         coefficients, denominator ≈ 1).
      2. `adams.inc` of order k, AB_k/AM_k coefficients for k = 4, 6, 8
         (tabulated or from the generating functions); an input
         `adams_order` (default 4, keeping parity), on the .dat line
         with the mesh parameters.
      3. The bootstrap from k: k start points (the series, with RK4 on
         the model equation for relativistic points, whose steps are
         far below al); the inward re-integration (12 points now) and
         the exponential inward start of k points; a matching window
         long enough for k-point histories.
      4. Check that the model potential's start error stays below the
         AM_k truncation at high order; extend the series radius or the
         re-integration if not.
      5. One binary, amesh 1.003 / 1.006 / 1.012 × order 4 / 6 / 8:
         eigenvalues, kernels and V_loc (not projectors) should converge
         as al⁴, al⁶, al⁸; then decide whether a coarser default mesh is
         safe.
      6. No step adaptivity: a fixed mesh and a fixed operation count
         keep the solution C∞ in E and V.
      7. With it, the audit's §6c item 2: the AE quadrature (radint,
         vout's eeel and its AM4 Hartree potential) by `lm_cumint`, whose
         head must then take a real power (2γ−2 for relativistic s
         densities; an integer one undoes C-3). Alone it does not pay,
         measured 2026-10-22 on 14_Si, 07_N, 74_W (sr) and 80_Hg (fr),
         amesh 1.024 … 1.00075 against 1.00075: the AE total energy
         converges as al⁴ in both builds, the quadrature is 5–6 % of
         the error at every mesh, and the ODE is the rest.

         | atom  | build     | 1.024   | 1.012   | 1.006   | 1.003  |
         |-------|-----------|---------|---------|---------|--------|
         | 80_Hg | AM4 + 4th | 1.19e-3 | 7.31e-5 | 3.52e-6 | 6.6e-7 |
         | 80_Hg | lm_cumint | 1.13e-3 | 6.94e-5 | 3.29e-6 | 6.5e-7 |
         | 74_W  | AM4 + 4th | 8.92e-4 | 5.80e-5 | 3.73e-6 | 7.6e-7 |
         | 74_W  | lm_cumint | 8.43e-4 | 5.49e-5 | 3.53e-6 | 7.5e-7 |

         14_Si and 07_N reach about 1e-7, where the differences stop
         falling smoothly, at amesh 1.012; there the quadrature changes
         etot by 1.4e-6 at 1.024 and by less than 1e-8 at 1.006.

#### Numerics audit (2026-10-22)

`doc/numerics_audit.md` (a peer session's audit of the branch at
9c78c3a) ranks 18 findings.
C-1 to C-11 are resolved here, each checked against the code before
changing it:

- [x] C-1, psatom's occupation ramp (fa0b95e). fat1 + 1·(fat2 − fat1)
      need not round to fat2, and the convergence test compared the two
      exactly: a test configuration 2 → 0.1 was dropped. Confirmed by a
      run; the occupations are fat2 itself from the 50th iteration; test
      `occupation_ramp_sr_14_Si`.
- [x] C-2, aeatom's ierr not tested for the reference configuration
      (fa0b95e): now an error. The audit's other half, ae_converged's
      stop on a residual stalled below 1e-9, is kept: it is a
      convergence test at the floor of the residual, with a warning,
      not an iteration limit.
- [x] C-3, the head below rr(1) of the AE energy integrals (260d2af):
      the density's power through its first two nodes added to the
      nominal one (`lm_head`). Hg's kinetic energy no longer depends on
      rr(1) (22675.6866 → 22675.7192 at the default mesh, 22675.7193 at
      rr1 → 0); the energy terms of heavy atoms move by up to 0.049 Ha,
      totals by ≤3.3e-7 Ha, kernels not at all.
- [x] C-4, the Dirac solvers' start: r^γ alone at five points, and
      PECEC. `dirac_m` (shared by ldiracfb, ldiracfs) starts from the
      series at the origin for −Z/r + d, integrates PECECE and the first
      points again inward, as `sr_flux_m`; ldiracfb's inward integration
      is PECECE too. Hg's total energy at its default mesh was 8.4e-6 Ha
      from its rr1 → 0 value, now 6e-7; fr outputs only: totals ≤7.9e-6
      Ha, AE eigenvalues ≤1.3e-6 relative 1.2e-8, kernels ≤1e-7.
- [x] C-5, the free tail's width: fixed in d603e40
- [x] C-6, modcore_budget's silent limits and roots at a jump. Partly
      as the audit says: the fill-type switch is already in the list of
      discrete events above. `zeroin`, `fmin` and `gauss_legendre` stop
      with an error on their iteration limits, as does modcore3's
      Newton refinement (its Nelder–Mead start keeps its warning: its
      point only seeds Newton). A Brent root whose bracket closes on a
      jump (|f| above 1e-6 of the bracket's values at both ends, a
      change of fill type) is reported with the jump, as G5 asks. No
      output changes; no reference input or model-core scan reports one.
- [x] C-7, the ghost test's barrier basis without the ionic switch:
      confirmed (`vlocpw` alone, while lschvkbb's reference energies
      use vp, switched to the ionic tail beyond rcmax). The basis now
      takes vp's difference from `vlocpw` on the nodes, interpolated, as
      `outwf_setup` does. Basis eigenvalues move by ≤9.8e-6 Ha (73_Ta_r),
      most by 1e-7 to 1e-6, below the error the wall at 3 rc and the
      basis cutoff put in the margins (on the matched bound rows 12
      agree better, 29 slightly worse)
- [x] C-8, ghost_r's hard entry at 0.05 Ha: states now enter from 0.025
      Ha above the reference, their ⟨r⟩/rc divided by a C∞ weight that
      reaches 1 at 0.05 Ha, in the margin and the GHOST(+) test alike
      (a first try, the weight on 0.05–0.10 Ha, hid 56_Ba's GHOST(+) at
      0.0588 Ha from the margin while the flag still showed it). States
      at or above 0.05 Ha are treated as before, so no margin or
      GHOST(+) of the references changes. The audit's other cases are
      discrete events: the pairing after a GHOST(−) skip is the ghost
      event itself, the π·anint alignment of run_phsft can only change
      where the phases differ by π/2 (a failed log-derivative match),
      and the basis cut where lschpsbar fails is reported
- [x] C-9, der2exc's occupation step 0.15: checked, kept. Its
      truncation error is 3.5e-6 of d2Exc (14_Si; 3.3e-7 for 83_Bi,
      steps 0.15 to 0.001), and icmod 7's r_s and mismatch do not move
      at printed precision for any step from 0.15 to 0.001, nor for the
      slope's difference step from 1e-2 to 2.5e-3. Negative densities
      do occur, in icmod 7's probes (14_Si_hm, 328 points inside rcut),
      not icmod 1–4 as the audit has it; they persist at a step of 0.001
      (221 points), so they are where the density is negligible, not
      made by the step. A step of 0.01 was tried: icmod 4's optimum
      moved (74_W_r V_loc by 6e-5) and 40_Zr's Newton refinement no
      longer converged, so icmod 4's objective, as modcore3 samples it,
      is not smooth at that scale.
- [x] icmod 4: why modcore3's optimum depends on der2exc's step at 1e-5
      (the clipped near-zero densities are the first suspect); then
      analytic f_xc (libxc can give it) in place of the differences.
      Found (2026-10-29): PZ's step at rs = 1 (see the 14_Si residual
      above), not the clipped densities. A point whose density crosses
      rs = 1 within the differences takes a spike of the step over h,
      narrower but not smaller as h falls, which the mesh samples
      unevenly. 40_Zr's prefactors (amplitude / scale) by the step h:
      PZ 2.3788 / 1.5525 at 0.15, 2.3832 / 1.5514 at 0.05, Newton failing
      from 0.01; PZ with its exact coefficients 2.3777 to 2.3802 from
      0.15 to 0.001; PW92 (−001012) 2.3757 / 1.5528 from 0.003 to 0.15
      (and Newton failing at 0.001, where the differences' round-off,
      ε/h, passes its stopping test; 38_Sr_lxc's PBE from 0.01). The
      AE matrix's asymmetry measures it: 1e-6 for PZ at 0.15, 5e-8 with
      the exact coefficients, 1.5e-8 for PW92. Done: der2exc takes
      dV_xc/dn_k = f_xc ρ_k from `fxc_dispatch` (analytic for the
      built-in LDAs and Libxc's LDAs and GGAs, the GGA's response on the
      mesh with derivs; differences at 0.15 for the built-in PBE or a
      Libxc functional without fxc). 40_Zr: 2.3801 / 1.5522, PW92 and
      38_Sr_lxc unchanged to the printed digits, the matrices symmetric.
      Unit test `fxc`: each f_xc against 4th-order differences of the
      potential, to 6e-11 of max |dV| (LDAs) and the matrix element to
      3e-13 (all, the GGA included). PZ's exact coefficients were not
      taken: they change the functional (19 references, total energies
      by up to 7e-4 Ha), and QE's and libxc's PZ keep the rounded ones
- [x] C-10, eresid's q spacing qmax/(nq − 1), nq an integer: confirmed.
      The grid is now q_i = qmax − (nq − i) dq with dq = 2π/(150 rc)
      exactly, from below dq up to qmax: the spacing is continuous and
      a node that appears enters at q = 0, where the integrands vanish
      as q⁴. Its effect is at the quadrature's floor: the Ecut crossings
      move by ≤1.3e-11 relative, as much as doubling the q points moves
      either grid
- [x] C-11, V_l near the origin by cancellation: confirmed. Inside rc,
      pspot formed ½(l(l+1)u/r² − u″) from two terms that go as r^(l−1)
      and cancel to r^(l+1); for u = r j_l(qr) it is ½q²u exactly
      (Bessel's equation), now used (and in chix). Only V_l at the first
      nodes changes: the psp8 value at r = 0 (written from rr(1)) of
      07_N_dl's V_loc (l = 2) was 6.4e-7 Ha off the even continuation
      of V_loc from r ≥ 0.01, now 5.6e-9; elsewhere ≤6e-12

### Decisions

These change behaviour beyond numerics. Decided on review (2026-10-20).

- **D-1 wellstate overrides** (R-12). Revised 2026-10-21: kept by default,
  behind an input flag. Today a bound state of the requested n whose
  outermost peak is inside 0.75 rc replaces the input `ep` with its own
  energy, one inside 1.5 rc forces `ep = 0.05`, and a negative `ep`
  becomes 0.25. The first decision was to remove them and always use the
  input `ep`; on the references that misbehaved: the "moderately
  localized" 6p of 74_W_r, 80_Hg_r, 79_Au_lxc, the 2p of 03_Li and the 5d
  of 56_Ba then had to be raised to `ep` by a wall 0.17 bohr past rc
  (74_W_r's solution was inconsistent at rc), and 56_Ba's `ep = 0.0` (a
  placeholder of the old convention) gave a state at exactly E = 0, on
  which the diagnostics fail. So the overrides stay as the default, and
  `ep_as_given = true` in `[pseudopotentials]` (or on the `.dat` lmax
  line, follow-ups) selects the other behaviour: `ep` is always used, a
  bound state at or above `ep` is an input error, and one below `ep`
  inside 1.5 rc is raised to it by a wall as close as that takes.
- **D-2 fr psp8 truncation** (O-2). All terms are written. The files
  then hold up to 2·nproj terms per l where they hold fewer now.
- **D-3 Basis wave vectors** (Q-5). The method stays: the roots are found
  robustly and the crossing of the lowest root through 0 is a discrete
  event.
- **D-4 The vr window** (R-8, R-11). As PseudoPotentialGenerators.jl
  (`relativistic_taper`): subtract f·Δ[u_AE]/2, where Δ = u″ − [l(l+1)/r²
  + 2(V − ε)]u is the difference of the all-electron u″ from the
  Schrödinger one, instead of f·v_rel·φ̃ with v_rel = Δ/2u. It never
  divides by u, so it needs neither `vrel`'s threshold nor its
  neighbour average at nodes; it needs no knowledge of the equation, so
  in fr it holds the spin-orbit term too; and it is evaluated at the true
  rc. It differs from today by f·v_rel·(φ̃ − u_AE) = O((0.05 rc)^ncon).
  The window f = (1 − x²)², x = 20(r/rc − 1), stays: it is C¹ at 0.95 rc,
  a fixed kink in r that moves continuously with rc, so the quadratures
  split there and the dependence on the inputs is smooth. Its shape is
  arbitrary (PseudoPotentialGenerators.jl TODO C4: the narrow window
  gives χ high-q content and B an asymmetry of up to 3e-3 of max|B|), but
  changing it changes every sr and fr pseudopotential and is not a
  smoothness question; it is left out of this plan.
- **D-5 icmod 6/7 matching order.** An option, `match_order` in
  `[model_core_charge]` (2 or 4, default 2, so the results do not
  change). With 4 the fill matches the core charge to C⁴ at r_s and the
  3rd and 4th derivative columns of psp8 are continuous there. The
  minimum tail keeps 4 of its 9 coefficients free; the Martins fill
  becomes exp(a₀ + a₁r² + … + a₄r⁸), whose monotonicity on [0, r_s] is
  decided by the roots of its derivative. C⁴ has shown no practical
  benefit in the tests so far; the option is for smoothness. Deferred
  (2026-10-27): with S-10 dropped, nothing needs it now.

- **D-6 Exchange-correlation** (F-1, F-2). The xc implementations are
  not changed; a variant would be a new iexc, not needed so far. The
  published PZ81 coefficients
  (libxc's XC_LDA_C_PZ, and the LDA of the plane-wave codes) leave a jump
  in V_xc at r_s = 1; the exactly continuous ones (C and D from the
  continuity of ε_c and dε_c/dr_s, libxc's XC_LDA_C_PZ_MOD) remove it
  but change the functional: on the 30 iexc 3 references, AE total
  energies by up to 9e-4 Ha and kernels by up to 3.6e-4. The jumps stay.
  Planned, not for now: integrate around them. A registry of the known
  breakpoints of each functional (a jump or kink in ε_xc or V_xc at a
  density, such as PZ81's r_s = 1, or at a reduced gradient), for the
  common LDAs and GGAs of solid-state work and for every one with a local
  implementation here (iexc 1–7); the radii where the density crosses a
  registered value are breakpoints of V, so the solvers' integrations,
  the quadratures and the interpolants split there (as at rc) and the
  jump of V_xc enters the Taylor data of `ode_rc_der` and the matching
  as a known step. For a functional with no entry (libxc in general) the
  fallback is today's treatment, integrating through, with its noise
  noted. Not needed yet: the measured effect is at the solvers'
  discretization error near the jump (1e-8 Ha in a 3p matched at it).

- **D-7 The tail of the well states** (phase 8). Decided 2026-10-22,
  after the validation: (a) behind `well_tail`, "well" the default, as
  D-1 (the optimizer sets "free"); (b) the windowed growing solution;
  (c) the eigenvalue of the j = l+1/2 well. The well states with the
  free tail are normalized to 1 inside rc (the free solution has no
  finite norm, and any other choice depends on the window), so their
  convergence tables are on that scale: 1 to 20 Ha higher at 1e-4 than
  with the well's normalization, 119 Ha in 73_Ta.
  - (a) Default: the free tail by default, or behind an input
    (`well_tail`) with today's well tail as the default, as D-1 was
    decided. The free tail removes an arbitrary choice and the node
    jumps, but it moves the kernels of every well state by up to the
    size of that arbitrariness (several % in 14_Si l=1).
  - (b) Negative ep: the windowed growing solution, with L from the
    bound in phase 8 (no branch at ep = 0), or the well kept for
    ep ≤ 0. Keeping the well puts a policy branch (g) at ep = 0.
  - (c) The fr code's j = l−1/2 energy: kept as the eigenvalue of the
    j = l+1/2 well (the kernels then depend on the well through it,
    smoothly), or the input ep for both j (a further change to the
    method).
- **D-8 The confined well** (`well_tail = "confined"`). Decided
  2026-10-23: a third tail, beside "well" and "free", for well states
  that stay confined. Hamann's well rises to ep + 0.5 Ha and holds only
  the states below that; as ep rises its wall steepens until no state of
  that n exists, and the node loop moves the state to n + 1, a jump in
  its tail and so in the kernels (14_Si's l = 0 at debl 0.84: 1.5e-2, 60
  times the median step; Pb's p state flips 9/10/9/10 between debl 4.7
  and 5.1). The confined wall rises without bound, x³ e^(−0.01/x) Ha, and
  wellstate keeps the n it is given: E(rwell) runs monotonically from
  the atom's own state of that n to infinity, so every ep above it has
  one width, and nothing switches. The n of each l's first well state is
  an input (`well_n`, its default the n the well states start from); a
  tight wall adds content above qcut, which the well's-share line
  reports (Pb's p at 4.57 Ha: kernel 1.3 from the free tail's at n = 7,
  7e-3 at n = 9). The states are solved on the continued mesh
  (well_mesh_m) and normalized inside rc, as with the free tail.

### Bugs found on the way

To fix in the phase that touches the code, or in phase 7.

- icmod 7 at a low budget (2026-10-23, found when the default budget
  became that of the ground-state waves): 14_Si_hm at budget 11.5 Ha,
  slack 0.69 to 0.71, takes r_s = icmod 6's smallest fitting radius, but
  the cutoff the model core needs there (7.934 Ha) is below the one it
  must fit (7.941 Ha), so r_s is not the budget's root but some other
  limit of the search, and its steps along slack vary by up to 1e-3
  (smooth in the middle of the span). The slack scan of 14_Si_hm keeps
  the budget of before ("max_all", 33.9 Ha) until this is understood.

- `vpinteg` skips mesh points 1–3 (harmless, about 1e-12; removed in
  phase 7, no caller left).
- `sbf8` with nm = 3 or 4 leaves `sb_out(2)` (j_1) unset (latent: `eresid`
  and `sbf_basis_con` read only `sb_out(nm)`; fixed in phase 3, Q-10).
- `pspot` kinetic-energy origin term has the wrong power (diagnostic).
  Fixed in phase 7: ro^(2l+3)/(2l+3); the r-space kinetic energy moves
  by ≤1.4e-7 relative.
- `sbf_basis_con` origin term divides by rr(1)^l instead of rr(1)^(l+1)
  (negligible; gone in phase 3, R-5: exact overlaps).
- `const_basis` leaves components below its threshold unscaled instead
  of zeroing them (fixed in phase 3, Q-8: an error instead).
- `eresid` leaves `eresidmat(:,:,1)` unset if qcut sits exactly halfway
  between grid points (gone in phase 3, Q-2).
- `modcore` a0 bisection: `aco` stale when the loop runs out; degenerate
  bracket when the first a0 is already monotone (gone in phase 5, M-3).
- `modcore2`, `modcore3`: `rhomod(1,3:5)` set to 0 by `dpnint` (rr(1) ≠ 0;
  gone in phase 5, M-5).
- `modcore3`: `xx` uninitialised if every grid value is ≥ 10 Ha; dead
  `rcross` loops with a stale argument (gone in phase 5, M-7). The first
  fixed in phase 7: the grid's minimum starts at `huge`.
- `gg1cc` series: sign of the x⁶ term (about 1e-17; gone in phase 5,
  M-5: `modcore_m`).
- `doc/core_correction.md` describes a blending (`fcrossover.f90`) that
  the code does not do. Fixed in phase 7: it describes the taper, and
  Newton after Nelder-Mead.
- `lschfb`, `lschvkbb`, `ldiracfb`: `up`, `qq`, `gup`, `fup` beyond `nin`
  left unscaled (not latent: `up(2281)/(al r)` = 0.76 for hydrogen 1s on
  the default mesh, where the state is 1e-16). Fixed on oncvpsp-opt
  (752a6e5), with `lschvkbbe` and `lschpb`.
- `lschfb`, `lschfs`, `ldiracfb`, `ldiracfs`, `boundstate`,
  `scatterstate` document `up` as d(uu)/dr; it is d(uu)/dn = al·r·d(uu)/dr.
  Fixed in phase 7, with the other solvers, `wellstate` and `vkboutwf`.
- `lschpsbar` 204 uses `cf(mbar)` for `cf(ii)`. Fixed in phase 7: the
  ghost test's basis eigenvalues move by ≤1e-6 (32_Ge, 4.6e-7
  relative), its margins by ≤3e-7.
- `wellstate` 157 tests `ierr<0`, which lschfb never returns, so box
  states (`ierr=10`) can enter the bracketing. Fixed in phase 7: they
  still steer it (in 56_Ba, 14_Si_r, 74_W_r one trial each is a box
  state, at an energy in line with its neighbours'), but the state taken
  must have converged. No output changes.
- `psatom` can declare convergence during the occupation ramp (fixed in
  phase 2, S-5), and its ramp need not end on the target occupations
  (fixed by audit C-1).
- `lschpb` has no caller (removed in phase 7, with `lschps`, `lschkb`).
- The `oncvpsp_main` comment on the rc rounding says "nearest" and "5
  decimals"; it is "next" and 6 (gone in phase 3, R-1: rc as given).
