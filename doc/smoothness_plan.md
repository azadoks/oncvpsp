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

`scripts/smooth_scan.py` steps one input through 41 values over two
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

For a value at an arbitrary radius the generator today uses `dpnint`, an
8-point Lagrange interpolant whose stencil jumps as the point crosses a
mesh node: its value is continuous there but its derivative is not.
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
  input keeps its node count and scales its nodes. This covers the q-grid
  of the residual kinetic energy on [qcut, qmax], the Gauss panels of
  modcore_budget, and the r-grid of eresid. This removes integer
  truncations like `nq = qmax/dq + 1`.

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
| r_s | model core, icmod 6/7 | C² by default, C⁴ with `match_order = 4` (D-5) |
| 3/b | Teter cutoff, icmod 2–4 | C¹ today; removed in phase 5 |
| r_w | wellstate well onset | C² today; C-infinity after phase 3 |
| ra, rb | ionic tail switch | C-infinity (`ionic_switch`) |
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
that ends on its iteration limit is an error, not a result.

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
| R-4 | `wf_rc_der`: 7-point stencils applied 4 times at `irc` | a, e | W | ODE-recursion derivatives at rc |
| R-5 | `fpovlp`, `sbf_basis_con`: overlaps to `irc` with the (23,28,9)/24 end rule | c | W | `lm_integral` for all-electron overlaps; `sbf_overlap` (Lommel) for Bessel ones |
| R-6 | `pspot`: inside/outside by mesh index; V_l joins V_AE at irc+1 | a, d | W | done in phase 3: `pieces_m` (a function of up to 4 segments, each piece held on the whole mesh as its smooth extension, integrals per segment by `lm_integral`; unit test `pieces`); `pspot` returns the inside piece and its seed (ε − T)φ̃ on the nodes to past rc; `irc` is the last node ≤ rc, so the node values are by side of rc |
| R-7 | `vploc`: v′, v″, v‴ at irc(5) by three nested 5-point stencils, polynomial for ii ≤ irc(5); round-off grows as al⁻³ (Fe: V_loc near 0 moves 3.6e-6, then 1.8e-5, as the mesh is halved twice) | a, e | W | done in phase 3 (`vploc_m`, unit test `vploc`): V, V′, V″, V‴ at rc(5) from `lm_derivs` (smooth in rc(5)); the polynomial evaluated analytically at any radius, its extension past rc(5) a piece of the projectors' integrals. The bump-weighted least-squares fit was tried: of degree 8 on the nodes within 0.25 bohr its V‴ was good to 2e-6 for V = −2/r + e^−r at 1.7 bohr (the truncation of a fit of fixed width for a potential with the derivatives of −Z/r), `lm_derivs` to 1e-9; its al⁻³ round-off (2e-11 on the default mesh) is a fixed pattern on the nodes, smooth in rc(5) through the interpolant |
| R-8 | `run_vkb`: support cut at irc; window `20(r/rr(irc) − 1)`, ff = (1−x²)², C¹; ff·vr·φ̃ with vr = Δ/2u | a, d | W | done in phase 3: the projectors are `piecewise_t` with segments split at 0.95 rc, rc(5) (and 0.95 rc(5) for lloc<4) when inside, and rc; on each segment χ = seed − V_loc^piece φ̃ − f Δ/2 (D-4, `ae_delta`), every term its smooth extension; the B matrix, normalizations and overlaps by `pw_integral`/`pw_product_integral`; the orthonormal combinations act on the pieces; `vkbpw` returned beside the node values. On the T3 inputs the fr kernels move by ≤1.1e-4 (the Δ form), sr/nr by ≤1e-6 |
| R-9 | `vpinteg` to `irc`; points 1–3 skipped | c | W | done in phase 3: `run_vkb` (R-8), `sr_so_r` and `run_ghosts` on the pieces (`pw_product_integral`, `pw_integral`), `vkboutwf` per segment on its grid (R-10); `vpinteg` remains only in `der2exc` (phase 5) |
| R-10 | `vkboutwf`, `lschvkbb`: rc found as the last nonzero vkb. Since e809a29 the overlaps end at irc and the Adams steps past irc use the source-free u″ at their back points (no O(al²) from the kink), both assuming the support ends on a node | a, e | W | done in phase 3 (`outwf_m`, unit test `outwf`): no integration on the nodes is smooth, since a step whose stencil spans a kink has an error that depends on where in the cell the kink falls (a restart on the nodes at the first node past the kink still jumps by one local truncation error, 1e-10, when the kink crosses a node); so each segment between kinks is integrated on its own grid, uniform in ln r with the mesh step, starting exactly at the segment's lower bound and running 32 points past its upper one: the first on the nodes from the power series at the first node, each later one from (u, u′) interpolated from the grid before; every segment primed by 3 Runge-Kutta steps in 4 sub-steps (so the Adams scheme starts from data exact to its own error), V and the projectors interpolated onto the grid from their pieces, the overlaps per segment by the cumulative 8th-order rule to the bound, the solution interpolated back to the nodes. Every operation then has the same form wherever the kinks sit. The method stays 4th order (the log derivative at r = 3 of the unit test's model falls 15.6, 15.4, 14.7 times on three halvings; truncation 2e-7 at the default step). The potential's pieces (`run_vkb`) and the projectors' are threaded to every solver; `lschvkbb`'s window stays beyond the last kink of either. The weights of the stencils and the cell rule are computed in quad once (they were recomputed per integral, 60 % of the runtime); the first segment's overlaps, most of the work, are dot products with the weights of `lm_integral` to its bound (`lm_integral_weights`, computed once per setup), in a compensated sum of four lanes or, with `ONCVPSP_BLAS_DOT`, one `dgemm`; the interpolants' weights are shared by u and u′ and by V and the sources. Against e3a587f the runs take 74_W_r 8.7 → 10.2 s (9.5 with BLAS), 14_Si_r 2.8 → 3.0 s, 32_Ge 2.5 → 3.1 s, 08_O_nr 1.0 → 1.2 s. More accurate: in the nr references the reference configuration's pseudo-atom eigenvalues meet the all-electron ones to 1.1e-9 Ha (03_Li, was 2.6e-8), 8.3e-10 (14_Si, was 2.0e-9), and 32_Ge's 3d to 4.1e-7 (was 1.1e-6); none got worse |
| R-11 | `vrel` window 0.5–1.5 rr(irc); 1e-10 threshold on u; vr = 0 points replaced by the average of their neighbours | a, g | W | done in phase 3: `ae_delta` (unit test `ae_delta`), Δ = u″ − [l(l+1)/r² + 2(V−ε)]u at the nodes about the window from the first-order system (u, u′ = up/(al r), V′ by `lm_derivs`), no division by u; `vrel` removed. The local potential of lloc<4 keeps Δ/(2u) for its own taper (S-1), with an error if u has a node within reach of rc |
| R-12 | `wellstate`: well onset at rr(irc+5), x³ onset (C²), inside the ±12-point stencils of `wf_rc_der` (uord 4, 5 of well states O(al²), O(al) off); thresholds 0.75, 1.5 rr(irc) override `ep`; negative `ep` reset to 0.25 | a, d, g | W | done in phase 3: onset at rc·amesh⁵ (continuous), then at rc·`well_onset`, an input of default 1.006⁵ that does not depend on the mesh (T4 at the end of phase 3); the well cwell x³/(1+x³) e^(−0.01/x), the well of before times a C-infinity factor within 5 % of 1 from x = 0.2 on, so that the states found are those of before (with e^(−1/x³), a different shape, 14_Si's kernels moved 7 % and 11 inputs found another n at the wall's lower limit; with e^(−0.1/x), still 11; with e^(−0.01/x) only 40_Zr's l=1, whose wall sat at the limit). The `rc[0]` scan 8.9e-6 → 4.4e-12. The overrides stay by default; `ep_as_given` (D-1, revised) |
| R-13 | `ionic_window`: ra = rcmax (snapped); rcore = first mesh point where `|vion+zion/r|` < 1e-10; rb = max(rcore, 2ra) | a, g, h | W | done in phase 3: ra the largest given rc; rcore the root of the interpolant of (V_ion + zion/r)² − tol² between the mesh points where it crosses (`lm_root`; the square, not the log, as V_ion + zion/r may cross zero there); rb the C-infinity maximum of rcore and 2ra within 0.05 bohr |
| R-14 | `run_diag` `mchf = max(irc)+5`; `run_phsft` rxpsh → irphs; `run_ghosts` mbar; `run_plot` n1..n4 | a | – | phase 6 |
| R-15 | `upfout` `nrlproj`, `cutoff_radius_index`; `linout` rc+2·drl format switch | a, g | W | β exactly 0 beyond rc (phase 4); indices then harmless |

#### Wave vectors, residual kinetic energy, optimisation (phase 3)

| id | where | cat | W | decision |
|---|---|---|---|---|
| Q-1 | `eresid` 268–274: qcut → nearest node of a q-grid with dq = 2π/(150 rc) | a | W | done in phase 3: E_r at exactly qcut, the value at the next node of the q grid plus the partial cell by 8-point Gauss-Legendre with the transforms at its nodes |
| Q-2 | `eresid`: nq = qmax/dq + 1, top 3 intervals dropped; irclin, mmaxlin integer, tail cut at \|u\| < 1e-10 | i, c, g | W | done in phase 3: the inside piece by `sbf_overlap` at each q (exact); the outside piece on a linear mesh from rc of spacing 2π/(75 qmax), continuous in rc and qmax, to the last nonzero point of u (the node count changes where u is at machine precision; a spacing set by the end point moved every node, and the 1e-12 interpolation error of u onto them jumped together), u by `lm_interp`, integrated by the 8th-order rule (`ugrid_m`); the q grid of continuous spacing qmax/(nq−1) with 150 points per 2π/rc, accumulated from qmax by the 8th-order cumulative rule (Kahan), from the small values up |
| Q-3 | `run_optimize` 158–164: qmax = dq·int(2·qroot(nbas)/dq), max with 20, branch on qcut | i, g | W | done in phase 3: qmax = smax(smax(2 qroot(nbas), 20), 2 qcut), smax the C-infinity maximum (a+b+√((a−b)²+1))/2, within 0.5 bohr⁻¹ of the larger. The truncation at qmax stays, as a smooth bias. Considered and not taken: Parseval, E_r = T − (1/π)∫₀^q q'⁴F² (PseudoPotentialGenerators.jl `plane_wave.jl`), which needs no qmax and F below q only, but leaves E_r as the difference of two O(1 Ha) numbers, a round-off floor of ε·T/E_r ≈ 1e-11 at E_r = 1e-5 Ha (PPG's test notes a 1e-9 Ha floor); the direct sum from qmax down has the floor of E_r itself. An analytic tail beyond qmax (the transform's asymptotic form from the C^(ncon−1) matching) would remove the bias; phase 7 if wanted |
| Q-4 | `run_optimize` cvgplt/ecurve: linear interpolation assuming spacing dqout while eresid moved the q's | e, h | W (icmod 6/7) | done in phase 3: the curve of the optimized wave on the q grid (its coefficients returned by `optimize`), the crossings by Brent on the FH interpolant of ln E in the grid's index (`ug_cross_down`), in `run_optimize` and in `modcore_budget_table` (which takes every nstep-th node of the grid as its table). `dqout` is no longer used |
| Q-5 | `qroots`: scan dq = 0.02, pole filter \|f·f_last\| < 1, bisection to \|f\| < 1e-12 | h, f | W | done in phase 2: the bisection's root Newton-polished on g(q) to round-off; the scan and pole filter stay (D-3) |
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
| S-10 | `exc_libxc`: ρ′, ρ″ of ρ_v + ρ̃_c by mesh stencils (`derivs`) across the kinks of the model core (r_s C², r* C⁴, the Teter cutoff C¹) and of ρ_v at rc | d, e | W (GGA) | ρ′, ρ″ from the pieces: the model core from its derivative columns, ρ_v from u, u′, u″ of the pseudo-waves |

#### Exchange-correlation

| id | where | cat | W | decision |
|---|---|---|---|---|
| F-1 | `excpzca`: PZ81's two branches meet at r_s = 1 with the published, rounded C = 0.0020 and D = −0.0116, so ε_xc and V_xc jump there (V_xc by about 1e-5 Ha); the radius where r_s = 1 moves with the density. The solvers integrate across the jump: the all-electron eigenvalue of a state matched near it depends on the matching point by up to 1e-8 Ha (3p of Ti, matched 1 to 3 points inside its turning point) | d | W | kept (D-6): noted, not changed |
| F-2 | `excpzca`, `excvwn`, `excwig`, `exchdl`, `exc_relx`: V_xc and ε_xc set to 0 below ρ = 1e-20 (1e-15 for Wigner), a jump of up to 1e-7 Ha at the radius where ρ crosses it; `exc_off_pbe` below 1e-18, about 1e-6 Ha; libxc its own thresholds | g | W (tails) | kept (D-6): noted, not changed. Measured: computing the LDA functionals at every ρ > 0 moves the kernels by at most 3.3e-9 and the AE terms by 1e-10 (`exchdl` would need ln(1+1/x) and its ε_c by a series at large x, where the closed form cancels); PBE and libxc overflow as ρ → 0 and would need a smooth damping of f(n, ∇n) with its potential |

#### Model core (phase 5)

| id | where | cat | W | decision |
|---|---|---|---|---|
| M-1 | `modcore` ircc, `modcore2` ircc, `modcore3` ircc: crossover rhoc = fcfact·rhotps at a mesh point | a | W | `lm_root` of ln rhoc − ln(fcfact·rhotps); match at r* |
| M-2 | `modcore`: derivatives at ircc by stencils | a, e | W | `lm_derivs` at r* |
| M-3 | `modcore`: a0 by bisection to 1e-7, monotonicity checked only at mesh points | f, h | W | closed form: p' is affine in a0; critical points by companion matrix |
| M-4 | `modcore2`: x by bisection to 1e-7 | f | W | Newton on the monotone x·g'/g |
| M-5 | `gg1cc`, `gpp1cc`: Teter function cut to 0 at x = 3 (C¹); columns 4, 5 by stencils across it; `dpnint` patches near the origin (zeroing rhomod(1,3:5)) | d, e | W | no cutoff (it decays as x⁻¹⁰); analytic derivatives with series near 0; delete the patches |
| M-6 | `modcore3`: grid argmin then Nelder–Mead to an f-spread of 1e-4; returns a vertex | h, f | W | keep grid and NM to choose the basin; then Newton (trust region) on ∇f = 0 to round-off |
| M-7 | `modcore3` ircross, dead rcross loops | a | W | remove (analytic derivatives everywhere) |
| M-8 | `modcore_budget` ifrom: first q of a 0.1 grid above the budget; the only way slack enters | a, b | W | tested cutoffs start at q_b exactly; max over q by grid plus Brent |
| M-9 | `modcore_budget` r_s by integer bisection; irtop, ilo from snapped rc | h, a | W | bracket on the mesh, Brent on g(r_s) = max_q err − errmax |
| M-10 | `modcore_budget` icmod 7 candidates snapped with nint, chosen by sorting | b, h | W | continuous 1-D minimisation of the mismatch on [r_s6, r_top] |
| M-11 | `modcore_budget`: panel count ceiling(qhi−qlo), nxg = 20 + ceiling(qhi·r_s) | i | W | fixed counts |
| M-12 | `modcore_budget` outer integral from mesh index `is`; forward transform a plain sum across the C² joint | c, d | W | split at r_s, `lm_integral` |
| M-13 | `modcore_budget_table`: linear interpolation in q with an assumed spacing, last crossing | e, h | W | root of E_i(q) = errmax on exact E_r (Q-4) |
| M-14 | `oncvpsp_main` irct (PSML matching radius) | a | W (metadata) | write r* or r_s |

#### Output (phase 4)

| id | where | cat | W | decision |
|---|---|---|---|---|
| O-1 | `linout`, `upfout`: `dpnint` (switching 8-point Lagrange) of kinked arrays onto the linear mesh; ringing beyond rc of about 1e-4 relative | d, e | W | each output point from its piece by FH on that piece's extension; β ≡ 0 beyond rc |
| O-2 | `oncvpsp_main` 882–891 + `linout`: fr psp8 drops sr/so terms below 1e-3 of the first | g | W | write all; see D-2 |
| O-3 | UPF numbers e20.10, PP_R f10.4 | b | W | 17 digits (done, phase 0) |
| O-4 | `upfout` nrl extended to a tail-norm radius (1e-6) | g, i | W | length from rlmax only |
| O-5 | `oncvpsp_main` 997: vpuns = −zion/r where rho == 0 | g | W (~1e-10) | switch at a continuous radius with `ionic_switch` |
| O-6 | printed diagnostics at f12.6–f14.8 and d10.2 | b | – | machine-readable block at 17 digits (phase 6) |

#### Diagnostics (phase 6)

| id | where | cat | decision |
|---|---|---|---|
| X-1 | `run_diag`, `run_diag_sr_so_r`: matching at mchf = max(irc)+5; renorm to the last nonzero point | a, c | continuous r_m = max(rc_l, rc_loc)·amesh⁵; (u, q) by FH at r_m |
| X-2 | `run_phsft`: rxpsh → irphs; npsh = int(range/depsh); 2π unwrapping on any decrease; `vkbphsft` heuristic thresholds and nearest-grid alignment | a, i, g, h | exact rxpsh; fixed point count; branch-free angle (doubled-angle form, sign-invariant); offset at ep exactly |
| X-3 | `run_ghosts`: barrier at mesh index mbar; basis = states below ebcut (count jumps) | a, i | continuous barrier radius; fixed basis with a smooth energy cutoff; print the ghost margins |
| X-4 | `run_plot`: index ranges and thinning | a, h | fixed output grid (plots only) |
| X-5 | `run_config`/`psatom`: see S-5; occupation ramp may exit early | f, g | explicit criterion; ramp completes |
| X-6 | `coulomb_tail`: α_z and its cut changes by the rectangle rule on the nodes; the 1e-10 radius the node past the last one above it, which sits at the noise floor of the unscreening (R-10 moved it 4.71 → 4.58 bohr in 14_Si nr as max \|V_loc − V_ion\| beyond rcmax fell 2.0e-9 → 5.8e-10) | c, a, g | `lm_integral`, cut at the radii exactly; the radius by `lm_root` as in R-13, reported with the noise floor beside it |

#### Truly discrete (kept)

- Node counts in every solver, and so the principal quantum number of a
  reference.
- A bound state exists or not; `boundstate`'s box-state fallback
  (flagged, `ierr=10`).
- The bound/scattering branch of the diagnostics at epa = 0.
- A ghost exists or not (the margin is printed, G5).
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
when it is green. All of these carry the ctest label `smooth` until the
plan is complete, so `ctest -LE smooth` keeps the suite green on the
branch in between; the default test preset (`ctest --preset default`)
leaves them out until phase 7 drops the label.

- **T1 Unit tests of the primitives** (`tests/unit/test_logmesh.f90`,
  written, red against the stub `src/logmesh_m.f90` and
  `src/sbf_overlap.f90`). They check:
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
  `smooth_<program>_<input>_<key>`, written, red). Every row of the
  baseline table, and with them:
  - each rc of 14_Si;
  - `rcloc` with lpopt 1 and 5;
  - `qcut`;
  - `ep` of a bound and a well state;
  - `debl`, `dvloc0`;
  - `fcfact` for icmod 1, 2, 3;
  - `rcfact` for icmod 3;
  - the icmod 4 grid;
  - `slack` for icmod 6 and 7;
  - rc for nr 08_O and fr 52_Te;
  - `rxpsh` and the test eigenvalues, added in phase 6 with the
    full-precision diagnostics they need.

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
- [ ] T2 rows `dvloc0`, `debl`, `ep` (well) at G2: `dvloc0` is at
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
      (TOML; the text input takes the default), the onset as a multiple
      of rc, at least 1, default 1.006⁵. 1 is the natural limit: the
      well is flat to every order at its onset, so it can start at rc.
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
      The origin keeps `dpnint`'s limit (see the follow-ups). On the log
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
- [ ] O-2 (D-2): every term of the fr decomposition written (no 1e-3
      cut), so the diagnostics of `run_diag_sr_so_r` use them all too
- [ ] O-4: reopened. The UPF mesh from rlmax only would cut the
      pseudo-wave functions (PP_CHI) there: of their norm, 14_Si_UPF's 3p
      would lose 5.8e-2, 32_Ge's 4p 1.7e-1, 52_Te_r's 5p 1.0e-2 (the
      extension to where 1e-6 is left is what keeps them). The extension
      changes only how many points are written, not a value at any
      radius; to be decided before the step
- [ ] O-5: the switch of V_loc to −zion/r where ρ = 0 by `ionic_switch`
      at continuous radii. It reaches no written function inside rlmax;
      it moves the pseudo-atom tests
- [ ] References patched

#### Phase 5: model core

- [ ] Unit tests: crossover root, closed-form a0, Teter derivatives
      against series and finite differences, written first
- [ ] M-1…M-14; S-10 GGA gradients from the pieces
- [ ] `match_order` (D-5): input key, minimum tail and Martins fill to
      C⁴, T2 row `slack` with `match_order = 4`
- [ ] T2 rows fcfact, rcfact, icmod 4 grid, slack at G1
- [ ] References patched

#### Phase 6: diagnostics

- [ ] X-1…X-5; machine-readable diagnostics block (O-6), parsed by
      pyoncvpsp
- [ ] T2 rows rxpsh and the test eigenvalues at G1
- [ ] References patched

#### Phase 7: wrap-up

- [ ] The numerics above moved into the permanent docs
      (`program_overview.md`); this page kept as the record
- [ ] Dead code removed (`lschpb`, `lschps`, `lschkb`); the bugs below
      fixed where the phases did not already
- [ ] `smooth` label dropped; T2 part of the default suite

#### Phase 8: the kernels free of the well (follow-up)

A follow-up, started once phases 4 to 7 are done. It changes the method,
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
- [ ] D-7 decided
- [ ] Unit test `eresid_tail`; the independence, scan and mesh tests
      (red)
- [ ] `free_tail` and the window in `eresid`; normalization inside rc
- [ ] Validation on the references recorded
- [ ] References patched

#### Follow-ups (once everything above is done)

- [ ] The `.dat` input brought as far as reasonably possible to parity
      with the TOML input: the keys added on this branch (`ep_as_given`,
      `well_onset`, `upf_mesh`, the log mesh, …) readable from `.dat`
      too, or a documented reason why not
- [ ] The value at the origin of the psp8 and linear-mesh UPF files:
      `outmesh_m` keeps `dpnint`'s limit, the polynomial of degree 7
      through the next 8 output points (0.01–0.08 bohr). Not yet trusted:
      check it against the functions' known behaviour at r = 0 (r β → 0,
      ρ'(0) = 0, V_loc finite) and against a fit on the log mesh

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
  `ep_as_given = true` in `[pseudopotentials]` (TOML only) selects the
  other behaviour: `ep` is always used, a bound state at or above `ep` is
  an input error, and one below `ep` inside 1.5 rc is raised to it by a
  wall as close as that takes.
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
  benefit in the tests so far; the option is for smoothness.

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

- **D-7 The tail of the well states** (phase 8). Open; decide at the
  start of the phase.
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

### Bugs found on the way

To fix in the phase that touches the code, or in phase 7.

- `vpinteg` skips mesh points 1–3 (harmless, about 1e-12).
- `sbf8` with nm = 3 or 4 leaves `sb_out(2)` (j_1) unset (latent: `eresid`
  and `sbf_basis_con` read only `sb_out(nm)`; fixed in phase 3, Q-10).
- `pspot` kinetic-energy origin term has the wrong power (diagnostic).
- `sbf_basis_con` origin term divides by rr(1)^l instead of rr(1)^(l+1)
  (negligible).
- `const_basis` leaves components below its threshold unscaled instead
  of zeroing them.
- `eresid` leaves `eresidmat(:,:,1)` unset if qcut sits exactly halfway
  between grid points.
- `modcore` a0 bisection: `aco` stale when the loop runs out; degenerate
  bracket when the first a0 is already monotone.
- `modcore2`, `modcore3`: `rhomod(1,3:5)` set to 0 by `dpnint` (rr(1) ≠ 0).
- `modcore3`: `xx` uninitialised if every grid value is ≥ 10 Ha; dead
  `rcross` loops with a stale argument.
- `gg1cc` series: sign of the x⁶ term (about 1e-17).
- `doc/core_correction.md` describes a blending (`fcrossover.f90`) that
  the code does not do.
- `lschfb`, `lschvkbb`, `ldiracfb`: `up`, `qq`, `gup`, `fup` beyond `nin`
  left unscaled (not latent: `up(2281)/(al r)` = 0.76 for hydrogen 1s on
  the default mesh, where the state is 1e-16). Fixed on oncvpsp-opt
  (752a6e5), with `lschvkbbe` and `lschpb`.
- `lschfb`, `lschfs`, `ldiracfb`, `ldiracfs`, `boundstate`,
  `scatterstate` document `up` as d(uu)/dr; it is d(uu)/dn = al·r·d(uu)/dr.
- `lschpsbar` 204 uses `cf(mbar)` for `cf(ii)`.
- `wellstate` 157 tests `ierr<0`, which lschfb never returns, so box
  states (`ierr=10`) can enter the bracketing.
- `psatom` can declare convergence during the occupation ramp (fixed in
  phase 2, S-5).
- `lschpb` has no caller.
- The solvers' doc comments call `up` d(uu)/dr; it is d(uu)/dn = al·r·du/dr.
- The `oncvpsp_main` comment on the rc rounding says "nearest" and "5
  decimals"; it is "next" and 6.
