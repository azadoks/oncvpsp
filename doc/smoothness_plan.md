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
| R-1 | `oncvpsp_main` 462–465: rc truncated to 6 decimals | b | W | remove |
| R-2 | `oncvpsp_main` 575–588: rc → first mesh point beyond | a | W | remove; `irc` becomes only the last mesh point ≤ rc, for bookkeeping |
| R-3 | `rcmax` = max of snapped rc | a | W | max of continuous rc |
| R-4 | `wf_rc_der`: 7-point stencils applied 4 times at `irc` | a, e | W | ODE-recursion derivatives at rc |
| R-5 | `fpovlp`, `sbf_basis_con`: overlaps to `irc` with the (23,28,9)/24 end rule | c | W | `lm_integral` for all-electron overlaps; `sbf_overlap` (Lommel) for Bessel ones |
| R-6 | `pspot`: inside/outside by mesh index; V_l joins V_AE at irc+1 | a, d | W | pieces with extensions; values by side of rc |
| R-7 | `vploc`: v′, v″, v‴ at irc(5) by three nested 5-point stencils, polynomial for ii ≤ irc(5); round-off grows as al⁻³ (Fe: V_loc near 0 moves 3.6e-6, then 1.8e-5, as the mesh is halved twice) | a, e | W | derivatives of V_AE at rc(5) from a least-squares fit weighted by a C-infinity bump of fixed physical width (points enter and leave with zero weight, so the fit is smooth in rc(5) and does not amplify the noise as the mesh refines); polynomial evaluated analytically |
| R-8 | `run_vkb`: support cut at irc; window `20(r/rr(irc) − 1)`, ff = (1−x²)², C¹; ff·vr·φ̃ with vr = Δ/2u | a, d | W | true rc; subtract ff·Δ[u_AE]/2 (D-4); same window shape; quadratures split at 0.95 rc, rc(5), rc |
| R-9 | `vpinteg` to `irc`; points 1–3 skipped | c | W | piecewise `lm_integral` |
| R-10 | `vkboutwf`, `lschvkbb`: rc found as the last nonzero vkb. Since e809a29 the overlaps end at irc and the Adams steps past irc use the source-free u″ at their back points (no O(al²) from the kink), both assuming the support ends on a node | a, e | W | rc passed in; one fixed [0, rc] rule shared with run_vkb, with the partial cell [rr(irc), rc]; the integration stops at rc exactly (partial step) and restarts beyond it from (u, u′) at rc, so no step straddles the kink |
| R-11 | `vrel` window 0.5–1.5 rr(irc); 1e-10 threshold on u; vr = 0 points replaced by the average of their neighbours | a, g | W | replaced by Δ = u″ − u″_Schrödinger of the all-electron wave on [0.95 rc, rc], from the ODE (D-4): no division by u, no threshold |
| R-12 | `wellstate`: well onset at rr(irc+5), x³ onset (C²), inside the ±12-point stencils of `wf_rc_der` (uord 4, 5 of well states O(al²), O(al) off); thresholds 0.75, 1.5 rr(irc) override `ep`; negative `ep` reset to 0.25 | a, d, g | W | onset at rc·amesh⁵ (continuous), C-infinity onset; the derivatives at rc come from the ODE (R-4), which sees V at rc only, so the onset beyond rc cannot reach them; overrides removed, the input `ep` is always used (D-1) |
| R-13 | `ionic_window`: ra = rcmax (snapped); rcore = first mesh point where `|vion+zion/r|` < 1e-10; rb = max(rcore, 2ra) | a, g, h | W | continuous ra; rcore by `lm_root` of ln\|vion+zion/r\| − ln tol; smooth max |
| R-14 | `run_diag` `mchf = max(irc)+5`; `run_phsft` rxpsh → irphs; `run_ghosts` mbar; `run_plot` n1..n4 | a | – | phase 6 |
| R-15 | `upfout` `nrlproj`, `cutoff_radius_index`; `linout` rc+2·drl format switch | a, g | W | β exactly 0 beyond rc (phase 4); indices then harmless |

#### Wave vectors, residual kinetic energy, optimisation (phase 3)

| id | where | cat | W | decision |
|---|---|---|---|---|
| Q-1 | `eresid` 268–274: qcut → nearest node of a q-grid with dq = 2π/(150 rc) | a | W | E_r at exactly qcut |
| Q-2 | `eresid`: nq = qmax/dq + 1, top 3 intervals dropped; irclin, mmaxlin integer, tail cut at \|u\| < 1e-10 | i, c, g | W | fixed-count quadrature on [qcut, qmax]; inside piece by Lommel at each q; outside piece on a fixed grid to a fixed radius |
| Q-3 | `run_optimize` 158–164: qmax = dq·int(2·qroot(nbas)/dq), max with 20, branch on qcut | i, g | W | qmax continuous (smooth max), or a fixed multiple plus an analytic tail |
| Q-4 | `run_optimize` cvgplt/ecurve: linear interpolation assuming spacing dqout while eresid moved the q's | e, h | W (icmod 6/7) | E_r at exact q; crossings by root-finding on E_r(q), whose q-derivative is known |
| Q-5 | `qroots`: scan dq = 0.02, pole filter \|f·f_last\| < 1, bisection to \|f\| < 1e-12 | h, f | W | done in phase 2: the bisection's root Newton-polished on g(q) to round-off; the scan and pole filter stay (D-3) |
| Q-6 | `optimize`: bisection to a 1e-10 bracket | f | W | done in phase 2: Newton polish of the norm constraint after the bisection |
| Q-7 | `optimize`: hard-case sign of x₁ | h | W | intrinsic (two global minima); keep, warn |
| Q-8 | `const_basis`: absolute 1e-8 singular-value threshold, values below left unscaled | g | W | relative threshold; error if rank-deficient |
| Q-9 | `wellstate`: rwell by ×/÷0.8 then bisection to \|E − ep\| < 1e-8; ep := E | f, h | W | done in phase 2: Illinois regula falsi on E(rwell) − ep once bracketed, to round-off; a bracket closed on a residual above 1e-10 is an error |
| Q-10 | `sbf8`: l = 2, 3 closed forms lose eps/x^l at small x | – | W | series below x ≈ 0.5 for l ≥ 2 |

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
| S-8 | matching point `mch` at the turning point index; `lschfb`'s mch and `inward_start`'s nin are discrete in E (about 1e-12 in the all-electron eigenvalues) | h | W | below 1e-13 at amesh 1.006 (measured, emulation); keep. The all-electron atom does not depend on the scanned inputs, so its noise is a constant offset in a scan |
| S-9 | `sr_flux_m`: series/RK4 switch at rs, RK4 step count by ceiling | g, i | W (~1e-13) | done: rs = a/(8 max(1, \|b\|)), constant in E for all bound and scattering states (\|b\| ≤ 1), so the series/RK4 branch and the step count are too |
| S-10 | `exc_libxc`: ρ′, ρ″ of ρ_v + ρ̃_c by mesh stencils (`derivs`) across the kinks of the model core (r_s C², r* C⁴, the Teter cutoff C¹) and of ρ_v at rc | d, e | W (GGA) | ρ′, ρ″ from the pieces: the model core from its derivative columns, ρ_v from u, u′, u″ of the pseudo-waves |

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
- [ ] R-1, R-2: rc continuous; `irc` the last point ≤ rc
- [ ] R-4 matching data at rc from the ODE; R-5 overlaps; R-6 pieces
- [ ] Q-1…Q-4 residual kinetic energy at exact qcut, fixed quadratures
- [ ] R-7 vploc at continuous rc(5); R-8…R-10 run_vkb, vpinteg,
      vkboutwf on a shared piecewise rule; R-11 vrel; S-1
- [ ] R-12 wellstate onset, overrides removed (D-1); R-3, R-13 ionic
      window
- [ ] Q-8, Q-10
- [ ] T2 rows rc, rcloc, qcut at G1; T3, T4 recorded
- [ ] References patched

#### Phase 4: output

- [ ] O-1 output meshes from the pieces; β ≡ 0 beyond rc
- [ ] O-2, O-4, O-5 (after D-2)
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

### Decisions

These change behaviour beyond numerics. Decided on review (2026-10-20).

- **D-1 wellstate overrides** (R-12). Removed. Today a bound state of
  the requested n whose outermost peak is inside 0.75 rc replaces the
  input `ep` with its own energy, one inside 1.5 rc forces `ep = 0.05`,
  and a negative `ep` becomes 0.25. The input `ep` is now always used,
  and the well is always built. The well only raises the energy, so an
  `ep` below the free bound state of that n has no well state: that is
  an input error, reported with the bound energy. 13 of the 51
  references take an override today and change: 03_Li (nr, sr, fr),
  20_Ca (sr, fr), 38_Sr_lxc (sr, fr), 56_Ba, 57_La, 74_W_r, 79_Au_lxc
  (sr, fr), 80_Hg_r. In all of them the input `ep` is above the bound
  state (57_La 4f: −0.10 against −0.111), so none becomes an error.
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

### Bugs found on the way

To fix in the phase that touches the code, or in phase 7.

- `vpinteg` skips mesh points 1–3 (harmless, about 1e-12).
- `sbf8` with nm = 3 or 4 leaves `sb_out(2)` (j_1) unset (latent: `eresid`
  and `sbf_basis_con` read only `sb_out(nm)`).
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
