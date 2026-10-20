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
- [Decisions for review](#decisions-for-review)
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
interpolants with weights that fall off as |x − x_i|^−(d+1). Summing only the nodes
whose weight relative to the nearest is above machine epsilon (about 60
cells each side for d = 6) keeps it local and cheap without a visible
switch. **d = 6** (order 7). At al = ln 1.006 this is round-off for every
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
  comparable. Phase 1 measures both on 14_Si (PBE) and hydrogen
  (analytic), and the plan keeps the ODE route unless the measurement
  says otherwise. The finite-difference route is still the cross-check
  (T1).

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

  For α near β the first is evaluated by its Taylor series in (α − β),
  so there is no cancellation. The same formula gives the Fourier–Bessel
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
| 0.95·rc_l | the vr window of run_vkb | C¹ today; C-infinity after phase 3 |
| r* | model core, icmod 1 | C⁴ |
| r_s | model core, icmod 6/7 | C² (C⁴ possible, see D-5) |
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
| R-7 | `vploc`: 5-point stencils at irc(5)±6, polynomial for ii ≤ irc(5) | a, e | W | `lm_derivs` of V_AE at rc(5); polynomial evaluated analytically |
| R-8 | `run_vkb`: support cut at irc; window `20(r/rr(irc) − 1)`, ff = (1−x²)², C¹ | a, d | W | true rc; C-infinity bump window; quadratures split at 0.95 rc, rc(5), rc |
| R-9 | `vpinteg` to `irc`; points 1–3 skipped | c | W | piecewise `lm_integral` |
| R-10 | `vkboutwf`, `lschvkbb`: rc found as the last nonzero vkb; overlaps to `mch` (rule changes with the turning point) | a, e | W | rc passed in; one fixed [0, rc] rule shared with run_vkb |
| R-11 | `vrel` window 0.5–1.5 rr(irc); 1e-10 threshold on u | a, g | W | evaluate only on [0.95 rc, rc] from the ODE; no threshold (u ≠ 0 there unless rc sits on a node, which pspot already rejects) |
| R-12 | `wellstate`: well onset at rr(irc+5), x³ onset (C²); thresholds 0.75, 1.5 rr(irc) | a, d, g | W | onset at rc·amesh⁵ (continuous), C-infinity onset; thresholds: D-1 |
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
| Q-5 | `qroots`: scan dq = 0.02, pole filter \|f·f_last\| < 1, bisection to \|f\| < 1e-12 | h, f | W | roots of the pole-free g(q) = ulgd·j_l(q rc) − q j_l'(q rc), Newton-polished; see D-3 |
| Q-6 | `optimize`: bisection to a 1e-10 bracket | f | W | Newton polish of the secular equation |
| Q-7 | `optimize`: hard-case sign of x₁ | h | W | intrinsic (two global minima); keep, warn |
| Q-8 | `const_basis`: absolute 1e-8 singular-value threshold, values below left unscaled | g | W | relative threshold; error if rank-deficient |
| Q-9 | `wellstate`: rwell by ×/÷0.8 then bisection to \|E − ep\| < 1e-8; ep := E | f, h | W | solve E(rwell) = ep by Brent to round-off |
| Q-10 | `sbf8`: l = 2, 3 closed forms lose eps/x^l at small x | – | W | series below x ≈ 0.5 for l ≥ 2 |

#### Solvers (phase 2)

| id | where | cat | W | decision |
|---|---|---|---|---|
| S-1 | sr: V_l inside rc reproduces V+vr, outside V: V_l jumps by vr(rc) at rc | d | W | the vr window (R-8) is what makes χ vanish at rc; also make the written V_l continuous for lloc < 4 (smooth switch of vr), phase 3 |
| S-2 | `lschvkbb`, `lschpsbar`: eps 1e-8; result is the guess before the last correction | f | W (rho, unscreening) | apply the correction and one more step; eps 1e-13 relative |
| S-3 | `lschvkbbe`, `ldiracfb`, `ldiracfs`: eps 1e-10, same | f | W | same |
| S-4 | `lschfb`: eps 1e-12, same | f | W | same |
| S-5 | `psatom` 204: SCF converged when the eigenvalue guess comes back unchanged | f | – | explicit criterion on eigenvalues and potential residual |
| S-6 | `ionic_screening`: fixed point to 1e-10, 50 iterations | f | W | converge to round-off (Anderson), then a fixed polish |
| S-7 | `lschvkbb` nin = mch + 2.3/al, truncation of the density tail | a, b | W (~0) | `inward_start` budget as in lschfb |
| S-8 | matching point `mch` at the turning point index | h | W | below 1e-13 at amesh 1.006 (measured, emulation); keep |
| S-9 | `sr_flux_m`: series/RK4 switch at rs, RK4 step count by ceiling | g, i | W (~1e-13) | RK4 from r0 always; fixed step count |

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
| O-1 | `linout`, `upfout`: `dpnint` (switching 8-point Lagrange) of kinked arrays onto the linear mesh; ringing beyond rc | d, e | W | each output point from its piece by FH on that piece's extension; β ≡ 0 beyond rc |
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
branch in between.

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
- [ ] Test comparisons through pyoncvpsp and upf-tools instead of
      fldiff (branch `pyoncvpsp-tests`, in progress); `smooth_scan.py`
      then reads outputs through them too

#### Phase 1: primitives

- [ ] `logmesh_m`: `lm_interp`, `lm_derivs`, `lm_cumint`,
      `lm_integral`, `lm_root` (FH d = 6, truncated support)
- [ ] `sbf_overlap` (Lommel, with the near-degenerate series)
- [ ] ODE-recursion derivatives of (u, q) for irel 0, 1, 2; measure
      against the stencil route on hydrogen and 14_Si (PBE)
- [ ] T1 green

#### Phase 2: solvers to round-off

- [ ] T5 written
- [ ] S-2…S-4: final correction applied plus one step, tolerances
      1e-13, in lschvkbb, lschpsbar, lschvkbbe, lschfb, ldiracfb,
      ldiracfs
- [ ] S-5 psatom criterion; S-6 ionic_screening; S-7; S-9
- [ ] Q-6 optimize, Q-5 qroots polish, Q-9 wellstate width by Brent
- [ ] T2 rows `dvloc0`, `debl`, `ep` (well) at G2
- [ ] References patched

#### Phase 3: continuous radii, qcut, projectors

- [ ] T6 written; T3 and T4 baselines recorded
- [ ] R-1, R-2: rc continuous; `irc` the last point ≤ rc
- [ ] R-4 matching data at rc from the ODE; R-5 overlaps; R-6 pieces
- [ ] Q-1…Q-4 residual kinetic energy at exact qcut, fixed quadratures
- [ ] R-7 vploc at continuous rc(5); R-8…R-10 run_vkb, vpinteg,
      vkboutwf on a shared piecewise rule; R-11 vrel; S-1
- [ ] R-12 wellstate onset; R-3, R-13 ionic window
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
- [ ] M-1…M-14
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

### Decisions for review

These change behaviour beyond numerics, so they need a decision before
their phase.

- **D-1 wellstate localisation thresholds** (R-12). Today a bound state
  whose outermost peak is inside 0.75 rc replaces the input `ep`, and one
  inside 1.5 rc forces `ep = 0.05`. Both thresholds are on a continuous
  ratio. Proposal: decide by discrete facts only (a bound state of that n
  exists, and the sign of the input `ep`), and keep the input `ep`
  otherwise. Alternative: keep them as they are, as documented discrete
  events.
- **D-2 fr psp8 truncation** (O-2). Today scalar-relativistic and
  spin-orbit terms below 1e-3 of the first are dropped. Proposal: write
  them all. The files then hold up to 2·nproj terms per l where they hold
  fewer now.
- **D-3 Basis wave vectors** (Q-5). They come from the all-electron log
  derivative at rc; when it crosses l/rc the lowest root passes through
  0 and the basis changes discretely. Proposal: keep the method, find
  the roots robustly, and treat the crossing as a discrete event.
  Alternative: q_n = z_{l,n}/rc, independent of the all-electron
  function, which is exactly smooth but changes every pseudopotential.
- **D-4 The vr window** (R-8). It makes χ vanish at rc in sr and fr. A
  C-infinity bump instead of (1 − x²)² changes the projectors slightly.
  Proposal: the bump.
- **D-5 icmod 6/7 matching order.** It is C² at r_s, so the 3rd and 4th
  derivative columns of psp8 jump at r_s. Proposal: match to C⁴ (5
  constraints of the 9 coefficients); Martins cannot, so it would need
  extending or would keep C².

### Bugs found on the way

To fix in the phase that touches the code, or in phase 7.

- `vpinteg` skips mesh points 1–3 (harmless, about 1e-12).
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
  left unscaled (latent).
- `lschpsbar` 204 uses `cf(mbar)` for `cf(ii)`.
- `wellstate` 157 tests `ierr<0`, which lschfb never returns, so box
  states (`ierr=10`) can enter the bracketing.
- `psatom` can declare convergence during the occupation ramp.
- The solvers' doc comments call `up` d(uu)/dr; it is d(uu)/dn = al·r·du/dr.
- The `oncvpsp_main` comment on the rc rounding says "nearest" and "5
  decimals"; it is "next" and 6.
