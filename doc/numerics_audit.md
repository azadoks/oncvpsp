# Numerics audit of ONCVPSP, branch `smooooth` (2026-10-22)

The audit covers the working tree at 9c78c3a plus the uncommitted changes: free_tail_m, `eresid`'s new tail argument, `wellprj` in the diagnostics, and `smooth_scan.py`. Another session was editing `free_tail_m.f90` and `doc/smoothness_plan.md` while the audit ran, so line numbers in those files are approximate.

**Status.** C-1 to C-11 are resolved; the status column of §6b gives the commits, and "Numerics audit" in doc/smoothness_plan.md how each was checked and what changed. C-12 to C-18 and §6c are open; §6c item 2 is deferred to the configurable Adams–Moulton order in doc/smoothness_plan.md, since with the AM4 ODE it removes only 5–6 % of the AE mesh error (measured there).

**Note (2026-10-27).** Two items are closed without a code change, by the user's decision:
- §6c item 4 (S-10, GGA ρ′ and ρ″ from the pieces) is dropped: RadialAtomicSolvers.jl showed that the 5-point stencils of `derivs` across the model-core joints do not matter. The question of the C² joint stays with D-5.
- R-16 in doc/smoothness_plan.md, the kink of `outwf_m`'s per-segment restart where rcloc crosses 0.95 rc(l) (1.7e-11 Ha in the 3s of 20_Ca, at most 1e-10 in the test configurations), is resolved as within tolerance: it is below G1 and is kept as is.

**Notation**
- Log mesh: r_i = rr1·amesh^(i−1), x = ln r, h = al = ln amesh. Default: amesh 1.006, rr1 = min(5e-4/Z, 5e-5) (nr and sr) or 5e-4/Z (fr), rmax 45, mmax = int(ln(rmax/rr1)/al).
- "Order p" means a global error of O(h^p).
- ✓ means I checked the claim in the code myself. Claims without ✓ come from reading the code during the sub-audits and were not re-checked line by line.

---

## 1. The shared primitives

Everything new on this branch goes through these modules, so their properties set the accuracy of most rows in the tables below.

| primitive | method | order | parameters |
|---|---|---|---|
| `lm_interp` (logmesh_m) ✓ | Windowed Floater–Hormann barycentric rational interpolant in x. The window is a C∞ bump, 1 within 10 steps and 0 beyond 20. Returns node values exactly. | 7 (d = 6) | Reach ±27 nodes. No poles. C∞ in x. |
| `lm_derivs` ✓ | Derivative fields d^k/dx^k (k ≤ 4) on the nodes from 11-point Fornberg stencils, with weights computed in quad precision. The fields are FH-interpolated to the radius, then converted to d/dr with Stirling numbers of the first kind (I checked the matrix). | Centred: 10, 10, 8, 8 for k = 1–4. Shifted one-sided near the ends: 10, 9, 8, **7**. | Round-off ~ eps/h^k. |
| `lm_cumint` / `lm_integral` ✓ | Cumulative integral of f·r dx. Each cell integrates the degree-7 Lagrange interpolant on 8 nodes (centred i−3..i+4, shifted at the ends). Kahan summation. Below rr(1), f ∝ r^m with m given. ∫_a^b is FH(C)(b) − FH(C)(a). | 8 at the nodes, 7 at limits between nodes | The cumulative array is recomputed on every call (cost only). |
| `lm_integral_weights` ✓ | The same linear functional written as weights, for repeated integrals to one radius. | as above | I checked the tail-sum bookkeeping. |
| `lm_root`, `ug_root` ✓ | Brent (zeroin) in ln r or in the grid index, with no tolerance of its own (stops at a few ulps). | superlinear | **The 200-iteration limit returns silently.** |
| ugrid_m ✓ | The same 8-point cell rule and windowed FH on uniform grids: `ug_weights`, `ug_cumint_down` (accumulated from the small end, Kahan), `ug_cross_down` (Brent on FH of ln f, last crossing). | 8 / 7 | Needs nn ≥ 8. |
| pieces_m ✓ | Piecewise functions (≤ 4 segments). Each piece is held as its smooth extension ≥ 64 nodes past its bounds. Integrals are `lm_integral` per segment; values come from the segment's piece. | as lm_* | No stencil or interpolant spans a breakpoint. |
| match_m ✓ | Continuous turning point (FH root). Bump-weighted energy correction over ±4 points. C∞ blend of the inward and outward solutions. | — | `match_smax` = (a+b+√((a−b)²+1))/2. |
| `cinf_taper` ✓ | f(1−x)/(f(1−x)+f(x)), f = e^(−1/t) | C∞ | Used for the ionic window, the error sphere, the vpuns tail and the ghost weights. |
| qlinalg_m `qjacobi` ✓ | Cyclic Jacobi (Rutishauser rotations) in quad precision | quadratic convergence | 60 sweeps. Stops when off == 0 exactly, which underflow reaches. |
| sbf_m `sbf_trig` ✓ | j_0..j_3 from sin and cos, switching to the power series below x = 0.5, 1.5, 2.5, 3.5 | ~3 eps | Hard switch, invisible to round-off. |
| `sbf8` ✓ | x < 1e-3: two-term series. nm ≤ 4: `sbf_trig`. Otherwise Miller downward recurrence from nm + 1.36x + 15, normalized by Σ(2l+1)j_l² = 1. | ~0.6 eps·x | — |
| sbf_q_m ✓ | Quad precision. Values: series below x = 1/4, else downward recurrence **normalized by j_0 = sin x / x**. Overlaps: Lommel closed form; a fixed 48×48-term double power series when βrc ≤ 1.5; a Taylor series in Δq from the Bessel ODE when Δq·rc ≤ 1e-13·αrc (I re-derived the recurrence). Derivatives at rc: the recursion (1/z d/dz)^m. | exact to ~1e-30 | l ≤ 4 for the derivatives. |
| outmesh_m ✓ | Output values from `lm_interp` per piece. r = 0 by degree-7 Lagrange extrapolation through output points 2..9. | 7 | The extrapolation amplifies noise by up to Σ\|c\| = 255. |
| `derivs` ✓ (legacy) | 5-point centred stencils in the index, one-sided 5/6-point at the ends | **4** | Used for the GGA ρ′ and the divergence form. |
| `dpnint` ✓ (legacy) | Switching 8-point Lagrange interpolation | 8 | **No callers left.** |

---

## 2. Differentiation

| where | what, why | method | order, parameters | notes |
|---|---|---|---|---|
| `ode_rc_der` | u, u′…u⁗ at rc for the matching constraints | Taylor recursion of the first-order system y′ = A(r)y, y^(n+1) = Σ C(n,j) A^(j) y^(n−j). y comes from FH; V′, V″, V‴ from `lm_derivs`. | Accuracy of u itself, ~1e-10. The 4th derivative carries V‴ round-off eps/h³ ~ 1e-9. | Cross-checked against stencils on u: the program stops if they disagree by more than (1e-6 + 100h⁴)·scale. Smooth in rc. |
| `sbf_rc_der`, `sbfq_rc_der` | Derivatives of j_l(q r) at rc (constraint rows) | Analytic recursion | exact | Small-x cancellation of 1/x^k terms (double precision version). |
| `vploc_fit` | V_loc polynomial matched to V, V′, V″, V‴ at rc(5) | `lm_derivs`, nder = 3 | 8 | Closed-form coefficients for lpopt 1–5, checked numerically to 1e-14. |
| `ae_delta` | Δ = u″ − [l(l+1)/r² + 2(V−ε)]u at the projector window | Pointwise from (u, up/(al r)), V′ from `lm_derivs` | 8 | No division by u. |
| `pspot` | AE u″ beyond rc (projector seed, kinetic integrand) | 7-point centred first derivative in the index, /(60 al r), applied twice | 6 | Uses nodes only beyond rc. |
| `pspot` inside rc | (ε − T)φ and V_l = ε − Tφ/φ | Analytic Bessel sums via `sbf_rc_der` | exact, but… | For l ≥ 1, work4 = ½(l(l+1)u/r² − u″) cancels two O(r^(l−1)) terms. The absolute error is ~eps·l(l+1)/r², about 1e-7 Ha at rr(1). |
| `modcore_budget` set_model | ρ_c, ρ_c′, ρ_c″ at r_s (C² match of the fill) | `lm_derivs` | 8 | |
| modcore_m, modcore_budget | Psp8 model-core columns 2–5 | Analytic: Teter jets (Leibniz rule, C∞-switch jets), Legendre-basis fill derivatives, Faà di Bruno for t = (r/r_s)², Martins recursion | exact | Beyond the joint, ρ_c columns come from `lm_derivs` out to mmax, against its "away from the ends" precondition (harmless, ρ_c ≈ 0 there). |
| `teter` sinc | h, h′…h⁗ | 16-term series for \|y\| ≤ 2, closed form beyond | ~eps | Hard switch. |
| `modcore_budget` slope_at | Slope of the icmod 7 mismatch in ln r_s | Central difference, step 1e-2 | 2 | h² bias on the minimiser. It differentiates `der2exc`, which reads a C² mesh-sampled core. |
| `der2exc` ✓ | dV_xc/dn_k for the icmod 1–4/7 mismatch | 4-point central difference in occupation, (1,−8,8,−1)/12h, **h = 0.15 e** ✓ | 4 | The perturbed ρ − 2h ρ_k can go negative where state k is thinly occupied, and the xc routines then clip. |
| `derivs` → excggc, exc_libxc, functionals | GGA ρ′ (and ρ″, ∇²ρ); divergence of ∂f/∂∇ρ | 5-point stencils in the index | **4** | Runs across the model-core joint (C² at r_s, C⁴ at r*) and across rc (S-10, open). **Bug ✓: the end formulas for dppn at mmax−1 and mmax negate the forward coefficients, which gives −f″** (`derivs.f90:84,89`). An even derivative keeps the sign of the forward stencil. Only ρ″ and ∇²ρ are affected, they are passed only to meta-GGA, which is unreachable, and ρ ≈ 0 there, so nothing written changes today. |
| run_diag, phsft_m | Log derivatives | Solver's own up = al·r·u′, divided by u, interpolated (phsft_m) or at the snapped mchf (run_diag) | solver order (4) | No differencing of interpolants. |
| `free_tail_m` | t′ = u′w + u w′ | Analytic erfc derivative | exact | |

---

## 3. Integration

### 3a. Radial quadrature on the log mesh

| where | integral | rule | order | ends |
|---|---|---|---|---|
| `fpovlp` | AE overlaps/norm inside rc | `lm_integral` | 8 (7 off-node) | Power law below rr1. Exactly to rc. |
| run_vkb, sr_so_r, run_ghosts | B matrix, projector norms, overlaps | `pw_integral` / `pw_product_integral` per segment | 8 | Split at 0.95 rc, rc(5), rc. |
| `outwf_m` | Overlaps ⟨β\|u⟩ in the nonlocal solvers | Segment 1: `lm_integral_weights` (4-lane Kahan or dgemm). Later segments: `ug_cumint_down` on the segment grid. | 8 | Per segment. |
| `run_optimize` (free-tail states) | ⟨ψ\|T\|ψ⟩ inside rc | `lm_integral` + `tail_kinetic` (`ug_weights`) | 8 | |
| `pspot` ske (other states) | Kinetic energy | Midpoint in x from rr(5)/√a to mmax−9, plus a power law below | **2**, across the kink at rc | Diagnostic only. |
| `der2exc` | d²E_xc matrix | `lm_integral` to rcut | 8 | |
| `modcore_budget` xc_error | xc error of the model core | `lm_integral`, split at r_s; inner part on its own log mesh ending at r_s; cinf_taper 6→10 bohr | 8 | |
| lschfb, lschfs, ldiracfb/fs, lschvkbb/be/bs, lschpsbar, renorm_r | Normalization ∫u² | Weight 1 in x, Gregory end (23,28,9)/24, head ro^(2γ+1)/(2γ+1) at ro = r1/√a | 4 at the right end, 2 at the head (negligible) | Head coefficient implicit (1), except renorm_r, which scales it by u1²/r1^(2γ). lschvkbbe and lschvkbs normalize only to mch. |
| `vout` | Hartree r·V_H | Inward cumulative Adams–Moulton 4 used as quadrature | **4** | rv = Q at the last 3 points (2 cells dropped). Q is imposed, not ∫ρ. |
| vout, `ae_energies` (radint) ✓ | E_ee, E_xc, E_nuc, E_scr | Gregory start (9,28,23)/24, weight 1 to mmax, head f(1)/pw | 4 at the start, open end | **Head exponent ✓ is non-relativistic (pw = 2 for E_nuc, 3 for the others).** For sr/fr s states the density goes as r^(2γ−2), γ = √(1−(Z/c)²), so the exact head is f(1)/(2γ) or f(1)/(2γ+1). The sub-audit estimates ~0.16 Ha moved between E_kin and E_nuc for U (it cancels in etot) and ~1e-6 Ha left in etot through the f/3 terms. Not run. |
| `ae_residual`, Anderson metrics | SCF residual, mixing inner products | Rectangle sums | 1–2 | Measures only. The metric differs between ae_residual (ρr³) and Anderson (r²). |
| `inward_start` | ∫κ dr to pick nin | Trapezoid in r | 2 | Only selects an index. |
| `coulomb_tail` | α_Z, its change at the cuts | Rectangle in x; cuts at 6 and 10 bohr **snapped to the mesh** | 2 | Diagnostic (X-6, open). |
| `upfout` uurcut | Tail norm → UPF nrl | Rectangle, inward | 1 | Picks a point count only (O-4, kept). |

### 3b. Bessel, Fourier and q-space integrals

| where | integral | rule | parameters |
|---|---|---|---|
| sbf_basis_con, eresid (inside rc) | ∫₀^rc r² j_l j_l, and the inside transform | Lommel closed form in quad precision; series where it cancels | exact |
| `eresid` (outside rc) | ∫_rc u(r) r j_l(qr) dr | 8th-order `ug_weights` on the linear mesh rc + (i−1)hh, hh = 2π/(75 qmax) | sin/cos rotated from one q to the next and re-seeded every 32 steps (~32 eps). The tail is cut at nr (u at machine precision, or the window's end). |
| `eresid` | E_r(q) = (1/π)∫_q^qmax q⁴\|F\|² | `ug_cumint_down` (8th order, from qmax down), plus an 8-point Gauss–Legendre partial cell to qcut exactly | dq = qmax/(nq−1), nq = 1 + **nint**(150 qmax rc/2π). qmax = smax(smax(2q_nbas, 20), 2qcut). Truncation at qmax is a smooth bias (Q-3). |
| `modcore_budget` | Forward and inverse core transforms, window ∫q²\|n(q)\|² | Private Gauss–Legendre (Newton on P_n, tol 1e-15, 100 iterations, silent at the limit). 400 nodes for the fill on [0, r_s]; 8 nodes per mesh interval (in x) beyond r_s; 80×8 window panels; 5 nodes per 0.1 of q for the inverse transform. | Fixed counts (M-11). The guard qmax·r_top ≤ 700 keeps the 400-node rule exact. |

---

## 4. Interpolation and extrapolation

| where | what | method | notes |
|---|---|---|---|
| everywhere on the branch | Values at arbitrary r | `lm_interp`, or `pw_value` per piece | Derivatives are never taken of an interpolant (by design). |
| `ode_rc_der`, phsft_m | u and up at rc, rxpsh | FH | |
| `mesh_tail` / `eresid` | u → linear tail mesh | FH of u to its last nonzero node | |
| `free_tail_m` | V at the RK4 stages | FH inside the mesh; **V(rmax)·rmax/r beyond** (Coulomb extrapolation) | Valid for the ionic V_AE. |
| `well_mesh_m` | V_AE on the log mesh continued past rmax, for the well states | **V(rmax)·rmax/r beyond**, as `free_tail_m` | The mesh ends where the state's WKB decay beyond the well's turning point reaches machine precision (C-19). |
| `outwf_m` | V and sources onto the segment grids, and the solution back to the nodes | FH with offset | Restart per segment by interpolation. |
| run_ghosts | V_loc, β onto the scaled mesh s·rr (wall at 3rc exactly) | `pw_value`, `pw_resample` | See concern C-7. |
| outmesh_m (psp8, UPF) | All written arrays on the linear mesh | FH per piece; 0 beyond the last bound | `om_values` interpolates ρ_v (C^(ncon−1) at rc) as one smooth array. That is the only remaining interpolant across a kink in the output (small). |
| outmesh `origin` | r = 0 | Degree-7 Lagrange extrapolation through output points 2..9 | V_loc is overridden by vpuns(1) in psp8 and in sr linear UPF, but **fr UPF keeps the polynomial**. PSML uses linear extrapolation through nodes 2 and 3 (`add_zero_r`). Three conventions for one quantity. |
| upfout `upf_mesh` | Log-strided UPF mesh | Nodes exactly, rab = istr·al·r | No interpolation. |
| ug_cross_down | Ecut at error levels | Brent on FH of ln E in the q index | Last crossing (a discrete event at shoulders). |
| vploc_m | V_loc inside rc(5) | Analytic polynomial | |
| `tfapot` | Thomas–Fermi start potential | Latter fit | Start guess only. |
| `dpnint` | — | 8-point switching Lagrange | **Dead.** |

---

## 5. Other algorithms

### ODE integration
| where | method | order | notes |
|---|---|---|---|
| `sr_flux_m` (lschfb, lschfs) | Koelling–Harmon flux form (no dV/dr). Frobenius/KH series at the origin (I re-derived the recurrences); RK4 in x (h ≤ 2e-3) from the series radius a/(8 max(1,\|b\|)) when rr1 is beyond it. AB4/AM4 PECECE outward, then re-integrated inward over points 1..11. | 4 | The series has silent 80/200-term limits. |
| lschfb inward | AB4/AM4 PECECE from `inward_start` (WKB decay to eps), exponential start | 4 | |
| ldiracfb, ldiracfs ✓ | AB4/AM4, two corrector passes, **no derivative re-evaluation after the last (PECEC)** ✓. **Start: leading term g = r^γ, f = cof·g at points 1..5 only** ✓ | 4 | Unlike the sr path: no series, no screening constant, no inward re-integration. The sr comment records a ~1e-5 Ha contact error from exactly this shortcut (Na). |
| lschvkbb, lschpsbar, lschvkbbe | Adams PC (adams.inc) on the nodes | 4 | lschpsbar still matches at the discrete mch, and its stencils cross the V_loc kinks at 0.95 rc_loc and rc_loc (R-14, X-3 caveat). |
| outwf_m | Per-segment uniform-x grids, RK4 priming (3 steps × 4 substeps), AB4 + 2 AM4 corrections, restart per segment | 4 | Measured ratio 15.4 per halving. |
| free_tail_m | RK4 on the flux/Dirac system: hh/2 on the linear mesh; ≤ 2e-3 bohr per step on the log mesh, `ns = ceiling(...)` | 4 | Two meshes for the same tail (eresid's and the stored uua) agree only to truncation (~1e-8 estimated). |

### Roots, eigenvalues, optimisation
| where | method | tolerance / limit |
|---|---|---|
| lschfb, ldiracfb, lschvkbb/be, lschpsbar | Node-count bisection plus perturbative (Newton) correction | 1e-12·max(\|E\|, 0.2), then polished to 8 ulp or until the step stops halving; 60 passes, then ierr. |
| boxstate | Doubling bracket plus bisection on the node count | 1e-12; the 60/200-step limits are silent. |
| qroots | Scan with step 0.02, pole filter \|f·f_last\| < 1, bisection to 1e-12, Newton polish on g = ulgd·j − q j′ | A failed bisection gives a **warning and the root is skipped**. |
| optimize | dsyev of E_r; secular equation by bisection (1e-10) + Newton (≤ 6 steps, not clipped to the bracket) | Hard-case sign branch (Q-7, intrinsic). |
| wellstate | rwell ×/÷ 0.8, then Illinois regula falsi | 4 eps; a closed bracket with a residual above 1e-10 is an error ✓ (plan). |
| const_basis | Unit-norm rows, quad Jacobi on AᵀA (squares the condition number; quad absorbs it); rank test (1e-10)² | |
| run_vkb | B symmetrized, dsyev; S by canonical V Λ^(±½); sort by \|evkb\| | |
| sr_so_r | Same, on 2·nproj | **No guard on S eigenvalues** (√λ, 1/√λ). |
| run_ghosts | dsyev of diag(e) + Σ evkb (w g)(w g)ᵀ, C∞ weights in energy | |
| modcore_m | lm_root crossover; closed-form a0; cubic by companion matrix + 2 Newton steps; Teter Newton-in-bracket | |
| modcore_budget | Bisection on the mesh then Brent (r_s); Brent `fmin` (1e-8, 100 iterations, silent); 0.1 q grid + fmin for max err; companion-matrix monotonicity test; dgglse constrained LSQ in a Legendre basis | The zeroin/fmin limits are silent. |
| modcore3 | Grid argmin → Nelder–Mead (100 steps, warns) → Newton on central-difference Hessian eigen-directions (30 iterations, warns) | |
| phsft_m | Continuation of the doubled angle, halving the step while \|Δ\| > π/2, depth 40 (accepted silently at the limit) | |

### SCF and fixed points
| where | method | criterion |
|---|---|---|
| sratom (`ae_mix`) | Anderson over 8 iterations, dgesv with 1e-10 regularization; Kerker-type filter in x for irel ≠ 0 (Thomas solve); 100 iterations | \|ΔE\| ≤ 1e-11 rel and res ≤ 1e-11, **or res ≤ 1e-9 after a 5-iteration stall (warning only)**. |
| relatom | **One-step Anderson, no filter, no sd = 0 guard** | Same criterion. |
| psatom | Linear 0.5, then one-step Anderson (no sd guard); occupations ramped over 50 iterations | res, ΔE ≤ 1e-12 and `all(fa == fat2)`. |
| ionic_screening | Picard iteration | 1e-13 Ha; 50 iterations, then error stop. |

---

## 6. Assessment

### 6a. Consistency

**Good.**
- The new path is one coherent system. It works in x = ln r throughout and uses one interpolant (FH, d = 6, order 7), one family of derivative stencils (11-point, ≥ 8th order) and one quadrature (the 8-point cell rule, order 8) on both the log mesh and the uniform grids.
- Breakpoints are handled by pieces, never crossed.
- Closed forms (Lommel, Teter jets, Legendre fills, the vploc polynomial) are used wherever the integrand is analytic, and quad precision where the basis is ill-conditioned.
- This is what doc/smoothness_plan.md specifies, and the code matches it closely.

**The inconsistency is between that system and the legacy core it sits on.** The all-electron/SCF layer and several solvers still run 4th-order numerics with their own conventions:

| quantity | new path | legacy path |
|---|---|---|
| ∫ on the log mesh | `lm_integral`, order 8, exact limits | Gregory/(23,28,9) sums, order 4, two different end conventions (normalization: midpoint head + Gregory tail; energies: Gregory head + open tail); vout's Hartree by AM4 |
| d/dr | `lm_derivs`, order 8 | `derivs`, order 4 (GGA); 7-point twice in pspot (order 6) |
| head below rr1 | power law with a given exponent | normalization uses γ; the energies use non-relativistic 2 and 3 |
| `al` | `mesh_step` from the ends | `0.01·ln(rr(101)/rr(1))` in ~6 places (needs mmax ≥ 101) |
| ODE | sr: series + RK4 + PECECE + inward re-integration | fr: leading term only, PECEC |
| SCF mixing | ae_mix (Anderson-8 + Kerker) | relatom and psatom: one-step Anderson |
| r = 0 output | — | three conventions (psp8/sr-UPF constant, fr-UPF degree-7, PSML linear) |

None of these breaks smoothness in the scanned inputs: they sit in the all-electron atom, which does not depend on rc, or they are C² in the moving radius, as S-10 argues. They do set the accuracy floor, and the order-4 rules dominate the mesh-convergence rate of everything downstream of the AE atom.

### 6b. Correctness: findings, ranked

| # | severity | finding | status |
|---|---|---|---|
| C-1 | **high** | `psatom.f90:272–274,147`: occupations ramp as fa = fat1 + wt·(fat2 − fat1), and convergence needs `fa == fat2` exactly. In IEEE arithmetic 2 + (0.1 − 2) = 0.10000000000000009 ✓, so such a test configuration never converges and is dropped as ierr = 101 ("no output for configuration"). The legacy code assigned fat2 exactly for it > 50. Fix: `if(it>=50) fa=fat(:,2)`. | **fixed in fa0b95e** |
| C-2 | high (policy) | `oncvpsp_main.F90:488`: `ierr` from `aeatom` is never tested, so a non-converged reference SCF only prints a warning. ae_converged also accepts a 1e-9 stall as converged. Both contradict the plan's rule that an iteration ending on its limit is an error. | **fixed in fa0b95e** (the aeatom check; the 1e-9 floor kept, see the plan) |
| C-3 | medium | AE energy heads below rr1 use non-relativistic exponents (`ae_scf.f90:109,115,121`, `vout.f90:108,188,192`). It mostly cancels in etot; E_kin and E_nuc are individually off for heavy atoms. | **fixed in 260d2af** |
| C-4 | medium | The fully-relativistic solvers start from r^γ only, and are PECEC rather than PECECE. The sr path had a ~1e-5 Ha error from the same start. | **fixed in 34befb5** |
| C-5 | ~~medium~~ fixed | `free_tail_width` (free_tail_m.f90:148) had `2*wbound` under the root, so it solved (ΔqL)²/4 − 9κL = 60 instead of the documented 30, making L √2 too long for κ = 0 (+41 % cost). | **fixed in d603e40** (now `wbound*dq**2`; free-tail leakage tests unchanged) |
| C-6 | medium | `modcore_budget`: the minimum-tail → Martins/none switch is a strict `p′ > 0` test, a policy branch on a continuous quantity that the inventory doesn't list. Brent `zeroin` then accepts a "root" at that jump with no residual check, and its 200-iteration limit returns silently (as do fmin, gauss_legendre, Nelder–Mead and the modcore3 Newton). | **fixed in 86f8e58** (the fill switch was already listed) |
| C-7 | medium | run_ghosts builds the barrier basis from `vlocpw` (no ionic-screening switch beyond rcmax), while its reference energies come from lschvkbb, which also uses the vp that `ionic_screening` modified. The GHOST(−) margin compares two different local potentials. | **fixed in 123954b** |
| C-8 | medium | Diagnostics meant for finite differences reselect discretely: `ghost_r`'s argmin, the jj index after ghost skips, the π·anint phase alignment, and the basis truncation when lschpsbar fails. | **fixed in 123954b** (ghost_r; the rest are discrete events) |
| C-9 | medium | `der2exc` uses an occupation step of 0.15 e ✓ in a 4-point stencil. Densities can go negative and get clipped (icmod 1–4 mostly). icmod 7's slope is then a 1e-2 central difference of that. | **checked, kept** (23d54d2; see the plan) |
| C-10 | low–medium | `eresid` dq = qmax/(nq−1) with nq = 1 + nint(...) ✓ jumps by ~1/nq when nint flips (a quadrature-error-sized step, called "continuous" in the comments). The free-tail resolution worry is answered: dq = 2π/(150 rc) ≈ 0.015 bohr⁻¹ at rc = 2.7, and the tail's spectral peak has width ~1/L ≥ 0.3, so it spans more than 20 q nodes. | **fixed in 23d54d2** |
| C-11 | low–medium | pspot/run_vkb V_l and V_loc (lloc < 4) near the origin for l ≥ 1: cancellation gives ~eps·l(l+1)/r², about 1e-7 Ha at rr1. Fix: −u″ + l(l+1)u/r² = Σ c_j q_j² r j_l(q_j r), exact. | **fixed in bbae0de** |
| C-12 | low | `derivs.f90:84,89`: the 2nd-derivative end formulas have the wrong sign. Latent (only meta-GGA reads them). | ✓ verified |
| C-13 | low | `sbfq_values` normalizes by j_0 = sin x / x, which loses relative accuracy as qrc approaches a zero of j_0 (nπ). Normalizing by the larger of j_0, j_1, or by Σ(2l+1)j_l² as sbf8 does, removes it. | ✓ code |
| C-14 | low | sr_so_r has no guard on small S eigenvalues (light elements, j = l ± ½ nearly equal). | sub-audit |
| C-15 | low | qroots skips a root on a bisection failure (warning; done 2026-10-23: the bisection on the pole-free g cannot fail, smoothness_plan Q-5). The Brent/fmin/gauss/phase-depth limits are silent. | sub-audit |
| C-16 | low | Origin values differ by format; fr PSML potentials lack the vpuns tail taper; UPF `nrlproj` is not clamped to nrl, and an odd nrl made even can pass rr(mmax). | sub-audit |
| C-17 | low | Order claims: `lm_derivs`' comment says "8th order or better", but its one-sided end stencils give order 7 for k = 4. `lm_integral` at limits between nodes is order 7, not 8. The plan says the derivative stencils are "centred 8th order". | ✓ |
| C-18 | info | Dead code: `dpnint`, `wf_rc_der` (only the test cross-check), the double-precision `sbf_overlap` (only a unit test), `rpk`, the negative-ierr branches in config_ae_m, assorted unused locals. | ✓ for the first three |
| C-19 | high | The end of the log mesh was a box for the well states. A well state at a small energy (a wide well) reaches out to and past rmax, where the solvers took u = 0: its energy, its tail in E_resid and so the kernel depended on rmax. 14_Si's d well state at 0.02 Ha (half-point radius 75 bohr): kernel 2e-2, its 3e-5 Ha cutoff 95 % between rmax 45 and 60; Pb's f well state at 0.05 Ha 8.7e-4; in the test inputs 74_W's f state 1.7e-7. Everything else moves with rmax only by the SCF's round-off floor (~1e-9, at most 1e-7 in 07_N's log derivatives), measured over every test input at rmax 45, 60, 80. The default well rises only 0.5 Ha above the state, so its WKB decay reaches machine precision about 36 bohr past the wall: most well states of the test inputs met inward_start's fallback, u = 0 at rr(mmax − 4), at e^-30 or so. | **fixed** (`well_mesh_m`: the well states solved on the mesh continued as far as they need, 56 to 70 bohr in the test inputs; their kernels move by at most 3e-8, 74_W's by 1.5e-6; test `well_rmax_*`) |

Documentation drift found along the way:
- program_overview says the SCF converges at 1e-8; the code uses 1e-11.
- program_overview says the series radius is "a quarter"; the code uses a/8.
- The plan's "two passes below 1e-12" is now an 8-ulp polish.
- core_correction.md's icmod 6/7 description (201 points, nearest grid point, T/F column) is out of date.
- The plan mentions `match_order = 4`, which does not exist.

### 6c. Suggested improvements, roughly by payoff

1. **Fix C-1 and C-2 now.** Both are one-liners and both are real failures.
2. **Retire the 4th-order legacy quadrature in the AE layer** (deferred: measured, the ODE limits the rate; see the status note). Replace radint, vout's Hartree (AM4) and the normalizations with `lm_cumint`/`lm_integral`, with the correct head exponent (2γ−2 for relativistic s states) passed in. V_H then follows from two cumulative integrals, which the 8th-order rule gives directly. This removes C-3 and the two Gregory conventions, and should lift the AE mesh-convergence rate from h⁴ toward h⁸ wherever the ODE is not the limit. Measure with T4.
3. **Bring the Dirac solvers to the sr standard:** the series start (the Dirac analogue of `origin_series`, with the screening constant), PECECE, and inward re-integration of the first points (C-4). Then make relatom use `ae_mix`.
4. **S-10 for the GGA** (dropped 2026-10-27, see the note at the top). ρ′ from the pieces and the 8th-order stencils (`lm_derivs`-style fields) instead of `derivs`. At the same time, decide whether the C² joint should stay C² (D-5).
5. **Make iteration limits errors everywhere:** zeroin, fmin, gauss_legendre, Nelder–Mead and Newton in modcore3, boxstate, the series term limits, phsft depth, qroots skip. Add a residual check after Brent, as wellstate does, so that a "root" at a jump is caught (C-6, C-15).
6. **Remove the remaining policy branches on continuous quantities,** or list them in the plan as discrete events with a printed margin (G5): the fill-type switch (C-6), the diagnostic reselections (C-8), and `eresid`'s nint (make dq a continuous function and keep nq fixed per input, as Q-2 intended) (C-10).
7. ~~Resolve C-5~~ (done, d603e40); the q-grid resolution for long free tails (C-10) is adequate.
8. **Smaller items:**
   - Use the exact Bessel identity for V_l near the origin (C-11).
   - Make the sbfq normalization robust (C-13).
   - Guard sr_so_r's S (C-14).
   - Give `der2exc` a smaller step, or analytic f_xc (libxc can supply it) (C-9).
   - Use one r = 0 convention for every format (C-16).
   - Fix the `derivs` end signs (C-12).
   - Replace the hard-coded `al = 0.01·ln(rr(101)/rr(1))` with `mesh_step`.
9. **Housekeeping:** delete `dpnint` and the dead branches; keep `wf_rc_der` and `sbf_overlap` only if the tests still want them as cross-checks; fix the doc drift listed above.
