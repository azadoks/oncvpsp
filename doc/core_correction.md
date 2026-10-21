## MODEL CORE CHARGE

Release 3.2.3 changed modcore2.f90 and modcore3.f90.
See the release notes for details.

Releases prior to 3.2.1 had the option of including a model core charge
to give a "non-linear core correction" in the spirit of Louie, Froyen, and
Cohen.[1]
The functional form was actually a polynomial model allowing
greater continuity at or near the core-valence charge "crossover" following,
more or less, Fuchs and Scheffler.[2]

Model core charges often produce convergence problems, especially  with
phonon calculations and violations of the acoustic sum rule at Gamma.  An
alternative approach to the LFC non-linearity argument was introduced by
Teter.[3]  He proposed that the 2nd derivatives of the exchange-correlation
energy with respect to the occupations of the valence states, known in the
chemical literature as the "hardness," was an important metric.  He
stated that a model core charge should be adjusted to optimize the rms
agreement of these 2nd derivatives calculated for pseudopotential atoms
and all-electron atoms.  In addition, he proposed a functional form for
such a model core charge designed to have good convergence properties
in Fourier space.  This function was used in some pseudopotentials
packaged with the ABINIT code [4}, and some tabulated on that web site.
It  does not appear that the "Teter metric" was actually used in the
generation of these potentials [5].  Routines gg1cc.f90, gp1cc.f90, and
gpp1cc.f90 to calculate  this function and its first two derivatives have
been adapted from ABINIT.

Since ONCVPSP can include shallow core states as valence while
maintaining acceptable connvergence, the model-core corrections are
often unnecessary.  In addition, if inner cores are modeled, the shallow-
core charge at the crossover may be significantly different from that
of the all-electron atom, so the 'non-linearity" correction may not be
the only issue.  Teter 2nd-derivative calculations of Exc are now
included, with 5 options specified by the input variable icmod.

To my knowledge there is no systematic evidence that "Teter hardness
metric" accuracy correlates with improved agreement between all-electron
and pseudopotential results. Its use to tune oncvpsp potentials should be
regarded as experimental.

Model core options:

icmod==0:

    No core correction as usual.

icmod==1:

    The previous polynomial continuation of the all-electron core density
    to the origin, performed at a radius at which the core charge is
    the fcfact times the valence pseudocharge.  fcfact is typically chosen
    to be less than 1, eg. 0.25-0.5.  Exc 2nd derivatives are calculated
    with and without this charge, and the rms discripancies reported in the
    *.out file, but play no active role.  This option is executed by
    modcore.f90.

icmod==2:

    Amplitude and scale of the Teter function are adjusted to fit the
    value and slope of the all-electron rhoc at the radius determined
    as by fcfact as above.  The Teter function is the model at every
    radius; it is not blended into the all-electron rhoc (the
    fcrossover.f90 that this note once described is not in the code).
    It is tapered to 0 by a C-infinity switch between dimensionless
    arguments 3 and 4, so that it and all its derivatives are
    continuous (doc/smoothness_plan.md, M-5).  This is in the spirit of
    LFC [1] but significantly improves convergence.
    The "before and after" Exc 2nd derivatives are computed but not
    used.  Amplitude and scale parameters consistent with options 3 and 4
    below are given in the out file.  This option is executed by modcore2.f90

icmod==3:

    The amlitude and scale parameters for the Teter function are taken
    from the input file, with the fcfact input re-interpreted and the
    rcfact argument added as follows:

    ```
        # MODEL CORE CHARGE
        # icmod, fcfact, rcfact
            3    5.0  1.3
    ```

    These parameters are defined consisently with option 4 below.  They
    are also given the the out  files of options 2 and 4 as " amplitude
    prefactor, scale prefactor."  Note that "fcfact" in this context
    has a DIFFERENT MEANING than in options 1 and 2.  The Teter function
    is tapered as in option 2 above.  Exc 2nd derivatives
    are computed but not used.  This option is executed by modcore3.f90.

icmod==4:

    The Teter function parameters are optimized to reduce Exc
    2nd-derivative discrepancies in modcore3.f90.  The optimization is
    initiated by finding the radius rmatch at which rhoc=pseudo-rho, and
    this charge value is stored as rhocmatch.  The core model is then

         rhomod(rr)=(fcfact*rhocmatch)*F_Teter(rr/(rcfact*rmatch))

    tapered as in option 2.  Exc 2nd derivatives are then computed on a
    coarse grid of plausible fcfact and rcfact values, the rms error
    is printed as a matrix in the out file.  The minimum on this grid
    is then refined using the Nelder-Mead simplex algorithm, and Newton's
    method finds the stationary point of that basin to round-off
    (doc/smoothness_plan.md, M-6).  This option
    is now used by tests/data/40_Zr.dat, and an annotated version of the
    relevant section of tests/refs/40_Zr.out follows this text.  A 1-2
    order-of-magntude reduction of the rms error is often obtained.
    However, the minimization may converge to yield a wildly implausalbe
    model, or may not converge.  An local minimum on the initial grid
    may be selected and used as input in option 3 above.

    fcfact is ignored (0.0 is fine). The psp8 file gives fchrg = 1, which
    only tells Abinit that a model core follows; the echoes repeat the
    input fcfact.

icmod==6:

    The all-electron core charge wherever the plane-wave cutoff affords
    it. Inside the radius r_s the model is a fill that matches rhoc and its
    first two derivatives at r_s; beyond r_s it is rhoc. fcfact is the
    plane-wave cutoff ecut (Ha) of the calculations the pseudopotential is
    for, rcfact (optional, default 1e-5) an exchange-correlation energy
    error (Ha), and slack (optional, 0 < slack <= 1, default 1) the share of
    ecut the model must fit:

    ```
        # MODEL CORE CHARGE
        # icmod, fcfact (ecut, Ha), rcfact (xc energy error, Ha), slack
            6   25.0   1.0e-5   1.0
    ```

    The fill is the minimum tail, the even polynomial

        sum_{k=0..8} c_k (r/r_s)^(2k)

    whose six free coefficients minimize the Fourier tail of the whole
    model, int_{qlo}^{qhi} q^2 |rhomod(q)|^2 dq, with qlo = 2*sqrt(2*ecut),
    the density cutoff of ecut, and qhi = min(3*qlo, 120 bohr^-1). Where
    that polynomial rises somewhere on [0, r_s] (checked at 201 points),
    the fill is Martins' exp(a + b r^2 + c r^4) instead, if that does not
    rise (b <= 0 and b + 2c r_s^2 <= 0); where neither is monotone, the
    radius has no model. The out file says which fill was used, and gives
    its coefficients.

    r_s is the smallest grid radius, no larger than the smallest rc (or
    the last radius with core charge, if that is smaller), whose model
    fits slack*ecut. A model fits a cutoff e if, for every density cutoff
    q from 2*sqrt(2*e) on (steps of 0.1 bohr^-1 up to 120 bohr^-1, 1800
    Ha), the core charge without its Fourier components beyond q changes
    the exchange-correlation energy by less than rcfact:
    |Exc[n + rhomod_q] - Exc[n + rhomod] - Vxc[n_bg] int (rhomod_q -
    rhomod)| within r <= 10 bohr, with n the pseudo valence charge on a
    uniform background n_bg = 5e-3 bohr^-3 that stands for the interstitial
    electrons of the neighbours. The run stops if even the largest r_s
    does not fit. The out file gives r_s and the cutoff the model needs.

    With slack < 1 the model is still made for the qlo of ecut, but must
    fit slack*ecut: a core made for its own, lower cutoff converges worse
    in plane waves. This is the BudgetedCore of PseudoPotentialGenerators.jl
    (PPG), executed by modcore_budget.f90; tests/data/14_Si_bc.dat uses it.
    On the same core and valence charges it gives PPG's r_s within the grid
    spacing, and the same fill and core cutoff (Si, Li, Na, K, Zn, Ga, N,
    I, Pb, Cs; TM, PBE).

    With fcfact <= 0, ecut is instead the cutoff the pseudo-waves need for
    the same error rcfact per electron, times -fcfact (times 1 for
    fcfact = 0): a calculation runs at that cutoff anyway, so the core
    costs nothing more.
    The psp8 file then gives fchrg = 1 (fchrg > 0 tells Abinit that a model
    core follows; otherwise it is fcfact).

    ```
        # MODEL CORE CHARGE
        # icmod, fcfact (0, or minus a multiple of the pseudo-wave cutoff), rcfact
            6   0   1.0e-5
    ```

    The pseudo-wave cutoff is the largest over the optimized pseudo wave
    functions (every l and projector, and both J fully relativistic) of
    the cutoff at which the energy error per electron falls to rcfact:
    the largest Ecut of the "Energy error per electron / Cutoff" tables of
    the out file at that error, interpolated as they are, so the cutoff
    one would run at. The out file gives it on the line "pseudo-wave
    cutoff, largest Ecut of the tables at ...", and modcore_budget_table
    (in modcore_budget.f90) computes it; a wave whose table does not reach
    rcfact does not count. tests/data/14_Si_bca.dat uses it.

    PPG's default budget (ecut_for_error) is instead the root of the total
    residual kinetic energy sum_i w_i E_i(ecut) = rcfact, with w_i the
    electrons that fill wave i (2(2l+1) for an empty or scattering one),
    an error per atom rather than per electron. On ONCVPSP's inputs it is
    up to 2.2 times the cutoff they are run at (ecut_recommended + 5 Ha in
    the QE equation-of-state runs), and the cores made for it at slack 1
    were truncated by ecutrho = 4 ecutwfc; it was the meaning of
    fcfact <= 0 here until the table maximum replaced it.

icmod==7:

    The core of icmod=6 (the same fill, inputs and cutoff test) at the
    radius whose model acts on the valence most as the all-electron core
    does: the HardnessMatchedCore of PPG, and the one to use in plane
    waves. The radii run from the largest r_s of icmod=6 down by factors of
    1.025 to the r_s of icmod=6 over 1.25, that r_s included (each on the
    nearest grid point). For each, the measure is the rms over i, j of

        d2Exc/dn_i dn_j [rho_v + rhomod] - d2Exc/dn_i dn_j [rho_v + rhoc]

    (der2exc, as icmod=3 and 4, but with the pseudo valence charge rho_v
    in both), for the probe densities: the eigenstates of the screened
    pseudo-Hamiltonian with iprj-1 nodes for every projector iprj of every
    l, at unit occupation, if bound. The models are tested against
    slack*ecut in the order of this measure, and the first that fits is
    taken; ties go to the larger radius.

    ```
        # MODEL CORE CHARGE
        # icmod, fcfact, rcfact, slack
            7   0   1.0e-5   0.7
    ```

    The out file adds the r_s of icmod=6, the measure of the chosen model,
    and a table of the radii: r_s, the fill (MT for the minimum tail, M for
    Martins, - for none), the measure, and whether it fits (T, F, or - if
    it was not tested), ended by a blank line. tests/data/14_Si_hm.dat
    uses it. The search costs 1.4-3 times that of icmod=6.

1) S. G. Louie, S. Froyen, and M. L. Cohen, Phys. Rev. B 43, 1993 (1991).

2) M. Fuchs and M. Scheffler, Comput. Phys. Commun. 119, 67 (1999).

3) M. Teter, Phys. Rev. B 48, 5031 (1993).

4) http://www.abinit.org

5_ D. C. Allan, private communication.


From tests/refs/40_Zr.out

```
Model core correction analysis
  based on d2Exc/dn_idn_j

# rows and columns represent the 4 valence states (4s, 4p, 4d, 5s)
# small asymmetries occur because one derivative is analytic and the
# second nnumeric
# d2ex** is an energy given in Ha

d2excae - all-electron derivatives

   -3.104084D-02   -2.792440D-02   -1.363085D-02   -2.288830D-03
   -2.792588D-02   -2.670777D-02   -1.542246D-02   -2.895952D-03
   -1.363104D-02   -1.542158D-02   -1.434836D-02   -4.959607D-03
   -2.288860D-03   -2.895925D-03   -4.959767D-03   -3.447872D-03

d2excps - pseudofunction derivatives with no core correction

   -4.100487D-02   -3.381006D-02   -1.500148D-02   -3.267544D-03
   -3.381063D-02   -3.131882D-02   -1.690321D-02   -3.414016D-03
   -1.500227D-02   -1.690403D-02   -1.502511D-02   -5.051477D-03
   -3.267810D-03   -3.414031D-03   -5.051637D-03   -3.541159D-03

# this serves as a reference to evaluate the improvement that can be
# obtained with the model core charge
# convergence of the Nelder-Mead refinement is based on 10E-4 times this
# value

rms 2nd derivative error    3.543488E-03

# these set the range and amplitude overall scales

rmatch, rhocmatch    0.8317    8.0367

Coarse scan for minimum error
  matrix elements: rms 2nd-derivative errors (mHa)
  column index : amplitude prefactor to rhocmatch
  row index : scale prefactor to rmatch

         1.500  2.000  2.500  3.000  3.500  4.000  4.500  5.000  5.500  6.000

  1.0    2.674  2.497  2.344  2.209  2.087  1.977  1.877  1.784  1.697  1.617
  1.1    2.430  2.201  2.002  1.827  1.670  1.528  1.398  1.279  1.169  1.068
  1.2    2.151  1.861  1.610  1.390  1.193  1.017  0.858  0.716  0.591  0.486
  1.3    1.836  1.479  1.171  0.903  0.668  0.469  0.322  0.271  0.335  0.454
  1.4    1.488  1.058  0.692  0.387  0.212  0.325  0.539  0.755  0.961  1.156
  1.5    1.112  0.609  0.221  0.307  0.631  0.940  1.226  1.491  1.736  1.964
  1.6    0.716  0.196  0.450  0.885  1.283  1.644  1.972  2.273  2.551  2.809
  1.7    0.340  0.461  1.017  1.520  1.971  2.374  2.740  3.075  3.382  3.665
  1.8    0.317  0.977  1.613  2.173  2.668  3.113  3.511  3.874  4.207  4.511
  1.9    0.714  1.514  2.217  2.827  3.364  3.840  4.268  4.654  5.006  5.328

Nelder-Mead iteration
 converged in   9 steps
amplitude prefactor, scale prefactor
    2.3184    1.5543

# these two numbers can be used as input to a option-3 calculation as
    # icmod, fcfact, rcfact
        3   2.3184    1.5543
# this will reproduce this model and can be used as a starting point
# for hand-tuning
# an option-2 calculation will produce different prefactors for the  same
# purpose


Optimized Teter model core charge

d2excps - pseudofunction derivatives with core correction

   -3.119948D-02   -2.772889D-02   -1.343294D-02   -2.463641D-03
   -2.772866D-02   -2.694824D-02   -1.561979D-02   -2.963984D-03
   -1.343341D-02   -1.562047D-02   -1.457252D-02   -4.954341D-03
   -2.463659D-03   -2.963951D-03   -4.954498D-03   -3.469214D-03

rms 2nd derivative error    1.654494E-04

# this is a factor of 20 reduction of the error
# For options 1-3, the coarse-grid and Nelder-Mead results are absent
# in the *.out file, but the 2nd-derivative matrices and rms values
# appear as above.
```
