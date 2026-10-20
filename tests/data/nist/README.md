# NIST atomic reference data

`LDA.dat`, `ScRLDA.dat` and `RLDA.dat` hold the NIST atomic reference data
for electronic-structure calculations (S. Kotochigova, Z. H. Levine, E. L.
Shirley, M. D. Stiles and C. W. Clark,
<https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations>):
the total energy, its terms and the eigenvalues (Ha) of the neutral atoms H to U
and their first cations, in the ground configurations NIST lists. Each block
starts with a line `Z symbol neutral|cation configuration`, followed by the
lines of the NIST file of that atom as NIST gives them. The data are a work of
the U.S. government.

| set | equation | exchange-correlation | ONCVPSP |
|---|---|---|---|
| LDA | Schrödinger | Slater, VWN | `oncvpspnr_ae.x`, iexc 5 |
| ScRLDA | scalar-relativistic (Koelling-Harmon) | MacDonald-Vosko relativistic Slater, VWN | `oncvpsp_ae.x`, iexc 6 |
| RLDA | Dirac | MacDonald-Vosko relativistic Slater, VWN | `oncvpspr_ae.x`, iexc 6 |

The relativistic sets use the speed of light 137.0359895 (CODATA 1986), in the
equation and in the exchange; ONCVPSP's default is 137.036, so the inputs set it
(`c_light`, the optional last value of the first input line).

NIST states the accuracy of the data as "numerical convergence to the target
accuracy of 1 microHartree in total energy" (the Procedure page; the codes
that produced it differ by at most 0.5 μHa), and "all individual orbital
eigenvalues agreed to within 2 microHartrees" (Reliability). It states none
for the energy terms. `scripts/nist_test.py` runs ONCVPSP's all-electron atom
on the data with exactly these tolerances, 1e-6 Ha for the total energy and
2e-6 Ha for every eigenvalue, and prints the largest difference of the terms
without testing it; ctest runs it on a selection (label `nist`).

The files were made from the copy of the data in RadialAtomicSolvers.jl
(`test/data/dftdata`), which agrees with NIST's pages (checked for U), and
keep NIST's values and labels as they are. The ScRLDA file of neutral Pm
(Z = 61) labels its 3d eigenvalue `32`; `scripts/nist_test.py` reads it as
3d (`LABEL_ERRATA`).
