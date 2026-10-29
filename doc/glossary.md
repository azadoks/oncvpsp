# Glossary of argument names

The dummy arguments that mean the same thing wherever they occur, with
the one `@param` text each has. A routine's `!! @param` line for one of
these names starts with the text below, word for word; it may go on
after a comma or a semicolon to say what the routine does with it (`,
output`, `; zero on return`); `a / b` gives two texts, either of which
may start it (the singular and the plural). `scripts/check_conventions.py` checks it
(rule `glossary`); the table is read from this file, so a change here is
a change of the rule.

Units are bohr for lengths, Ha for energies and potentials and bohr⁻¹
for wave numbers. A radial function is stored times r (`uu` is r R).
A charge density is stored as the sum of f (u/r)² over the states, 4π
times the density, so that the integral of r² rho is the number of
electrons.

## Fixed meanings

<!-- glossary: name | text -->
| name | text |
|---|---|
| `iout` | unit of the text output |
| `ifail` | 0, or the kind of failure (failure_m), whose message is written |
| `at` | where the run was (failure_m): its stage, l and ikap |
| `warns` | the warnings, to which it adds its own (failure_m) |
| `c_light` | speed of light (Hartree atomic units) |
| `rr` | log radial mesh (bohr) |
| `mmax` | size of the log radial mesh |
| `al` | log of the mesh ratio |
| `zz` | atomic number |
| `zion` | ion charge |
| `irel` | 0 non-relativistic, 1 scalar-relativistic |
| `iexc` | exchange-correlation functional |
| `ll` | angular momentum |
| `kap` | Dirac kappa, l or -(l+1) for j = l -/+ 1/2 |
| `nkap` | number of kappas per l: 1, or 2 for irel=2 |
| `lmax` | maximum angular momentum |
| `lloc` | l of the local potential, 4 for the polynomial |
| `nc` | number of core states |
| `nv` | number of valence states |
| `na` | principal quantum numbers |
| `la` | angular momenta |
| `fa` | occupancies |
| `ncnf` | number of test configurations (not counting the reference) |
| `nproj` | number of projectors for each l |
| `evkb` | coefficients of the projectors |
| `vkbpw` | the projectors as pieces |
| `vlocpw` | the screened local potential as pieces |
| `up` | d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr) |
| `rc` | core radius (bohr) / core radii (bohr) |
| `rcmax` | largest core radius (bohr) |
| `irc` | last mesh point at or inside rc |
| `rhoc` | core charge density |
| `rhomod` | model core charge |
| `rhotae` | all-electron valence charge density |
| `rhotps` | pseudo valence charge density |
| `vfull` | all-electron potential (Ha) |
| `vxc` | exchange-correlation potential (Ha) |
| `exc` | exchange-correlation energy per electron (Ha) |
| `eterm` | kinetic, electron-nucleus, Hartree and exchange-correlation energies (Ha) |
| `nrl` | size of the output mesh |
| `rl` | output mesh (bohr) |
| `nder` | highest derivative wanted, 0<=nder<=4 |
| `atsym` | atomic symbol |
| `inp` | the input values (input_m) |
| `cnf` | the reference and test configurations (input_m) |
| `epsh1` | lower limit of the log derivative analysis energy range (Ha) |
| `epsh2` | upper limit of the log derivative analysis energy range (Ha) |
| `depsh` | energy step of the log derivative analysis (Ha) |

## Names with several meanings

These names mean different things in different routines, so they have
no fixed text. Their `@param` line says what the name holds in that
routine, with its unit when it is an energy, a potential or a length.

| name | meanings |
|---|---|
| `ee` | an energy (Ha): a state's eigenvalue, a scattering energy, a trial energy |
| `uu` | a radial wave function (times r) |
| `vv` | a potential (Ha): all-electron, local or screened |
| `ff` | a function on the mesh or on the nodes of a grid |
| `xx` | the argument of a function |
| `nn` | a number of nodes; a principal quantum number |
| `ss` | an integral; an overlap; step weights |
| `mch` | a matching point of the inward and outward integrations; the last point integrated to is `mlast` |
| `rho` | a charge density: total, valence, or pseudo valence |
| `pw` | a function as pieces (pieces_m) |
| `vp` | the semi-local pseudopotentials, screened or not |
| `ierr` | 0 on success; the non-zero values differ by routine |
| `etot` | a total energy (Ha): all-electron or pseudo-atom |
| `fcfact`, `rcfact` | the model core's amplitude and scale factors, whose meaning depends on icmod |
