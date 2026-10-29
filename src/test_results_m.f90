!
! Copyright (c) 2026 by Austin Zadoks
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!

!> The results of the tests of a pseudopotential, one type per test, as
!> oncvpsp_run returns them and diagout writes them (the DIAGNOSTICS
!> block). ikap is 1 for j=l+1/2 (and for nkap=1), 2 for j=l-1/2; the
!> arrays per projector, l and ikap are (mxprj,mxl+1,nkap).
 module test_results_m

 use constants_m, only: dp, mxl, mxstate, mxkap, mxcnf, nerrlv
 implicit none
 private
 public :: tests_t, logder_t, ghosts_t, configs_t, convergence_t, &
&           weighted_cut_t, vloc_cut_t, mxghost

!> the GHOST(-) and GHOST(+) states kept per l and ikap
 integer, parameter :: mxghost=20

!> the log derivatives, atan(r u'/u) (run_phsft)
 type :: logder_t
!> the number of energies, 0 with --no-log-derivatives
   integer :: ne=0
!> the energies (Ha), from epsh2 down by depsh
   real(dp), allocatable :: ee(:)
!> the radius of each l (bohr)
   real(dp) :: rl(4)=0.0d0
!> of the all-electron wave, per energy, l and ikap
   real(dp), allocatable :: ae(:,:,:)
!> of the pseudo wave, per energy, l and ikap
   real(dp), allocatable :: ps(:,:,:)
 end type logder_t

!> the ghost tests (run_ghosts), per l and ikap
 type :: ghosts_t
!> the number of bound states of the nonlocal Schroedinger equation
   integer :: nrow(mxl+1,mxkap)=0
!> their n
   integer :: nn(mxstate,mxl+1,mxkap)=0
!> their energy and that of the barrier-basis eigenstate matched to it
!> (Ha; e_ref 0 for the vacuum level when l has none); e_basis below
!> e_ref is a GHOST(-)
   real(dp) :: ee(2,mxstate,mxl+1,mxkap)=0.0d0
!> of the basis eigenstates more than 0.05 Ha above the reference, the
!> one of least <r>/rc: its energy (Ha) and <r>/rc, a GHOST(+) below 1;
!> 0 if none
   real(dp) :: least(2,mxl+1,mxkap)=0.0d0
!> the number of GHOST(-) states, basis eigenvalues below a bound state
!> of the nonlocal Schroedinger equation (or the vacuum level)
   integer :: nminus(mxl+1,mxkap)=0
!> their energies (Ha), the first mxghost
   real(dp) :: minus(mxghost,mxl+1,mxkap)=0.0d0
!> the number of GHOST(+) states, basis states above the reference whose
!> weighted <r>/rc (as least) is below 1
   integer :: nplus(mxl+1,mxkap)=0
!> their energy (Ha) and weighted <r>/rc, the first mxghost
   real(dp) :: plus(2,mxghost,mxl+1,mxkap)=0.0d0
!> the energy (Ha) the basis goes to, the cutoff of the test's lines:
!> 1.25 times the largest kinetic energy of the sbf basis, or lower where
!> the basis was cut
   real(dp) :: ecut(mxl+1,mxkap)=0.0d0
!> the basis: 0 complete to its cutoff, 1 cut at a barrier state below the
!> cutoff that did not converge (the test goes to ecut), 2 full below the
!> cutoff (its largest size reached)
   integer :: basis(mxl+1,mxkap)=0
 end type ghosts_t

!> the test configurations, the reference first (run_config_ps)
 type :: configs_t
!> their number
   integer :: ncfg=0
!> the rows of each one's eigenvalue table, -1 if either atom failed
   integer :: nrow(mxcnf)=0
!> n, l, ikap of each row
   integer :: state(3,2*mxstate,mxcnf)=0
!> f, the all-electron and the pseudo eigenvalue (Ha) of each row
   real(dp) :: ee(3,2*mxstate,mxcnf)=0.0d0
!> the all-electron and pseudo total-energy differences from the
!> reference configuration and the excitation error de_ae - de_ps (Ha)
   real(dp) :: de(3,mxcnf)=0.0d0
 end type configs_t

!> the convergence of the optimized pseudo-waves with the cutoff
!> (run_optimize), per projector, l and ikap
 type :: convergence_t
!> (Ecut, err) where E_resid per electron crosses each level err, from
!> 1e-2 down by half decades, Ecut 0 for a level not reached
!> (2,nerrdg,mxprj,mxl+1,nkap)
   real(dp), allocatable :: cross(:,:,:,:,:)
!> the curves, (Ecut, E_resid per electron) on the q grid, zero beyond
!> their ends (2,:,mxprj,mxl+1,nkap)
   real(dp), allocatable :: curve(:,:,:,:,:)
!> (1) E_resid from q_infinity to three times that, the curve's tail;
!> (2) the r-space less the Fourier kinetic energy to three times
!> q_infinity (2,mxprj,mxl+1,nkap)
   real(dp), allocatable :: tail(:,:,:,:)
!> E_resid at qcut, the regularization's term, and the smallest
!> eigenvalue e1 and |f1| of the minimization at qcut (4,mxprj,mxl+1,nkap)
   real(dp), allocatable :: opt(:,:,:,:)
!> of a well state confined by its well, E_resid at qcut with the well's
!> tail and with the free tail; (3) 1 for these, -1 where the free tail
!> is undefined, 0 for the other waves (3,mxprj,mxl+1,nkap)
   real(dp), allocatable :: wshare(:,:,:,:)
 end type convergence_t

!> the weighted cutoff of the model core's budget, the weighted mean of
!> the pseudo-waves' convergence curves (modcore_budget_table)
 type :: weighted_cut_t
!> whether the budget took it (icmod 6, 7 with fcfact <= 0 and
!> budget_cutoff "weighted")
   logical :: used=.false.
!> its error per electron (Ha)
   real(dp) :: errw=0.0d0
!> the cutoff (Ha)
   real(dp) :: ecutw=0.0d0
!> the weight of each pseudo-wave in the mean (mxprj,mxl+1,nkap)
   real(dp), allocatable :: weight(:,:,:)
 end type weighted_cut_t

!> the convergence of the local energy of the reference pseudo valence
!> density (vloc_cutoff)
 type :: vloc_cut_t
!> (Ecut, err) per level, from 1e-2 down by half decades, Ecut 0 for a
!> level not reached
   real(dp) :: cut(2,nerrlv)=0.0d0
!> the local energy per electron from the q grid less that from r space
!> (Ha), the error of the transforms
   real(dp) :: check=0.0d0
 end type vloc_cut_t

!> the results of all the tests
 type :: tests_t
!> the log derivatives
   type(logder_t) :: logder
!> the ghost tests
   type(ghosts_t) :: ghosts
!> the test configurations
   type(configs_t) :: configs
!> the convergence of the optimized pseudo-waves with the cutoff
   type(convergence_t) :: convergence
!> the weighted cutoff of the model core's budget
   type(weighted_cut_t) :: weighted_cut
!> the convergence of the local energy
   type(vloc_cut_t) :: vloc_cut
!> the eigenvalues of the projectors' symmetrized B matrices (run_vkb),
!> their condition (mxprj,mxl+1,nkap)
   real(dp), allocatable :: bbeig(:,:,:)
 end type tests_t

 end module test_results_m
