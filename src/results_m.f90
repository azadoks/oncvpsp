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

!> What oncvpsp_run returns: the pseudopotential, what the writers of the
!> psp8, UPF and PSML files take besides, and the results of the tests
 module results_m

 use constants_m, only: dp, mxlin, mxstate, mxkap
 use pseudo_m, only: pseudo_t
 use test_results_m, only: tests_t
 implicit none
 private
 public :: results_t, model_core_t, core_budget_t

!> the input values as read, before oncvpsp_run changes them, for the
!> echoes of the input in the files
 type :: echo_t
!> the core radii (bohr)
   real(dp) :: rc_in(mxlin+1)=0.0d0
!> the model core's amplitude factor fcfact
   real(dp) :: fcfact_in=0.0d0
!> the icmod=4 search grid of fcfact (min, max, step)
   real(dp) :: fcgrid(3)=0.0d0
!> the icmod=4 search grid of rcfact (min, max, step)
   real(dp) :: rcgrid(3)=0.0d0
!> the log mesh line (amesh, rr1, rmax and the Adams order; 0 when the
!> input does not set it)
   real(dp) :: log_mesh(4)=0.0d0
 end type echo_t

!> the budgeted model core (modcore_budget, icmod 6 and 7), as its report
 type :: core_budget_t
!> the cutoff budget (Ha)
   real(dp) :: ecut=0.0d0
!> the exchange-correlation energy error (Ha) that defines the needed
!> cutoff
   real(dp) :: error=0.0d0
!> the share of the budget the model must fit
   real(dp) :: slack=0.0d0
!> the pseudo-waves' cutoff (Ha) the budget is factor times, when the
!> input's fcfact <= 0; else 0
   real(dp) :: ecutw=0.0d0
!> that multiple, 0 when the input's fcfact > 0
   real(dp) :: factor=0.0d0
!> the model radius r_s (bohr)
   real(dp) :: r_s=0.0d0
!> the smallest rc (bohr)
   real(dp) :: rc_min=0.0d0
!> the last radius with core charge (bohr) when it is below rc_min and so
!> bounds r_s; else 0
   real(dp) :: r_core=0.0d0
!> the fill inside r_s: 1 the minimum tail, 2 Martins'
   integer :: fill=0
!> the minimum tail's coefficients c0..c8 of (r/r_s)**(2k)
   real(dp) :: cpoly(0:8)=0.0d0
!> Martins' fill exp(a + b r**2 + c r**4): a, b, c
   real(dp) :: martins(3)=0.0d0
!> the cutoff the model core needs (Ha)
   real(dp) :: ecut_needed=0.0d0
!> icmod 7: the radius of icmod 6, the smallest that fits (bohr)
   real(dp) :: r_s6=0.0d0
!> icmod 7: the rms hardness mismatch of the model (Ha)
   real(dp) :: mismatch=0.0d0
!> icmod 7: the radii of the basin (bohr)
   real(dp), allocatable :: crad(:)
!> their fills, 0 for none
   integer, allocatable :: cfill(:)
!> their hardness mismatches (Ha)
   real(dp), allocatable :: cmis(:)
 end type core_budget_t

!> the model core's parameters
 type :: model_core_t
!> icmod 3, 4 (modcore3): the Teter function's amplitude prefactor fcfact
!> and scale prefactor rcfact (those of the input, or optimized), the
!> core-valence crossover radius rmatch (bohr) and the core charge there
!> rhocmatch
   real(dp) :: teter(4)=0.0d0
!> icmod 6, 7: the budgeted core
   type(core_budget_t) :: budget
 end type model_core_t

!> the results of a run
 type :: results_t
!> the pseudopotential
   type(pseudo_t) :: pot
!> the size of the linear mesh of the files
   integer :: nrl=0
!> the pseudo-atom's total energy (Ha)
   real(dp) :: epstot=0.0d0
!> the all-electron eigenvalues (Ha) of the reference configuration
!> (mxstate,mxkap), the UPF's orbital energies
   real(dp) :: ea(mxstate,mxkap)=0.0d0
!> all-electron valence charge density
   real(dp), allocatable :: rhotae(:)
!> core charge density, that of the psp8 file
   real(dp), allocatable :: rhoc(:)
!> the PSML file's radius beyond which the model core is the core
   real(dp) :: rcore=0.0d0
!> the PSML file's semi-local potentials (mmax,mxl+2,nkap)
   real(dp), allocatable :: vpsml(:,:,:)
!> the pseudo-atomic orbitals (mmax,nkap,nv), the UPF's
   real(dp), allocatable :: uupsa(:,:,:)
!> the model core's parameters
   type(model_core_t) :: core
!> the input as read
   type(echo_t) :: echo
!> the tests
   type(tests_t) :: tests
 end type results_t

 end module results_m
