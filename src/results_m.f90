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
 public :: results_t, echo_t

!> the input values as read, before oncvpsp_run changes them, for the
!> echoes of the input in the files
 type :: echo_t
!> the core radii (bohr)
   real(dp) :: rc0(mxlin+1)=0.0d0
!> the model core's amplitude factor fcfact
   real(dp) :: fcfact0=0.0d0
!> the icmod=4 search grid of fcfact (min, max, step)
   real(dp) :: fcgrid(3)=0.0d0
!> the icmod=4 search grid of rcfact (min, max, step)
   real(dp) :: rcgrid(3)=0.0d0
!> the log mesh (amesh, rr1, rmax; 0 when the input does not set it)
   real(dp) :: log_mesh(3)=0.0d0
 end type echo_t

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
!> the input as read
   type(echo_t) :: echo
!> the tests
   type(tests_t) :: tests
 end type results_t

 end module results_m
