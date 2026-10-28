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

!> A pseudopotential: what oncvpsp_run makes (results_m) and the
!> pseudopotential data file holds (ps_file_m). ikap is 1 for j=l+1/2
!> (and for nkap=1), 2 for j=l-1/2.
 module pseudo_m

 use constants_m, only: dp, mxlin, nrhomod
 use pieces_m, only: piecewise_t
 implicit none
 private
 public :: pseudo_t

!> the pseudopotential
 type :: pseudo_t
!> atomic symbol
   character(len=2) :: atsym='  '
!> atomic number
   real(dp) :: zz=0.0d0
!> ion charge
   real(dp) :: zion=0.0d0
!> speed of light (Hartree atomic units) of the run (relativity_m)
   real(dp) :: c_light=0.0d0
!> 0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
   integer :: irel=0
!> number of kappas per l: 1, or 2 when fully relativistic
   integer :: nkap=0
!> exchange-correlation functional
   integer :: iexc=0
!> size of the log radial mesh
   integer :: mmax=0
!> log radial mesh (bohr)
   real(dp), allocatable :: rr(:)
!> the largest l of the projectors
   integer :: lmax=0
!> l of the local potential, 4 for the polynomial
   integer :: lloc=0
!> number of projectors for each l
   integer :: nproj(mxlin+1)=0
!> largest core radius (bohr)
   real(dp) :: rcmax=0.0d0
!> coefficients of the projectors (Ha) (mxprj,mxl+1,nkap)
   real(dp), allocatable :: evkb(:,:,:)
!> the projectors on the mesh (mmax,mxprj,mxl+1,nkap; not in the data
!> file)
   real(dp), allocatable :: vkb(:,:,:,:)
!> the projectors as pieces (mxprj,mxl+1,nkap)
   type(piecewise_t), allocatable :: vkbpw(:,:,:)
!> fully relativistic (made for any irel, written for nkap=2): the
!> scalar-relativistic projectors' coefficients (2*mxprj,mxl+1; not in
!> the data file)
   real(dp), allocatable :: esr(:,:)
!> the spin-orbit projectors' coefficients, as esr
   real(dp), allocatable :: eso(:,:)
!> the scalar-relativistic projectors as pieces (2*mxprj,mxl+1; not in
!> the data file)
   type(piecewise_t), allocatable :: vsrpw(:,:)
!> the spin-orbit projectors as pieces, as vsrpw
   type(piecewise_t), allocatable :: vsopw(:,:)
!> the screened local potential's pieces
   type(piecewise_t) :: vlocpw
!> the unscreened semi-local potentials (mmax,mxl+2), vpuns(:,lloc+1)
!> the local one (the only one in the data file)
   real(dp), allocatable :: vpuns(:,:)
!> the valence pseudo-density of the reference configuration
   real(dp), allocatable :: rho(:)
!> the model core charge and its first four derivatives (mmax,nrhomod;
!> only the charge in the data file)
   real(dp), allocatable :: rhomod(:,:)
!> the pieces of rhomod (not in the data file)
   type(piecewise_t) :: corepw(nrhomod)
 end type pseudo_t

 end module pseudo_m
