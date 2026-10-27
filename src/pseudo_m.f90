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
!> the atom: its symbol, atomic number and ion charge, the equation (irel
!> 0 non-, 1 scalar-, 2 fully-relativistic), the number of kappas per l
!> (2 when fully relativistic), the exchange-correlation functional and
!> the speed of light (relativity_m)
   character(len=2) :: atsym='  '
   real(dp) :: zz=0.0d0,zion=0.0d0,c_light=0.0d0
   integer :: irel=0,nkap=0,iexc=0
!> the log radial mesh and its size
   integer :: mmax=0
   real(dp), allocatable :: rr(:)
!> the largest l of the projectors, l of the local potential (4 for the
!> polynomial), the number of projectors of each l, and the largest core
!> radius
   integer :: lmax=0,lloc=0,nproj(mxlin+1)=0
   real(dp) :: rcmax=0.0d0
!> the projectors: their coefficients (mxprj,mxl+1,nkap), on the mesh
!> (mmax,mxprj,mxl+1,nkap; not in the data file) and as pieces
!> (mxprj,mxl+1,nkap)
   real(dp), allocatable :: evkb(:,:,:),vkb(:,:,:,:)
   type(piecewise_t), allocatable :: vkbpw(:,:,:)
!> fully relativistic (made for any irel, written for nkap=2): the
!> scalar-relativistic and spin-orbit projectors' coefficients and pieces
!> (2*mxprj,mxl+1; not in the data file)
   real(dp), allocatable :: esr(:,:),eso(:,:)
   type(piecewise_t), allocatable :: vsrpw(:,:),vsopw(:,:)
!> the screened local potential's pieces, and the unscreened semi-local
!> potentials (mmax,mxl+2), vpuns(:,lloc+1) the local one (the only one
!> in the data file)
   type(piecewise_t) :: vlocpw
   real(dp), allocatable :: vpuns(:,:)
!> the valence pseudo-density of the reference configuration, the model
!> core charge and its first four derivatives (mmax,nrhomod; only the
!> charge in the data file), and their pieces (not in the data file)
   real(dp), allocatable :: rho(:),rhomod(:,:)
   type(piecewise_t) :: corepw(nrhomod)
 end type pseudo_t

 end module pseudo_m
