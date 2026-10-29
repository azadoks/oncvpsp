!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
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

 module exchdl_m

 implicit none
 private
 public :: exchdl, fxchdl

 contains

 subroutine exchdl(rho,vxc,exc,mmax)

!! calculates Hedin-Lundquist exchange-correlation potential and energy
!! density
!!
!! @param rho  charge density times 4 pi
!! @param vxc  exchange-correlation potential (Ha)
!! @param exc  exchange-correlation energy per electron (Ha)
!! @param mmax  size of the log radial mesh

 use constants_m, only: dp, pi

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rho(mmax)

!Output variables
 real(dp), intent(out) :: vxc(mmax),exc(mmax)

!Local variables
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4
 real(dp), parameter :: thrd=1.0d0/3.0d0
 integer :: ii
 real(dp) :: conrs,convx,conex
 real(dp) :: rs,ecp,aln,xx,rh

 conrs = (3.0d0/(4.0d0*pi))**thrd
 convx = (1.5d0/pi)**(2.0d0/3.0d0)
 conex = 0.75d0*convx

! Hedin-Lundqvist correlation

 do ii=1,mmax
   if(rho(ii)>1.0d-20) then
     rh=rho(ii)*pi4i
     rs=conrs/rh**thrd
     xx=rs/21.0d0
     aln=log(1.0d0 + 1.0d0/xx)
     ecp = aln+(xx**3*aln-xx*xx)+0.5d0*xx-thrd
     exc(ii)=-conex/rs - 0.0225d0*ecp
     vxc(ii)=-convx/rs - 0.0225d0*aln
   else
     vxc(ii)=0.0d0 ; exc(ii)=0.0d0
   end if
 end do

 return
 end subroutine exchdl

 subroutine fxchdl(rho,fxc,mmax)

!! the derivative by rho of exchdl's potential, zero where exchdl's is
!!
!! @param rho  charge density
!! @param fxc  derivative of the exchange-correlation potential by rho
!!             (Ha bohr^3)
!! @param mmax  size of the log radial mesh

 use constants_m, only: dp, pi

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rho(mmax)

!Output variables
 real(dp), intent(out) :: fxc(mmax)

!Local variables
 real(dp), parameter :: pi4i=1.0d0/(4.0d0*pi)
 real(dp), parameter :: thrd=1.0d0/3.0d0
 real(dp) :: rs,xx,conrs,convx,dvdrs
 integer :: ii

 conrs = (3.0d0/(4.0d0*pi))**thrd
 convx = (1.5d0/pi)**(2.0d0/3.0d0)

! dv/drs, and drs/drho = -rs/(3 rho)
 do ii=1,mmax
   if(rho(ii)>1.0d-20) then
     rs=conrs/(rho(ii)*pi4i)**thrd
     xx=rs/21.0d0
     dvdrs=convx/rs**2 + 0.0225d0/(21.0d0*xx*(xx + 1.0d0))
     fxc(ii)=-dvdrs*rs*thrd/rho(ii)
   else
     fxc(ii)=0.0d0
   end if
 end do
 return
 end subroutine fxchdl

 end module exchdl_m
