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

 module excwig_m

 implicit none
 private
 public :: excwig

 contains

 subroutine excwig(rho,vxc,exc,mmax)

!! Wigner ingerpolation formula for exchange-correlation potential and
!! energy density
!!
!! @param rho  charge density
!! @param vxc  exchange-correlation potential
!! @param exc  exchange-correlation energy density
!! @param mmax  dimension of log radial grid

 use constants_m, only: dp, pi
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) ::  rho(mmax)

!Output variables
 real(dp), intent(out) :: vxc(mmax),exc(mmax)

!Local variables
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: pi4i=1.0d0/pi4
 real(dp), parameter :: thrd=1.0d0/3.0d0
 real(dp) ::  rh,rs
 real(dp) :: conrs,convx,conex
 integer :: ii

 conrs = (3.d0/(4.d0*pi))**thrd
 convx = (1.5d0/pi)**(2.0d0/3.0d0)
 conex = 0.75d0*convx

! Wigner interpolation formula, E. Wigner, Phys. Rev. 46, 1002 (1934).

 do ii=1,mmax
   if(rho(ii)>=1.0d-15) then
     rh=rho(ii)*pi4i
     rs=conrs*rh**(-thrd)
     vxc(ii)=-convx/rs - 0.44d0*(4.0d0*thrd*rs+7.8d0)/(rs+7.8d0)**2
     exc(ii)=-conex/rs - 0.44d0/(rs+7.8d0)
   else
     vxc(ii)=0.0d0 ; exc(ii)=0.0d0
   end if
 end do
 return
 end subroutine excwig

 end module excwig_m
