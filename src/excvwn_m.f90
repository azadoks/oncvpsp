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

 module excvwn_m

 implicit none
 private
 public :: excvwn

 contains

 subroutine excvwn(rho,vxc,exc,mmax)

!! calculates the Slater exchange and Vosko-Wilk-Nusair correlation
!! (the fit to the Ceperley-Alder results, "VWN5" of Libxc's lda_c_vwn)
!! potential and energy density, S. H. Vosko, L. Wilk and M. Nusair,
!! Can. J. Phys. 58, 1200 (1980): the LDA of the NIST atomic reference data
!!
!! @param rho  charge density
!! @param vxc  exchange-correlation potential
!! @param exc  exchange-correlation energy density
!! @param mmax  dimension of log radial grid

 use constants_m, only: dp, pi

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rho(mmax)

!Output variables
 real(dp), intent(out) :: vxc(mmax),exc(mmax)

!Local variables
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: thrd=1.0d0/3.0d0
 real(dp), parameter :: conrs = (3.0d0/(4.0d0*pi))**thrd
 real(dp), parameter :: convx = (1.5d0/pi)**(2.0d0/3.0d0)
 real(dp), parameter :: conex = 0.75d0*convx

!paramagnetic VWN5 parameters
 real(dp), parameter :: aa=0.0310907d0
 real(dp), parameter :: bb=3.72744d0
 real(dp), parameter :: cc=12.9352d0
 real(dp), parameter :: x0=-0.10498d0

 integer :: ii
 real(dp) :: qq,xx0,rs,xx,xxx,at,tq,ec,decdx

 qq=sqrt(4.0d0*cc - bb**2)
 xx0=x0**2 + bb*x0 + cc

 do ii=1,mmax
   if(rho(ii)>1.0d-20) then
     rs=conrs*(rho(ii)/pi4)**(-thrd)
     xx=sqrt(rs)
     xxx=xx**2 + bb*xx + cc
     at=atan(qq/(2.0d0*xx + bb))
     tq=(2.0d0*xx + bb)**2 + qq**2
     ec=aa*(log(xx**2/xxx) + 2.0d0*bb/qq*at &
&      - bb*x0/xx0*(log((xx - x0)**2/xxx) + 2.0d0*(bb + 2.0d0*x0)/qq*at))
     decdx=aa*(2.0d0/xx - (2.0d0*xx + bb)/xxx - 4.0d0*bb/tq &
&      - bb*x0/xx0*(2.0d0/(xx - x0) - (2.0d0*xx + bb)/xxx &
&      - 4.0d0*(bb + 2.0d0*x0)/tq))
! v = e - (rs/3) de/drs, with de/drs = (de/dx)/(2x)
     exc(ii)=-conex/rs + ec
     vxc(ii)=-convx/rs + ec - xx*decdx/6.0d0
   else
     vxc(ii)=0.0d0 ; exc(ii)=0.0d0
   end if
 end do

 return
 end subroutine excvwn

 end module excvwn_m
