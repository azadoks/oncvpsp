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

 module excpzca_m

 implicit none
 private
 public :: excpzca, fxcpzca

 contains

 subroutine excpzca(rho,vxc,exc,mmax)

!! calculates Perdew-Zunger-Ceperly_Alder exchange-correlation potential
!! and energy density J. P. Perdew and A. Zunger, Phys. Rev. B23, 5048 (1981)
!!
!! @param rho  charge density
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
 real(dp), parameter :: thrd=1.0d0/3.0d0
 integer :: ii
 real(dp) :: rs,rh,sqrs,rsl,den

 real(dp), parameter :: conrs = (3.0d0/(4.0d0*pi))**thrd
 real(dp), parameter :: convx = (1.5d0/pi)**(2.0d0/3.0d0)
 real(dp), parameter ::conex = 0.75d0*convx


!rs>1, Ceperley-Alder fit to QMC
 real(dp), parameter :: gam = -0.1423d0
 real(dp), parameter :: bt1 =  1.0529d0
 real(dp), parameter :: bt2 =  0.3334d0

!rs<1, Perdew-Zunger leading terms of RPA high-density expansion
!AA and BB known from Gell-Mann & Bruckner
!PZ adjust CC and DD for contiunity of exc and vxc at rs=1
 real(dp), parameter :: AA =   0.0311d0
 real(dp), parameter :: BB =  -0.0480d0
 real(dp), parameter :: CC =   0.0020d0
 real(dp), parameter :: DD =  -0.0116d0

!libxc uses PZ's modified coefficients CC=0.0020191519406228 and
!DD=-0.0116320663789130

 do ii=1,mmax
   if(rho(ii)>1.0d-20) then
     rh=rho(ii)/pi4
     rs=conrs*rh**(-thrd)
     if(rs > 1.0d0) then
       sqrs=sqrt(rs)
       den=1.0d0 + bt1*sqrs +bt2*rs
       exc(ii)=-conex/rs + gam/den
       vxc(ii)=-convx/rs + gam*(1.0d0 + (3.5d0*bt1*sqrs &
&              + 4.0d0*bt2*rs)*thrd)/den**2
     else
       rsl=log(rs)
       exc(ii)=-conex/rs + AA*rsl + BB + CC*rs*rsl + DD*rs
       vxc(ii)=-convx/rs + AA*rsl + (BB-thrd*AA) + 2.0d0*thrd*CC*rs*rsl &
&              + thrd*(2.0d0*DD-CC)*rs
     end if
   else
     vxc(ii)=0.0d0 ; exc(ii)=0.0d0
   end if
 end do

 return
 end subroutine excpzca

 subroutine fxcpzca(rho,fxc,mmax)

!! the derivative by rho of excpzca's potential, zero where excpzca's is,
!! on the branch of rs that excpzca takes: the potential steps at rs = 1
!! (its coefficients CC and DD are rounded), whose delta is left out
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
 real(dp), parameter :: pi4=4.0d0*pi
 real(dp), parameter :: thrd=1.0d0/3.0d0
 real(dp), parameter :: conrs = (3.0d0/(4.0d0*pi))**thrd
 real(dp), parameter :: convx = (1.5d0/pi)**(2.0d0/3.0d0)
 real(dp), parameter :: gam = -0.1423d0
 real(dp), parameter :: bt1 =  1.0529d0
 real(dp), parameter :: bt2 =  0.3334d0
 real(dp), parameter :: AA =   0.0311d0
 real(dp), parameter :: CC =   0.0020d0
 real(dp), parameter :: DD =  -0.0116d0
 real(dp) :: rs,sqrs,den,dden,num,dnum,dvdrs
 integer :: ii

! dv/drs, and drs/drho = -rs/(3 rho)
 do ii=1,mmax
   if(rho(ii)>1.0d-20) then
     rs=conrs*(rho(ii)/pi4)**(-thrd)
     if(rs > 1.0d0) then
       sqrs=sqrt(rs)
       den=1.0d0 + bt1*sqrs + bt2*rs
       dden=0.5d0*bt1/sqrs + bt2
       num=1.0d0 + (3.5d0*bt1*sqrs + 4.0d0*bt2*rs)*thrd
       dnum=(1.75d0*bt1/sqrs + 4.0d0*bt2)*thrd
       dvdrs=convx/rs**2 + gam*(dnum*den - 2.0d0*num*dden)/den**3
     else
       dvdrs=convx/rs**2 + AA/rs + 2.0d0*thrd*CC*(log(rs) + 1.0d0) &
&            + thrd*(2.0d0*DD - CC)
     end if
     fxc(ii)=-dvdrs*rs*thrd/rho(ii)
   else
     fxc(ii)=0.0d0
   end if
 end do
 return
 end subroutine fxcpzca

 end module excpzca_m
