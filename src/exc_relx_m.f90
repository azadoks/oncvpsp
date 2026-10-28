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

 module exc_relx_m

 implicit none
 private
 public :: exc_relx

 contains

 subroutine exc_relx(rho,vxc,exc,mmax,c_light)

!! adds the relativistic correction of the LDA (Slater) exchange of
!! A. H. MacDonald and S. H. Vosko, J. Phys. C 12, 2977 (1979) to an LDA
!! exchange-correlation potential and energy density: with
!! beta = kF/c_light and mu = sqrt(1 + beta**2), the exchange energy density
!! is multiplied by 1 - 1.5*((beta*mu - asinh(beta))/beta**2)**2 and the
!! potential by 1.5*asinh(beta)/(beta*mu) - 0.5.  This is the exchange of
!! the NIST atomic reference data (ScRLDA, RLDA; iexc=7, after excvwn) and
!! of Libxc's lda_x_rel, whose speed of light is fixed.  rho, vxc and exc
!! are as in excpzca
!!
!! @param rho  charge density
!! @param vxc  exchange-correlation potential (Ha), corrected on output
!! @param exc  exchange-correlation energy per electron (Ha), corrected on
!!              output
!! @param mmax  size of the log radial mesh
!! @param c_light  speed of light (Hartree atomic units)

 use constants_m, only: dp, pi

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: rho(mmax)

!Input/output variables
 real(dp), intent(inout) :: vxc(mmax),exc(mmax)

!Local variables
 real(dp), parameter :: pi4=4.0d0*pi
 integer :: ii
 real(dp) :: akf,bb,bmu,gg,ex,vx

 do ii=1,mmax
   if(rho(ii)>1.0d-20) then
! Fermi wave vector, and the Slater exchange energy density and potential
     akf=(3.0d0*pi**2*rho(ii)/pi4)**(1.0d0/3.0d0)
     ex=-0.75d0*akf/pi
     vx=-akf/pi
     bb=akf/c_light
     bmu=sqrt(1.0d0 + bb**2)
! the corrections R - 1 and S - 1 cancel for small beta: series
     if(bb<1.0d-3) then
       gg=bb*(2.0d0/3.0d0 - bb**2*(0.2d0 - 3.0d0/28.0d0*bb**2))
       vxc(ii)=vxc(ii) + vx*bb**2*(-1.0d0 + 0.8d0*bb**2)
     else
       gg=(bb*bmu - asinh(bb))/bb**2
       vxc(ii)=vxc(ii) + vx*1.5d0*(asinh(bb)/(bb*bmu) - 1.0d0)
     end if
     exc(ii)=exc(ii) - ex*1.5d0*gg**2
   end if
 end do

 return
 end subroutine exc_relx

 end module exc_relx_m
