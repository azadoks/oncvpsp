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
! calculates all-electron wave function overlap intergrals

 subroutine fpovlp(gg,hh,nn,ll,kap,zz,ss,rr,irel)

!! product of all-electron wave functions gg*hh goes like rr**(2*gamma) as
!! rr -> 0, with the exponent of the equation that gave them: l+1, that of
!! the Koelling-Harmon equation (sr_gamma) or, for the large components of
!! the Dirac equation, sqrt(kap**2-(zz/c)**2), as in renorm_r
!! integral on usual log mesh from rr=0 to rr(nn)
!!
!! @param gg  wave function one
!! @param hh  wave function two
!! @param nn  outer integration limit
!! @param ll  angular momentum
!! @param kap  Dirac kappa (irel=2 only)
!! @param zz  nuclear charge
!! @param ss  overlap output
!! @param rr  log radial mesh
!! @param irel  0 non-relativistic, 1 scalar-, 2 fully-relativistic

 use constants_m, only: dp
 use relativity_m, only: c_light
 use sr_flux_m, only: sr_gamma
 implicit none

!Input variables
 integer, intent(in) :: nn,ll,kap
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: gg(nn),hh(nn),rr(nn)
 integer, intent(in) :: irel

!Output variables
 real(dp), intent(out) :: ss

!Local variables
 real(dp) :: r0,amesh,al,fss,gamma
 integer :: ii

 al = 0.01d0 * log(rr(101)/rr(1))
 amesh = exp(al)

 if(irel/=0) then
   fss=(1.0d0/c_light)**2
 else
   fss=0.0d0
 end if

 if(irel==2) then
   gamma=sqrt(kap**2-fss*zz**2)
 else
   gamma=sr_gamma(ll,zz,fss)
 end if

 r0=rr(1)/sqrt(amesh)
 ss=r0**(2.0d0*gamma+1.0d0)/(2.d0*gamma+1.0d0)
 ss=(ss*gg(1)*hh(1))/rr(1)**(2.d0*gamma)

 do ii = 1, nn - 3
   ss =  ss + al*gg(ii)*hh(ii)*rr(ii)
 end do

 ss=ss + al*(23.d0*rr(nn-2)*gg(nn-2)*hh(nn-2) &
&        + 28.d0*rr(nn-1)*gg(nn-1)*hh(nn-1) &
&        +  9.d0*rr(nn  )*gg(nn  )*hh(nn  ))/24.d0


 return
 end subroutine fpovlp
