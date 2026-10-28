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

 module fpovlp_m

 implicit none
 private
 public :: fpovlp

 contains

! calculates all-electron wave function overlap intergrals

 subroutine fpovlp(uu1,uu2,mmax,rc,ll,kap,zz,ss,rr,irel,c_light)

!! product of all-electron wave functions uu1*uu2 goes like rr**(2*gamma) as
!! rr -> 0, with the exponent of the equation that gave them: l+1, that of
!! the Koelling-Harmon equation (sr_gamma) or, for the large components of
!! the Dirac equation, sqrt(kap**2-(zz/c)**2), as in renorm_r
!! integral from rr=0 to rc: the power law below rr(1), the interpolated
!! cumulative integral of lm_integral from there, so smooth in rc
!! (doc/smoothness_plan.md, R-5)
!!
!! @param uu1  radial wave function one (times r)
!! @param uu2  radial wave function two (times r)
!! @param mmax  size of the log radial mesh
!! @param rc  core radius (bohr), the upper limit of the integral
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2 (irel=2 only)
!! @param zz  atomic number
!! @param ss  overlap output
!! @param rr  log radial mesh (bohr)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)

 use constants_m, only: dp
 use sr_flux_m, only: sr_gamma
 use logmesh_m, only: lm_integral

!Input variables
 integer, intent(in) :: mmax,ll,kap
 real(dp), intent(in) :: zz,rc
 real(dp), intent(in) :: uu1(mmax),uu2(mmax),rr(mmax)
 integer, intent(in) :: irel
 real(dp), intent(in) :: c_light

!Output variables
 real(dp), intent(out) :: ss

!Local variables
 real(dp) :: fss,gamma,si

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

! the piece below rr(1) from the power law, the rest by the cumulative
! integral (whose own start cancels between its two limits)
 ss=uu1(1)*uu2(1)*rr(1)/(2.0d0*gamma+1.0d0)
 call lm_integral(rr,uu1*uu2,mmax,2*ll+2,rr(1),rc,si)
 ss=ss+si

 return
 end subroutine fpovlp

 end module fpovlp_m
