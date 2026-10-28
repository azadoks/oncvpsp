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

 module renorm_r_m

 implicit none
 private
 public :: renorm_r

 contains

 subroutine renorm_r(uu,rr,ll,kap,c_light,zz,mmax,cnorm,rnorm,iout,ifail)

!! renormalize Dirac wave function so that large component is normalized
!!
!! @param uu  Dirac wave function
!! @param rr  log radial mesh (bohr)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2
!! @param c_light  speed of light (Hartree atomic units)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param cnorm  renormalization coefficient
!! @param rnorm  the radius the large component is normalized to (a
!!               scattering state's matching radius); absent, its last
!!               nonzero point
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use quad_m, only: quad_log_norm, quad_log_norm_to
 use logmesh_m, only: lm_step

!Input variables
 integer, intent(in) :: ll,kap,mmax
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: c_light
 real(dp), intent(in), optional :: rnorm
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2)

!Output variables
 real(dp), intent(out) :: cnorm
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: sn
 real(dp) :: al,amesh
 real(dp) :: gam,cc,cci
 integer :: ii,nin

 ifail=0
 al=lm_step(rr,mmax)
 amesh = exp(al)

 cc=c_light
 cci=1.0d0/cc

 gam=sqrt(kap**2-(zz*cci)**2)

 nin=0
 do ii=mmax,1,-1
   if(abs(uu(ii,1))>0.0d0) then
     nin=ii
     exit
   end if
 end do
 if(nin<4) then
   write(iout,'(/a,i0,a,i0)') 'renorm_r: ERROR wave function is zero, ll=',ll, &
&        ', kap=',kap
   ifail=1
   return
 end if

 if(present(rnorm)) then
   sn=quad_log_norm_to(rr,uu(1:nin,1)**2,nin,ll,gam,rnorm)
 else
   sn=quad_log_norm(rr,uu(1:nin,1)**2,nin,ll,gam)
 end if


 cnorm=sqrt(1.0d0/sn)

 uu(:,:)=cnorm*uu(:,:)

 return
 end subroutine renorm_r

 end module renorm_r_m
