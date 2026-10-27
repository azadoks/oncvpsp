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

 subroutine renorm_r(uu,rr,ll,kap,c_light,zz,mmax,cnorm,iout,ifail)

!! renormalize Dirac wave function so that large component is normalized
!!
!! @param uu  Dirac wave function
!! @param rr  log radial grid
!! @param ll  angular momentum
!! @param kap  Dirac kappa
!! @param c_light  speed of light (Hartree atomic units)
!! @param zz  atomic number
!! @param mmax  size of radial grid
!! @param cnorm  renormalization coefficient
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use quad_m, only: log_norm
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: ll,kap,mmax
 real(dp), intent(in) :: c_light

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2)

!Output variables
 real(dp), intent(out) :: cnorm
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: sn
 real(dp) :: al,amesh
 real(dp) :: gam,cc,cci
 integer :: ii,nin

 ifail=0
 al = 0.01d0 * log(rr(101)/rr(1))
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

 sn=log_norm(rr,uu(1:nin,1)**2,nin,ll,gam)


 cnorm=sqrt(1.0d0/sn)

 uu(:,:)=cnorm*uu(:,:)

 return
 end subroutine renorm_r

 end module renorm_r_m
