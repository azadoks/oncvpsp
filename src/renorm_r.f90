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
 subroutine renorm_r(uu,rr,ll,kap,zz,mmax,cnorm)

!! renormalize Dirac wave function so that large component is normalized
!!
!! @param uu  Dirac wave function
!! @param rr  log radial grid
!! @param ll  angular momentum
!! @param kap  Dirac kappa
!! @param zz  atomic number
!! @param mmax  size of radial grid
!! @param cnorm  renormalization coefficient

 use constants_m, only: dp
 use relativity_m, only: c_light
 use quad_m, only: log_norm_sum
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: ll,kap,mmax

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2)

!Output variables
 real(dp), intent(out) :: cnorm

!Local variables
 real(dp) :: sn
 real(dp) :: al,amesh
 real(dp) :: r0,gam,cc,cci
 integer :: ii,nin

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
   write(6,'(/a,i0,a,i0)') 'renorm_r: ERROR wave function is zero, ll=',ll, &
&        ', kap=',kap
   error stop 1
 end if

 r0=rr(1)/sqrt(amesh)
 sn=r0**(2.0d0*gam+1.0d0)/(2.d0*gam+1.0d0)
 sn=sn*(uu(1,1)**2)/rr(1)**(2.d0*gam)

 sn=log_norm_sum(sn,al,rr,uu(:,1),nin)


 cnorm=sqrt(1.0d0/sn)

 uu(:,:)=cnorm*uu(:,:)

 return
 end subroutine renorm_r
