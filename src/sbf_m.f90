!
! Copyright (c) 2026 by Austin Zadoks
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
!> Spherical Bessel functions of low order in double precision from the
!> sine and cosine of the argument (doc/smoothness_plan.md, Q-10). The
!> closed form of j_l is a sum of terms of order x/x**(l+1) that cancels
!> to x**l at small x: from exact sin and cos it loses about 3 eps/x**2,
!> 45 eps/x**4, 1575 eps/x**6 relative for l=1, 2, 3, and from sin and cos
!> with an absolute error d (eresid rotates them from one q to the next)
!> d (2l-1)!!(2l+1)!!/x**(2l+1). Below a threshold per l the power series
!> takes over, which does not use them; at the thresholds no term of the
!> series exceeds its first and the sum is at least 0.48 of it (l=3), so
!> both forms are good to a few eps there and the switch is invisible.

 module sbf_m
 use constants_m, only: dp
 implicit none
 private
 public :: sbf_trig, sbf_trig_rising

!> below these arguments j_0..j_3 come from the power series
 real(dp), parameter :: xser(0:3)=[0.5d0,1.5d0,2.5d0,3.5d0]

 contains

 elemental function sbf_trig(ll,xx,ss,cc) result(jl)

!! j_l(xx) for l=0..3 from ss=sin(xx) and cc=cos(xx): the power series
!! below xser(l), else the closed form
!!
!! @param ll  angular momentum l, 0..3
!! @param xx  argument, >=0
!! @param ss,cc  sin(xx), cos(xx)

!Input variables
 integer, intent(in) :: ll
 real(dp), intent(in) :: xx,ss,cc

!Output variables
 real(dp) :: jl

 if(xx<xser(ll)) then
   jl=sbf_series(ll,xx)
   return
 end if
 select case(ll)
 case(0)
   jl=ss/xx
 case(1)
   jl=(ss-xx*cc)/xx**2
 case(2)
   jl=((3.0d0-xx**2)*ss - 3.0d0*xx*cc)/xx**3
 case default
   jl=((15.0d0-6.0d0*xx**2)*ss - (15.0d0-xx**2)*xx*cc)/xx**4
 end select
 end function sbf_trig

 subroutine sbf_trig_rising(ll,nn,qq,rlin,ss,cc,rjl)

!! rlin*j_l(qq*rlin) by sbf_trig on nondecreasing radii, each form in a loop
!! of its own (the closed forms vectorize): the points below xser(l) come
!! first. The radial transforms' integrand, without a pass of its own
!!
!! @param ll  angular momentum l, 0..3
!! @param nn  number of points
!! @param qq  wave vector, >=0
!! @param rlin  radii of the linear mesh (bohr), >=0, nondecreasing
!! @param ss,cc  sin(qq*rlin), cos(qq*rlin)
!! @param rjl  rlin*j_l(qq*rlin)

!Input variables
 integer, intent(in) :: ll,nn
 real(dp), intent(in) :: qq,rlin(nn),ss(nn),cc(nn)

!Output variables
 real(dp), intent(out) :: rjl(nn)

!Local variables
 real(dp) :: xx
 integer :: ii,i0

 i0=nn+1
 do ii=1,nn
   xx=qq*rlin(ii)
   if(xx>=xser(ll)) then
     i0=ii
     exit
   end if
   rjl(ii)=rlin(ii)*sbf_series(ll,xx)
 end do

 select case(ll)
 case(0)
   do ii=i0,nn
     xx=qq*rlin(ii)
     rjl(ii)=rlin(ii)*(ss(ii)/xx)
   end do
 case(1)
   do ii=i0,nn
     xx=qq*rlin(ii)
     rjl(ii)=rlin(ii)*((ss(ii)-xx*cc(ii))/xx**2)
   end do
 case(2)
   do ii=i0,nn
     xx=qq*rlin(ii)
     rjl(ii)=rlin(ii)*(((3.0d0-xx**2)*ss(ii) - 3.0d0*xx*cc(ii))/xx**3)
   end do
 case default
   do ii=i0,nn
     xx=qq*rlin(ii)
     rjl(ii)=rlin(ii)*(((15.0d0-6.0d0*xx**2)*ss(ii) &
&           - (15.0d0-xx**2)*xx*cc(ii))/xx**4)
   end do
 end select
 end subroutine sbf_trig_rising

 elemental function sbf_series(ll,xx) result(jl)

!! j_l(xx) = xx**l/(2l+1)!! sum_k (-xx**2/2)**k/(k!(2l+3)(2l+5)...(2l+2k+1)),
!! summed until a term is below eps of the sum
!!
!! @param ll  angular momentum l, >=0
!! @param xx  argument, >=0

!Input variables
 integer, intent(in) :: ll
 real(dp), intent(in) :: xx

!Output variables
 real(dp) :: jl

!Local variables
 real(dp) :: xn,xs,term,sum_k
 integer :: nn,kk

 xn=1.0d0
 do nn=1,ll
   xn=xx*xn/dble(2*nn+1)
 end do
 xs=0.5d0*xx**2
 term=1.0d0
 sum_k=1.0d0
 do kk=1,40
   term=-term*xs/(dble(kk)*dble(2*ll+2*kk+1))
   sum_k=sum_k+term
   if(abs(term)<=epsilon(1.0d0)*abs(sum_k)) exit
 end do
 jl=xn*sum_k
 end function sbf_series

 end module sbf_m
