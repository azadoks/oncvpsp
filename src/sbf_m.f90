!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
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
!> Every order from sbf_values (the downward recurrence beyond l=3) and the
!> derivatives at a radius from sbf_rc_der: the double-precision j_l of
!> the programs; sbf_q_m holds those in quad precision.

 module sbf_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: sbf_trig, sbf_trig_rising, sbf_values, sbf_rc_der

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

 subroutine sbf_values(nm,xx,sb_out)

!! the spherical Bessel functions j_l(xx) for l=0,...,nm-1: a series at
!! small xx, the closed forms for nm<=4, else downward recursion
!!
!! @param nm  maximum angular momentum wanted plus 1
!! @param xx  argument of spherical bessel function
!! @param sb_out  output of sbf_l(xx) for l=0,...,nm-1

!Input variables
 integer, intent(in) :: nm
 real(dp), intent(in) :: xx

!Output variables
 real(dp), intent(out) :: sb_out(nm)

!Local variables
 integer :: nlim,nn
 real(dp) :: cc,fn,sn,ss,xi,xn,xs
 real(dp),allocatable :: sb(:)


 if(xx<= 1.0d-36) then
!  zero argument section
   sb_out(:)=0.0d0
   sb_out(1)=1.0d0
 else if(xx<1.0d-3) then
!  small argument section
   xn=1.0d0
   xs=0.5d0*xx**2
   do nn=1,nm
     sb_out(nn)=xn*(1.0d0 - xs*(1.0d0 - xs/(4*nn+6))/(2*nn+1))
     xn=xx*xn/(2*nn+1)
   end do
!  trigonometric section (small l values), every order from sin and cos
!  (sbf_trig: the series where the closed forms cancel)
 else if(nm<=4) then
   ss=sin(xx)
   cc=cos(xx)
   do nn=1,nm
     sb_out(nn)=sbf_trig(nn-1,xx,ss,cc)
   end do
 else
!  recursion method (accurate for large l, slow for large arguments),
!  started at order nm+1.36 xx+15 at every xx (the start's error is then
!  below round-off; at small xx a start at nm+1 had cost the top orders
!  1e-9 and jumped where it moved), normalized by
!  sum_l (2l+1) j_l**2 = 1. Its round-off grows with the xx oscillatory
!  orders it passes, about 0.6 eps xx
   nlim=nm+int(1.36d0*xx)+15
   allocate(sb(nlim+1))
   xi=1.0d0/xx
   sb(nlim+1)=0.0d0
   sb(nlim)=1.0d-18
   do nn=nlim-1,1,-1
     sb(nn)=dble(2*nn+1)*xi*sb(nn+1) - sb(nn+2)
     if(abs(sb(nn))>1.0d150) then
       sb(nn:nlim+1)=sb(nn:nlim+1)*1.0d-150
     end if
   end do
   sn=0.0d0
   do nn=1,nlim
     sn=sn + dble(2*nn-1)*sb(nn)*sb(nn)
   end do
   fn=1.0d0/sqrt(sn)
   sb_out(:)=fn*sb(1:nm)
   deallocate(sb)
 end if
 return
 end subroutine sbf_values

 subroutine sbf_rc_der(llin,qq,rc,sbfder)

!! calculate spherical Bessel function values and first four derivatives
!!
!! @param llin  angular momentum l
!! @param qq  wave vector
!! @param rc  core radius (bohr)
!! @param sbfder  values and first 4 derivatives of j_l(qq*rc)

!Input variables
 integer, intent(in) :: llin
 real(dp), intent(in) :: qq,rc

!Output variables
 real(dp), intent(out) :: sbfder(5)

!Local variables
 real(dp) :: xx
 real(dp) :: sb_out(10),sbfad(5,10)
 integer :: ii,ll

 if(llin>4) then
   write(error_unit,'(a,i4,a)') 'sbf_rc_der: argument ERROR, llin = ',llin,&
&   ' >5'
   error stop 1
 end if

! calculate sbf for needed l values
 xx=qq*rc
 call sbf_values(llin+5,xx,sb_out)

! derivatives based on derivative formula applied recursively
!
! (1/z d/dz)^m [z^-n j_n(z)]= (-1)^m z^(-n-m) j_(n+m)(z)
!

! store values (orders llin..llin+4, those the derivatives use)
 do ll=llin,llin+4
   sbfad(1,ll+1)=sb_out(ll+1)
 end do

! first derivatives
 do ll=llin,llin+3
   sbfad(2,ll+1)=ll*sbfad(1,ll+1)/xx &
&                 - sbfad(1,ll+2)
 end do

! second derivatives
 do ll=llin,llin+2
   sbfad(3,ll+1)=-ll*sbfad(1,ll+1)/xx**2 &
&               + ll*sbfad(2,ll+1)/xx&
&                   -sbfad(2,ll+2)
 end do

! third derivatives
 do ll=llin,llin+1
   sbfad(4,ll+1)= 2*ll*sbfad(1,ll+1)/xx**3 &
&                -2*ll*sbfad(2,ll+1)/xx**2 &
&                  +ll*sbfad(3,ll+1)/xx &
&                     -sbfad(3,ll+2)
 end do

! fourth derivatives
 do ll=llin,llin
   sbfad(5,ll+1)=-6*ll*sbfad(1,ll+1)/xx**4 &
&                +6*ll*sbfad(2,ll+1)/xx**3 &
&                -3*ll*sbfad(3,ll+1)/xx**2 &
&                  +ll*sbfad(4,ll+1)/xx &
&                     -sbfad(4,ll+2)
 end do

 do ii=1,5
   sbfder(ii)=sbfad(ii,llin+1)*qq**(ii-1)
 end do

 return
 end subroutine sbf_rc_der

 end module sbf_m
