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

 module sbf_overlap_m

 implicit none
 private
 public :: sbf_overlap

 contains

 subroutine sbf_overlap(ll,q1,q2,rc,ss)

!! overlap of two spherical Bessel functions of one l inside a radius, in
!! closed form (Lommel's integrals), so exact and smooth in q1, q2 and rc:
!! int_0^rc r**2 j_l(q1 r) j_l(q2 r) dr. Where Lommel's form cancels, a
!! series takes over: the power series of the integral when both q rc are
!! small, the Taylor series in q2-q1 when q1 is near q2 (and q1=q2) (see
!! doc/smoothness_plan.md)
!!
!! @param ll  angular momentum l
!! @param q1,q2  wave vectors (bohr^-1), >0
!! @param rc  radius (bohr), >0
!! @param ss  the overlap

 use constants_m, only: dp
 use sbf8_m, only: sbf8
 implicit none

!Input variables
 integer, intent(in) :: ll
 real(dp), intent(in) :: q1,q2,rc

!Output variables
 real(dp), intent(out) :: ss

!Local variables
 real(dp), allocatable :: cc(:)
 real(dp) :: alpha,beta,de,xa,xb,tt,ja,jda,jb,jdb,term,sum_n,tpow,xl
 integer :: kk,nn
 integer, parameter :: nterm=60,nser=20

! alpha<=beta=alpha+de; x=alpha rc, and t=de rc the distance to beta rc
 alpha=min(q1,q2)
 beta=max(q1,q2)
 de=beta-alpha
 xa=alpha*rc
 xb=beta*rc
 tt=de*rc

! both small: the leading powers of the two Lommel terms cancel exactly, so
! integrate the power series, j_l(z) = z**l sum_k (-z**2/2)**k/(k!(2l+2k+1)!!),
! term by term with a fixed number of terms (a polynomial, so smooth)
 if(xb<=1.5d0) then
   allocate(cc(0:nser))
   cc(0)=1.0d0
   do kk=1,2*ll+1,2
     cc(0)=cc(0)/dble(kk)
   end do
   do kk=1,nser
     cc(kk)=-cc(kk-1)/(2.0d0*dble(kk)*dble(2*ll+2*kk+1))
   end do
   ss=0.0d0
   do kk=nser,0,-1
     do nn=nser,0,-1
       ss=ss + cc(kk)*cc(nn)*xa**(2*kk)*xb**(2*nn)/dble(2*ll+3+2*kk+2*nn)
     end do
   end do
   ss=ss*rc**3*(xa*xb)**ll
   deallocate(cc)
   return
 end if

 call sbf_value_slope(ll,xa,ja,jda)

! Lommel: rc**2 [beta j(xa) j'(xb) - alpha j'(xa) j(xb)]/(alpha**2-beta**2),
! where the numerator cancels to O(de) relative; used while the
! cancellation costs at most a few digits
 if(tt>min(1.0d0,0.25d0*xa)) then
   call sbf_value_slope(ll,xb,jb,jdb)
   ss=rc**2*(beta*ja*jdb - alpha*jda*jb)/((alpha-beta)*(alpha+beta))
   return
 end if

! Otherwise expand j(xa+t) = sum_n c_n t**n about xa, with the Taylor
! coefficients from the Bessel equation
!   x**2 j'' + 2x j' + (x**2 - l(l+1)) j = 0
! by the recurrence
!   c_k+2 = -[2x(k+1)**2 c_k+1 + (k(k+1) + x**2 - l(l+1)) c_k
!            + 2x c_k-1 + c_k-2] / (x**2 (k+1)(k+2)).
! The numerator divided by de is then
!   c_0 c_1 + sum_n>=1 rc t**(n-1) [alpha((n+1) c_0 c_n+1 - c_1 c_n)
!                                   + de (n+1) c_0 c_n+1],
! with no cancellation (t<=x/4 keeps the recurrence's error small), and
! the overlap is -rc**2 (numerator/de)/(2 alpha + de).
 allocate(cc(-2:nterm+1))
 xl=dble(ll*(ll+1))
 cc(-2:-1)=0.0d0
 cc(0)=ja
 cc(1)=jda
 do kk=0,nterm-1
   cc(kk+2)=-(2.0d0*xa*dble((kk+1)**2)*cc(kk+1) &
&             + (dble(kk*(kk+1)) + xa**2 - xl)*cc(kk) &
&             + 2.0d0*xa*cc(kk-1) + cc(kk-2)) &
&           /(xa**2*dble((kk+1)*(kk+2)))
 end do
 sum_n=cc(0)*cc(1)
 tpow=1.0d0
 do nn=1,nterm
   term=rc*tpow*(alpha*(dble(nn+1)*cc(0)*cc(nn+1) - cc(1)*cc(nn)) &
&               + de*dble(nn+1)*cc(0)*cc(nn+1))
   sum_n=sum_n + term
   if(abs(term)<=0.25d0*epsilon(1.0d0)*abs(sum_n) .and. nn>2) exit
   tpow=tpow*tt
 end do
 deallocate(cc)
 ss=-rc**2*sum_n/(2.0d0*alpha + de)
 return

 contains

 subroutine sbf_value_slope(lq,xx,jv,jd)

!! j_l(x) and j_l'(x)=(l j_l-1 - (l+1) j_l+1)/(2l+1) (-j_1 for l=0), from
!! sbf8's downward recurrence
!!
!! @param lq  angular momentum l
!! @param xx  argument, >0
!! @param jv  j_l(x)
!! @param jd  j_l'(x)

!Input variables
 integer, intent(in) :: lq
 real(dp), intent(in) :: xx

!Output variables
 real(dp), intent(out) :: jv,jd

!Local variables
 real(dp) :: sb(max(lq+2,5))

! sbf8 is complete for every order only with 5 or more
 call sbf8(max(lq+2,5),xx,sb)
 jv=sb(lq+1)
 if(lq==0) then
   jd=-sb(2)
 else
   jd=(dble(lq)*sb(lq) - dble(lq+1)*sb(lq+2))/dble(2*lq+1)
 end if
 return
 end subroutine sbf_value_slope

 end subroutine sbf_overlap

 end module sbf_overlap_m
