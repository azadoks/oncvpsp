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
!> Spherical Bessel functions, their derivatives at a radius and their
!> overlaps inside a radius in quad precision: the data of the pseudo
!> wave's basis inside rc (doc/smoothness_plan.md, R-5). The basis is
!> nearly dependent, so a quantity of one of its orthonormal functions is
!> a sum of terms up to 1e4 times its size: in double precision such a sum
!> is good to 1e-12, in quad to 1e-30. sbf_overlap and sbf_rc_der are the
!> double-precision versions, for the sums with no such cancellation.

 module sbf_q_m
 use constants_m, only: qp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: sbfq_values, sbfq_value_slope, sbfq_overlap, sbfq_overlap_vs
 public :: sbfq_rc_der

!> the round-off of quad precision
 real(qp), parameter :: epsq=epsilon(1.0_qp)

 contains

 subroutine sbfq_values(nm,xx,sb)

!! j_l(xx) for l=0..nm-1: the power series below xx=1/4 (its terms fall
!! by (xx**2/2)/(k(2l+2k+1)) a term), else the downward recurrence of
!! sbf_values normalized by sum_l (2l+1) j_l**2 = 1, its sign that of
!! (j_0, j_1)
!!
!! @param nm  number of orders wanted
!! @param xx  argument, >=0
!! @param sb  j_l(xx) in sb(l+1)

!Input variables
 integer, intent(in) :: nm
 real(qp), intent(in) :: xx

!Output variables
 real(qp), intent(out) :: sb(nm)

!Local variables
 real(qp), allocatable :: tt(:)
 real(qp) :: xn,xs,term,sum_k,xi,scale,sn,j0,j1
 integer :: nn,kk,nlim

 if(xx==0.0_qp) then
   sb(:)=0.0_qp
   sb(1)=1.0_qp
   return
 end if

 if(xx<0.25_qp) then
   xn=1.0_qp
   xs=0.5_qp*xx**2
   do nn=1,nm
! j_l = xx**l/(2l+1)!! sum_k (-xs)**k/(k!(2l+3)(2l+5)...(2l+2k+1))
     term=1.0_qp
     sum_k=1.0_qp
     do kk=1,60
       term=-term*xs/(real(kk,qp)*real(2*nn-1+2*kk,qp))
       sum_k=sum_k+term
       if(abs(term)<=epsq*abs(sum_k)) exit
     end do
     sb(nn)=xn*sum_k
     xn=xx*xn/real(2*nn+1,qp)
   end do
   return
 end if

! downward recurrence from an order where j_l is negligible (the
! recurrence is stable downward; the start's error decays by about
! (xx/2n)**2 an order, so the margin of 50 orders is ample for 1e-34),
! normalized by sum_l (2l+1) j_l**2 = 1: j_0 = sin(xx)/xx alone lost its
! relative accuracy near its zeros (numerics audit C-13). The sign is
! that of the projection on (j_0, j_1), which are never both small
 nlim=nm+int(1.4_qp*xx)+50
 allocate(tt(nlim+2))
 xi=1.0_qp/xx
 tt(nlim+2)=0.0_qp
 tt(nlim+1)=1.0e-30_qp
 do nn=nlim,1,-1
   tt(nn)=real(2*nn+1,qp)*xi*tt(nn+1) - tt(nn+2)
   if(abs(tt(nn))>1.0e300_qp) then
     tt(nn:nlim+2)=tt(nn:nlim+2)*1.0e-300_qp
   end if
 end do
 sn=0.0_qp
 do nn=1,nlim+1
   sn=sn + real(2*nn-1,qp)*tt(nn)**2
 end do
 scale=1.0_qp/sqrt(sn)
 j0=sin(xx)*xi
 j1=(j0-cos(xx))*xi
 if(tt(1)*j0+tt(2)*j1<0.0_qp) scale=-scale
 do nn=1,nm
   sb(nn)=scale*tt(nn)
 end do
 deallocate(tt)
 return
 end subroutine sbfq_values

 subroutine sbfq_value_slope(lq,xx,jv,jd)

!! j_l(x) and j_l'(x)=(l j_l-1 - (l+1) j_l+1)/(2l+1) (-j_1 for l=0)
!!
!! @param lq  angular momentum l
!! @param xx  argument, >0
!! @param jv  j_l(x)
!! @param jd  j_l'(x)

!Input variables
 integer, intent(in) :: lq
 real(qp), intent(in) :: xx

!Output variables
 real(qp), intent(out) :: jv,jd

!Local variables
 real(qp) :: sb(lq+2)

 call sbfq_values(lq+2,xx,sb)
 jv=sb(lq+1)
 if(lq==0) then
   jd=-sb(2)
 else
   jd=(real(lq,qp)*sb(lq) - real(lq+1,qp)*sb(lq+2))/real(2*lq+1,qp)
 end if
 return
 end subroutine sbfq_value_slope

 subroutine sbfq_overlap(ll,q1,q2,rc,ss)

!! overlap of two spherical Bessel functions of one l inside a radius, in
!! closed form (Lommel's integrals), so exact and smooth in q1, q2 and rc:
!! int_0^rc r**2 j_l(q1 r) j_l(q2 r) dr. Where Lommel's form cancels, a
!! series takes over: the power series of the integral when both q rc are
!! small, the Taylor series in q2-q1 when q1 is near q2 (and q1=q2); the
!! method of sbf_overlap, with the terms that quad precision needs. With 34
!! digits, Lommel's form is kept until its cancellation has cost 13 of
!! them (q2-q1 above 1e-13 q1), so the series is for q1=q2 in practice
!!
!! @param ll  angular momentum l
!! @param q1,q2  wave vectors (bohr^-1), >0
!! @param rc  core radius (bohr), >0
!! @param ss  the overlap

!Input variables
 integer, intent(in) :: ll
 real(qp), intent(in) :: q1,q2,rc

!Output variables
 real(qp), intent(out) :: ss

!Local variables
 real(qp), allocatable :: cc(:)
 real(qp) :: alpha,beta,de,xa,xb,tt,ja,jda,jb,jdb,term,sum_n,tpow,xl
 integer :: kk,nn
 integer, parameter :: nterm=160,nser=48

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
 if(xb<=1.5_qp) then
   allocate(cc(0:nser))
   cc(0)=1.0_qp
   do kk=1,2*ll+1,2
     cc(0)=cc(0)/real(kk,qp)
   end do
   do kk=1,nser
     cc(kk)=-cc(kk-1)/(2.0_qp*real(kk,qp)*real(2*ll+2*kk+1,qp))
   end do
   ss=0.0_qp
   do kk=nser,0,-1
     do nn=nser,0,-1
       ss=ss + cc(kk)*cc(nn)*xa**(2*kk)*xb**(2*nn)/real(2*ll+3+2*kk+2*nn,qp)
     end do
   end do
   ss=ss*rc**3*(xa*xb)**ll
   deallocate(cc)
   return
 end if

 call sbfq_value_slope(ll,xa,ja,jda)

! Lommel: rc**2 [beta j(xa) j'(xb) - alpha j'(xa) j(xb)]/(alpha**2-beta**2),
! where the numerator cancels to O(de) relative; used while the
! cancellation costs at most 13 digits
 if(tt>1.0e-13_qp*xa) then
   call sbfq_value_slope(ll,xb,jb,jdb)
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
 xl=real(ll*(ll+1),qp)
 cc(-2:-1)=0.0_qp
 cc(0)=ja
 cc(1)=jda
 do kk=0,nterm-1
   cc(kk+2)=-(2.0_qp*xa*real((kk+1)**2,qp)*cc(kk+1) &
&             + (real(kk*(kk+1),qp) + xa**2 - xl)*cc(kk) &
&             + 2.0_qp*xa*cc(kk-1) + cc(kk-2)) &
&           /(xa**2*real((kk+1)*(kk+2),qp))
 end do
 sum_n=cc(0)*cc(1)
 tpow=1.0_qp
 do nn=1,nterm
   term=rc*tpow*(alpha*(real(nn+1,qp)*cc(0)*cc(nn+1) - cc(1)*cc(nn)) &
&               + de*real(nn+1,qp)*cc(0)*cc(nn+1))
   sum_n=sum_n + term
   if(abs(term)<=0.25_qp*epsq*abs(sum_n) .and. nn>2) exit
   tpow=tpow*tt
 end do
 deallocate(cc)
 ss=-rc**2*sum_n/(2.0_qp*alpha + de)
 return
 end subroutine sbfq_overlap

 subroutine sbfq_overlap_vs(ll,q1,q2,rc,j1,jd1,j2,jd2,ss)

!! sbfq_overlap given j_l and j_l' at q1 rc and q2 rc, for the many
!! overlaps of eresid that share them: Lommel's form where it applies,
!! else sbfq_overlap itself
!!
!! @param ll  angular momentum l
!! @param q1,q2  wave vectors (bohr^-1), >0
!! @param rc  core radius (bohr), >0
!! @param j1,jd1  j_l(q1 rc) and j_l'(q1 rc)
!! @param j2,jd2  the same at q2 rc
!! @param ss  the overlap

!Input variables
 integer, intent(in) :: ll
 real(qp), intent(in) :: q1,q2,rc,j1,jd1,j2,jd2

!Output variables
 real(qp), intent(out) :: ss

 if(max(q1,q2)*rc>1.5_qp .and. abs(q2-q1)>1.0e-13_qp*min(q1,q2)) then
   ss=(rc*rc)*(q2*j1*jd2 - q1*jd1*j2)/((q1-q2)*(q1+q2))
 else
   call sbfq_overlap(ll,q1,q2,rc,ss)
 end if
 return
 end subroutine sbfq_overlap_vs

 subroutine sbfq_rc_der(llin,qq,rc,sbfder,ifail)

!! spherical Bessel function values and first four derivatives in r of
!! j_l(qq r) at r=rc, by the derivative formula applied recursively,
!! (1/z d/dz)^m [z^-n j_n(z)] = (-1)^m z^(-n-m) j_(n+m)(z), as sbf_rc_der
!!
!! @param llin  angular momentum l, <=4
!! @param qq  wave vector
!! @param rc  core radius (bohr)
!! @param sbfder  values and first 4 derivatives of j_l(qq*rc)
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_internal

!Input variables
 integer, intent(in) :: llin
 real(qp), intent(in) :: qq,rc

!Output variables
 real(qp), intent(out) :: sbfder(5)
 integer, intent(out) :: ifail

!Local variables
 real(qp) :: xx
 real(qp) :: sb_out(10),sbfad(5,10)
 integer :: ii,ll

 ifail=0
 if(llin>4) then
   write(error_unit,'(a,i4,a)') &
&       'sbfq_rc_der: argument ERROR, llin = ',llin,' >4'
   sbfder(:)=0.0_qp
   ifail=fail_internal
   return
 end if

 xx=qq*rc
 call sbfq_values(llin+6,xx,sb_out)

 do ll=llin,llin+5
   sbfad(1,ll+1)=sb_out(ll+1)
 end do
 do ll=llin,llin+3
   sbfad(2,ll+1)=ll*sbfad(1,ll+1)/xx - sbfad(1,ll+2)
 end do
 do ll=llin,llin+2
   sbfad(3,ll+1)=-ll*sbfad(1,ll+1)/xx**2 + ll*sbfad(2,ll+1)/xx &
&                -sbfad(2,ll+2)
 end do
 do ll=llin,llin+1
   sbfad(4,ll+1)= 2*ll*sbfad(1,ll+1)/xx**3 - 2*ll*sbfad(2,ll+1)/xx**2 &
&                + ll*sbfad(3,ll+1)/xx - sbfad(3,ll+2)
 end do
 do ll=llin,llin
   sbfad(5,ll+1)=-6*ll*sbfad(1,ll+1)/xx**4 + 6*ll*sbfad(2,ll+1)/xx**3 &
&                -3*ll*sbfad(3,ll+1)/xx**2 + ll*sbfad(4,ll+1)/xx &
&                -sbfad(4,ll+2)
 end do

! derivatives with respect to r
 do ii=1,5
   sbfder(ii)=sbfad(ii,llin+1)*qq**(ii-1)
 end do
 return
 end subroutine sbfq_rc_der

 end module sbf_q_m
