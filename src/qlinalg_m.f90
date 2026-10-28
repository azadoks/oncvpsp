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
!> Quad-precision linear algebra for the small, ill-conditioned problems
!> of the pseudo wave inside rc (doc/smoothness_plan.md, R-5): the basis
!> of spherical Bessel functions is nearly dependent (its overlap matrix
!> has eigenvalues down to 5e-9), so its orthonormalization and the
!> matching constraints solved in double precision carry round-off at
!> eps times its condition number into the pseudo wave's derivatives.

 module qlinalg_m
 implicit none
 private
 public :: qjacobi

 contains

 subroutine qjacobi(nn,aa,ww,evec)

!! eigenvalues and eigenvectors of a symmetric matrix by cyclic Jacobi
!! rotations, to the round-off of quad precision
!!
!! @param nn  order
!! @param aa  the matrix (its upper triangle is read)
!! @param ww  eigenvalues, ascending
!! @param evec  eigenvectors in the columns, orthonormal

 use constants_m, only: qp

!Input variables
 integer, intent(in) :: nn
 real(qp), intent(in) :: aa(nn,nn)

!Output variables
 real(qp), intent(out) :: ww(nn),evec(nn,nn)

!Local variables
 real(qp) :: bb(nn,nn),off,theta,tt,cc,ss,tau,gg,hh,apq,tmp
 integer :: pp,qq,kk,sweep,isort
 integer, parameter :: nsweep=60

 bb(:,:)=0.0_qp
 do pp=1,nn
   do qq=pp,nn
     bb(pp,qq)=aa(pp,qq)
     bb(qq,pp)=aa(pp,qq)
   end do
 end do
 evec(:,:)=0.0_qp
 do pp=1,nn
   evec(pp,pp)=1.0_qp
 end do

 do sweep=1,nsweep
   off=0.0_qp
   do pp=1,nn-1
     do qq=pp+1,nn
       off=off+bb(pp,qq)**2
     end do
   end do
   if(off==0.0_qp) exit
   do pp=1,nn-1
     do qq=pp+1,nn
       apq=bb(pp,qq)
! an element below eps**2 of its diagonal moves nothing: zeroed rather
! than rotated on into underflow, where theta**2 overflows
       if(abs(apq)<=epsilon(1.0_qp)**2*(abs(bb(pp,pp))+abs(bb(qq,qq)))) &
&        then
         bb(pp,qq)=0.0_qp
         bb(qq,pp)=0.0_qp
         cycle
       end if
! the rotation that annihilates bb(pp,qq) (Rutishauser's formulas)
       theta=0.5_qp*(bb(qq,qq)-bb(pp,pp))/apq
       tt=sign(1.0_qp,theta)/(abs(theta)+sqrt(1.0_qp+theta**2))
       cc=1.0_qp/sqrt(1.0_qp+tt**2)
       ss=tt*cc
       tau=ss/(1.0_qp+cc)
       bb(pp,pp)=bb(pp,pp)-tt*apq
       bb(qq,qq)=bb(qq,qq)+tt*apq
       bb(pp,qq)=0.0_qp
       bb(qq,pp)=0.0_qp
       do kk=1,nn
         if(kk==pp .or. kk==qq) cycle
         gg=bb(kk,pp)
         hh=bb(kk,qq)
         bb(kk,pp)=gg-ss*(hh+gg*tau)
         bb(kk,qq)=hh+ss*(gg-hh*tau)
         bb(pp,kk)=bb(kk,pp)
         bb(qq,kk)=bb(kk,qq)
       end do
       do kk=1,nn
         gg=evec(kk,pp)
         hh=evec(kk,qq)
         evec(kk,pp)=gg-ss*(hh+gg*tau)
         evec(kk,qq)=hh+ss*(gg-hh*tau)
       end do
     end do
   end do
 end do

 do pp=1,nn
   ww(pp)=bb(pp,pp)
 end do

! sort ascending
 do pp=1,nn-1
   isort=pp
   do qq=pp+1,nn
     if(ww(qq)<ww(isort)) isort=qq
   end do
   if(isort/=pp) then
     tmp=ww(pp)
     ww(pp)=ww(isort)
     ww(isort)=tmp
     do kk=1,nn
       tmp=evec(kk,pp)
       evec(kk,pp)=evec(kk,isort)
       evec(kk,isort)=tmp
     end do
   end if
 end do
 return
 end subroutine qjacobi

 end module qlinalg_m
