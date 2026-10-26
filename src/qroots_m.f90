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

 module qroots_m

 implicit none
 private
 public :: qroots

 contains

 subroutine qroots(ll,rc,ulgd,nroot,dltq,qmax,qroot,ifail)

!! calculate q values for spherical Bessel function orthogonal basis inside r_c
!! all having specified log derivative ulgd at rc
!! 1st and 3rd output roots do not correspond to these conditions, rather
!! 0.5 times the first matching root, and the average of the 1st and second
!! matching roots
!!
!! @param ll  angular momentum
!! @param rc  core radius
!! @param ulgd  log derivaive of all-electron radial wave function at rc
!! @param nroot  number of spherical Bessel functions matgching ulgd at rc (form
!!               orthogonal basis)
!! @param dltq  step for search for qroots satisfying log derivative condition
!! @param qmax  maximum q for search
!! @param qroot  set of q's satisfying dj_l(q*r)/dr / jl(q*r) = ulgd
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use sbf_rc_der_m, only: sbf_rc_der
 implicit none

!Input variables
 integer, intent(in) :: ll,nroot
 real(dp), intent(in) :: rc,ulgd,dltq,qmax

!Output variables
 real(dp), intent(out) :: qroot(nroot)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: sbfd(5)
 real(dp) :: qq,gg,glast,glow,gt,qhi,qlow,qt
 integer :: ii,jj,iroot,nq

 ifail=0
 nq=int(qmax/dltq)+1

! the roots are those of the pole-free g(q) = ulgd*j_l(q rc) - q j_l'(q rc),
! which are those of the log-derivative mismatch where j_l(q rc) is not 0
! and g is never 0 where it is: every sign change of g on the grid brackets
! a root, so no pole has to be told from a root. The mismatch itself, with a
! pole beside each root for a large ulgd, made the grid's sign changes of
! roots look like poles: for a Pb l=3 well state (ulgd = 65) every root
! was skipped and the search ran out
 qroot(:)=0.0d0
 iroot=2
 glast=gfun(dltq)
 do ii=2,nq
   qq=dltq*ii
   gg=gfun(qq)
   if(gg*glast<0.0d0 .or. gg==0.0d0) then
! bisection to a bracket of a few ulps (or g = 0), then the Newton polish
     qlow=qq-dltq
     qhi=qq
     glow=glast
     do jj=1,200
       qt=0.5d0*(qlow+qhi)
       if(qhi-qlow<=4.0d0*spacing(qt)) exit
       gt=gfun(qt)
       if(gt==0.0d0) exit
       if(gt*glow<0.0d0) then
         qhi=qt
       else
         qlow=qt
         glow=gt
       end if
     end do
     if(gg==0.0d0) qt=qq
     call polish(qt,ifail)
     if(ifail/=0) return
     iroot=iroot+1
     qroot(iroot)=qt
     if(iroot==nroot) exit
! a root on the grid point itself is not counted again at the next one
     if(gg==0.0d0) then
       glast=gfun(qq+0.5d0*dltq)
       cycle
     end if
   end if
   glast=gg
 end do
 if(iroot<nroot) then
   write(6,'(a)') 'qroots: ERROR - failed to find nroot roots'
   ifail=1
   return
 end if

!extra q values for needed flexibility to satisfy constraints
 qroot(1)=0.5d0*qroot(3)
 qt=0.5d0*(qroot(3)+qroot(4))
 qroot(2)=qroot(3)
 qroot(3)=qt

 return

 contains

 real(dp) function gfun(qp)
!! g(q) = ulgd*j_l(q rc) - q j_l'(q rc), whose roots polish refines
 real(dp), intent(in) :: qp
 call sbf_rc_der(ll,qp,rc,sbfd)
 gfun=ulgd*sbfd(1) - sbfd(2)
 end function gfun

 subroutine polish(qp,ifail)

!! the root qp polished to round-off by Newton's method on the pole-free
!! g(q) = ulgd*j_l(q rc) - q j_l'(q rc), whose root it is where j_l(q rc)
!! is not 0: the bisection stops anywhere within a few ulps of it, the
!! polished root does not depend on where (doc/smoothness_plan.md, Q-5).
!! Converged when a step was within the round-off of g or a few ulps of q,
!! g then at the moved q
!!
!! @param qp  the root (in/out)
!! @param ifail  0, or 1 after a failure, whose message is written

 real(dp), intent(inout) :: qp
 real(dp) :: gg,dg,dq,dgold,gscale
 integer :: kk
 integer, intent(out) :: ifail

 ifail=0
 dq=0.0d0
 dgold=0.0d0
 do kk=1,6
   call sbf_rc_der(ll,qp,rc,sbfd)
! sbfd(k+1) = d^k/dr^k j_l(q r) at rc, so d/dq of sbfd(k+1) is
! (k sbfd(k+1) + rc sbfd(k+2))/q
   gg=ulgd*sbfd(1) - sbfd(2)
   dg=(ulgd*rc*sbfd(2) - sbfd(2) - rc*sbfd(3))/qp
   gscale=abs(ulgd*sbfd(1)) + abs(sbfd(2))
   if(kk>1 .and. (abs(dq)*dgold<=8.0d0*epsilon(1.0d0)*gscale .or. &
&                 abs(dq)<=4.0d0*spacing(qp))) return
   dq=-gg/dg
   dgold=abs(dg)
   qp=qp+dq
 end do
 write(6,'(/a,i0,a,es12.4,a,es10.2)') &
&     'qroots: ERROR Newton polish not converged, ll=',ll,', q=',qp, &
&     ', last step',dq
 ifail=1
 return
 end subroutine polish

 end subroutine qroots

 end module qroots_m
