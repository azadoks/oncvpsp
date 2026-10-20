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
 subroutine qroots(ll,rc,ulgd,nroot,dltq,qmax,qroot)

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

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: ll,nroot
 real(dp), intent(in) :: rc,ulgd,dltq,qmax

!Output variables
 real(dp), intent(out) :: qroot(nroot)

!Local variables
 real(dp), parameter :: eps=1.0d-12

 real(dp) :: sbfd(5)
 real(dp) :: al,qq,dlgd,dlgd_grid,dlgd_last,qhi,qlow,qt
 integer :: ii,jj,ll1,lmax,mmax,iroot,nq
 logical :: found_root

 nq=int(qmax/dltq)+1

 qroot(:)=0.0d0
 iroot=2
!the sign test reads dlgd_last at ii=1 too (.and. need not short-circuit)
 dlgd_last=0.0d0
 do ii=1,nq
   qq=dltq*ii
   call sbf_rc_der(ll,qq,rc,sbfd)
   dlgd=ulgd - sbfd(2)/sbfd(1)
!keep the grid value for the next sign test; bisection overwrites dlgd
   dlgd_grid=dlgd
   found_root=.false.
!interval halving to find root
   if(ii>1 .and. dlgd*dlgd_last < 0.0d0 .and. &
&          abs(dlgd*dlgd_last)<1.0d0) then
     if(dlgd>0.0d0) then
       qhi=qq
       qlow=dltq*(ii-1)
     else
       qhi=dltq*(ii-1)
       qlow=qq
     end if

     do jj=1,100
       qt=0.5d0*(qhi+qlow)
       call sbf_rc_der(ll,qt,rc,sbfd)
       dlgd=ulgd - sbfd(2)/sbfd(1)
       if(abs(dlgd)<eps) then
         iroot=iroot+1
         call polish(qt)
         qroot(iroot)=qt
         found_root=.true.
         exit
       end if
       if(dlgd>0.0d0) then
         qhi=qt
       else
         qlow=qt
       end if
     end do
     if(.not. found_root) then
       write(6,'(/a,i0,a,es12.4,a,es12.4)') &
&            'qroots: WARNING bisection failed, skipping', &
&            ', ll=',ll,', qlow=',qlow,', qhi=',qhi
     end if
   end if
   if(iroot == nroot) exit
   dlgd_last=dlgd_grid
 end do
 if(iroot<nroot) then
   write(6,'(a)') 'qroots: ERROR - failed to find nroot roots'
   error stop 1
 end if

!extra q values for needed flexibility to satisfy constraints
 qroot(1)=0.5d0*qroot(3)
 qt=0.5d0*(qroot(3)+qroot(4))
 qroot(2)=qroot(3)
 qroot(3)=qt

 return

 contains

 subroutine polish(qp)

!! the root qp polished to round-off by Newton's method on the pole-free
!! g(q) = ulgd*j_l(q rc) - q j_l'(q rc), whose root it is where j_l(q rc)
!! is not 0: the bisection stops anywhere within eps of it, the polished
!! root does not depend on where (doc/smoothness_plan.md, Q-5). Converged
!! when a step was within the round-off of g or a few ulps of q, g then at
!! the moved q
!!
!! @param qp  the root (in/out)

 real(dp), intent(inout) :: qp
 real(dp) :: gg,dg,dq,dgold,gscale
 integer :: kk

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
 error stop 1
 end subroutine polish

 end subroutine qroots
