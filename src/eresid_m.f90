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

 module eresid_m

 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: eresid

 contains

! calculates residual energy operator matrix elements as a function of q

 subroutine eresid(ll,rc,nnull,nbas,nq,dq,nqout,qout,qmax,qroot,hh,nr, &
&                  ulin,pswf0_sb,pswfnull_sb,e0g,edotg,ematg, &
&                  eresid0,eresiddot,eresidmat,iout,ifail)

!! residual kinetic energy matrix elements as functions of the lower
!! cutoff q: E0(q) = (1/pi) int_q^{qmax} q'**4 |ps0(q')|**2 dq', and the
!! same with the null-space functions in place of one or both ps0; on the
!! uniform q grid q_i = qmax - (nq-i) dq, i=1..nq, from below dq up to
!! qmax (e0g, edotg, ematg), and at the wave vectors qout (eresid0,
!! eresiddot, eresidmat). The grid's spacing is a continuous function of
!! the inputs and its nodes move continuously with qmax; a node that
!! appears enters at q = 0, where the integrands vanish as q**4 (the
!! spacing qmax/(nq-1) with an integer nq jumped where nq did;
!! doc/numerics_audit.md, C-10). The
!! transform ps(q) = int_0^inf (r psi(r)) (r j_l(q r)) dr is the exact
!! transform of the Bessel sum inside rc (sbf_overlap, Lommel's integrals)
!! plus that of the all-electron tail outside rc on a linear mesh from rc,
!! of spacing hh = 2 pi/(75 qmax) (75 points per period of qmax): the
!! nodes depend on rc and qmax only, so that the tail on them (free_tail_m:
!! the wave function interpolated, 1e-12 of it far out, or the windowed
!! free tail of a well state) is a smooth function of the inputs;
!! integrated by the 8th-order rule of ugrid_m. For l<=3, r j_l(q r) is
!! a sum of q**-k r**(1-k) sin or cos(q r), k=1..l+1, so the tail's
!! transforms on the grid are those sums over the mesh, at every node at
!! once by the chirp-z transform (czt_m: O(L log L), L>=nr+nq, for the
!! O(nr nq) of the sums; their round-off is 1e-16 of max|ps|, that of the
!! sums of sines being some 4e-16); at the lowest nodes, where those
!! terms cancel (their sizes q**-k sum|u r**(1-k)| above a quarter of
!! sum|u r|), and for l>3, they are the sums of rjl, as at the partial
!! cells' points below. The
!! integrals over q accumulate from qmax by the same rule, from the small
!! values up, so their round-off is that of their own size (Parseval,
!! E_resid as the kinetic energy less the part below q, would leave the
!! round-off of the kinetic energy, 1e-11 of an E_resid of 1e-5 Ha); at a
!! q between the nodes, E is the value at the next node plus the integral
!! over the partial cell by 8-point Gauss-Legendre, the transforms at its
!! nodes evaluated directly. Everything is smooth in rc, qmax, qout and
!! the functions (doc/smoothness_plan.md, Q-1 to Q-3)
!!
!! @param ll  angular momentum
!! @param rc  core radius (bohr)
!! @param nnull  number of null-space functions
!! @param nbas  number of basis functions
!! @param nq  number of nodes of the q grid
!! @param dq  its spacing
!! @param nqout  number of wave vectors at which E is wanted
!! @param qout  those wave vectors, 0<=qout<qmax (bohr^-1)
!! @param qmax  largest q, taken as infinity (bohr^-1)
!! @param qroot  wave vectors of the basis functions
!! @param hh  spacing of the linear mesh of the tail, 2 pi/(75 qmax)
!! @param nr  number of its nodes, at least 2 lm_ncell
!! @param ulin  the all-electron tail (*r) at rc + (i-1) hh
!! @param pswf0_sb  basis coefficients of the constrained pseudo wave (quad)
!! @param pswfnull_sb  basis coefficients of the null-space functions (quad):
!!                     the basis is nearly dependent, so their transforms
!!                     are sums with 1e4 cancellation, formed in quad
!! @param e0g  <ps0| E_resid |ps0> at each node of the grid
!! @param edotg  <psnull| E_resid |ps0> at each node
!! @param ematg  <psnull| E_resid |psnull'> at each node
!! @param eresid0  <ps0| E_resid |ps0> at each qout
!! @param eresiddot  <psnull| E_resid |ps0> at each qout
!! @param eresidmat  <psnull| E_resid |psnull'> at each qout
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, qp, pi
 use logmesh_m, only: lm_ncell
 use ugrid_m, only: ug_weights, ug_cumint_down
 use sbf_q_m, only: sbfq_overlap_vs, sbfq_value_slope
 use sbf_m, only: sbf_trig_rising, sbf_values
 use czt_m, only: czt_plan_t, czt_plan, czt_apply, czt_free

!Input variables
 integer, intent(in) :: ll,nnull,nbas,nq,nqout,nr
 real(dp), intent(in) :: rc,qmax,hh
 real(dp), intent(in) :: ulin(nr),qout(nqout),dq
 real(qp), intent(in) :: pswf0_sb(nbas),pswfnull_sb(nbas,nnull)
 real(dp), intent(in) :: qroot(nbas)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: e0g(nq),edotg(nq,nnull),ematg(nq,nnull,nnull)
 real(dp), intent(out) :: eresid0(nqout),eresiddot(nnull,nqout)
 real(dp), intent(out) :: eresidmat(nnull,nnull,nqout)
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: rlin(:),wr(:),sbf(:)
 real(dp), allocatable :: sqr(:),cqr(:),sdqr(:),cdqr(:)
 real(dp), allocatable :: g0(:),gdot(:,:),gmat(:,:,:),work(:)
 real(dp), allocatable :: aw(:),tailg(:)
 integer :: nqd
 real(dp) :: sb_out(5)
 real(qp) :: ovl(nbas),ps0in,psnin(nnull),rcq,jq(nbas),jdq(nbas),jb,jdb
 real(dp) :: q0
 integer :: ii,jj,ll1
!steps between exact sin/cos(qq*rlin) evaluations in the q loop
 integer, parameter :: nqseed=32
!fewest nodes for the chirp-z transform
 integer, parameter :: nczt=16
!|c_k| of the terms q**-k r**(1-k) sin or cos(q r) of r j_l(q r), l=0..3
 real(dp), parameter :: ccoef(4,0:3)=reshape([1.0d0,0.0d0,0.0d0,0.0d0, &
&  1.0d0,1.0d0,0.0d0,0.0d0, 1.0d0,3.0d0,3.0d0,0.0d0, &
&  1.0d0,6.0d0,15.0d0,15.0d0],[4,4])
!8-point Gauss-Legendre on [-1,1], for the partial cells
 integer, parameter :: ngl=8
 real(dp), parameter :: xgl(ngl)=[-0.9602898564975363d0, &
&  -0.7966664774136267d0,-0.5255324099163290d0,-0.1834346424956498d0, &
&   0.1834346424956498d0, 0.5255324099163290d0, 0.7966664774136267d0, &
&   0.9602898564975363d0]
 real(dp), parameter :: wgl(ngl)=[0.1012285362903763d0, &
&   0.2223810344533745d0,0.3137066458778873d0,0.3626837833783620d0, &
&   0.3626837833783620d0,0.3137066458778873d0,0.2223810344533745d0, &
&   0.1012285362903763d0]

 ifail=0
 ll1=ll+1
 q0=qmax-dble(nq-1)*dq
 if(q0<0.0d0) then
   write(error_unit,'(/a,i0,a,es12.4,a,es12.4)') &
&       'eresid: ERROR q grid below 0, nq=', &
&        nq,', dq=',dq,', qmax=',qmax
   error stop 1
 end if

!the basis functions' values and slopes at rc, for their overlaps with
!j_l(q r) at each q (quad precision: the overlaps enter sums with 1e4
!cancellation)
 rcq=real(rc,qp)
 do jj=1,nbas
   call sbfq_value_slope(ll,real(qroot(jj),qp)*rcq,jq(jj),jdq(jj))
 end do

!the linear mesh of the tail
 if(nr<2*lm_ncell) then
   write(iout,'(/a,i0)') 'eresid: ERROR tail too short, nr=',nr
   ifail=1
   return
 end if

 allocate(rlin(nr),wr(nr),sbf(nr))
 allocate(sqr(nr),cqr(nr),sdqr(nr),cdqr(nr))
 allocate(g0(nq),gdot(nq,nnull),gmat(nq,nnull,nnull),work(nq))

 do ii=1,nr
   rlin(ii)=rc+dble(ii-1)*hh
 end do
 call ug_weights(nr,hh,wr)

!the tail's transforms at the nodes (above): the weighted tail, and the
!nodes below the first where the terms of the closed forms stop cancelling
 allocate(aw(nr),tailg(nq))
 aw(:)=wr(:)*ulin(:)
 tailg(:)=0.0d0
 call czt_split()

!the nodes below, by the sums
 call tail_by_sums()

!the nodes above, by the chirp-z transform
 if(nqd<nq) call tail_by_czt()

!the integrands at each node of the q grid, their integrals from each node
!to qmax, and the values at each qout
 call node_integrals()
 call at_qout()

 deallocate(rlin,wr,sbf)
 deallocate(sqr,cqr,sdqr,cdqr)
 deallocate(g0,gdot,gmat,work)
 deallocate(aw,tailg)

 return

 contains

 subroutine czt_split()

!! the first node nqd+1 from which the tail's transforms are those of the
!! chirp-z transform: above the last where the closed forms' terms cancel
!! (nqd=nq for l>3, or for fewer than nczt nodes above)

!Local variables
 real(dp) :: pk(0:4),est,qq
 integer :: ii,iq,kk

 nqd=nq
 if(ll<=3) then
   pk(:)=0.0d0
   do ii=1,nr
     if(rlin(ii)>0.0d0) then
       do kk=0,ll+1
         pk(kk)=pk(kk)+abs(aw(ii))*rlin(ii)**(1-kk)
       end do
     end if
   end do
   do iq=1,nq
     qq=q0+dble(iq-1)*dq
     if(qq==0.0d0) cycle
     est=0.0d0
     do kk=1,ll+1
       est=est+ccoef(kk,ll)*pk(kk)/qq**kk
     end do
     if(est<=0.25d0*pk(0)) then
       nqd=iq-1
       exit
     end if
   end do
   if(nq-nqd<nczt) nqd=nq
 end if
 return
 end subroutine czt_split

 subroutine tail_by_sums()

!! the tail's transforms at the nodes 1..nqd by the sums: sin and cos of
!! dq*rlin step sin and cos of qq*rlin from one qq to the next by rotation

!Local variables
 real(dp) :: qq,tn
 integer :: ii,iq

 do ii=1,nr
   sdqr(ii)=sin(dq*rlin(ii))
 end do
 do ii=1,nr
   cdqr(ii)=cos(dq*rlin(ii))
 end do
 do iq=1,nqd
   qq=q0+dble(iq-1)*dq

!sin and cos of qq*rlin, evaluated directly every nqseed steps and
!otherwise rotated by dq*rlin from the previous qq
!(separate loops so sin and cos each vectorize)
   if(mod(iq-1,nqseed)==0) then
     do ii=1,nr
       sqr(ii)=sin(qq*rlin(ii))
     end do
     do ii=1,nr
       cqr(ii)=cos(qq*rlin(ii))
     end do
   else
     do ii=1,nr
       tn=sqr(ii)*cdqr(ii)+cqr(ii)*sdqr(ii)
       cqr(ii)=cqr(ii)*cdqr(ii)-sqr(ii)*sdqr(ii)
       sqr(ii)=tn
     end do
   end if
   if(qq==0.0d0) cycle
   tailg(iq)=tail_sum(qq)
 end do
 return
 end subroutine tail_by_sums

 subroutine tail_by_czt()

!! the tail's transforms at the nodes nqd+1..nq by the chirp-z transform:
!! F_k(q) = sum aw r**(1-k) exp(i q r) on q = qlo + m dq, r = rlin(1) +
!! n hh, is exp(i q rlin(1)) times the transform of aw r**(1-k) exp(i qlo
!! n hh) at the angle dq hh

!Local variables
 complex(dp), allocatable :: xin(:),fk(:,:)
 type(czt_plan_t) :: plan
 real(dp) :: qlo,qq
 integer :: ii,iq,kk,mq

 mq=nq-nqd
 qlo=q0+dble(nqd)*dq
 allocate(xin(0:nr-1),fk(0:mq-1,ll+1))
 call czt_plan(nr,mq,dq*hh,plan)
 do kk=1,ll+1
   do ii=1,nr
     xin(ii-1)=(0.0d0,0.0d0)
     if(rlin(ii)>0.0d0) then
       xin(ii-1)=aw(ii)*rlin(ii)**(1-kk) &
&                 *cmplx(cos(qlo*dble(ii-1)*hh),sin(qlo*dble(ii-1)*hh),dp)
     end if
   end do
   call czt_apply(plan,xin,fk(:,kk))
 end do
 call czt_free(plan)
 do iq=nqd+1,nq
   qq=q0+dble(iq-1)*dq
   fk(iq-nqd-1,:)=fk(iq-nqd-1,:)*cmplx(cos(qq*rlin(1)),sin(qq*rlin(1)),dp)
   tailg(iq)=czt_combine(qq,fk(iq-nqd-1,:))
 end do
 deallocate(xin,fk)
 return
 end subroutine tail_by_czt

 subroutine node_integrals()

!! the integrands at each node of the q grid and their integrals from
!! each node to qmax: e0g, edotg, ematg

!Local variables
 real(dp) :: qq,qq4,ps0q,psnullq(nnull)
 integer :: iq,jj,kk

 do iq=1,nq
   qq=q0+dble(iq-1)*dq

!q**4 makes the integrands zero at q=0 whatever the transforms
   if(qq==0.0d0) then
     g0(iq)=0.0d0
     gdot(iq,:)=0.0d0
     gmat(iq,:,:)=0.0d0
     cycle
   end if

   call transforms(qq,tailg(iq),ps0q,psnullq)

!integrands of the E_resid q integrals
   qq4=qq**4
   g0(iq)=qq4*ps0q**2
   gdot(iq,:)=qq4*ps0q*psnullq(:)
   do jj=1,nnull
     gmat(iq,:,jj)=qq4*psnullq(:)*psnullq(jj)
   end do

 end do !iq

!integrals from each node to qmax
!1/pi factor from implicit sbf*sbf integral
 call ug_cumint_down(g0,nq,dq,e0g)
 e0g(:)=e0g(:)/pi
 do jj=1,nnull
   call ug_cumint_down(gdot(:,jj),nq,dq,work)
   edotg(:,jj)=work(:)/pi
   do kk=1,nnull
     call ug_cumint_down(gmat(:,kk,jj),nq,dq,work)
     ematg(:,kk,jj)=work(:)/pi
   end do
 end do
 return
 end subroutine node_integrals

 subroutine at_qout()

!! E at each qout: the value at the next node plus the partial cell by
!! Gauss-Legendre, with the transforms at its nodes

!Local variables
 real(dp) :: qq,qq4,ps0q,psnullq(nnull),qa,qb,wg
 integer :: ii,jj,kk,ig,icell

 do ii=1,nqout
   if(qout(ii)<0.0d0 .or. qout(ii)>=qmax) then
     write(error_unit,'(/a,es12.4,a,es12.4)') &
&         'eresid: ERROR q off the grid, q=', &
&        qout(ii),' qmax=',qmax
     error stop 1
   end if
! icell -1 for a qout below the first node
   icell=floor((qout(ii)-q0)/dq)
   qa=qout(ii)
   qb=q0+dble(icell+1)*dq
   if(qa==q0+dble(icell)*dq) then
! on a node
     eresid0(ii)=e0g(icell+1)
     eresiddot(:,ii)=edotg(icell+1,:)
     do jj=1,nnull
       eresidmat(:,jj,ii)=ematg(icell+1,:,jj)
     end do
     cycle
   end if
   eresid0(ii)=e0g(icell+2)
   eresiddot(:,ii)=edotg(icell+2,:)
   do jj=1,nnull
     eresidmat(:,jj,ii)=ematg(icell+2,:,jj)
   end do
   do ig=1,ngl
     qq=0.5d0*(qa+qb)+0.5d0*(qb-qa)*xgl(ig)
     wg=0.5d0*(qb-qa)*wgl(ig)/pi
     do kk=1,nr
       sqr(kk)=sin(qq*rlin(kk))
     end do
     do kk=1,nr
       cqr(kk)=cos(qq*rlin(kk))
     end do
     call transforms(qq,tail_sum(qq),ps0q,psnullq)
     qq4=wg*qq**4
     eresid0(ii)=eresid0(ii)+qq4*ps0q**2
     eresiddot(:,ii)=eresiddot(:,ii)+qq4*ps0q*psnullq(:)
     do jj=1,nnull
       eresidmat(:,jj,ii)=eresidmat(:,jj,ii)+qq4*psnullq(:)*psnullq(jj)
     end do
   end do
 end do
 return
 end subroutine at_qout

 real(dp) function tail_sum(qq)

!! the tail's transform at qq by its sum, from sin and cos of qq*rlin in
!! sqr and cqr

!Input variables
 real(dp), intent(in) :: qq

!Local variables
 integer :: ii

 call rsbf_lin(ll1,qq,rlin,sqr,cqr,nr,sbf)
 tail_sum=0.0d0
 do ii=1,nr
   tail_sum=tail_sum+aw(ii)*sbf(ii)
 end do
 return
 end function tail_sum

 real(dp) function czt_combine(qq,ff)

!! the tail's transform at qq from F_k = sum aw r**(1-k) exp(i qq r):
!! the closed forms of r j_l(qq r), their sines from Im F_k and cosines
!! from Re F_k

!Input variables
 real(dp), intent(in) :: qq
 complex(dp), intent(in) :: ff(ll+1)

 select case(ll)
 case(0)
   czt_combine=aimag(ff(1))/qq
 case(1)
   czt_combine=aimag(ff(2))/qq**2-real(ff(1),dp)/qq
 case(2)
   czt_combine=3.0d0*aimag(ff(3))/qq**3-aimag(ff(1))/qq &
&             -3.0d0*real(ff(2),dp)/qq**2
 case default
   czt_combine=15.0d0*aimag(ff(4))/qq**4-6.0d0*aimag(ff(2))/qq**2 &
&             -15.0d0*real(ff(3),dp)/qq**3+real(ff(1),dp)/qq
 end select
 return
 end function czt_combine

 subroutine transforms(qq,tailq,ps0q,psnullq)

!! the transforms at qq of ps0 (the Bessel sum inside rc, the all-electron
!! tail outside, whose transform is tailq) and of the null functions
!! (inside rc only)

!Input variables
 real(dp), intent(in) :: qq,tailq

!Output variables
 real(dp), intent(out) :: ps0q,psnullq(nnull)

!Local variables
 integer :: jj

 ps0q=tailq

!transforms inside rc: the basis functions' overlaps with j_l(qq r), in
!closed form, combined in quad precision
 call sbfq_value_slope(ll,real(qq,qp)*rcq,jb,jdb)
 do jj=1,nbas
   call sbfq_overlap_vs(ll,real(qroot(jj),qp),real(qq,qp),rcq,jq(jj),jdq(jj), &
&                       jb,jdb,ovl(jj))
 end do
 ps0in=0.0_qp
 psnin(:)=0.0_qp
 do jj=1,nbas
   ps0in=ps0in+pswf0_sb(jj)*ovl(jj)
   psnin(:)=psnin(:)+pswfnull_sb(jj,:)*ovl(jj)
 end do
 ps0q=ps0q+real(ps0in,dp)
 psnullq(:)=real(psnin(:),dp)
 return
 end subroutine transforms

 subroutine rsbf_lin(ll1,qq,rlin,sqr,cqr,nr,rsbf)

!! rlin * j_l(qq rlin) on the linear mesh from sin and cos of qq rlin, by
!! sbf_trig_rising for l<=3 (the closed forms, the series where they
!! cancel) and sbf_values beyond

!Input variables
 integer, intent(in) :: ll1,nr
 real(dp), intent(in) :: qq,rlin(nr),sqr(nr),cqr(nr)

!Output variables
 real(dp), intent(out) :: rsbf(nr)

!Local variables
 integer :: ii

 if(ll1>4) then
   do ii=1,nr
     call sbf_values(ll1,qq*rlin(ii),sb_out)
     rsbf(ii)=rlin(ii)*sb_out(ll1)
   end do
   return
 end if

 call sbf_trig_rising(ll1-1,nr,qq,rlin,sqr,cqr,rsbf)

 return
 end subroutine rsbf_lin

 end subroutine  eresid

 end module eresid_m
