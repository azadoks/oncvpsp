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
! calculates residual energy operator matrix elements as a function of q

 subroutine eresid(ll,rc,nnull,nbas,nq,dq,nqout,qout,qmax,qroot,hh,nr, &
&                  ulin,pswf0_sb,pswfnull_sb,e0g,edotg,ematg, &
&                  eresid0,eresiddot,eresidmat)

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
!! integrated by the 8th-order rule of ugrid_m. The
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

 use constants_m, only: dp, pi
 use, intrinsic :: iso_fortran_env, only: qp => real128
 use logmesh_m, only: lm_ncell
 use ugrid_m, only: ug_weights, ug_cumint_down
 use sbf_q_m, only: sbfq_overlap_vs, sbfq_value_slope
 use sbf_m, only: sbf_trig_rising
 implicit none
 external :: sbf8

!Input variables
 integer, intent(in) :: ll,nnull,nbas,nq,nqout,nr
 real(dp), intent(in) :: rc,qmax,hh
 real(dp), intent(in) :: ulin(nr),qout(nqout),dq
 real(qp), intent(in) :: pswf0_sb(nbas),pswfnull_sb(nbas,nnull)
 real(dp), intent(in) :: qroot(nbas)

!Output variables
 real(dp), intent(out) :: e0g(nq),edotg(nq,nnull),ematg(nq,nnull,nnull)
 real(dp), intent(out) :: eresid0(nqout),eresiddot(nnull,nqout)
 real(dp), intent(out) :: eresidmat(nnull,nnull,nqout)

!Local variables
 real(dp), allocatable :: rlin(:),wr(:),sbf(:)
 real(dp), allocatable :: sqr(:),cqr(:),sdqr(:),cdqr(:)
 real(dp), allocatable :: g0(:),gdot(:,:),gmat(:,:,:),work(:)
 real(dp) :: sb_out(5),psnullq(nnull)
 real(qp) :: ovl(nbas),ps0in,psnin(nnull),rcq,jq(nbas),jdq(nbas),jb,jdb
 real(dp) :: q0,qq,qq4,ps0q,tn,qa,qb,wg
 integer :: ii,iq,jj,kk,ll1,ig,icell
!steps between exact sin/cos(qq*rlin) evaluations in the q loop
 integer, parameter :: nqseed=32
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

 ll1=ll+1
 q0=qmax-dble(nq-1)*dq
 if(q0<0.0d0) then
   write(6,'(/a,i0,a,es12.4,a,es12.4)') 'eresid: ERROR q grid below 0, nq=', &
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
   write(6,'(/a,i0)') 'eresid: ERROR tail too short, nr=',nr
   error stop 1
 end if

 allocate(rlin(nr),wr(nr),sbf(nr))
 allocate(sqr(nr),cqr(nr),sdqr(nr),cdqr(nr))
 allocate(g0(nq),gdot(nq,nnull),gmat(nq,nnull,nnull),work(nq))

 do ii=1,nr
   rlin(ii)=rc+dble(ii-1)*hh
 end do
 call ug_weights(nr,hh,wr)

!sin and cos of dq*rlin, which step sin and cos of qq*rlin from one qq to
!the next by rotation
 do ii=1,nr
   sdqr(ii)=sin(dq*rlin(ii))
 end do
 do ii=1,nr
   cdqr(ii)=cos(dq*rlin(ii))
 end do

!the integrands at each node of the q grid
 do iq=1,nq
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

!q**4 makes the integrands zero at q=0 whatever the transforms
   if(qq==0.0d0) then
     g0(iq)=0.0d0
     gdot(iq,:)=0.0d0
     gmat(iq,:,:)=0.0d0
     cycle
   end if

   call transforms(qq,ps0q,psnullq)

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

!at each qout: the value at the next node plus the partial cell by
!Gauss-Legendre, with the transforms at its nodes
 do ii=1,nqout
   if(qout(ii)<0.0d0 .or. qout(ii)>=qmax) then
     write(6,'(a,es12.4,a,es12.4)') 'eresid: ERROR q off the grid, q=', &
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
     call transforms(qq,ps0q,psnullq)
     qq4=wg*qq**4
     eresid0(ii)=eresid0(ii)+qq4*ps0q**2
     eresiddot(:,ii)=eresiddot(:,ii)+qq4*ps0q*psnullq(:)
     do jj=1,nnull
       eresidmat(:,jj,ii)=eresidmat(:,jj,ii)+qq4*psnullq(:)*psnullq(jj)
     end do
   end do
 end do

 deallocate(rlin,wr,sbf)
 deallocate(sqr,cqr,sdqr,cdqr)
 deallocate(g0,gdot,gmat,work)

 return

 contains

 subroutine transforms(qq,ps0q,psnullq)

!! the transforms at qq of ps0 (the Bessel sum inside rc, the all-electron
!! tail outside) and of the null functions (inside rc only), from sin and
!! cos of qq*rlin in sqr and cqr

!Input variables
 real(dp), intent(in) :: qq

!Output variables
 real(dp), intent(out) :: ps0q,psnullq(nnull)

!Local variables
 integer :: ii,jj

!transform of the all-electron tail: rlin * j_l(qq rlin) on the linear
!mesh, times uu, integrated
 call rsbf_lin(ll1,qq,rlin,sqr,cqr,nr,sbf)
 ps0q=0.0d0
 do ii=1,nr
   ps0q=ps0q+wr(ii)*ulin(ii)*sbf(ii)
 end do

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
!! cancel) and sbf8 beyond

!Input variables
 integer, intent(in) :: ll1,nr
 real(dp), intent(in) :: qq,rlin(nr),sqr(nr),cqr(nr)

!Output variables
 real(dp), intent(out) :: rsbf(nr)

!Local variables
 integer :: ii

 if(ll1>4) then
   do ii=1,nr
     call sbf8(ll1,qq*rlin(ii),sb_out)
     rsbf(ii)=rlin(ii)*sb_out(ll1)
   end do
   return
 end if

 call sbf_trig_rising(ll1-1,nr,qq,rlin,sqr,cqr,rsbf)

 return
 end subroutine rsbf_lin

 end subroutine  eresid
