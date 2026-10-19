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

 subroutine eresid(ll,irc,nnull,nbas,mmax,rr,drdum,dqdum,qmax,qroot, &
&                  uu,pswf0_sb,pswfnull_sb, nqout,qout, &
&                  eresid0,eresiddot,eresidmat)

!ll  anglar momentum
!irc  log mesh index of rc
!nnull  number of unconstrained basis vectors for residual minimization
!nbas  number of sbf basis functions
!mmax  index of largest radius for non-zero value of wave function
!rr(mmax)  log radial mesh
!drdum, dqdum  (not used)
!qmax  largest wave vector for E_res integration
!qroot(nbas)  q values for sbfs in basis
!uu(mmax)  r*all-electron wave function
!pswf0_sb(nbas)  sbf basis  coefficients for constrained part of pseudo
!                    wave function
!pswfnull_sb(nbas,nnull) sbf coefficients for unconstrained basis of pseudo
!                       wave function
!nqout  number of wave vector lower cutoffs
!qout(nqout)  wave vector lower cutoffs for E_res for which output is wanted,
!       rounded on output to nearest multiple of dq
!nqout  number of q values in [0,qmax] for which results are to be saved
!qout(nqout)  set of these q values
!eresid0(nqout)  set of <pswf0| E_resid |pswf0> matrix elements
!eresiddot(nnull,nqout) set of <pswfnull| E_resid |pswf0> matrix elements
!eresidevec(nnull,null,nqout)  set of <pswfnull| E_res |pswfnull'> matrices

 use constants_m, only: dp, pi
 implicit none

!Input variables
 integer, intent(in) :: irc,ll,nnull,nbas,mmax,nqout
 real(dp), intent(in) :: rr(mmax),uu(mmax),pswf0_sb(nbas),pswfnull_sb(nbas,nnull)
 real(dp), intent(in) :: qroot(nbas)
!real(dp) :: dr,dq,qmax
 real(dp), intent(in) :: drdum,dqdum,qmax

!Input/output variables
 real(dp), intent(inout) :: qout(nqout)
 real(dp), intent(inout) :: eresid0(nqout),eresiddot(nnull,nqout)
 real(dp), intent(inout) :: eresidmat(nnull,nnull,nqout)

!Local variables
 real(dp), allocatable :: rlin(:),uulin(:),sbf(:),ps0(:)
 real(dp), allocatable :: erdot(:),ereval(:),erevec(:,:)
 real(dp), allocatable :: psnull(:,:),tpsnull(:),psnullq(:)
 real(dp), allocatable :: dac0(:),dacdot(:,:),dacmat(:,:,:)
 real(dp), allocatable :: acdot(:),acmat(:,:)
 real(dp), allocatable :: work(:),wmat(:,:),wev(:),wvec(:)
 real(dp), allocatable :: sqr(:),cqr(:),sdqr(:),cdqr(:)
 real(dp) :: sb_out(5)
 real(dp) :: ac0,rc,ps0q,tps0,xx,qq,qq4,sn,tn,td
 real(dp) :: dr,dq

 integer :: nrlin,nq,irclin,mmaxlin
 integer :: ii,info,iq,jj,kk,ll1
!steps between exact sin/cos(qq*rlin) evaluations in the q loop
 integer, parameter :: nqseed=32


! new calculations of dq, dr
! larger denominators yield better convergence at the cost of more time
! dq=2.0d0*pi/(300.0d0*rr(irc))
  dq=2.0d0*pi/(150.0d0*rr(irc))
! dq=2.0d0*pi/(75.0d0*rr(irc))
! dr=2.0d0*pi/(150.0d0*qmax)
  dr=2.0d0*pi/(75.0d0*qmax)

!set up linear radial mesh
 rc=rr(irc)
 irclin=rc/dr + 1.5d0
 dr=rc/(irclin-1)

 mmaxlin=0
 do ii=mmax,1,-1
   if(abs(uu(ii))>1.0d-10) then
     mmaxlin=rr(ii)/dr - 2.0d0
     exit
    end if
 end do
!mmaxlin=rr(mmax)/dr - 1.5d0 !leave non-zero room for cubic interpolation
 if(mmaxlin<irclin+6) then
   write(6,'(/2a,i0,a,i0)') 'eresid: ERROR wave function vanishes too close', &
&       ' to rc, mmaxlin=',mmaxlin,', irclin=',irclin
   error stop 1
 end if

 allocate(rlin(mmaxlin),uulin(mmaxlin),ps0(mmaxlin))
 allocate(psnull(irclin,nnull),tpsnull(nnull),psnullq(nnull))
 allocate(sbf(mmaxlin))
 allocate(sqr(mmaxlin),cqr(mmaxlin),sdqr(mmaxlin),cdqr(mmaxlin))
 allocate(dac0(4),dacdot(nnull,4),dacmat(nnull,nnull,4))
 allocate(acdot(nnull),acmat(nnull,nnull))


 do ii=1,mmaxlin
  rlin(ii)=(ii-1)*dr
 end do

!interpolate all-electron wave function onto linear mesh, even though we
!will only use the tail region

 call dpnint(rr,uu,mmax,rlin,uulin,mmaxlin)

 do ii=irclin,mmaxlin
  ps0(ii)=uulin(ii)
 end do

!set up constrained basis functions on linear grid inside rc (multiplied by r)
 ll1=ll+1
 do ii=1,irclin
  tps0=0.0d0
  tpsnull(:)=0.0d0
  do jj=1,nbas
   xx=rlin(ii)*qroot(jj)
   call sbf8(ll1,xx,sb_out)
   tps0      =tps0      +pswf0_sb(jj)     *sb_out(ll1)
   tpsnull(:)=tpsnull(:)+pswfnull_sb(jj,:)*sb_out(ll1)
  end do
  ps0(ii)=rlin(ii)*tps0
  psnull(ii,:)=rlin(ii)*tpsnull(:)
 end do

!big q integration loop

!q integrals are performed on all matrix elements simultaneously starting
!inward from qmax (considered infinity) to qq=0, saving snapshots along the
!way

 nq=qmax/dq+1

!sin and cos of dq*rlin, which step sin and cos of qq*rlin from one qq to
!the next by rotation
 do ii=1,mmaxlin
  sdqr(ii)=sin(dq*rlin(ii))
 end do
 do ii=1,mmaxlin
  cdqr(ii)=cos(dq*rlin(ii))
 end do

!null accumulators and running registers

 ac0=0.0d0
 acdot(:)=0.0d0
 acmat(:,:)=0.0d0

 dac0(:)=0.0d0
 dacdot(:,:)=0.0d0
 dacmat(:,:,:)=0.0d0

 do iq=nq,1,-1
  qq=dq*(iq-1)

!sin and cos of qq*rlin, evaluated directly every nqseed steps and
!otherwise rotated inward by dq*rlin from the previous qq
!(separate loops so sin and cos each vectorize)
  if(mod(nq-iq,nqseed)==0) then
   do ii=1,mmaxlin
    sqr(ii)=sin(qq*rlin(ii))
   end do
   do ii=1,mmaxlin
    cqr(ii)=cos(qq*rlin(ii))
   end do
  else
   do ii=1,mmaxlin
    tn=sqr(ii)*cdqr(ii)-cqr(ii)*sdqr(ii)
    cqr(ii)=cqr(ii)*cdqr(ii)+sqr(ii)*sdqr(ii)
    sqr(ii)=tn
   end do
  end if

!sbfs on full rlin mesh for this qq
!multiply by rlin for integration
  call rsbf_lin(ll1,qq,rlin,sqr,cqr,mmaxlin,sbf)

!Fourier transforms of basis functions
!null functions exist for r<rc

  ps0q=( 9.0d0*ps0(1       )*sbf(1       ) &
&      +28.0d0*ps0(2       )*sbf(2       ) &
&      +23.0d0*ps0(3       )*sbf(3    )    &
&      +23.0d0*ps0(irclin-2)*sbf(irclin-2) &
&      +28.0d0*ps0(irclin-1)*sbf(irclin-1) &
&      + 9.0d0*ps0(irclin  )*sbf(irclin  ))/24.0d0

  psnullq(:)=( 9.0d0*psnull(1       ,:)*sbf(1       ) &
&            +28.0d0*psnull(2       ,:)*sbf(2       ) &
&            +23.0d0*psnull(3       ,:)*sbf(3    )    &
&            +23.0d0*psnull(irclin-2,:)*sbf(irclin-2) &
&            +28.0d0*psnull(irclin-1,:)*sbf(irclin-1) &
&            + 9.0d0*psnull(irclin  ,:)*sbf(irclin  ))/24.0d0

  do ii=4,irclin-3
   ps0q=ps0q+ps0(ii)*sbf(ii)
  end do
  do jj=1,nnull
   do ii=4,irclin-3
    psnullq(jj)=psnullq(jj)+psnull(ii,jj)*sbf(ii)
   end do
  end do

  ps0q=ps0q +( 9.0d0*ps0(irclin   )*sbf(irclin   ) &
&            +28.0d0*ps0(irclin+1 )*sbf(irclin+1 ) &
&            +23.0d0*ps0(irclin+2 )*sbf(irclin+2 )    &
&            +23.0d0*ps0(mmaxlin-2)*sbf(mmaxlin-2) &
&            +28.0d0*ps0(mmaxlin-1)*sbf(mmaxlin-1) &
&            + 9.0d0*ps0(mmaxlin  )*sbf(mmaxlin  ))/24.0d0

  do ii=irclin+3,mmaxlin-3
   ps0q=ps0q+ps0(ii)*sbf(ii)
  end do

  ps0q=dr*ps0q
  psnullq(:)=dr*psnullq(:)

!stuff integrands for E_resid q integration
  qq4=dq*qq**4
  dac0(1)=qq4*ps0q**2
  dacdot(:,1)=qq4*ps0q*psnullq(:)
  do jj=1,nnull
    dacmat(:,jj,1)=qq4*psnullq(:)*psnullq(jj)
  end do

!inward integration as soon as dac arrays ("derivative accumulator arrays")
!are full
  if(iq<nq-3) then
   ac0=ac0+(9.0d0*dac0(1)+19.0d0*dac0(2) &
&      -5.0d0*dac0(3)+dac0(4))/24.0d0

   acdot(:)=acdot(:)+(9.0d0*dacdot(:,1)+19.0d0*dacdot(:,2) &
&      -5.0d0*dacdot(:,3)+dacdot(:,4))/24.0d0

   acmat(:,:)=acmat(:,:)+(9.0d0*dacmat(:,:,1)+19.0d0*dacmat(:,:,2) &
&      -5.0d0*dacmat(:,:,3)+dacmat(:,:,4))/24.0d0
  end if

!shift dac* arrays one step inward

  do ii=4,2,-1
   dac0(ii)=dac0(ii-1)
   dacdot(:,ii)=dacdot(:,ii-1)
   dacmat(:,:,ii)=dacmat(:,:,ii-1)
  end do
  dac0(1)=0.0d0
  dacdot(:,1)=0.0d0
  dacmat(:,:,1)=0.0d0

!test if present qq value is a "save" point
!1/pi factor from implicit sbf*sbf integral
  do ii=1,nqout
   if(abs(qout(ii)-qq)<0.5d0*dq) then
    eresid0(ii)=ac0/pi
    eresiddot(:,ii)=acdot(:)/pi
    eresidmat(:,:,ii)=acmat(:,:)/pi
    qout(ii)=qq
   end if
  end do

 end do !big iq loop


 deallocate(rlin,uulin,ps0)
 deallocate(psnull,tpsnull,psnullq)
 deallocate(sbf)
 deallocate(sqr,cqr,sdqr,cdqr)
 deallocate(dac0,dacdot,dacmat)
 deallocate(acdot,acmat)

 return

 contains

 subroutine rsbf_lin(ll1,qq,rlin,sqr,cqr,nr,rsbf)
! rlin(ii)*j_l(qq*rlin(ii)) for l=ll1-1 on the linear mesh rlin(nr), given
! sqr=sin(qq*rlin) and cqr=cos(qq*rlin).  Same branches and formulas as
! sbf8 for l<=3, written as loops over the mesh that vectorize; this is
! the inner loop of the E_resid q integration.

!Input variables
 integer, intent(in) :: ll1,nr
 real(dp), intent(in) :: qq,rlin(nr),sqr(nr),cqr(nr)

!Output variables
 real(dp), intent(out) :: rsbf(nr)

!Local variables
 real(dp) :: xx,xn,xs
 integer :: ii,i0,nn

 if(ll1>4) then
  do ii=1,nr
   call sbf8(ll1,qq*rlin(ii),sb_out)
   rsbf(ii)=rlin(ii)*sb_out(ll1)
  end do
  return
 end if

!rlin is increasing, so the small-argument points come first
 i0=nr+1
 do ii=1,nr
  xx=qq*rlin(ii)
  if(xx>=1.0d-3) then
   i0=ii
   exit
  end if
  if(xx<=1.0d-36) then
   if(ll1==1) then
    rsbf(ii)=rlin(ii)
   else
    rsbf(ii)=0.0d0
   end if
  else
   xn=1.0d0
   xs=0.5d0*xx**2
   do nn=1,ll1-1
    xn=xx*xn/(2*nn+1)
   end do
   rsbf(ii)=rlin(ii)*xn*(1.0d0 - xs*(1.0d0 - xs/(4*ll1+6))/(2*ll1+1))
  end if
 end do

 select case(ll1)
 case(1)
  do ii=i0,nr
   xx=qq*rlin(ii)
   rsbf(ii)=rlin(ii)*(sqr(ii)/xx)
  end do
 case(2)
  do ii=i0,nr
   xx=qq*rlin(ii)
   rsbf(ii)=rlin(ii)*((sqr(ii)-xx*cqr(ii))/xx**2)
  end do
 case(3)
  do ii=i0,nr
   xx=qq*rlin(ii)
   rsbf(ii)=rlin(ii)*(((3.0d0-xx**2)*sqr(ii) - 3.0d0*xx*cqr(ii))/xx**3)
  end do
 case(4)
  do ii=i0,nr
   xx=qq*rlin(ii)
   rsbf(ii)=rlin(ii)*((15.d0*sqr(ii)-15.d0*xx*cqr(ii) - 6.d0*xx**2*sqr(ii) &
&           + xx**3*cqr(ii))/xx**4)
  end do
 case default
! ll1>4 returned above
 end select

 return
 end subroutine rsbf_lin

 end subroutine  eresid
