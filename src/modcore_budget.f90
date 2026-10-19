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
! Creates the model core charge of icmod=6 and 7: the all-electron core
! charge wherever the plane-wave cutoff affords it.  Inside the radius r_s
! the model is a fill that matches the core charge and its first two
! derivatives at r_s; beyond r_s it is the core charge.
!
! The fill is the minimum tail: the even polynomial
!
!   sum_{k=0..8} c_k (r/r_s)**(2k)
!
! whose six free coefficients minimize the Fourier tail of the whole model,
! int_{qlo}^{qhi} q**2 |n(q)|**2 dq, with qlo = 2*sqrt(2*ecut) and
! qhi = min(3*qlo, qmax).  Where that polynomial rises somewhere on
! [0, r_s], the fill is Martins' exp(a + b r**2 + c r**4) instead, if
! that does not rise (b <= 0 and b + 2c r_s**2 <= 0); where neither is
! monotone, the radius has no model.
!
! ecut is fcfact (Ha) if fcfact > 0.  Otherwise it is the cutoff the
! pseudo-waves need for the same error (ecutw, from modcore_budget_ecut),
! times -fcfact if fcfact < 0: a calculation needs that cutoff anyway.
! The model must fit slack*ecut (0 < slack <= 1), and is still made for
! the qlo of ecut.
!
! The error of a density cutoff q (a plane-wave cutoff q**2/8) is that of
! the core charge without its Fourier components beyond q, in a solid:
!
!   |Exc[n + nc_q] - Exc[n + nc] - Vxc[n_bg] * int (nc_q - nc)|
!
! within the sphere r <= 10 bohr, with n the pseudo valence charge on a
! uniform background n_bg = 5e-3 bohr**-3 that stands for the interstitial
! electrons of the neighbours.  A model fits the cutoff e if the error
! stays below errmax for every density cutoff from 2*sqrt(2*e) on, in
! steps of 0.1 bohr**-1 up to 120 bohr**-1 (1800 Ha).
!
! icmod=6 (BudgetedCore of PseudoPotentialGenerators.jl): r_s is the
! smallest grid radius, no larger than the smallest rc (or the last radius
! with core charge), whose model fits slack*ecut.  The error almost always
! peaks at the first of the cutoffs, so r_s is bisected for with the first
! nquick of them only, and then confirmed with all of them (and bisected
! for again with all of them if that fails).
!
! icmod=7 (HardnessMatchedCore): of the radii from that largest one down by
! factors of 1.025 to the r_s of icmod=6 over 1.25 (that r_s included, on
! the nearest grid points), the one whose model is closest to the core
! charge in the exchange-correlation hardness of the pseudo valence: the
! rms over i,j of d2Exc/dn_i dn_j at rhotps + model against that at
! rhotps + rhoc (der2exc), for the probe densities of modcore_probes.
! The models are checked against slack*ecut in that order, the first that
! fits is taken; ties go to the larger radius.

 subroutine modcore_budget(icmod,rhotps,rhoc,rhomod,fcfact,ecutw,errmax, &
&                          slack,irmax,mmax,rr,iexc,zion,probes,lprobe, &
&                          nprobe)

!! @param icmod  6 for the smallest radius, 7 for the hardness-matched one
!! @param rhotps  total pseudocharge density
!! @param rhoc  core-charge density
!! @param rhomod  model core density and 4 derivatives
!! @param fcfact  plane-wave cutoff (Ha) the model is made for, if > 0; else 0
!!                for ecutw, or minus the multiple of ecutw
!! @param ecutw  cutoff (Ha) the pseudo-waves need for errmax (used if fcfact <=
!!               0)
!! @param errmax  exchange-correlation energy error (Ha) that defines the needed
!!                cutoff
!! @param slack  share of the cutoff the model must fit, 0 < slack <= 1
!! @param irmax  rr index of the largest r_s (smallest rc)
!! @param mmax  dimension of log grid
!! @param rr  log radial grid
!! @param iexc  exchange-correlation function to be used
!! @param zion  ion charge (for der2exc)
!! @param probes  probe densities of the hardness (icmod=7), one electron each
!! @param lprobe  angular momenta of the probes
!! @param nprobe  number of probes (0 for icmod=6)

 use constants_m, only: dp, pi
 implicit none

!Input variables
 integer, intent(in) :: icmod,mmax,irmax,iexc,nprobe
 integer, intent(in) :: lprobe(nprobe)
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: probes(mmax,nprobe)
 real(dp), intent(in) :: fcfact,ecutw,errmax,slack,zion

!Output variables
 real(dp), intent(out) :: rhomod(mmax,5)

!Local variables
 real(dp), parameter :: pi4=4.0_dp*pi
!background density (bohr**-3) and radius (bohr) of the error,
!largest density cutoff and its step (bohr**-1),
!panels of the band-limited transform per step,
!cutoffs of the bisection before the confirmation with all of them
 real(dp), parameter :: background=5.0d-3, radius=10.0d0
 real(dp), parameter :: qmax=120.0d0, qstep=0.1d0
 integer, parameter :: npanel=4, nquick=10
!highest power (r/r_s)**(2*npow) of the minimum tail, Gauss-Legendre nodes
!per 1 bohr**-1 panel of its window, points of its monotonicity check
 integer, parameter :: npow=8, ngq=8, nmono=201
!icmod=7: step and span of the radii
 real(dp), parameter :: ratio=1.025d0, span=1.25d0
 real(dp) :: al,dq,vbg,ecore,aco,bco,cco,e0,ecut,budget,factor,qlo,qhi
 real(dp) :: cpoly(0:npow),mismatch
 real(dp), allocatable :: dc(:,:),core(:),nval(:),sn(:),cs(:),sh(:),ch(:)
 real(dp), allocatable :: rq(:),work(:),eps(:),vxc(:),err(:),rqs(:,:)
 real(dp), allocatable :: qn(:),wn(:)
 real(dp), allocatable :: d2ae(:,:),d2mod(:,:),d2zero(:,:),cmis(:)
 integer, allocatable :: cidx(:),cfill(:),corder(:)
 character(len=1), allocatable :: cfit(:)
 integer :: ii,jj,nq,nqn,nin,nx,ilo,irs,irs6,ifrom,iquick,ibad,irtop
 integer :: ifill,ncand,kk
 real(dp) :: rad
 character(len=2), parameter :: fillname(0:2)=['- ','MT','M ']

 if(fcfact>0.0d0) then
   ecut=fcfact
 else
   factor=1.0d0
   if(fcfact<0.0d0) factor=-fcfact
   ecut=factor*ecutw
 end if
 budget=slack*ecut

 al = 0.01d0 * log(rr(101) / rr(1))
 dq = qstep/npanel
 nq = nint(qmax/qstep)

 qlo=2.0d0*sqrt(2.0d0*ecut)
 if(qlo >= qmax) then
   write(6,'(/a,f10.2,a)') 'modcore_budget: ERROR stop - ecut beyond', &
&        qmax**2/8.0d0,' Ha'
   error stop 1
 end if
 qhi=min(3.0d0*qlo,qmax)

 allocate(dc(mmax,5),core(mmax),sn(mmax),cs(mmax),sh(mmax),ch(mmax))
 allocate(err(nq))

! core charge derivatives,
! 7-point numerical first derivatives applied successively
 dc(:,:)=0.0d0
 dc(:,1)=rhoc(:)
 do jj=2,5
   do ii=4,mmax-3
     dc(ii,jj)=(-dc(ii-3,jj-1)+ 9.d0*dc(ii-2,jj-1)&
&     -45.d0*dc(ii-1,jj-1)+45.d0*dc(ii+1,jj-1)&
&     -9.d0*dc(ii+2,jj-1)+dc(ii+3,jj-1))&
&     /(60.d0*al*rr(ii))
   end do
 end do

! the window of the minimum tail: Gauss-Legendre panels of at most
! 1 bohr**-1 on [qlo, qhi], the weights times q**2
 kk=max(1,ceiling(qhi-qlo))
 nqn=ngq*kk
 allocate(qn(nqn),wn(nqn))
 do ii=1,kk
   call gauss_legendre(qlo+(ii-1)*(qhi-qlo)/kk,qlo+ii*(qhi-qlo)/kk,ngq, &
&                      qn((ii-1)*ngq+1:ii*ngq),wn((ii-1)*ngq+1:ii*ngq))
 end do
 wn(:)=wn(:)*qn(:)**2

! the sphere of the error, and a few more points for the gradients
 nin=0
 do ii=1,mmax
   if(rr(ii)>radius) exit
   nin=ii
 end do
 nx=min(mmax,nin+8)
 allocate(nval(nx),rq(nx),work(nx),eps(nx),vxc(nx))

 nval(:)=pi4*background + rhotps(1:nx)

! xc potential of the background alone
 work(:)=pi4*background
 call xc_eps(work,eps,vxc)
 vbg=vxc(nx/2)

! the first density cutoff that decides whether the model fits budget:
! the needed cutoff is (q+qstep)**2/8 for the last q above errmax
 ifrom=nq
 do ii=1,nq
   if(((ii+1)*qstep)**2/8.0d0 > budget) then
     ifrom=ii
     exit
   end if
 end do
 iquick=min(nq,ifrom+nquick-1)

! the largest r_s: the smallest rc, or the last radius at which the core
! charge and the derivatives matched there are not zero
 irtop=0
 do ii=1,min(irmax,mmax-6)
   if(.not. all(rhoc(ii:ii+6)>0.0d0)) exit
   irtop=ii
 end do
 if(irtop<4) then
   write(6,'(/a)') 'modcore_budget: ERROR stop - no core charge'
   error stop 1
 end if

! smallest r_s on the grid, by bisection in the index between rr(irtop)/50
! and rr(irtop) (the needed cutoff falls as the model softens).
! A model that fails at the first cutoffs fails at all of them
 ilo=irtop
 do ii=irtop,4,-1
   if(rr(ii)<rr(irtop)/50.0d0) exit
   ilo=ii
 end do

 if(.not. fits(irtop,iquick)) call no_model()
 irs=bisect(ilo,irtop,iquick)
 if(.not. fits(irs,nq)) then
   if(irs==irtop) call no_model()
   if(.not. fits(irtop,nq)) call no_model()
   irs=bisect(irs,irtop,nq)
 end if
 irs6=irs

! icmod=7: the candidate radii, largest first, their hardness mismatch,
! and the first in the order of the mismatch that fits
 if(icmod==7) then
   kk=floor(log(rr(irtop)*span/rr(irs6))/log(ratio))
   allocate(cidx(kk+2),cfill(kk+2),cmis(kk+2),cfit(kk+2),corder(kk+2))
   ncand=0
   do ii=0,kk+1
     if(ii<=kk) then
       rad=rr(irtop)/ratio**ii
       jj=nint(log(rad/rr(1))/al)+1
       jj=max(4,min(irtop,jj))
     else
       jj=irs6
     end if
     if(any(cidx(1:ncand)==jj)) cycle
     ncand=ncand+1
     cidx(ncand)=jj
   end do
! largest first (irs6 is among them)
   do ii=2,ncand
     jj=cidx(ii)
     kk=ii-1
     do while(kk>=1)
       if(cidx(kk)>=jj) exit
       cidx(kk+1)=cidx(kk)
       kk=kk-1
     end do
     cidx(kk+1)=jj
   end do

   allocate(d2ae(nprobe,nprobe),d2mod(nprobe,nprobe),d2zero(nprobe,nprobe))
   d2zero(:,:)=0.0d0
   call der2exc(rhotps,rhoc,probes,rr,d2ae,d2zero,mismatch, &
&               zion,iexc,0,nprobe,lprobe,mmax,mmax)
   do ii=1,ncand
     call model(cidx(ii))
     cfill(ii)=ifill
     cfit(ii)='-'
     cmis(ii)=huge(1.0d0)
     if(ifill==0) cycle
     call der2exc(rhotps,core,probes,rr,d2mod,d2ae,cmis(ii), &
&                 zion,iexc,0,nprobe,lprobe,mmax,mmax)
   end do

! the order of the mismatch, stable (ties to the larger radius)
   do ii=1,ncand
     corder(ii)=ii
   end do
   do ii=2,ncand
     jj=corder(ii)
     kk=ii-1
     do while(kk>=1)
       if(cmis(corder(kk))<=cmis(jj)) exit
       corder(kk+1)=corder(kk)
       kk=kk-1
     end do
     corder(kk+1)=jj
   end do

   irs=0
   do ii=1,ncand
     jj=corder(ii)
     if(cfill(jj)==0) exit
     if(fits(cidx(jj),nq)) then
       cfit(jj)='T'
       irs=cidx(jj)
       mismatch=cmis(jj)
       exit
     end if
     cfit(jj)='F'
   end do
   if(irs==0) call no_model()
 end if

! the cutoff the model needs: (q+qstep)**2/8 for the last density cutoff q
! above errmax, which is below 2*sqrt(2*budget)
 call model(irs)
 ibad=0
 if(ifrom>1) then
   allocate(rqs(nx,ifrom-1))
   call scan_errors(1,0,huge(1.0d0),ii)
   do ii=ifrom-1,1,-1
     err(ii)=xc_error(rqs(:,ii))
     if(err(ii)>errmax) then
       ibad=ii
       exit
     end if
   end do
   deallocate(rqs)
 end if
 if(ibad==0) then
   ecore=0.0d0
 else
   ecore=((ibad+1)*qstep)**2/8.0d0
 end if

! fill in model and derivatives, all-electron beyond r_s
 rhomod(:,:)=dc(:,:)
 do ii=1,irs-1
   if(ifill==1) then
     call tail_derivatives(rr(ii),rr(irs),rhomod(ii,:))
   else
     call martins_derivatives(rr(ii),rhomod(ii,:))
   end if
 end do

 if(icmod==7) then
   write(6,'(/a)') 'Model core correction: hardness-matched core (icmod 7)'
 else
   write(6,'(/a)') 'Model core correction: budgeted core (icmod 6)'
 end if
 if(fcfact<=0.0d0) then
   write(6,'(a,f12.4,a,f8.4)') '  pseudo-wave cutoff',ecutw, &
&        ' Ha, times',factor
 end if
 write(6,'(a,f12.4,a,1p,e10.2,a)') '  cutoff budget',ecut, &
&      ' Ha for an xc energy error of',errmax,' Ha'
 write(6,'(a,f8.4,a,f12.4,a)') '  slack',slack, &
&      ', the model core must fit',budget,' Ha'
 write(6,'(a,f12.6,a,f10.6)') '  model radius r_s',rr(irs), &
&      '  smallest rc',rr(irmax)
 if(irtop<irmax) then
   write(6,'(a,f12.6)') '  last radius with core charge',rr(irtop)
 end if
 if(ifill==1) then
   write(6,'(a)') '  fill minimum tail'
   write(6,'(a,1p,9e16.8)') '  coefficients c0..c8',cpoly(:)
 else
   write(6,'(a)') '  fill Martins'
   write(6,'(a,1p,3e16.8)') '  a, b, c',aco,bco,cco
 end if
 write(6,'(a,f12.4,a)') '  cutoff the model core needs',ecore,' Ha'
 if(icmod==7) then
   write(6,'(a,f12.6)') '  smallest fitting radius (icmod 6)',rr(irs6)
   write(6,'(a,1p,e12.4,a)') '  hardness mismatch',mismatch,' Ha'
   write(6,'(a)') '  candidates: r_s, fill, hardness mismatch (Ha), fits'
   do ii=1,ncand
     if(cfill(ii)==0) then
       write(6,'(f12.6,2x,a2,a12,2x,a1)') rr(cidx(ii)),fillname(0), &
&            '-',cfit(ii)
     else
       write(6,'(f12.6,2x,a2,1p,e12.4,2x,a1)') rr(cidx(ii)), &
&            fillname(cfill(ii)),cmis(ii),cfit(ii)
     end if
   end do
   write(6,'(a)') ''
   deallocate(cidx,cfill,cmis,cfit,corder,d2ae,d2mod,d2zero)
 end if

 deallocate(dc,core,sn,cs,sh,ch,err,qn,wn)
 deallocate(nval,rq,work,eps,vxc)

 return

 contains

 logical function fits(is,ilast)

!! whether the model of radius rr(is) exists and has errors below errmax
!! at the density cutoffs ifrom..ilast
!!
!! @param is  rr index of the model radius
!! @param ilast  last density cutoff index tested

!Input variables
 integer, intent(in) :: is,ilast

!Local variables
 integer :: ifail
 call model(is)
 fits=.false.
 if(ifill==0) return
 call scan_errors(ifrom,ilast,errmax,ifail)
 fits = ifail==0
 return
 end function fits

 integer function bisect(ilo,ihi,ilast)

!! the smallest index above ilo, which fails, that fits up to ilast;
!! ihi fits
!!
!! @param ilo  lower rr index of the bracket, a model radius that fails
!! @param ihi  upper rr index of the bracket, a model radius that fits
!! @param ilast  last density cutoff index tested

!Input variables
 integer, intent(in) :: ilo,ihi,ilast

!Local variables
 integer :: ia,ib,imid
 ia=ilo
 ib=ihi
 if(fits(ia,ilast)) then
   bisect=ia
   return
 end if
 do while(ib-ia>1)
   imid=(ia+ib)/2
   if(fits(imid,ilast)) then
     ib=imid
   else
     ia=imid
   end if
 end do
 bisect=ib
 return
 end function bisect

 subroutine no_model()
!! reports that no model radius fits and stops
 write(6,'(/a,f10.6,a,f10.4,a,1p,e10.2,a)') &
&      'modcore_budget: ERROR stop - no model core inside rc =',rr(irtop), &
&      ' needs less than',budget,' Ha for an error of',errmax,' Ha'
 error stop 1
 end subroutine no_model

 subroutine model(is)

!! the fill of radius rr(is) (ifill 1 for the minimum tail in cpoly, 2 for
!! Martins' in aco, bco, cco, 0 for none), and in core the model inside
!! rr(is) and the core charge beyond
!!
!! @param is  rr index of the model radius

!Input variables
 integer, intent(in) :: is

!Local variables
 real(dp) :: rs,l1,l2,tt,vv,vlast
 integer :: kk
 rs=rr(is)
 core(:)=rhoc(:)

 call min_tail(is)
 ifill=1
 vlast=huge(1.0d0)
 do kk=0,nmono-1
   tt=(dble(kk)/dble(nmono-1))**2
   vv=poly(tt)
   if(vv>vlast) then
     ifill=0
     exit
   end if
   vlast=vv
 end do
 if(ifill==1) then
   do kk=1,is-1
     core(kk)=poly((rr(kk)/rs)**2)
   end do
   return
 end if

 l1=dc(is,2)/dc(is,1)
 l2=dc(is,3)/dc(is,1) - l1**2
 cco=(rs*l2 - l1)/(8.0d0*rs**3)
 bco=(l1 - 4.0d0*cco*rs**3)/(2.0d0*rs)
 aco=log(dc(is,1)) - bco*rs**2 - cco*rs**4
 if(bco>0.0d0 .or. bco + 2.0d0*cco*rs**2>0.0d0) return
 ifill=2
 do kk=1,is-1
   core(kk)=exp(aco + bco*rr(kk)**2 + cco*rr(kk)**4)
 end do
 return
 end subroutine model

 real(dp) function poly(tt)

!! the minimum tail at (r/r_s)**2 = tt
!!
!! @param tt  (r/r_s)**2

!Input variables
 real(dp), intent(in) :: tt

!Local variables
 integer :: kk
 poly=cpoly(npow)
 do kk=npow-1,0,-1
   poly=poly*tt + cpoly(kk)
 end do
 return
 end function poly

 subroutine min_tail(is)

!! the coefficients cpoly of the minimum tail of radius rr(is): the least
!! squares of W (A c + b) with C c = d (LAPACK dgglse), where A c is the
!! transform of the polynomial at the window nodes (Gauss-Legendre on
!! [0, r_s]), b that of the core charge beyond r_s, W the square root of
!! the weights, and C c = d the match of the core charge and its first two
!! derivatives at r_s.  b is Gauss-Legendre on each grid interval, of the
!! 4-point Lagrange interpolant of r**3 rhoc in the grid index: the fill
!! depends on b to its last digits (the trapezoidal rule of the grid, off
!! by 3e-4 at q = 17 bohr**-1, moves r_s by 10% on Si)
!!
!! @param is  rr index of the model radius

!Input variables
 integer, intent(in) :: is

!Local variables
 integer, parameter :: ngi=4
 real(dp), allocatable :: xx(:),wx(:),amat(:,:),cvec(:),wk(:)
 real(dp), allocatable :: rt(:),gt(:)
 real(dp) :: bmat(3,0:npow),dvec(3),rs,tt,qx,jx,ww,tg(ngi),wg(ngi),fl(4)
 integer :: nxg,ii,jj,kk,ll,info,lwork,iend,nt
 rs=rr(is)
 nxg=20+ceiling(qhi*rs)
 lwork=nqn+npow+1+64*(nqn+npow+1)
 allocate(xx(nxg),wx(nxg),amat(nqn,0:npow),cvec(nqn),wk(lwork))
 call gauss_legendre(0.0d0,rs,nxg,xx,wx)

! the nodes beyond r_s, to the last grid interval where r**3 rhoc is above
! 1e-14 of its largest value: r**3 rhoc interpolated in the grid index, the
! weights al dx
 tt=maxval(rr(is:mmax)**3*rhoc(is:mmax))
 iend=is
 do ii=is,mmax-2
   if(rr(ii)**3*rhoc(ii)>1.0d-14*tt) iend=ii+1
 end do
 call gauss_legendre(0.0d0,1.0d0,ngi,tg,wg)
 allocate(rt(ngi*(iend-is)),gt(ngi*(iend-is)))
 nt=0
 do ii=is,iend-1
   do kk=1,ngi
     tt=tg(kk)
! Lagrange weights of the points ii-1..ii+2 at ii+tt
     fl(1)=-tt*(tt-1.0d0)*(tt-2.0d0)/6.0d0
     fl(2)=(tt+1.0d0)*(tt-1.0d0)*(tt-2.0d0)/2.0d0
     fl(3)=-(tt+1.0d0)*tt*(tt-2.0d0)/2.0d0
     fl(4)=(tt+1.0d0)*tt*(tt-1.0d0)/6.0d0
     nt=nt+1
     rt(nt)=rr(ii)*exp(al*tt)
     gt(nt)=0.0d0
     do ll=1,4
       gt(nt)=gt(nt) + fl(ll)*rr(ii+ll-2)**3*rhoc(ii+ll-2)
     end do
     gt(nt)=al*wg(kk)*gt(nt)
   end do
 end do

 amat(:,:)=0.0d0
 do ii=1,nqn
   do jj=1,nxg
     qx=qn(ii)*xx(jj)
     jx=wx(jj)*xx(jj)**2*j0(qx)
     tt=(xx(jj)/rs)**2
     do kk=0,npow
       amat(ii,kk)=amat(ii,kk) + jx
       jx=jx*tt
     end do
   end do
! the core charge beyond r_s
   ww=0.0d0
   do jj=1,nt
     ww=ww + gt(jj)*j0(qn(ii)*rt(jj))
   end do
   cvec(ii)=-sqrt(wn(ii))*ww
   amat(ii,:)=sqrt(wn(ii))*amat(ii,:)
 end do

! the match, row m scaled by rs**m: (2k)!/(2k-m)!
 do kk=0,npow
   bmat(1,kk)=1.0d0
   bmat(2,kk)=dble(2*kk)
   bmat(3,kk)=dble(2*kk*(2*kk-1))
 end do
 dvec(:)=[dc(is,1),dc(is,2)*rs,dc(is,3)*rs**2]

 call dgglse(nqn,npow+1,3,amat,nqn,bmat,3,cvec,dvec,cpoly,wk,lwork,info)
 if(info/=0) then
   write(6,'(/a,i4)') 'modcore_budget: ERROR stop - dgglse info =',info
   error stop 1
 end if
 deallocate(xx,wx,amat,cvec,wk,rt,gt)
 return
 end subroutine min_tail

 real(dp) function j0(xx)

!! the spherical Bessel function j0, by its series near 0
!!
!! @param xx  argument

!Input variables
 real(dp), intent(in) :: xx
 if(abs(xx)<1.0d-3) then
   j0=1.0d0 - xx**2/6.0d0 + xx**4/120.0d0
 else
   j0=sin(xx)/xx
 end if
 return
 end function j0

 subroutine tail_derivatives(rad,rs,dd)

!! the minimum tail and its first 4 derivatives at rad <= rs:
!! d**m/dr**m (r/rs)**(2k) = (2k)!/(2k-m)! r**(2k-m)/rs**(2k)
!!
!! @param rad  radius (bohr)
!! @param rs  model radius r_s (bohr)
!! @param dd  the model and its first 4 derivatives at rad

!Input variables
 real(dp), intent(in) :: rad,rs

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 real(dp) :: ff
 integer :: kk,mm,jj
 dd(:)=0.0d0
 do kk=0,npow
   do mm=0,min(4,2*kk)
     ff=1.0d0
     do jj=0,mm-1
       ff=ff*dble(2*kk-jj)
     end do
     dd(mm+1)=dd(mm+1) + cpoly(kk)*ff*rad**(2*kk-mm)/rs**(2*kk)
   end do
 end do
 return
 end subroutine tail_derivatives

 subroutine martins_derivatives(rad,dd)

!! the model and its first 4 derivatives at rad, from n' = p' n with
!! p = a + b r**2 + c r**4: n(k) = sum_j C(k-1,j) p(j+1) n(k-1-j)
!!
!! @param rad  radius (bohr)
!! @param dd  the model and its first 4 derivatives at rad

!Input variables
 real(dp), intent(in) :: rad

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 real(dp) :: dp1(4)
 integer :: kk,jj
 integer, parameter :: binom(0:3,0:3)=reshape([1,1,1,1, 0,1,2,3, &
&                                              0,0,1,3, 0,0,0,1],[4,4])
 dp1=[2.0d0*bco*rad + 4.0d0*cco*rad**3, 2.0d0*bco + 12.0d0*cco*rad**2, &
&     24.0d0*cco*rad, 24.0d0*cco]
 dd(1)=exp(aco + bco*rad**2 + cco*rad**4)
 do kk=1,4
   dd(kk+1)=0.0d0
   do jj=0,kk-1
     dd(kk+1)=dd(kk+1) + binom(kk-1,jj)*dp1(jj+1)*dd(kk-jj)
   end do
 end do
 return
 end subroutine martins_derivatives

 subroutine scan_errors(ifirst,ilast,elimit,ifail)

!! the errors err(iq) at the density cutoffs iq*qstep, iq=ifirst..ilast, of
!! the model in core, up to the first above elimit (its index in ifail,
!! else 0).  With ilast=0 instead the band-limited models of the cutoffs
!! 1..ifrom-1 are kept in rqs.
!! The band-limited model is built up panel by panel (midpoints, width dq),
!! the sines and cosines of the midpoints by rotation
!!
!! @param ifirst  first density cutoff index whose error is tested
!! @param ilast  last density cutoff index, or 0 to keep the band-limited models
!! @param elimit  error limit (Ha)
!! @param ifail  index of the first cutoff with an error above elimit, else 0

!Input variables
 integer, intent(in) :: ifirst,ilast
 real(dp), intent(in) :: elimit

!Output variables
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: qq,ww,tt
 integer :: iq,kk,ll,iend

 ifail=0
 err(:)=-1.0d0
 e0=exc_sphere(core(1:nx))
 iend=ilast
 if(ilast==0) iend=ifrom-1
 sn(:)=sin(0.5d0*dq*rr(:))
 cs(:)=cos(0.5d0*dq*rr(:))
 sh(:)=sin(dq*rr(:))
 ch(:)=cos(dq*rr(:))
 rq(:)=0.0d0
 do iq=1,iend
   do kk=1,npanel
     qq=((iq-1)*npanel + kk - 0.5d0)*dq
! weight dq (1/2pi**2) int 4pi r**2 nc j0(qr) dr of the panel
     ww=0.0d0
     do ll=1,mmax
       ww=ww + rr(ll)**2*core(ll)*sn(ll)
     end do
     ww=dq*al*ww*2.0d0/(pi*qq)
! q**2 j0(qr) = q sin(qr)/r
     do ll=1,nx
       rq(ll)=rq(ll) + ww*qq*sn(ll)/rr(ll)
     end do
     do ll=1,mmax
       tt=sn(ll)*ch(ll) + cs(ll)*sh(ll)
       cs(ll)=cs(ll)*ch(ll) - sn(ll)*sh(ll)
       sn(ll)=tt
     end do
   end do
   if(ilast==0) then
     rqs(:,iq)=rq(:)
     cycle
   end if
   if(iq<ifirst) cycle
   err(iq)=xc_error(rq)
   if(.not. err(iq)<=elimit) then
     ifail=iq
     return
   end if
 end do
 return
 end subroutine scan_errors

 real(dp) function xc_error(bq)

!! the error of the band-limited model bq of the model in core
!!
!! @param bq  band-limited model core charge on rr(1:nx)

!Input variables
 real(dp), intent(in) :: bq(nx)

!Local variables
 real(dp) :: tt
 integer :: ll
 tt=0.0d0
 do ll=1,nin
   tt=tt + rr(ll)**3*(bq(ll) - core(ll))
 end do
 xc_error=abs(exc_sphere(bq) - e0 - vbg*al*tt)
 return
 end function xc_error

 real(dp) function exc_sphere(cc)

!! the exchange-correlation energy within the sphere of the valence charge
!! on the background plus the core charge cc
!!
!! @param cc  core charge on rr(1:nx)

!Input variables
 real(dp), intent(in) :: cc(nx)

!Local variables
 integer :: ll
 work(:)=nval(:) + cc(:)
 call xc_eps(work,eps,vxc)
 exc_sphere=0.0d0
 do ll=1,nin
   exc_sphere=exc_sphere + rr(ll)**3*work(ll)*eps(ll)
 end do
 exc_sphere=al*exc_sphere
 return
 end function exc_sphere

 subroutine xc_eps(rho,exc,vx)

!! the xc energy per electron and potential of the charge rho on rr(1:nx)
!!
!! @param rho  charge density on rr(1:nx)
!! @param exc  exchange-correlation energy per electron
!! @param vx  exchange-correlation potential

!Input/output variables
 real(dp), intent(inout) :: rho(nx)

!Output variables
 real(dp), intent(out) :: exc(nx),vx(nx)
 if(iexc == 1) then
   call excwig(rho,vx,exc,nx)
 else if(iexc == 2) then
   call exchdl(rho,vx,exc,nx)
 else if(iexc == 3) then
   call excpzca(rho,vx,exc,nx)
 else if(iexc == 4) then
   call excggc(rho,vx,exc,rr,nx)
 else if(iexc < 0) then
   call exc_libxc(iexc,al,rho,vx,exc,rr,nx)
 else
   write(6,'(/a,i4)') 'modcore_budget: ERROR bad input iexc =',iexc
   error stop 1
 end if
 return
 end subroutine xc_eps

 subroutine gauss_legendre(xa,xb,nn,xx,ww)

!! the nn-point Gauss-Legendre nodes and weights on [xa, xb], the roots of
!! P_nn by Newton's method
!!
!! @param xa  lower end
!! @param xb  upper end
!! @param nn  number of nodes
!! @param xx  nodes
!! @param ww  weights

!Input variables
 real(dp), intent(in) :: xa,xb
 integer, intent(in) :: nn

!Output variables
 real(dp), intent(out) :: xx(nn),ww(nn)

!Local variables
 real(dp) :: zz,z1,p1,p2,p3,pp
 integer :: ii,jj,it
 do ii=1,(nn+1)/2
   zz=cos(pi*(ii-0.25d0)/(nn+0.5d0))
   do it=1,100
     p1=1.0d0
     p2=0.0d0
     do jj=1,nn
       p3=p2
       p2=p1
       p1=((2*jj-1)*zz*p2 - (jj-1)*p3)/jj
     end do
     pp=nn*(zz*p1 - p2)/(zz**2 - 1.0d0)
     z1=zz
     zz=z1 - p1/pp
     if(abs(zz-z1)<1.0d-15) exit
   end do
   xx(ii)=0.5d0*(xa+xb) - 0.5d0*(xb-xa)*zz
   xx(nn+1-ii)=0.5d0*(xa+xb) + 0.5d0*(xb-xa)*zz
   ww(ii)=(xb-xa)/((1.0d0 - zz**2)*pp**2)
   ww(nn+1-ii)=ww(ii)
 end do
 return
 end subroutine gauss_legendre

 end subroutine modcore_budget

 subroutine modcore_probes(lmax,lloc,nproj,epa,vp,vkb,evkb,rr,mmax,mxprj, &
&                         nkap,probes,lprobe,nprobe)

!! the probe densities of the hardness of icmod=7: for every projector
!! iprj of every l, the eigenstate with iprj-1 nodes of the screened
!! pseudo-Hamiltonian, as for run_diag, at unit occupation (the average
!! over J = l +/- 1/2 for nkap=2), if it is bound.  The occupied valence
!! states always are; an unoccupied one counts if bound, so that an atom
!! with one valence state still has a matrix to match
!!
!! @param lmax  maximum angular momentum
!! @param lloc  index of local potential
!! @param nproj  number of projectors for each l
!! @param epa  reference energies of the projectors
!! @param vp  screened semi-local and local pseudopotentials
!! @param vkb  Vanderbilt-Kleinman-Bylander projectors
!! @param evkb  coefficients of the projectors
!! @param rr  log radial grid
!! @param mmax  dimension of log grid
!! @param mxprj  dimension of number of projector
!! @param nkap  1 for scalar-relativistic or non-relativistic, 2 for
!!              relativistic
!! @param probes  probe densities, (u/r)**2 normalized to one electron
!! @param lprobe  angular momenta of the probes
!! @param nprobe  number of probes

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,mxprj,nkap
 integer, intent(in) :: nproj(6)
 real(dp), intent(in) :: rr(mmax),vp(mmax,5,nkap)
 real(dp), intent(in) :: vkb(mmax,mxprj,4,nkap)
 real(dp), intent(in) :: epa(mxprj,6,nkap),evkb(mxprj,4,nkap)

!Output variables
 real(dp), intent(out) :: probes(mmax,mxprj*(lmax+1))
 integer, intent(out) :: lprobe(mxprj*(lmax+1))
 integer, intent(out) :: nprobe

!Local variables
 real(dp), allocatable :: uu(:),up(:),dens(:)
 real(dp) :: et,emin,emax,fj
 integer :: l1,ll,iprj,ikap,mkap,npr,ierr,mch
 logical :: bound

 allocate(uu(mmax),up(mmax),dens(mmax))
 nprobe=0
 do l1=1,lmax+1
   ll=l1-1
   mkap=1
   if(nkap==2 .and. l1>1) mkap=2
   npr=nproj(l1)
   if(l1==lloc+1) npr=0
   do iprj=1,max(1,nproj(l1))
     dens(:)=0.0d0
     bound=.true.
     do ikap=1,mkap
       fj=1.0d0
       if(mkap==2 .and. ikap==1) fj=dble(ll+1)/dble(2*ll+1)
       if(mkap==2 .and. ikap==2) fj=dble(ll)/dble(2*ll+1)
       et=epa(iprj,l1,ikap)
       if(epa(iprj,l1,1)<0.0d0) then
         emax=min(et+0.1d0*abs(et)+0.05d0,0.0d0)
         emin=et-0.1d0*abs(et)-0.05d0
       else
! a scattering reference: its state, if bound, below zero
         emax=0.0d0
         emin=minval(vp(:,lloc+1,ikap))
         et=0.5d0*emin
       end if
       call lschvkbb(ll+iprj,ll,npr,ierr,et,emin,emax, &
&                    rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap), &
&                    evkb(1,l1,ikap),uu,up,mmax,mch)
       if(ierr/=0 .or. et>=0.0d0) then
         bound=.false.
         exit
       end if
       dens(:)=dens(:) + fj*(uu(:)/rr(:))**2
     end do
     if(.not. bound) cycle
     nprobe=nprobe+1
     probes(:,nprobe)=dens(:)
     lprobe(nprobe)=ll
   end do
 end do
 deallocate(uu,up,dens)
 return
 end subroutine modcore_probes

 subroutine modcore_budget_ecut(ecurve,wcurve,ncurve,nprof,errmax,ecut)

!! the cutoff ecut (Ha) the pseudo-waves need for the error errmax (Ha per
!! atom): the root of the total residual kinetic energy sum_i w_i E_i(ecut),
!! with E_i the error per electron of pseudo-wave i and w_i the electrons
!! that fill it, the ecut_for_error of PseudoPotentialGenerators.jl.
!! ecurve(:,:,i) holds the cutoffs and E_i of run_optimize on its mesh,
!! zero beyond its end; wcurve(i)=0 skips it.  log E is interpolated
!! linearly in q between the mesh points, as for the convergence profiles
!!
!! @param ecurve  residual kinetic energy curves: ecurve(1,:,i) cutoffs (Ha),
!!                ecurve(2,:,i) error per electron E_i (Ha) of pseudo-wave i
!! @param wcurve  electrons filling each pseudo-wave (0 skips it)
!! @param ncurve  number of points per curve
!! @param nprof  number of curves (pseudo-waves)
!! @param errmax  error (Ha per atom) that defines the needed cutoff
!! @param ecut  cutoff (Ha) the pseudo-waves need for errmax

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: ncurve,nprof
 real(dp), intent(in) :: ecurve(2,ncurve,nprof),wcurve(nprof),errmax

!Output variables
 real(dp), intent(out) :: ecut

!Local variables
 real(dp), allocatable :: total(:),qq(:)
 real(dp) :: qt
 integer :: ip,kk,nend,klast

 allocate(total(ncurve),qq(ncurve))
 total(:)=0.0d0
 qq(:)=0.0d0
 nend=0
 do ip=1,nprof
   if(wcurve(ip)<=0.0d0) cycle
   kk=ncurve
   do while(kk>1 .and. ecurve(1,kk,ip)<=0.0d0)
     kk=kk-1
   end do
   if(nend==0 .or. kk<nend) nend=kk
   total(:)=total(:)+wcurve(ip)*ecurve(2,:,ip)
   qq(:)=sqrt(2.0d0*ecurve(1,:,ip))
 end do
 if(nend<2) then
   write(6,'(/a)') 'modcore_budget_ecut: ERROR stop - no residual energy curves'
   error stop 1
 end if
 if(total(nend)>errmax) then
   write(6,'(/a,1p,e10.2,a,0p,f10.2,a)') &
&        'modcore_budget_ecut: ERROR stop - residual kinetic energy above', &
&        errmax,' Ha at',0.5d0*qq(nend)**2,' Ha, the end of its tables'
   error stop 1
 end if

 klast=0
 do kk=1,nend-1
   if(total(kk)>errmax) klast=kk
 end do
 if(klast==0) then
   ecut=0.0d0
 else if(total(klast+1)>0.0d0) then
   qt=qq(klast)+(qq(klast+1)-qq(klast)) &
&    *log(total(klast)/errmax)/log(total(klast)/total(klast+1))
   ecut=0.5d0*qt**2
 else
   qt=qq(klast)+(qq(klast+1)-qq(klast)) &
&    *(total(klast)-errmax)/(total(klast)-total(klast+1))
   ecut=0.5d0*qt**2
 end if

 deallocate(total,qq)
 return
 end subroutine modcore_budget_ecut
