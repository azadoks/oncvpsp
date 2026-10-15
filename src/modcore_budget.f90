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
! Creates the model core charge of icmod=6: the all-electron core charge
! wherever the plane-wave cutoff affords it.  Inside the radius r_s the
! model is Martins' exp(a + b r**2 + c r**4), matching the core charge and
! its first two derivatives at r_s; beyond r_s it is the core charge.
! r_s is the smallest grid radius, no larger than the smallest rc, whose
! model needs no more than the cutoff ecut for an exchange-correlation
! energy error below errmax.
!
! The error of a density cutoff q (a plane-wave cutoff q**2/8) is that of
! the core charge without its Fourier components beyond q, in a solid:
!
!   |Exc[n + nc_q] - Exc[n + nc] - Vxc[n_bg] * int (nc_q - nc)|
!
! within the sphere r <= 10 bohr, with n the pseudo valence charge on a
! uniform background n_bg = 5e-3 bohr**-3 that stands for the interstitial
! electrons of the neighbours.  The model needs the cutoff ecut if the
! error stays below errmax for every density cutoff from 2*sqrt(2*ecut) on,
! in steps of 0.1 bohr**-1 up to 40 bohr**-1 (200 Ha).  This is the
! BudgetedCore of PseudoPotentialGenerators.jl with its Martins candidate.

 subroutine modcore_budget(rhotps,rhoc,rhomod,ecut,errmax,irmax,mmax,rr,iexc)

!rhotps  total pseudocharge density
!rhoc  core-charge density
!rhomod  model core density and 4 derivatives
!ecut  plane-wave cutoff (Ha) the model must not need more than
!errmax  exchange-correlation energy error (Ha) that defines the needed cutoff
!irmax  rr index of the largest r_s (smallest rc)
!mmax  dimension of log grid
!rr log radial grid
!iexc  exchange-correlation function to be used

 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: pi4=4.0_dp*pi

!Input variables
 integer :: mmax,irmax,iexc
 real(dp) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp) :: ecut,errmax

!Output variables
 real(dp) :: rhomod(mmax,5)

!background density (bohr**-3) and radius (bohr) of the error,
!largest density cutoff and its step (bohr**-1),
!panels of the band-limited transform per step
 real(dp), parameter :: background=5.0d-3, radius=10.0d0
 real(dp), parameter :: qmax=40.0d0, qstep=0.1d0
 integer, parameter :: npanel=4

!Local variables
 real(dp) :: al,dq,vbg,ecore,aco,bco,cco
 real(dp), allocatable :: dc(:,:),core(:),nval(:),sn(:),cs(:),sh(:),ch(:)
 real(dp), allocatable :: rq(:),work(:),eps(:),vxc(:),err(:)
 integer :: ii,jj,nq,nin,nx,ilo,ihi,imid,irs,ifrom,ibad

 al = 0.01d0 * dlog(rr(101) / rr(1))
 dq = qstep/npanel
 nq = nint(qmax/qstep)

 if(2.0d0*sqrt(2.0d0*ecut) >= qmax) then
   write(6,'(/a,f10.2,a)') 'modcore_budget: ERROR stop - ecut beyond', &
&        qmax**2/8.0d0,' Ha'
   stop
 end if

 if(.not. rhoc(irmax)>0.0d0) then
   write(6,'(/a)') 'modcore_budget: ERROR stop - no core charge at rc'
   stop
 end if

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

! the first density cutoff that decides whether the model needs no more
! than ecut: the needed cutoff is (q+qstep)**2/8 for the last q above errmax
 ifrom=nq
 do ii=1,nq
   if(((ii+1)*qstep)**2/8.0d0 > ecut) then
     ifrom=ii
     exit
   end if
 end do

! smallest r_s on the grid, by bisection in the index between rc/50 and rc
! (the needed cutoff falls as the model softens)
 ihi=irmax
 ilo=irmax
 do ii=irmax,1,-1
   if(rr(ii)<rr(irmax)/50.0d0) exit
   ilo=ii
 end do

 if(.not. fits(ihi)) then
   write(6,'(/a,f10.6,a,f10.4,a,1p,e10.2,a)') &
&        'modcore_budget: ERROR stop - no model core inside rc =',rr(irmax), &
&        ' needs less than',ecut,' Ha for an error of',errmax,' Ha'
   stop
 end if
 if(fits(ilo)) then
   irs=ilo
 else
   do while(ihi-ilo>1)
     imid=(ilo+ihi)/2
     if(fits(imid)) then
       ihi=imid
     else
       ilo=imid
     end if
   end do
   irs=ihi
 end if

! the cutoff the model needs, from the full scan
 call model(irs)
 call scan_errors(1,huge(1.0d0),ibad)
 ibad=0
 do ii=1,nq
   if(err(ii)>errmax) ibad=ii
 end do
 if(ibad==0) then
   ecore=0.0d0
 else if(ibad==nq) then
   ecore=huge(1.0d0)
 else
   ecore=((ibad+1)*qstep)**2/8.0d0
 end if

! fill in model and derivatives, all-electron beyond r_s
 rhomod(:,:)=dc(:,:)
 do ii=1,irs-1
   call martins_derivatives(rr(ii),rhomod(ii,:))
 end do

 write(6,'(/a)') 'Model core correction: budgeted Martins core (icmod 6)'
 write(6,'(a,f12.4,a,1p,e10.2,a)') '  cutoff budget',ecut, &
&      ' Ha for an xc energy error of',errmax,' Ha'
 write(6,'(a,f12.6,a,f10.6)') '  model radius r_s',rr(irs), &
&      '  smallest rc',rr(irmax)
 write(6,'(a,1p,3e16.8)') '  a, b, c',aco,bco,cco
 if(ecore<huge(1.0d0)) then
   write(6,'(a,f12.4,a)') '  cutoff the model core needs',ecore,' Ha'
 end if

 deallocate(dc,core,sn,cs,sh,ch,err)
 deallocate(nval,rq,work,eps,vxc)

 return

 contains

! whether the model of radius rr(is) needs no more than ecut
 logical function fits(is)
   integer, intent(in) :: is
   integer :: ifail
   call model(is)
   call scan_errors(ifrom,errmax,ifail)
   fits = ifail==0
 end function fits

! the Martins coefficients at rr(is), and in core the model inside rr(is)
! and the core charge beyond
 subroutine model(is)
   integer, intent(in) :: is
   real(dp) :: rs,l1,l2
   integer :: kk
   rs=rr(is)
   l1=dc(is,2)/dc(is,1)
   l2=dc(is,3)/dc(is,1) - l1**2
   cco=(rs*l2 - l1)/(8.0d0*rs**3)
   bco=(l1 - 4.0d0*cco*rs**3)/(2.0d0*rs)
   aco=log(dc(is,1)) - bco*rs**2 - cco*rs**4
   core(:)=rhoc(:)
   do kk=1,is-1
     core(kk)=exp(aco + bco*rr(kk)**2 + cco*rr(kk)**4)
   end do
 end subroutine model

! the model and its first 4 derivatives at r, from n' = p' n with
! p = a + b r**2 + c r**4: n(k) = sum_j C(k-1,j) p(j+1) n(k-1-j)
 subroutine martins_derivatives(r,d)
   real(dp), intent(in) :: r
   real(dp), intent(out) :: d(5)
   real(dp) :: dp1(4)
   integer :: kk,jj
   integer, parameter :: binom(0:3,0:3)=reshape([1,1,1,1, 0,1,2,3, 0,0,1,3, 0,0,0,1],[4,4])
   dp1=[2.0d0*bco*r + 4.0d0*cco*r**3, 2.0d0*bco + 12.0d0*cco*r**2, &
&       24.0d0*cco*r, 24.0d0*cco]
   d(1)=exp(aco + bco*r**2 + cco*r**4)
   do kk=1,4
     d(kk+1)=0.0d0
     do jj=0,kk-1
       d(kk+1)=d(kk+1) + binom(kk-1,jj)*dp1(jj+1)*d(kk-jj)
     end do
   end do
 end subroutine martins_derivatives

! the errors err(iq) at the density cutoffs iq*qstep, iq=ifirst..nq, of the
! model in core, up to the first above elimit (its index in ifail, else 0).
! The band-limited model is built up panel by panel (midpoints, width dq),
! the sines and cosines of the midpoints by rotation
 subroutine scan_errors(ifirst,elimit,ifail)
   integer, intent(in) :: ifirst
   real(dp), intent(in) :: elimit
   integer, intent(out) :: ifail
   real(dp) :: e0,q,ww,tt
   integer :: iq,kk,ll

   ifail=0
   err(:)=-1.0d0
   e0=exc_sphere(core(1:nx))
   sn(:)=sin(0.5d0*dq*rr(:))
   cs(:)=cos(0.5d0*dq*rr(:))
   sh(:)=sin(dq*rr(:))
   ch(:)=cos(dq*rr(:))
   rq(:)=0.0d0
   do iq=1,nq
     do kk=1,npanel
       q=((iq-1)*npanel + kk - 0.5d0)*dq
! weight dq (1/2pi**2) int 4pi r**2 nc j0(qr) dr of the panel
       ww=0.0d0
       do ll=1,mmax
         ww=ww + rr(ll)**2*core(ll)*sn(ll)
       end do
       ww=dq*al*ww*2.0d0/(pi*q)
! q**2 j0(qr) = q sin(qr)/r
       do ll=1,nx
         rq(ll)=rq(ll) + ww*q*sn(ll)/rr(ll)
       end do
       do ll=1,mmax
         tt=sn(ll)*ch(ll) + cs(ll)*sh(ll)
         cs(ll)=cs(ll)*ch(ll) - sn(ll)*sh(ll)
         sn(ll)=tt
       end do
     end do
     if(iq<ifirst) cycle
     tt=0.0d0
     do ll=1,nin
       tt=tt + rr(ll)**3*(rq(ll) - core(ll))
     end do
     err(iq)=abs(exc_sphere(rq) - e0 - vbg*al*tt)
     if(.not. err(iq)<=elimit) then
       ifail=iq
       return
     end if
   end do
 end subroutine scan_errors

! the exchange-correlation energy within the sphere of the valence charge
! on the background plus the core charge cc
 real(dp) function exc_sphere(cc)
   real(dp), intent(in) :: cc(nx)
   integer :: ll
   work(:)=nval(:) + cc(:)
   call xc_eps(work,eps,vxc)
   exc_sphere=0.0d0
   do ll=1,nin
     exc_sphere=exc_sphere + rr(ll)**3*work(ll)*eps(ll)
   end do
   exc_sphere=al*exc_sphere
 end function exc_sphere

! the xc energy per electron and potential of the charge rho on rr(1:nx)
 subroutine xc_eps(rho,exc,vx)
   real(dp), intent(inout) :: rho(nx)
   real(dp), intent(out) :: exc(nx),vx(nx)
   if(iexc .eq. 1) then
     call excwig(rho,vx,exc,nx)
   else if(iexc .eq. 2) then
     call exchdl(rho,vx,exc,nx)
   else if(iexc .eq. 3) then
     call excpzca(rho,vx,exc,nx)
   else if(iexc .eq. 4) then
     call excggc(rho,vx,exc,rr,nx)
   else if (iexc < 0) then
     call exc_libxc(iexc,al,rho,vx,exc,rr,nx)
   else
     write(6,'(/a,i4)') 'modcore_budget: ERROR bad input iexc =',iexc
     stop
   end if
 end subroutine xc_eps

 end subroutine modcore_budget
