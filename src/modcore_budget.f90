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
! r_s is the smallest grid radius, no larger than the smallest rc (or the
! last radius with core charge), whose model needs no more than the cutoff
! ecut for an exchange-correlation energy error below errmax.
!
! ecut is fcfact (Ha) if fcfact > 0.  Otherwise it is the cutoff the
! pseudo-waves need for the same error (ecutw, from modcore_budget_ecut),
! times -fcfact if fcfact < 0: a calculation needs that cutoff anyway.
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
!
! The error almost always peaks at the first of these cutoffs, so r_s is
! bisected for with the first nquick of them only, and then confirmed with
! all of them (and bisected for again with all of them if that fails).

 subroutine modcore_budget(rhotps,rhoc,rhomod,fcfact,ecutw,errmax,irmax, &
&                           mmax,rr,iexc)

!rhotps  total pseudocharge density
!rhoc  core-charge density
!rhomod  model core density and 4 derivatives
!fcfact  plane-wave cutoff (Ha) the model must not need more than, if > 0;
!        else 0 for ecutw, or minus the multiple of ecutw
!ecutw  cutoff (Ha) the pseudo-waves need for errmax (used if fcfact <= 0)
!errmax  exchange-correlation energy error (Ha) that defines the needed cutoff
!irmax  rr index of the largest r_s (smallest rc)
!mmax  dimension of log grid
!rr log radial grid
!iexc  exchange-correlation function to be used

 use constants_m, only: dp, pi
 implicit none

!Input variables
 integer, intent(in) :: mmax,irmax,iexc
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: fcfact,ecutw,errmax

!Output variables
 real(dp), intent(out) :: rhomod(mmax,5)

!Local variables
 real(dp), parameter :: pi4=4.0_dp*pi
!background density (bohr**-3) and radius (bohr) of the error,
!largest density cutoff and its step (bohr**-1),
!panels of the band-limited transform per step,
!cutoffs of the bisection before the confirmation with all of them
 real(dp), parameter :: background=5.0d-3, radius=10.0d0
 real(dp), parameter :: qmax=40.0d0, qstep=0.1d0
 integer, parameter :: npanel=4, nquick=10
 real(dp) :: al,dq,vbg,ecore,aco,bco,cco,e0,ecut,factor
 real(dp), allocatable :: dc(:,:),core(:),nval(:),sn(:),cs(:),sh(:),ch(:)
 real(dp), allocatable :: rq(:),work(:),eps(:),vxc(:),err(:),rqs(:,:)
 integer :: ii,jj,nq,nin,nx,ilo,ihi,irs,ifrom,iquick,ibad,irtop

 if(fcfact>0.0d0) then
   ecut=fcfact
 else
   factor=1.0d0
   if(fcfact<0.0d0) factor=-fcfact
   ecut=factor*ecutw
 end if

 al = 0.01d0 * log(rr(101) / rr(1))
 dq = qstep/npanel
 nq = nint(qmax/qstep)

 if(2.0d0*sqrt(2.0d0*ecut) >= qmax) then
   write(6,'(/a,f10.2,a)') 'modcore_budget: ERROR stop - ecut beyond', &
&        qmax**2/8.0d0,' Ha'
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
   stop
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

! the cutoff the model needs: (q+qstep)**2/8 for the last density cutoff q
! above errmax, which is below 2*sqrt(2*ecut)
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
   call martins_derivatives(rr(ii),rhomod(ii,:))
 end do

 write(6,'(/a)') 'Model core correction: budgeted Martins core (icmod 6)'
 if(fcfact<=0.0d0) then
   write(6,'(a,f12.4,a,f8.4)') '  pseudo-wave cutoff',ecutw, &
&        ' Ha, times',factor
 end if
 write(6,'(a,f12.4,a,1p,e10.2,a)') '  cutoff budget',ecut, &
&      ' Ha for an xc energy error of',errmax,' Ha'
 write(6,'(a,f12.6,a,f10.6)') '  model radius r_s',rr(irs), &
&      '  smallest rc',rr(irmax)
 if(irtop<irmax) then
   write(6,'(a,f12.6)') '  last radius with core charge',rr(irtop)
 end if
 write(6,'(a,1p,3e16.8)') '  a, b, c',aco,bco,cco
 write(6,'(a,f12.4,a)') '  cutoff the model core needs',ecore,' Ha'

 deallocate(dc,core,sn,cs,sh,ch,err)
 deallocate(nval,rq,work,eps,vxc)

 return

 contains

! whether the model of radius rr(is) has errors below errmax at the
! density cutoffs ifrom..ilast
 logical function fits(is,ilast)

!is  rr index of the model radius
!ilast  last density cutoff index tested

!Input variables
 integer, intent(in) :: is,ilast

!Local variables
 integer :: ifail
 call model(is)
 call scan_errors(ifrom,ilast,errmax,ifail)
 fits = ifail==0
 end function fits

! the smallest index above ilo, which fails, that fits up to ilast;
! ihi fits
 integer function bisect(ilo,ihi,ilast)

!ilo  lower rr index of the bracket, a model radius that fails
!ihi  upper rr index of the bracket, a model radius that fits
!ilast  last density cutoff index tested

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
 end function bisect

! reports that no model radius fits and stops
 subroutine no_model()
 write(6,'(/a,f10.6,a,f10.4,a,1p,e10.2,a)') &
&      'modcore_budget: ERROR stop - no model core inside r =',rr(irtop), &
&      ' needs less than',ecut,' Ha for an error of',errmax,' Ha'
 stop
 end subroutine no_model

! the Martins coefficients at rr(is), and in core the model inside rr(is)
! and the core charge beyond
 subroutine model(is)

!is  rr index of the model radius

!Input variables
 integer, intent(in) :: is

!Local variables
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

! the model and its first 4 derivatives at rad, from n' = p' n with
! p = a + b r**2 + c r**4: n(k) = sum_j C(k-1,j) p(j+1) n(k-1-j)
 subroutine martins_derivatives(rad,dd)

!rad  radius (bohr)
!dd  the model and its first 4 derivatives at rad

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
 end subroutine martins_derivatives

! the errors err(iq) at the density cutoffs iq*qstep, iq=ifirst..ilast, of
! the model in core, up to the first above elimit (its index in ifail,
! else 0).  With ilast=0 instead the band-limited models of the cutoffs
! 1..ifrom-1 are kept in rqs.
! The band-limited model is built up panel by panel (midpoints, width dq),
! the sines and cosines of the midpoints by rotation
 subroutine scan_errors(ifirst,ilast,elimit,ifail)

!ifirst  first density cutoff index whose error is tested
!ilast  last density cutoff index, or 0 to keep the band-limited models
!elimit  error limit (Ha)
!ifail  index of the first cutoff with an error above elimit, else 0

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
 end subroutine scan_errors

! the error of the band-limited model bq of the model in core
 real(dp) function xc_error(bq)

!bq  band-limited model core charge on rr(1:nx)

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
 end function xc_error

! the exchange-correlation energy within the sphere of the valence charge
! on the background plus the core charge cc
 real(dp) function exc_sphere(cc)

!cc  core charge on rr(1:nx)

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
 end function exc_sphere

! the xc energy per electron and potential of the charge rho on rr(1:nx)
 subroutine xc_eps(rho,exc,vx)

!rho  charge density on rr(1:nx)
!exc  exchange-correlation energy per electron
!vx  exchange-correlation potential

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
 else if (iexc < 0) then
   call exc_libxc(iexc,al,rho,vx,exc,rr,nx)
 else
   write(6,'(/a,i4)') 'modcore_budget: ERROR bad input iexc =',iexc
   stop
 end if
 end subroutine xc_eps

 end subroutine modcore_budget

! the cutoff ecut (Ha) the pseudo-waves need for the error errmax (Ha per
! atom): the root of the total residual kinetic energy sum_i w_i E_i(ecut),
! with E_i the error per electron of pseudo-wave i and w_i the electrons
! that fill it, the ecut_for_error of PseudoPotentialGenerators.jl.
! ecurve(:,:,i) holds the cutoffs and E_i of run_optimize on its mesh,
! zero beyond its end; wcurve(i)=0 skips it.  log E is interpolated
! linearly in q between the mesh points, as for the convergence profiles
 subroutine modcore_budget_ecut(ecurve,wcurve,ncurve,nprof,errmax,ecut)

!ecurve  residual kinetic energy curves: ecurve(1,:,i) cutoffs (Ha),
!        ecurve(2,:,i) error per electron E_i (Ha) of pseudo-wave i
!wcurve  electrons filling each pseudo-wave (0 skips it)
!ncurve  number of points per curve
!nprof  number of curves (pseudo-waves)
!errmax  error (Ha per atom) that defines the needed cutoff
!ecut  cutoff (Ha) the pseudo-waves need for errmax

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
   stop
 end if
 if(total(nend)>errmax) then
   write(6,'(/a,1p,e10.2,a,0p,f10.2,a)') &
&        'modcore_budget_ecut: ERROR stop - residual kinetic energy above', &
&        errmax,' Ha at',0.5d0*qq(nend)**2,' Ha, the end of its tables'
   stop
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
