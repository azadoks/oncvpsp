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
! (found and evaluated as a series of Legendre polynomials of
! 2(r/r_s)**2 - 1: in powers, its coefficients reach 1e5 and alternate,
! and carry 1e-11 of noise into the model) whose six free coefficients
! minimize the Fourier tail of the whole model,
! int_{qlo}^{qhi} q**2 |n(q)|**2 dq, with qlo = 2*sqrt(2*ecut) and
! qhi = min(3*qlo, qmax).  Where that polynomial rises somewhere on
! [0, r_s] (p' > 0 at 0, at r_s or at a root of p''), the fill is Martins'
! exp(a + b r**2 + c r**4) instead, if that does not rise (b <= 0 and
! b + 2c r_s**2 <= 0); where neither is monotone, the radius has no model.
!
! ecut is fcfact (Ha) if fcfact > 0.  Otherwise it is the cutoff of the
! pseudo-waves at the same error per electron (ecutw, from
! modcore_budget_table: the largest Ecut of the "Energy error per electron"
! tables, or with budget_cutoff = "weighted" that of their mean weighted by
! the occupations), times -fcfact if fcfact < 0: a calculation runs at that
! cutoff anyway.
! The model must fit slack*ecut (0 < slack <= 1), and is still made for
! the qlo of ecut.
!
! The error of a density cutoff q (a plane-wave cutoff q**2/8) is that of
! the core charge without its Fourier components beyond q, in a solid:
!
!   |Exc[n + nc_q] - Exc[n + nc] - Vxc[n_bg] * int (nc_q - nc)|
!
! within the sphere r <= 10 bohr (its edge a C-infinity taper from 6 bohr,
! cinf_taper: a hard edge rippled the error in q with period pi/10, and
! the largest error hopped between ripples), with n the pseudo valence
! charge on a uniform background n_bg = 5e-3 bohr**-3 that stands for the
! interstitial electrons of the neighbours.  A model fits the cutoff e if
! the error stays below errmax for every density cutoff
! q >= q_b = 2*sqrt(2*e), up to 120 bohr**-1 (1800 Ha): the largest error
! over q, on a grid of 0.1 bohr**-1 from q_b and refined by Brent's method
! at its maxima.  The
! transforms are exact in r_s and q (M-8, M-12 of doc/smoothness_plan.md):
! the model's by Gauss-Legendre on [0, r_s] and on the core charge's
! interpolant (lm_interp) from r_s out, the band-limited model by
! Gauss-Legendre panels in q, the last one cut at q.
!
! icmod=6 (BudgetedCore of PseudoPotentialGenerators.jl): r_s is the
! smallest radius, no larger than the smallest rc (or the last radius with
! core charge), whose model fits slack*ecut: bracketed on the mesh, then
! the root of the largest error less errmax by Brent's method (M-9).  The
! error almost always peaks at the first of the cutoffs, so r_s is found
! with those up to q_b + 1 bohr**-1 only, and then confirmed with all of
! them (and found again with all of them if that fails).
!
! icmod=7 (HardnessMatchedCore): the radius in [r_s of icmod 6, its largest
! one] whose model is closest to the core charge in the exchange-
! correlation hardness of the pseudo valence: the rms over i,j of
! d2Exc/dn_i dn_j at rhotps + model against that at rhotps + rhoc
! (der2exc), for the probe densities of modcore_probes.  Radii from the
! largest one down by factors of 1.025, and r_s of icmod 6, choose the
! basin; the minimum is the root of the mismatch's derivative in ln r_s,
! by central differences, in it (M-10).  The radii above r_s of icmod 6
! fit, as the needed cutoff falls as the model softens.

 subroutine modcore_budget(icmod,rhotps,rhoc,rhomod,corepw,fcfact,ecutw, &
&                          weighted,errmax,slack,rcmin,mmax,rr,iexc,zion, &
&                          probes, &
&                          lprobe,nprobe)

!! @param icmod  6 for the smallest radius, 7 for the hardness-matched one
!! @param rhotps  total pseudocharge density
!! @param rhoc  core-charge density
!! @param rhomod  model core density and 4 derivatives
!! @param corepw  its pieces, split at r_s (modcore_m, core_pieces)
!! @param fcfact  plane-wave cutoff (Ha) the model is made for, if > 0; else 0
!!                for ecutw, or minus the multiple of ecutw
!! @param ecutw  cutoff (Ha) of the pseudo-waves at errmax per electron (used
!!               if fcfact <= 0)
!! @param weighted  whether ecutw is that of the weighted mean of the
!!                  pseudo-waves (modcore_budget_table), for the output
!! @param errmax  exchange-correlation energy error (Ha) that defines the needed
!!                cutoff
!! @param slack  share of the cutoff the model must fit, 0 < slack <= 1
!! @param rcmin  the smallest rc, the largest r_s (bohr)
!! @param mmax  dimension of log grid
!! @param rr  log radial grid
!! @param iexc  exchange-correlation function to be used
!! @param zion  ion charge (for der2exc)
!! @param probes  probe densities of the hardness (icmod=7), one electron each
!! @param lprobe  angular momenta of the probes
!! @param nprobe  number of probes (0 for icmod=6)

 use constants_m, only: dp, pi, nrhomod
 use modcore_m, only: core_pieces
 use ugrid_m, only: ug_brent, ug_fmin
 use quad_m, only: gauss_legendre
 use xc_dispatch_m, only: exc_dispatch
 use pieces_m, only: piecewise_t
 use logmesh_m, only: lm_derivs, lm_interp, lm_integral, lm_fh_reach, &
&                     lm_blend
 implicit none

!Input variables
 integer, intent(in) :: icmod,mmax,iexc,nprobe
 integer, intent(in) :: lprobe(nprobe)
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: probes(mmax,nprobe)
 real(dp), intent(in) :: fcfact,ecutw,errmax,slack,zion,rcmin
 logical, intent(in) :: weighted

!Output variables
 real(dp), intent(out) :: rhomod(mmax,nrhomod)
 type(piecewise_t), intent(out) :: corepw(nrhomod)

 abstract interface
   real(dp) function scalar_fn(xx)
   import :: dp
   implicit none
   real(dp), intent(in) :: xx
   end function scalar_fn
 end interface

!Local variables
 real(dp), parameter :: pi4=4.0_dp*pi
!background density (bohr**-3) and radius (bohr) of the error,
!largest density cutoff and its grid (bohr**-1), the cutoffs past q_b
!of the first search
 real(dp), parameter :: background=5.0d-3, radius=10.0d0, taper=0.6d0
 real(dp), parameter :: qmax=120.0d0, qstep=0.1d0, quick=1.0d0
!panels of the band-limited transform per step and their Gauss-Legendre
!nodes, nodes of the transform on [0, r_s], nodes per mesh interval beyond
!(each converged: twice as many move nothing printed)
 integer, parameter :: npanel=1, ngk=5, nxg=400, ngi=8
!highest power (r/r_s)**(2*npow) of the minimum tail, Gauss-Legendre nodes
!per panel of its window, panels of the window
 integer, parameter :: npow=8, ngq=8, nwin=80
!icmod=7: step of the radii, step in ln r_s of the central differences
 real(dp), parameter :: ratio=1.025d0, hmis=1.0d-2
 real(dp) :: al,dq,vbg,ecore,aco,bco,cco,ecut,budget,factor,qlo,qhi,qb
 real(dp) :: cpoly(0:npow),cleg(0:npow),mismatch,rs,rs6,rtop,rlo,qend
 real(dp) :: qneed,tt,ee,qpk,qext
 real(dp) :: tk(ngk),wk(ngk),tgi(ngi),wgi(ngi),xg(nxg),wxg(nxg)
 real(dp) :: pt(ngi),pg(ngi)
 real(dp), allocatable :: core(:),nval(:),work(:),eps(:),vxc(:),wf(:),wc(:)
 real(dp), allocatable :: fext(:),af(:),ac(:),swin(:),xr(:),ri(:),nvi(:)
 real(dp), allocatable :: wki(:),epi(:),vxi(:)
 real(dp), external :: cinf_taper
 real(dp), allocatable :: qn(:),wn(:),ot(:),og(:),sho(:),cho(:)
 real(dp), allocatable :: so0(:,:),co0(:,:),tab(:,:,:),tabw(:,:)
 real(dp), allocatable :: fx(:),fw(:),plg(:,:)
 real(dp), allocatable :: shr(:),chr(:),sr0(:,:),cr0(:,:),sr(:,:),cr(:,:)
 real(dp), allocatable :: ntk(:,:),rqs(:,:),rqp(:)
 real(dp), allocatable :: d2ae(:,:),d2mod(:,:),d2zero(:,:),cmis(:),crad(:)
 integer, allocatable :: cfill(:)
 integer :: ii,jj,kk,nq,npq,nqn,nin,nx,irtop,ifill,ncand,i0,iend,nout,js
 integer :: kin,ni,nt
 integer :: irmx,irow,ip,ig
 real(dp) :: sa,ca,sb
 integer :: npdone,msnap,ib
 logical :: full
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
 npq = nq*npanel

 qlo=2.0d0*sqrt(2.0d0*ecut)
 if(qlo >= qmax) then
   write(6,'(/a,f10.2,a)') 'modcore_budget: ERROR stop - ecut beyond', &
&        qmax**2/8.0d0,' Ha'
   error stop 1
 end if
 qhi=min(3.0d0*qlo,qmax)
! the first density cutoff that decides whether the model fits budget
 qb=2.0d0*sqrt(2.0d0*budget)

! the largest r_s: the smallest rc, or the last radius at which the core
! charge and the derivatives matched there are not zero
 irtop=0
 do ii=1,mmax-6
   if(.not. all(rhoc(ii:ii+6)>0.0d0)) exit
   irtop=ii
 end do
 if(irtop<4) then
   write(6,'(/a)') 'modcore_budget: ERROR stop - no core charge'
   error stop 1
 end if
 rtop=min(rcmin,rr(irtop))
 rlo=rtop/50.0d0
 if(qmax*rtop > 2.0d0*(nxg-50)) then
   write(6,'(/a,f10.4,a)') 'modcore_budget: ERROR stop - r_s up to', &
&        rtop,' bohr needs more nodes of the transform'
   error stop 1
 end if

 call gauss_legendre(0.0d0,1.0d0,ngk,tk,wk)
 call gauss_legendre(0.0d0,1.0d0,ngi,tgi,wgi)
 call gauss_legendre(0.0d0,1.0d0,nxg,xg,wxg)

! the window of the minimum tail: nwin Gauss-Legendre panels on
! [qlo, qhi], the weights times q**2
 nqn=ngq*nwin
 allocate(qn(nqn),wn(nqn))
 do ii=1,nwin
   call gauss_legendre(qlo+(ii-1)*(qhi-qlo)/nwin,qlo+ii*(qhi-qlo)/nwin, &
&                      ngq,qn((ii-1)*ngq+1:ii*ngq),wn((ii-1)*ngq+1:ii*ngq))
 end do
 wn(:)=wn(:)*qn(:)**2

! the core charge's transform beyond r_s: Gauss-Legendre in ln r on each
! mesh interval from below r_lo (and the steps of the central differences)
! to the last one where r**3 rhoc is above 1e-14 of its largest value,
! of the interpolant of rhoc; og the weights times r**2 rhoc, so that
! int r**2 rhoc j0(kr) dr = sum og sin(k ot)/k.  Their sums at the
! nodes of the band-limited transform and of the window, from each mesh
! point r_s can lie below on: tab(:,:,i), tabw(:,i) over the intervals
! from i on, the sines rotated from one panel to the next
 tt=maxval(rr(:)**3*rhoc(:))
 iend=1
 do ii=1,mmax-2
   if(rr(ii)**3*rhoc(ii)>1.0d-14*tt) iend=ii+1
 end do
 i0=max(1,floor(log(rlo/(1.05d0*rr(1)))/al))
 if(iend<=i0+1) then
   write(6,'(/a)') 'modcore_budget: ERROR stop - no core charge'
   error stop 1
 end if
 nout=ngi*(iend-i0)
 allocate(ot(nout),og(nout),sho(nout),cho(nout))
 allocate(so0(nout,ngk),co0(nout,ngk))
 do ii=i0,iend-1
   do kk=1,ngi
     jj=(ii-i0)*ngi+kk
     ot(jj)=rr(ii)*exp(al*tgi(kk))
     call lm_interp(rr,rhoc,mmax,ot(jj),tt)
     og(jj)=al*wgi(kk)*ot(jj)**2*tt
   end do
 end do
 sho(:)=sin(dq*ot(:))
 cho(:)=cos(dq*ot(:))
 do kk=1,ngk
   so0(:,kk)=sin(tk(kk)*dq*ot(:))
   co0(:,kk)=cos(tk(kk)*dq*ot(:))
 end do
 irmx=min(iend,floor(log(rtop*exp(2.0d0*hmis)/rr(1))/al)+3)
 allocate(tab(npq,ngk,i0:irmx),tabw(nqn,i0:irmx))
 allocate(rqp(npq*ngk))
 rqp(:)=0.0d0
 tabw(:,i0)=0.0d0
 do ii=iend-1,i0,-1
   if(ii+1<=irmx) then
     tab(:,:,ii+1)=reshape(rqp,[npq,ngk])
     tabw(:,ii+1)=tabw(:,i0)
   end if
   do kk=1,ngi
     jj=(ii-i0)*ngi+kk
     do ig=1,ngk
       sa=so0(jj,ig)
       ca=co0(jj,ig)
       do ip=1,npq
         rqp(ip+(ig-1)*npq)=rqp(ip+(ig-1)*npq) + og(jj)*sa
         sb=sa*cho(jj) + ca*sho(jj)
         ca=ca*cho(jj) - sa*sho(jj)
         sa=sb
       end do
     end do
     tabw(:,i0)=tabw(:,i0) + og(jj)*sin(qn(:)*ot(jj))
   end do
 end do
 tab(:,:,i0)=reshape(rqp,[npq,ngk])
 deallocate(rqp)
 allocate(fx(nxg),fw(nxg),plg(nxg,0:npow))
! the Legendre polynomials of the fill at the nodes on [0, r_s]
 do ii=1,nxg
   call legendre_derivs(2.0d0*xg(ii)**2-1.0d0,0,plg(ii,:))
 end do

! the sphere of the error, and the points its integrals' interpolant
! reads beyond it
 nin=0
 do ii=1,mmax
   if(rr(ii)>radius) exit
   nin=ii
 end do
 nx=min(mmax,nin+lm_fh_reach()+lm_blend+10)
! the mesh of the integral below r_s (set_model): ri(k) = r_s
! exp((k-kin) al), from at least rr(1) to the points its interpolant reads
! beyond r_s (r_s from r_lo, less the steps of the central differences).
! The band-limited models are taken at xr: rr(1:nx), then ri
 kin=floor((log(rlo/rr(1)) - 2.0d0*hmis)/al)+1
 ni=kin+lm_fh_reach()+lm_blend+10
 nt=nx+ni
 allocate(core(mmax),nval(nx),work(nx),eps(nx),vxc(nx),rqp(nt))
 allocate(wf(ni),wc(nx),fext(ni),af(ni),ac(nx),swin(nx))
 allocate(xr(nt),ri(ni),nvi(ni),wki(ni),epi(ni),vxi(ni))
 do ii=1,nx
   swin(ii)=cinf_taper(rr(ii),taper*radius,radius)
 end do
 allocate(shr(nt),chr(nt),sr0(nt,ngk),cr0(nt,ngk),sr(nt,ngk),cr(nt,ngk))
 allocate(ntk(npq,ngk),rqs(nt,0:nq))
 xr(1:nx)=rr(1:nx)
 shr(1:nx)=sin(dq*rr(1:nx))
 chr(1:nx)=cos(dq*rr(1:nx))
 do kk=1,ngk
   sr0(1:nx,kk)=sin(tk(kk)*dq*rr(1:nx))
   cr0(1:nx,kk)=cos(tk(kk)*dq*rr(1:nx))
 end do

 nval(:)=pi4*background + rhotps(1:nx)

! the xc energy density of the valence charge on the core charge
 work(:)=nval(:) + rhoc(1:nx)
 call exc_dispatch(iexc,al,work,vxc,eps,rr,nx,'modcore_budget')
 wc(:)=work(:)*eps(:)

! xc potential of the background alone
 work(:)=pi4*background
 call exc_dispatch(iexc,al,work,vxc,eps,rr,nx,'modcore_budget')
 vbg=vxc(nx/2)

! smallest r_s, first with the cutoffs to q_b + quick; if a larger one
! has a larger error at that r_s, again with those up to it + quick: each
! r_s found is a root of the largest error over all cutoffs less errmax
! once it is confirmed, whatever the range it was found with
 full=.false.
 qend=min(qmax,qb+quick)
 do
   rs=smallest_radius()
   call set_model(rs)
   if(qend>=qmax) exit
   if(err_max(qend,qmax,errmax)<=errmax) exit
   full=.true.
   qend=min(qmax,qpk+quick)
 end do
 rs6=rs
 qext=qend
 qend=qmax

! icmod=7: the radii from the largest one down to r_s of icmod 6, their
! hardness mismatch, and the minimum in the basin of the least of them
 if(icmod==7) then
   kk=max(0,floor(log(rtop/rs6)/log(ratio)))
   allocate(crad(kk+2),cfill(kk+2),cmis(kk+2))
   ncand=0
   do ii=0,kk
     tt=rtop/ratio**ii
     if(tt<=rs6) exit
     ncand=ncand+1
     crad(ncand)=tt
   end do
   ncand=ncand+1
   crad(ncand)=rs6

   allocate(d2ae(nprobe,nprobe),d2mod(nprobe,nprobe),d2zero(nprobe,nprobe))
   d2zero(:,:)=0.0d0
   call der2exc(rhotps,rhoc,probes,rr,d2ae,d2zero,mismatch, &
&               zion,iexc,0,nprobe,lprobe,rr(mmax),mmax)
   ib=0
   do ii=1,ncand
     cmis(ii)=mismatch_at(log(crad(ii)))
     cfill(ii)=ifill
     if(ifill==0) cycle
     if(ib==0) then
       ib=ii
     else if(cmis(ii)<cmis(ib)) then
       ib=ii
     end if
   end do
   if(ib==0) call no_model()
   rs=minimum_radius(ib)
   call set_model(rs)
   mismatch=mismatch_at(log(rs))
! r_s6 is a root of this error less errmax, so at r_s = r_s6 the two are
! equal to round-off; the margin is above that, far below a radius that
! does not fit
   if(err_max(qb,qmax,errmax)>errmax*(1.0d0+1.0d-7)) then
     write(6,'(/a,f12.6,a)') 'modcore_budget: WARNING the radius',rs, &
&          ' of the least mismatch does not fit the budget'
   end if
 end if

! the cutoff the model needs: q**2/8 for the last density cutoff q below
! q_b whose error is errmax
 call set_model(rs)
 call spectrum(qb)
 qneed=0.0d0
 tt=qb
 ee=err_at(qb)-errmax
 do kk=ceiling(qb/qstep)-1,1,-1
   if(ee>=0.0d0) then
     qneed=tt
     exit
   end if
   if(err_grid(kk)>errmax) then
     qneed=zeroin(excess,kk*qstep,tt)
     exit
   end if
   tt=kk*qstep
 end do
 ecore=qneed**2/8.0d0

! fill in model and derivatives, all-electron beyond r_s
 call core_pieces(rr,rhoc,mmax,rs,fill_derivatives,corepw,rhomod)

 if(icmod==7) then
   write(6,'(/a)') 'Model core correction: hardness-matched core (icmod 7)'
 else
   write(6,'(/a)') 'Model core correction: budgeted core (icmod 6)'
 end if
 if(fcfact<=0.0d0 .and. weighted) then
   write(6,'(a,1p,e10.2,a,0p,f12.4,a,f8.4)') &
&        '  pseudo-wave cutoff, weighted mean Ecut of the tables at',errmax, &
&        ' Ha/electron',ecutw,' Ha, times',factor
 else if(fcfact<=0.0d0) then
   write(6,'(a,1p,e10.2,a,0p,f12.4,a,f8.4)') &
&        '  pseudo-wave cutoff, largest Ecut of the tables at',errmax, &
&        ' Ha/electron',ecutw,' Ha, times',factor
 end if
 write(6,'(a,f12.4,a,1p,e10.2,a)') '  cutoff budget',ecut, &
&      ' Ha for an xc energy error of',errmax,' Ha'
 write(6,'(a,f8.4,a,f12.4,a)') '  slack',slack, &
&      ', the model core must fit',budget,' Ha'
 write(6,'(a,f12.6,a,f10.6)') '  model radius r_s',rs, &
&      '  smallest rc',rcmin
 if(rr(irtop)<rcmin) then
   write(6,'(a,f12.6)') '  last radius with core charge',rr(irtop)
 end if
 if(full) then
   write(6,'(a,f10.4,a)') '  (the error peaks beyond the first &
&cutoffs: found with those up to',qext,' bohr**-1)'
 end if
 if(ifill==1) then
   write(6,'(a)') '  fill minimum tail'
   call legendre_powers(.true.,cleg,cpoly)
   write(6,'(a,1p,9e16.8)') '  coefficients c0..c8',cpoly(:)
 else
   write(6,'(a)') '  fill Martins'
   write(6,'(a,1p,3e16.8)') '  a, b, c',aco,bco,cco
 end if
 write(6,'(a,f12.4,a)') '  cutoff the model core needs',ecore,' Ha'
 if(icmod==7) then
   write(6,'(a,f12.6)') '  smallest fitting radius (icmod 6)',rs6
   write(6,'(a,1p,e12.4,a)') '  hardness mismatch',mismatch,' Ha'
   write(6,'(a)') '  radii of the basin: r_s, fill, hardness mismatch (Ha)'
   do ii=1,ncand
     if(cfill(ii)==0) then
       write(6,'(f12.6,2x,a2,a12)') crad(ii),fillname(0),'-'
     else
       write(6,'(f12.6,2x,a2,1p,e12.4)') crad(ii),fillname(cfill(ii)), &
&            cmis(ii)
     end if
   end do
   write(6,'(a)') ''
   deallocate(crad,cfill,cmis,d2ae,d2mod,d2zero)
 end if

 deallocate(core,nval,work,eps,vxc,wf,wc,fext,af,ac,swin,rqp,qn,wn)
 deallocate(xr,ri,nvi,wki,epi,vxi)
 deallocate(ot,og,sho,cho)
 deallocate(so0,co0,tab,tabw,fx,fw,plg)
 deallocate(shr,chr,sr0,cr0,sr,cr,ntk,rqs)

 return

 contains

 real(dp) function smallest_radius()

!! the smallest r_s in [r_lo, r_top] whose model fits, with the density
!! cutoffs up to qend: bisection over the mesh points between them, then
!! Brent's method on the largest error less errmax in ln r_s (excess_rs)

!Local variables
 integer :: ia,ib,imid,ilo,ihi
 real(dp) :: ga
 if(excess_rs(log(rtop))>0.0d0) call no_model()
 smallest_radius=rlo
 ga=excess_rs(log(rlo))
 if(ga<=0.0d0) return
! the mesh points strictly between r_lo and r_top, ilo..ihi; ia fails,
! ib fits, ilo-1 standing for r_lo and ihi+1 for r_top
 ilo=1
 do while(rr(ilo)<=rlo)
   ilo=ilo+1
 end do
 ihi=ilo-1
 do while(rr(ihi+1)<rtop)
   ihi=ihi+1
 end do
 ia=ilo-1
 ib=ihi+1
 do while(ib-ia>1)
   imid=(ia+ib)/2
   if(excess_rs(log(rr(imid)))>0.0d0) then
     ia=imid
   else
     ib=imid
   end if
 end do
 smallest_radius=exp(zeroin(excess_rs,log(radius_of(ia)), &
&                           log(radius_of(ib))))
 return
 end function smallest_radius

 real(dp) function radius_of(ii)

!! rr(ii), r_lo below the points of the bisection and r_top above them
!!
!! @param ii  index

!Input variables
 integer, intent(in) :: ii
 if(rr(ii)<=rlo) then
   radius_of=rlo
 else if(rr(ii)>=rtop) then
   radius_of=rtop
 else
   radius_of=rr(ii)
 end if
 return
 end function radius_of

 real(dp) function excess_rs(uu)

!! the largest error over the density cutoffs q_b..qend of the model of
!! radius exp(uu), less errmax; errmax if the radius has no model
!!
!! @param uu  ln r_s

!Input variables
 real(dp), intent(in) :: uu
 call set_model(exp(uu))
 if(ifill==0) then
   excess_rs=errmax
   return
 end if
 excess_rs=err_max(qb,qend,errmax) - errmax
 return
 end function excess_rs

 real(dp) function excess(qq)

!! the error at the density cutoff qq less errmax
!!
!! @param qq  density cutoff (bohr**-1)

!Input variables
 real(dp), intent(in) :: qq
 excess=err_at(qq) - errmax
 return
 end function excess

 real(dp) function mismatch_at(uu)

!! the hardness mismatch of the model of radius exp(uu); huge if the
!! radius has no model
!!
!! @param uu  ln r_s

!Input variables
 real(dp), intent(in) :: uu

!Local variables
 real(dp) :: mm
 call set_model(exp(uu))
 mismatch_at=huge(1.0d0)
 if(ifill==0) return
 call der2exc(rhotps,core,probes,rr,d2mod,d2ae,mm, &
&             zion,iexc,0,nprobe,lprobe,rr(mmax),mmax)
 mismatch_at=mm
 return
 end function mismatch_at

 real(dp) function slope_at(uu)

!! the derivative of the hardness mismatch in ln r_s, by central
!! differences of step hmis
!!
!! @param uu  ln r_s

!Input variables
 real(dp), intent(in) :: uu
 slope_at=(mismatch_at(uu+hmis) - mismatch_at(uu-hmis))/(2.0d0*hmis)
 return
 end function slope_at

 real(dp) function minimum_radius(ib)

!! the radius of the least hardness mismatch in the basin of the radius
!! crad(ib): the root of its slope between the neighbours of crad(ib),
!! or an end of [r_s of icmod 6, r_top] where the slope points out of it
!!
!! @param ib  index of the radius of the least mismatch of the list

!Input variables
 integer, intent(in) :: ib

!Local variables
 real(dp) :: ua,ub,sa,sb
 integer :: ia,ic
 ia=min(ncand,ib+1)
 ic=max(1,ib-1)
 minimum_radius=crad(ib)
 if(ia==ib) then
   if(slope_at(log(crad(ib)))>=0.0d0) return
 end if
 if(ic==ib) then
   if(slope_at(log(crad(ib)))<=0.0d0) return
 end if
 ua=log(crad(ia))
 ub=log(crad(ic))
 sa=slope_at(ua)
 sb=slope_at(ub)
 if(sa>0.0d0 .or. sb<0.0d0) then
   write(6,'(/a,f12.6)') 'modcore_budget: WARNING no minimum of the &
&mismatch bracketed about',crad(ib)
   return
 end if
 minimum_radius=exp(zeroin(slope_at,ua,ub))
 return
 end function minimum_radius

 subroutine no_model()
!! reports that no model radius fits and stops
 write(6,'(/a,f10.6,a,f10.4,a,1p,e10.2,a)') &
&      'modcore_budget: ERROR stop - no model core inside rc =',rtop, &
&      ' needs less than',budget,' Ha for an error of',errmax,' Ha'
 error stop 1
 end subroutine no_model

 subroutine set_model(rad)

!! the model of radius rad: the fill (ifill 1 for the minimum tail in cpoly,
!! 2 for Martins' in aco, bco, cco, 0 for none), in core the model inside
!! rad and the core charge beyond, the nodes of its transform, the mesh
!! ri of the integral below rad, and the xc energy density of the valence
!! charge on the fill there
!!
!! @param rad  model radius r_s (bohr)

!Input variables
 real(dp), value :: rad

!Local variables
 real(dp) :: dm(3),l1,l2,xs,xe
 integer :: kk,is

 rs=rad
 npdone=0
 msnap=0
 core(:)=rhoc(:)
 call lm_derivs(rr,rhoc,mmax,rs,2,dm)

! the nodes of the transform beyond rs: the mesh intervals from is+1 on,
! and Gauss-Legendre on the part [rs, rr(is+1)] of interval is
 is=max(1,floor(log(rs/rr(1))/al)+1)
 do while(rr(is+1)<=rs)
   is=is+1
 end do
 js=(is+1-i0)*ngi+1
 irow=is+1
 if(irow>irmx) then
   write(6,'(/a)') 'modcore_budget: ERROR stop - r_s beyond the table'
   error stop 1
 end if
 xs=log(rs)
 xe=log(rr(is+1))
 do kk=1,ngi
   pt(kk)=exp(xs + (xe-xs)*tgi(kk))
   call lm_interp(rr,rhoc,mmax,pt(kk),l1)
   pg(kk)=(xe-xs)*wgi(kk)*pt(kk)**2*l1
 end do

 call min_tail(dm)
 ifill=0
 if(tail_monotone()) then
   ifill=1
   do kk=1,mmax
     if(rr(kk)>=rs) exit
     core(kk)=poly((rr(kk)/rs)**2)
   end do
 else
   l1=dm(2)/dm(1)
   l2=dm(3)/dm(1) - l1**2
   cco=(rs*l2 - l1)/(8.0d0*rs**3)
   bco=(l1 - 4.0d0*cco*rs**3)/(2.0d0*rs)
   aco=log(dm(1)) - bco*rs**2 - cco*rs**4
   if(bco>0.0d0 .or. bco + 2.0d0*cco*rs**2>0.0d0) return
   ifill=2
   do kk=1,mmax
     if(rr(kk)>=rs) exit
     core(kk)=exp(aco + bco*rr(kk)**2 + cco*rr(kk)**4)
   end do
 end if

! the fill's transform on [0, rs]: int r**2 n j0(kr) dr = sum fw sin(k fx)/k
 fx(:)=rs*xg(:)
 do kk=1,nxg
   fw(kk)=rs*wxg(kk)*fx(kk)*fill_value(fx(kk))
 end do

! the mesh ri, its points of the band-limited models, the valence charge
! there and the xc energy density of the valence charge on the fill,
! extended beyond rs as far as the interpolant of the integral to rs reads
 do kk=1,ni
   ri(kk)=rs*exp((kk-kin)*al)
   call lm_interp(rr,rhotps,mmax,ri(kk),nvi(kk))
   fext(kk)=fill_value(ri(kk))
 end do
 nvi(:)=pi4*background + nvi(:)
 xr(nx+1:nt)=ri(:)
 shr(nx+1:nt)=sin(dq*ri(:))
 chr(nx+1:nt)=cos(dq*ri(:))
 do kk=1,ngk
   sr0(nx+1:nt,kk)=sin(tk(kk)*dq*ri(:))
   cr0(nx+1:nt,kk)=cos(tk(kk)*dq*ri(:))
 end do
 wki(:)=nvi(:) + fext(:)
 call exc_dispatch(iexc,al,wki,vxi,epi,ri,ni,'modcore_budget')
 wf(:)=wki(:)*epi(:)
 return
 end subroutine set_model

 real(dp) function fill_value(rad)

!! the fill at rad <= rs
!!
!! @param rad  radius (bohr)

!Input variables
 real(dp), intent(in) :: rad
 if(ifill==1) then
   fill_value=poly((rad/rs)**2)
 else
   fill_value=exp(aco + bco*rad**2 + cco*rad**4)
 end if
 return
 end function fill_value

 logical function tail_monotone()

!! whether the minimum tail p, as a function of x = 2(r/rs)**2 - 1, does
!! not rise on [-1, 1]: p' <= 0 at -1, at 1, and at the roots of p'' in
!! between (companion matrix of p'' in powers of x, the roots polished by
!! Newton's method; a nearly real root counts by its real part, a point of
!! [-1, 1] all the same); p' from the Legendre series

!Local variables
 real(dp) :: cx(0:npow),d2(0:npow-2),hess(npow-2,npow-2)
 real(dp) :: wr(npow-2),wi(npow-2),zz(1,1),wkh(4*npow),tr,vv,vd
 integer :: kk,jj,nd,info,it
 external :: dhseqr
 tail_monotone=.false.
 if(leg_value(-1.0d0,1)>0.0d0) return
 if(leg_value(1.0d0,1)>0.0d0) return
 tail_monotone=.true.
 call legendre_powers(.false.,cleg,cx)
 do kk=2,npow
   d2(kk-2)=kk*(kk-1)*cx(kk)
 end do
 nd=npow-2
 do while(nd>0)
   if(d2(nd)/=0.0d0) exit
   nd=nd-1
 end do
 if(nd==0) return
 hess(:,:)=0.0d0
 do kk=1,nd
   hess(1,kk)=-d2(nd-kk)/d2(nd)
   if(kk>1) hess(kk,kk-1)=1.0d0
 end do
 call dhseqr('E','N',nd,1,nd,hess,npow-2,wr,wi,zz,1,wkh,size(wkh),info)
 if(info/=0) then
   write(6,'(/a,i4)') 'modcore_budget: ERROR dhseqr info =',info
   error stop 1
 end if
 do kk=1,nd
   tr=wr(kk)
   do it=1,2
     vv=d2(nd)
     vd=0.0d0
     do jj=nd-1,0,-1
       vd=vd*tr + vv
       vv=vv*tr + d2(jj)
     end do
     if(vd==0.0d0) exit
     tr=tr - vv/vd
   end do
   if(tr>-1.0d0 .and. tr<1.0d0) then
     if(leg_value(tr,1)>0.0d0) then
       tail_monotone=.false.
       return
     end if
   end if
 end do
 return
 end function tail_monotone

 real(dp) function leg_value(xx,mm)

!! the mm-th derivative in x of the Legendre series cleg at xx
!!
!! @param xx  argument, 2(r/rs)**2 - 1
!! @param mm  order of the derivative, 0..4

!Input variables
 real(dp), intent(in) :: xx
 integer, intent(in) :: mm

!Local variables
 real(dp) :: pl(0:npow,0:4)
 call legendre_derivs(xx,mm,pl(:,0:mm))
 leg_value=sum(cleg(:)*pl(:,mm))
 return
 end function leg_value

 subroutine legendre_derivs(xx,nder,pl)

!! the Legendre polynomials P_k(xx), k=0..npow, and their derivatives up
!! to order nder: P_(k+1) = ((2k+1) x P_k - k P_(k-1))/(k+1) and
!! P_(k+1)^(m) = P_(k-1)^(m) + (2k+1) P_k^(m-1)
!!
!! @param xx  argument
!! @param nder  highest derivative
!! @param pl  pl(k,m) = P_k^(m)(xx)

!Input variables
 real(dp), intent(in) :: xx
 integer, intent(in) :: nder

!Output variables
 real(dp), intent(out) :: pl(0:npow,0:nder)

!Local variables
 integer :: kk,mm
 pl(:,:)=0.0d0
 pl(0,0)=1.0d0
 pl(1,0)=xx
 do kk=1,npow-1
   pl(kk+1,0)=((2*kk+1)*xx*pl(kk,0) - kk*pl(kk-1,0))/(kk+1)
 end do
 do mm=1,nder
   pl(1,mm)=pl(0,mm-1)
   do kk=1,npow-1
     pl(kk+1,mm)=pl(kk-1,mm) + (2*kk+1)*pl(kk,mm-1)
   end do
 end do
 return
 end subroutine legendre_derivs

 subroutine legendre_powers(shifted,cl,cm)

!! the coefficients cm of the Legendre series cl in powers of x, or of t
!! for x = 2t - 1 if shifted
!!
!! @param shifted  powers of t = (x+1)/2 instead of x
!! @param cl  Legendre coefficients
!! @param cm  coefficients of the powers

!Input variables
 logical, intent(in) :: shifted
 real(dp), intent(in) :: cl(0:npow)

!Output variables
 real(dp), intent(out) :: cm(0:npow)

!Local variables
 real(dp) :: pk(0:npow,0:npow),xc(0:1)
 integer :: kk
! x as a polynomial, x = xc(0) + xc(1) t
 xc=[0.0d0,1.0d0]
 if(shifted) xc=[-1.0d0,2.0d0]
 pk(:,:)=0.0d0
 pk(0,0)=1.0d0
 pk(0:1,1)=xc(:)
 do kk=1,npow-1
   pk(:,kk+1)=-kk*pk(:,kk-1)
   pk(0:npow,kk+1)=pk(0:npow,kk+1) + (2*kk+1)*xc(0)*pk(0:npow,kk)
   pk(1:npow,kk+1)=pk(1:npow,kk+1) + (2*kk+1)*xc(1)*pk(0:npow-1,kk)
   pk(:,kk+1)=pk(:,kk+1)/(kk+1)
 end do
 cm(:)=matmul(pk,cl)
 return
 end subroutine legendre_powers

 real(dp) function poly(tt)

!! the minimum tail at (r/r_s)**2 = tt
!!
!! @param tt  (r/r_s)**2

!Input variables
 real(dp), intent(in) :: tt
 poly=leg_value(2.0d0*tt-1.0d0,0)
 return
 end function poly

 subroutine min_tail(dm)

!! the Legendre coefficients cleg of the minimum tail of radius rs: the least
!! squares of W (A c + b) with C c = d (LAPACK dgglse), where A c is the
!! transform of the polynomial at the window nodes (Gauss-Legendre on
!! [0, r_s]), b that of the core charge beyond r_s (tabw and the part of
!! the mesh interval at r_s), W
!! the square root of the weights, and C c = d the match of the core
!! charge and its first two derivatives at r_s
!!
!! @param dm  the core charge and its first two derivatives at r_s

!Input variables
 real(dp), intent(in) :: dm(3)

!Local variables
 real(dp), allocatable :: amat(:,:),cvec(:),wk(:)
 real(dp) :: bmat(3,0:npow),dvec(3),tt,qx,jx,xx,wx
 integer :: ii,jj,kk,info,lwork
 lwork=nqn+npow+1+64*(nqn+npow+1)
 allocate(amat(nqn,0:npow),cvec(nqn),wk(lwork))

 amat(:,:)=0.0d0
 do ii=1,nqn
   do jj=1,nxg
     xx=rs*xg(jj)
     wx=rs*wxg(jj)
     qx=qn(ii)*xx
     jx=wx*xx**2*j0(qx)
     do kk=0,npow
       amat(ii,kk)=amat(ii,kk) + jx*plg(jj,kk)
     end do
   end do
   tt=tabw(ii,irow)
   do jj=1,ngi
     tt=tt + pg(jj)*sin(qn(ii)*pt(jj))
   end do
   cvec(ii)=-sqrt(wn(ii))*tt/qn(ii)
   amat(ii,:)=sqrt(wn(ii))*amat(ii,:)
 end do

! the match, row m scaled by rs**m: with x = 2(r/rs)**2 - 1, rs d/dr =
! 4 d/dx and rs**2 d2/dr2 = 16 d2/dx2 + 4 d/dx at x = 1, where
! P_k' = k(k+1)/2 and P_k'' = (k-1)k(k+1)(k+2)/8
 do kk=0,npow
   bmat(1,kk)=1.0d0
   bmat(2,kk)=2.0d0*kk*(kk+1)
   bmat(3,kk)=2.0d0*(kk-1)*kk*(kk+1)*(kk+2) + 2.0d0*kk*(kk+1)
 end do
 dvec(:)=[dm(1),dm(2)*rs,dm(3)*rs**2]

 call dgglse(nqn,npow+1,3,amat,nqn,bmat,3,cvec,dvec,cleg,wk,lwork,info)
 if(info/=0) then
   write(6,'(/a,i4)') 'modcore_budget: ERROR stop - dgglse info =',info
   error stop 1
 end if
 deallocate(amat,cvec,wk)
 return
 end subroutine min_tail

 real(dp) function outer_transform(kk)

!! int_{r_s}^infinity r**2 rhoc j0(kk r) dr
!!
!! @param kk  wave number (bohr**-1), > 0

!Input variables
 real(dp), intent(in) :: kk

!Local variables
 real(dp) :: ss
 integer :: jj
 ss=0.0d0
 do jj=1,ngi
   ss=ss + pg(jj)*sin(kk*pt(jj))
 end do
 do jj=js,nout
   ss=ss + og(jj)*sin(kk*ot(jj))
 end do
 outer_transform=ss/kk
 return
 end function outer_transform

 real(dp) function model_transform(kk)

!! int_0^infinity r**2 n j0(kk r) dr of the model n
!!
!! @param kk  wave number (bohr**-1), > 0

!Input variables
 real(dp), intent(in) :: kk

!Local variables
 real(dp) :: ss
 integer :: jj
 ss=0.0d0
 do jj=1,nxg
   ss=ss + fw(jj)*sin(kk*fx(jj))
 end do
 model_transform=ss/kk + outer_transform(kk)
 return
 end function model_transform

 subroutine spectrum(qe)

!! the model's transform ntk at the nodes of the panels of width dq up to
!! qe (and the next grid point), and the band-limited models rqs(:,m) at
!! the density cutoffs m*qstep: n_q(r) = (2/pi) int_0^q k**2 n(k) j0(kr)
!! dk.  The sines at the nodes of each family (the same node of every
!! panel) rotated from one panel to the next, those of the partial
!! interval at rs computed
!!
!! @param qe  largest density cutoff needed (bohr**-1)

!Input variables
 real(dp), intent(in) :: qe

!Local variables
 real(dp) :: kq,ss,tt
 integer :: mlast,ip,ig,jj,ll
 mlast=min(nq,ceiling(qe/qstep)+1)
 if(msnap>=mlast) return
 if(npdone==0) then
   sr(:,:)=sr0(:,:)
   cr(:,:)=cr0(:,:)
   rqs(:,0)=0.0d0
   rqp(:)=0.0d0
 end if
 do ip=npdone+1,mlast*npanel
   do ig=1,ngk
     kq=(ip-1+tk(ig))*dq
     ss=0.0d0
     do jj=1,nxg
       ss=ss + fw(jj)*sin(kq*fx(jj))
     end do
     do jj=1,ngi
       ss=ss + pg(jj)*sin(kq*pt(jj))
     end do
     ntk(ip,ig)=(ss + tab(ip,ig,irow))/kq
! the band-limited model: (2/pi) w k n(k) sin(kr)/r
     tt=2.0d0/pi*wk(ig)*dq*kq*ntk(ip,ig)
     do ll=1,nt
       rqp(ll)=rqp(ll) + tt*sr(ll,ig)/xr(ll)
     end do
! the next panel
     do ll=1,nt
       ss=sr(ll,ig)*chr(ll) + cr(ll,ig)*shr(ll)
       cr(ll,ig)=cr(ll,ig)*chr(ll) - sr(ll,ig)*shr(ll)
       sr(ll,ig)=ss
     end do
   end do
   if(mod(ip,npanel)==0) rqs(:,ip/npanel)=rqp(:)
 end do
 npdone=mlast*npanel
 msnap=mlast
 return
 end subroutine spectrum

 real(dp) function err_at(qq)

!! the error of the density cutoff qq: the band-limited model from the
!! last grid point below qq, the whole panels after it and the part of
!! the last panel up to qq, at its own Gauss-Legendre nodes
!!
!! @param qq  density cutoff (bohr**-1), 0 <= qq <= qmax

!Input variables
 real(dp), intent(in) :: qq

!Local variables
 real(dp) :: bq(nt),kq,tt,qa,ww
 integer :: m0,ip,ig,pf
 call spectrum(qq)
 m0=min(msnap,max(0,floor(qq/qstep)))
 bq(:)=rqs(:,m0)
 pf=max(m0*npanel,floor(qq/dq))
 do ip=m0*npanel+1,pf
   do ig=1,ngk
     kq=(ip-1+tk(ig))*dq
     tt=2.0d0/pi*wk(ig)*dq*kq*ntk(ip,ig)
     bq(:)=bq(:) + tt*sin(kq*xr(:))/xr(:)
   end do
 end do
 qa=pf*dq
 ww=qq-qa
 if(ww>0.0d0) then
   do ig=1,ngk
     kq=qa + ww*tk(ig)
     tt=2.0d0/pi*wk(ig)*ww*kq*model_transform(kq)
     bq(:)=bq(:) + tt*sin(kq*xr(:))/xr(:)
   end do
 end if
 err_at=xc_error(bq)
 return
 end function err_at

 real(dp) function err_grid(mm)

!! the error of the density cutoff mm*qstep of the grid
!!
!! @param mm  index of the cutoff

!Input variables
 integer, intent(in) :: mm
 call spectrum(mm*qstep)
 err_grid=xc_error(rqs(:,mm))
 return
 end function err_grid

 real(dp) function err_max(qa,qe,elim)

!! the largest error over the density cutoffs qa..qe: on the grid qa,
!! m*qstep, qe; if that is above elim, the grid's largest (only whether
!! the error exceeds elim is wanted); otherwise refined by Brent's method
!! between the neighbours of the grid's largest, and of each other local
!! maximum of the grid whose parabola peaks within 1% of the result.
!! qpk the grid point of the largest
!!
!! @param qa  first density cutoff (bohr**-1)
!! @param qe  last density cutoff (bohr**-1)
!! @param elim  error limit (Ha)

!Input variables
 real(dp), intent(in) :: qa,qe,elim

!Local variables
 real(dp), allocatable :: tq(:),eq(:)
 integer :: nn,mm,ii,ilo,ihi,imx
 call spectrum(qe)
 allocate(tq(nint((qe-qa)/qstep)+3),eq(nint((qe-qa)/qstep)+3))
 nn=1
 tq(1)=qa
 eq(1)=err_at(qa)
 do mm=floor(qa/qstep)+1,nq
   if(mm*qstep>=qe) exit
   if(mm*qstep<=qa) cycle
   nn=nn+1
   tq(nn)=mm*qstep
   eq(nn)=err_grid(mm)
 end do
 if(qe>qa) then
   nn=nn+1
   tq(nn)=qe
   eq(nn)=err_at(qe)
 end if
 imx=maxloc(eq(1:nn),1)
 qpk=tq(imx)
 err_max=eq(imx)
 if(err_max>elim .or. nn==1) then
   deallocate(tq,eq)
   return
 end if
 err_max=max(err_max,-ug_fmin(neg_err,tq(max(1,imx-1)),tq(min(nn,imx+1))))
 do ii=1,nn
   if(ii==imx) cycle
   ilo=max(1,ii-1)
   ihi=min(nn,ii+1)
   if(eq(ii)<eq(ilo) .or. eq(ii)<eq(ihi)) cycle
   if(nn>=3) then
     if(parabola_peak(tq,eq,nn,ii)<0.99d0*err_max) cycle
   end if
   err_max=max(err_max,-ug_fmin(neg_err,tq(ilo),tq(ihi)))
 end do
 deallocate(tq,eq)
 return
 end function err_max

 real(dp) function parabola_peak(tq,eq,nn,ii)

!! the largest value on [tq(ii-1), tq(ii+1)] of the parabola through the
!! three points of tq, eq nearest tq(ii)
!!
!! @param tq  abscissae, rising
!! @param eq  values
!! @param nn  number of points, >= 3
!! @param ii  index of the point

!Input variables
 integer, intent(in) :: nn,ii
 real(dp), intent(in) :: tq(nn),eq(nn)

!Local variables
 real(dp) :: x0,x1,x2,y0,y1,y2,d1,d2,aa,bb,xv
 integer :: kc
 kc=min(nn-1,max(2,ii))
 x0=tq(kc-1)
 x1=tq(kc)
 x2=tq(kc+1)
 y0=eq(kc-1)
 y1=eq(kc)
 y2=eq(kc+1)
 parabola_peak=eq(ii)
 d1=(y1-y0)/(x1-x0)
 d2=(y2-y1)/(x2-x1)
 aa=(d2-d1)/(x2-x0)
 if(aa>=0.0d0) return
 bb=d1 - aa*(x0+x1)
 xv=-bb/(2.0d0*aa)
 if(xv<tq(max(1,ii-1)) .or. xv>tq(min(nn,ii+1))) return
 parabola_peak=max(eq(ii),y0 + (xv-x0)*(d1 + aa*(xv-x1)))
 return
 end function parabola_peak

 real(dp) function neg_err(qq)

!! minus the error of the density cutoff qq
!!
!! @param qq  density cutoff (bohr**-1)

!Input variables
 real(dp), intent(in) :: qq
 neg_err=-err_at(qq)
 return
 end function neg_err

 real(dp) function zeroin(fun,xa0,xb0)

!! a root of fun between xa0 and xb0, where fun changes sign, by Brent's
!! method with no tolerance of its own (ug_brent): it stops when the
!! bracket is a few rounding units wide (as lm_root). Where fun jumps
!! across zero (the fill type changes, a discrete event of
!! doc/smoothness_plan.md) that is the root, and the jump is reported; the
!! iteration limit is an error
!!
!! @param fun  the function
!! @param xa0,xb0  the bracket

!Input variables
 procedure(scalar_fn) :: fun
 real(dp), intent(in) :: xa0,xb0

!Local variables
 real(dp) :: fb,fc,fe(2)
 integer :: ierr
 logical :: conv

 call ug_brent(fun,xa0,xb0,zeroin,ierr,fq=fb,fo=fc,fends=fe,converged=conv)
 if(ierr/=0) then
   write(6,'(/a,1p,2e14.6)') 'modcore_budget: ERROR stop - no root &
&between',xa0,xb0
   error stop 1
 end if
 if(.not. conv) then
   write(6,'(/a,es12.4,a,es12.4)') &
&        'modcore_budget: ERROR zeroin not converged in 200 steps, x=', &
&        zeroin,', f=',fb
   error stop 1
 end if
! a continuous function is at its round-off on a bracket a few ulps wide
 if(fb/=0.0d0 .and. min(abs(fb),abs(fc))>1.0d-6*max(abs(fe(1)),abs(fe(2)))) &
& then
   write(6,'(/2a,es12.4,a,es12.4,a,es12.4)') 'modcore_budget: WARNING', &
&        ' root at a jump of the function, x=',zeroin,', f from ',fc,' to ',fb
 end if
 return
 end function zeroin

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

!! the minimum tail and its first 4 derivatives at rad (beyond rs, its
!! continuation, for the pieces of the model): those of the
!! Legendre series in x = 2t - 1, t = (r/rs)**2, by the chain rule with
!! t' = 2r/rs**2, t'' = 2/rs**2 (d/dt = 2 d/dx)
!!
!! @param rad  radius (bohr)
!! @param rs  model radius r_s (bohr)
!! @param dd  the model and its first 4 derivatives at rad

!Input variables
 real(dp), intent(in) :: rad,rs

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 real(dp) :: pl(0:npow,0:4),pt(0:4),t1,t2
 integer :: mm
 call legendre_derivs(2.0d0*(rad/rs)**2-1.0d0,4,pl)
 do mm=0,4
   pt(mm)=2.0d0**mm*sum(cleg(:)*pl(:,mm))
 end do
 t1=2.0d0*rad/rs**2
 t2=2.0d0/rs**2
 dd(1)=pt(0)
 dd(2)=pt(1)*t1
 dd(3)=pt(2)*t1**2 + pt(1)*t2
 dd(4)=pt(3)*t1**3 + 3.0d0*pt(2)*t1*t2
 dd(5)=pt(4)*t1**4 + 6.0d0*pt(3)*t1**2*t2 + 3.0d0*pt(2)*t2**2
 return
 end subroutine tail_derivatives

 subroutine fill_derivatives(rad,dd)

!! the fill of the model and its first 4 derivatives at rad
!!
!! @param rad  radius (bohr)
!! @param dd  the fill and its first 4 derivatives at rad

!Input variables
 real(dp), intent(in) :: rad

!Output variables
 real(dp), intent(out) :: dd(5)

 if(ifill==1) then
   call tail_derivatives(rad,rs,dd)
 else
   call martins_derivatives(rad,dd)
 end if
 return
 end subroutine fill_derivatives

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

 real(dp) function xc_error(bq)

!! the error of the band-limited model bq of the model: the integral
!! split at rs (lm_integral), of the fill's terms below and the core
!! charge's above, each taken point by point, so that no quadrature spans
!! the joint and the error is not the small difference of two sums over
!! the sphere
!!
!! The part below rs is taken on the mesh ri of step al ending at rs, so
!! that the fill, a polynomial in r/rs, is sampled at the same relative
!! points whatever rs: on the fixed mesh its quadrature error rippled as
!! rs slid over the mesh (4e-12 on an error of 1e-4, period ~0.003 bohr,
!! ~1e-9 in r_s; doc/smoothness_plan.md, M-12)
!!
!! @param bq  band-limited model core charge at xr: rr(1:nx), then ri

!Input variables
 real(dp), intent(in) :: bq(nt)

!Local variables
 real(dp) :: si,so1
 work(:)=nval(:) + bq(1:nx)
 call exc_dispatch(iexc,al,work,vxc,eps,rr,nx,'modcore_budget')
 work(:)=work(:)*eps(:)
 ac(:)=swin(:)*rr(1:nx)**2*(work(:) - wc(:) - vbg*(bq(1:nx) - rhoc(1:nx)))
 wki(:)=nvi(:) + bq(nx+1:nt)
 call exc_dispatch(iexc,al,wki,vxi,epi,ri,ni,'modcore_budget')
 wki(:)=wki(:)*epi(:)
 af(:)=ri(:)**2*(wki(:) - wf(:) - vbg*(bq(nx+1:nt) - fext(:)))
 call lm_integral(ri,af,ni,2,0.0d0,rs,si)
 call lm_integral(rr,ac,nx,2,rs,radius,so1)
 xc_error=abs(si + so1)
 return
 end function xc_error

 end subroutine modcore_budget


 subroutine modcore_probes(lmax,lloc,nproj,epa,wellprj,vp,vlocpw,vkbpw, &
&                         evkb,rr,mmax,mxprj,nkap,probes,lprobe,nprobe)

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
!! @param wellprj  whether each projector is a well state (its energy then
!!                 of either sign), else a bound state
!! @param vp  screened semi-local and local pseudopotentials
!! @param vlocpw  the screened local potential's pieces
!! @param vkbpw  Vanderbilt-Kleinman-Bylander projectors as functions of
!!               pieces
!! @param evkb  coefficients of the projectors
!! @param rr  log radial grid
!! @param mmax  dimension of log grid
!! @param mxprj  dimension of number of projector
!! @param nkap  1 for scalar-relativistic or non-relativistic, 2 for
!!              relativistic
!! @param probes  probe densities, (u/r)**2 normalized to one electron
!! @param lprobe  angular momenta of the probes
!! @param nprobe  number of probes

 use constants_m, only: dp, mxl, mxlin
 use pieces_m, only: piecewise_t
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,mxprj,nkap
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap)
 real(dp), intent(in) :: epa(mxprj,mxlin+1,nkap),evkb(mxprj,mxl+1,nkap)
 logical, intent(in) :: wellprj(mxprj,mxlin+1)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)

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
       if(.not. wellprj(iprj,l1)) then
         emax=min(et+0.1d0*abs(et)+0.05d0,0.0d0)
         emin=et-0.1d0*abs(et)-0.05d0
       else
! a scattering reference: its state, if bound, below zero
         emax=0.0d0
         emin=minval(vp(:,lloc+1,ikap))
         et=0.5d0*emin
       end if
       call lschvkbb(ll+iprj,ll,npr,ierr,et,emin,emax, &
&                    rr,vp(1,lloc+1,ikap),vlocpw,vkbpw(1:npr,l1,ikap), &
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

 subroutine modcore_budget_table(ecurve,ncurve,nprof,wcurve,weighted, &
&                                allwaves,errmax,ecut)

!! the cutoff ecut (Ha) of the pseudo-waves at the error errmax (Ha per
!! electron).
!! By default (weighted and allwaves false) the largest over the
!! ground-state pseudo-waves, those with a weight (the waves of occupied
!! states, and the first of each channel with none occupied), of the
!! cutoff at which the energy error per electron of each falls to errmax,
!! the "Energy error per electron / Cutoff" tables of run_optimize (the
!! same curve and the same root-finding on the interpolant of log E over
!! its points, ug_cross_down, from the last crossing); a pseudo-wave whose
!! table does not reach errmax does not count, as it has no line in its
!! table. The other waves, of well states at the energies debl adds,
!! carry no charge in the ground state; with allwaves (budget_cutoff =
!! "max_all") they count too, as before.
!! With weighted true (budget_cutoff = "weighted"), the cutoff at which
!! the weighted mean sum_i w_i E_i / sum_i w_i of the curves with a weight
!! falls to errmax: each curve's log E interpolated as in ug_cross_down
!! up to its end or its first value that is not positive, and held at its
!! last value beyond, from where all of them start to where the last ends;
!! the last crossing from
!! above on a grid of the finest spacing of the curves, and the root of
!! the log of the mean less log errmax in it by Brent's method (ug_brent).
!! Unlike the largest cutoff, the mean is smooth where two curves cross.
!! ecurve(:,:,i) holds the cutoffs and errors of pseudo-wave i on a uniform
!! grid in q from below its spacing (eresid), zero beyond its end (all zero
!! for no pseudo-wave)
!!
!! @param ecurve  residual kinetic energy curves: ecurve(1,:,i) cutoffs (Ha),
!!                ecurve(2,:,i) error per electron (Ha) of pseudo-wave i
!! @param ncurve  number of points per curve
!! @param nprof  number of curves (pseudo-waves)
!! @param wcurve  weight of each curve in the weighted mean (weighted only)
!! @param weighted  whether ecut is that of the weighted mean
!! @param allwaves  whether the largest is over all waves, not only those
!!                  with a weight (not weighted)
!! @param errmax  error per electron (Ha) that defines the cutoff
!! @param ecut  cutoff (Ha) of the pseudo-waves at errmax

 use constants_m, only: dp
 use logmesh_m, only: lm_fh_value
 use ugrid_m, only: ug_cross_down, ug_brent
 implicit none

!Input variables
 integer, intent(in) :: ncurve,nprof
 real(dp), intent(in) :: ecurve(2,ncurve,nprof),wcurve(nprof),errmax
 logical, intent(in) :: weighted,allwaves

!Output variables
 real(dp), intent(out) :: ecut

!Local variables
 real(dp) :: qerr,ecut1,qa,qb,hs,xa,xb,fa,fb
 real(dp) :: q1s(nprof),dqs(nprof)
 real(dp), allocatable :: lf(:,:)
 integer :: ip,kk,nk,ns,ierr
 integer :: nks(nprof)
 logical :: found,found1

 ecut=0.0d0
 found=.false.
 allocate(lf(ncurve,nprof))
! the curve's points, q=q1+(k-1)*dq, up to its end (and, for the
! weighted mean, up to its first value that is not positive)
 do ip=1,nprof
   nk=ncurve
   do kk=2,ncurve
     if(ecurve(1,kk,ip)<=0.0d0) then
       nk=kk-1
       exit
     end if
   end do
   nks(ip)=nk
   if(nk<2) cycle
   q1s(ip)=sqrt(2.0d0*ecurve(1,1,ip))
   dqs(ip)=sqrt(2.0d0*ecurve(1,2,ip))-q1s(ip)
   if(weighted) then
     do kk=1,nk
       if(.not. (ecurve(2,kk,ip)>0.0d0)) then
         nks(ip)=kk-1
         exit
       end if
       lf(kk,ip)=log(ecurve(2,kk,ip))
     end do
   end if
 end do

 if(.not. weighted) then
   do ip=1,nprof
     if(.not. allwaves .and. .not. (wcurve(ip)>0.0d0)) cycle
     ecut1=-1.0d0
     if(nks(ip)>=2) then
       call ug_cross_down(ecurve(2,1:nks(ip),ip),nks(ip),dqs(ip),q1s(ip), &
&                         errmax,qerr,found1)
       if(found1) ecut1=0.5d0*qerr**2
     end if
     if(ecut1>=0.0d0) then
       found=.true.
       ecut=max(ecut,ecut1)
     end if
   end do
   if(.not. found) then
     write(6,'(/a,1p,e10.2,a)') &
&          'modcore_budget_table: ERROR stop - no pseudo-wave reaches', &
&          errmax,' Ha per electron in its table'
     error stop 1
   end if
   return
 end if

! the weighted mean: from where all the curves with a weight start to
! where the last ends, each held at its last value past its end (the
! error falls with q, so the mean there is an upper bound); the s waves'
! curves end far below the p waves' at a small p rc
 qa=0.0d0
 qb=0.0d0
 hs=huge(1.0d0)
 do ip=1,nprof
   if(.not. (wcurve(ip)>0.0d0)) cycle
   if(nks(ip)<2) then
     write(6,'(/a,i0,a)') 'modcore_budget_table: ERROR stop - pseudo-wave ', &
&          ip,' has a weight and no convergence curve'
     error stop 1
   end if
   qa=max(qa,q1s(ip))
   qb=max(qb,q1s(ip)+dble(nks(ip)-1)*dqs(ip))
   hs=min(hs,dqs(ip))
 end do
 if(.not. (hs<huge(1.0d0)) .or. qa>=qb) then
   write(6,'(/a)') 'modcore_budget_table: ERROR stop - no weighted curves'
   error stop 1
 end if
 fb=excess_mean(qb)
 if(fb>0.0d0) then
   write(6,'(/a,1p,e10.2,a)') &
&        'modcore_budget_table: ERROR stop - the weighted mean does not'// &
&        ' reach',errmax,' Ha per electron'
   error stop 1
 end if
! the last crossing from above, from qb down
 ns=max(1,ceiling((qb-qa)/hs))
 xb=qb
 do kk=1,ns
   xa=max(qa,qb-dble(kk)*hs)
   fa=excess_mean(xa)
   if(fa>0.0d0) then
     call ug_brent(excess_mean,xa,xb,qerr,ierr)
     if(ierr/=0) then
       write(6,'(a)') 'modcore_budget_table: ERROR the bracket lost its'// &
&            ' sign change'
       error stop 1
     end if
     ecut=0.5d0*qerr**2
     return
   end if
   xb=xa
 end do
 write(6,'(/a,1p,e10.2,a)') &
&      'modcore_budget_table: ERROR stop - the weighted mean is below', &
&      errmax,' Ha per electron where the curves start'
 error stop 1

 contains

 real(dp) function excess_mean(qq)
!! log of the weighted mean of the curves at the wave vector qq, less
!! log errmax
 real(dp), intent(in) :: qq
 real(dp) :: tt,sw,se
 integer :: jp
 sw=0.0d0
 se=0.0d0
 do jp=1,nprof
   if(.not. (wcurve(jp)>0.0d0)) cycle
   tt=min(max((qq-q1s(jp))/dqs(jp),0.0d0),dble(nks(jp)-1))
   sw=sw+wcurve(jp)
   se=se+wcurve(jp)*exp(lm_fh_value(lf(1:nks(jp),jp),nks(jp),tt))
 end do
 excess_mean=log(se/sw)-log(errmax)
 end function excess_mean

 end subroutine modcore_budget_table
