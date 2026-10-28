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

!> The model cores of icmod 6 and 7 (modcore_budget), the probes of the
!> hardness of icmod 7 (modcore_probes), and the table of the pseudo-wave
!> cutoffs the budget starts from (modcore_budget_table).  The search of
!> modcore_budget runs in phases over its state, a budget_t.
 module modcore_budget_m

 use constants_m, only: dp, pi, nrhomod
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: modcore_budget,modcore_probes,modcore_budget_table

 real(dp), parameter :: pi4=4.0d0*pi !< 4 pi
!> background density of the error (bohr**-3)
 real(dp), parameter :: background=5.0d-3
!> radius of the sphere of the error (bohr)
 real(dp), parameter :: radius=10.0d0
!> where the taper of the sphere starts, as a share of its radius
 real(dp), parameter :: taper=0.6d0
!> largest density cutoff (bohr**-1)
 real(dp), parameter :: qmax=120.0d0
!> step of the grid of density cutoffs (bohr**-1)
 real(dp), parameter :: qstep=0.1d0
!> the cutoffs past q_b of the first search (bohr**-1)
 real(dp), parameter :: quick=1.0d0
!> panels of the band-limited transform per step
 integer, parameter :: npanel=1
!> Gauss-Legendre nodes of each panel of the band-limited transform
 integer, parameter :: ngk=5
!> nodes of the transform on [0, r_s]
 integer, parameter :: nxg=400
!> nodes per mesh interval beyond r_s (each of these converged: twice as
!> many move nothing printed)
 integer, parameter :: ngi=8
!> highest power (r/r_s)**(2*npow) of the minimum tail
 integer, parameter :: npow=8
!> Gauss-Legendre nodes per panel of the window of the minimum tail
 integer, parameter :: ngq=8
!> panels of the window of the minimum tail
 integer, parameter :: nwin=80
!> icmod=7: ratio of the radii of the basin
 real(dp), parameter :: ratio=1.025d0
!> icmod=7: step in ln r_s of the central differences
 real(dp), parameter :: hmis=1.0d-2
!> names of the fills: none, minimum tail, Martins
 character(len=2), parameter :: fillname(0:2)=['- ','MT','M ']

 abstract interface
   real(dp) function scalar_fn(xx)
   import :: dp
   implicit none
   real(dp), intent(in) :: xx
   end function scalar_fn
 end interface

!> The fill of a model core inside its radius r_s
 type :: fill_t
!> the model radius r_s (bohr)
   real(dp) :: rs=0.0d0
!> the fill's kind: 0 none, 1 the minimum tail, 2 Martins'
   integer :: ifill=0
!> the minimum tail's Legendre coefficients in x = 2(r/r_s)**2 - 1
   real(dp) :: cleg(0:npow)=0.0d0
!> the minimum tail's coefficients in powers of (r/r_s)**2 (for the
!> report)
   real(dp) :: cpoly(0:npow)=0.0d0
!> Martins' fill, exp(a + b r**2 + c r**4): a
   real(dp) :: aco=0.0d0
!> Martins' fill: b (bohr**-2)
   real(dp) :: bco=0.0d0
!> Martins' fill: c (bohr**-4)
   real(dp) :: cco=0.0d0
 end type fill_t

!> The state of the search of modcore_budget, which its phases and the
!> procedures of the model, its transform and its error share
 type :: budget_t
!> size of the log radial mesh
   integer :: mmax
!> exchange-correlation functional
   integer :: iexc
!> number of probes (0 for icmod=6)
   integer :: nprobe
!> speed of light (Hartree atomic units)
   real(dp) :: c_light
!> angular momenta of the probes
   integer, allocatable :: lprobe(:)
!> log radial mesh (bohr)
   real(dp), allocatable :: rr(:)
!> core charge density
   real(dp), allocatable :: rhoc(:)
!> pseudo valence charge density
   real(dp), allocatable :: rhotps(:)
!> probe densities of the hardness (icmod=7), one electron each
!> (mmax,nprobe)
   real(dp), allocatable :: probes(:,:)
!> exchange-correlation energy error (Ha) that defines the needed cutoff
   real(dp) :: errmax
!> ion charge (for der2exc)
   real(dp) :: zion
!> the budget's plane-wave cutoff ecut (Ha)
   real(dp) :: ecut
!> the share of ecut the model must fit, slack*ecut (Ha)
   real(dp) :: budget
!> the multiple of the pseudo-waves' cutoff ecutw that is ecut, when
!> fcfact <= 0
   real(dp) :: factor
!> the lower end of the window of the minimum tail, 2*sqrt(2*ecut)
!> (bohr**-1)
   real(dp) :: qlo
!> the upper end of the window, min(3*qlo, qmax) (bohr**-1)
   real(dp) :: qhi
!> q_b, the first density cutoff that decides whether the model fits
!> budget, 2*sqrt(2*budget) (bohr**-1)
   real(dp) :: qb
!> the last density cutoff of the search (bohr**-1)
   real(dp) :: qend
!> the qend r_s of icmod 6 was found with, for the report (bohr**-1)
   real(dp) :: qext
!> the grid's density cutoff of the largest error of err_max (bohr**-1)
   real(dp) :: qpk
!> the largest r_s, the smallest rc or the last radius with core charge
!> (bohr)
   real(dp) :: rtop
!> the smallest r_s, rtop/50 (bohr)
   real(dp) :: rlo
!> r_s of icmod 6, the smallest radius that fits (bohr)
   real(dp) :: rs6
!> the last mesh point with core charge
   integer :: irtop
!> log of the mesh ratio
   real(dp) :: al
!> the width of the panels of the band-limited transform (bohr**-1)
   real(dp) :: dq
!> the number of density cutoffs of the grid, qmax/qstep
   integer :: nq
!> the number of panels of the band-limited transform
   integer :: npq
!> the number of nodes of the window of the minimum tail
   integer :: nqn
!> the first mesh interval of the nodes beyond r_s
   integer :: i0
!> the mesh point after the last interval of the nodes beyond r_s
   integer :: iend
!> the number of nodes beyond r_s, ngi per interval
   integer :: nout
!> the last mesh point r_s can lie below on in the tables
   integer :: irmx
!> the last mesh point within the sphere of the error
   integer :: nin
!> the size of the error's part of the log mesh, rr(1:nx)
   integer :: nx
!> the index of r_s on the mesh ri
   integer :: kin
!> the size of the mesh ri
   integer :: ni
!> the size of xr, nx + ni
   integer :: nt
!> the model being tried: its radius and fill
   type(fill_t) :: fill
!> the row of the tables for the current r_s, the mesh point above it
   integer :: irow
!> the first node beyond r_s of the whole intervals
   integer :: js
!> the number of panels of ntk done for the current model
   integer :: npdone
!> the last grid cutoff of rqs done for the current model
   integer :: msnap
!> Gauss-Legendre nodes on [0, 1] of the panels in q
   real(dp) :: tk(ngk)
!> Gauss-Legendre weights of tk
   real(dp) :: wk(ngk)
!> Gauss-Legendre nodes on [0, 1] of each mesh interval in ln r
   real(dp) :: tgi(ngi)
!> Gauss-Legendre weights of tgi
   real(dp) :: wgi(ngi)
!> Gauss-Legendre nodes on [0, 1] of the transform on [0, r_s]
   real(dp) :: xg(nxg)
!> Gauss-Legendre weights of xg
   real(dp) :: wxg(nxg)
!> the nodes of the partial interval [r_s, the next mesh point] (bohr)
   real(dp) :: pt(ngi)
!> their weights times r**2 rhoc
   real(dp) :: pg(ngi)
!> the nodes of the window of the minimum tail (bohr**-1)
   real(dp), allocatable :: qn(:)
!> their weights times q**2
   real(dp), allocatable :: wn(:)
!> the nodes of the core charge's transform beyond r_s (bohr)
   real(dp), allocatable :: ot(:)
!> their weights times r**2 rhoc
   real(dp), allocatable :: og(:)
!> sin(dq ot), the rotation from one panel to the next
   real(dp), allocatable :: sho(:)
!> cos(dq ot)
   real(dp), allocatable :: cho(:)
!> sin(tk dq ot), the sines of the first panel (nout,ngk)
   real(dp), allocatable :: so0(:,:)
!> cos(tk dq ot) (nout,ngk)
   real(dp), allocatable :: co0(:,:)
!> the sums of og sin(k ot) at the nodes of the band-limited transform
!> over the intervals from each mesh point on (npq,ngk,i0:irmx)
   real(dp), allocatable :: tab(:,:,:)
!> the same at the nodes of the window (nqn,i0:irmx)
   real(dp), allocatable :: tabw(:,:)
!> the nodes of the fill's transform on [0, r_s] (bohr)
   real(dp), allocatable :: fx(:)
!> their weights times r times the fill
   real(dp), allocatable :: fw(:)
!> the Legendre polynomials of 2 xg**2 - 1 (nxg,0:npow)
   real(dp), allocatable :: plg(:,:)
!> the model inside r_s and the core charge beyond (mmax)
   real(dp), allocatable :: core(:)
!> the pseudo valence charge on the background, on rr(1:nx)
   real(dp), allocatable :: nval(:)
!> work space on rr(1:nx)
   real(dp), allocatable :: work(:)
!> the exchange-correlation energy per electron (Ha) on rr(1:nx)
   real(dp), allocatable :: eps(:)
!> exchange-correlation potential (Ha) on rr(1:nx)
   real(dp), allocatable :: vxc(:)
!> the xc energy density of the valence charge on the fill, on ri
   real(dp), allocatable :: wf(:)
!> the xc energy density of the valence charge on the core charge, on
!> rr(1:nx)
   real(dp), allocatable :: wc(:)
!> the fill on ri, extended beyond r_s
   real(dp), allocatable :: fext(:)
!> the integrand of the error below r_s, on ri
   real(dp), allocatable :: af(:)
!> the integrand of the error from r_s, on rr(1:nx)
   real(dp), allocatable :: ac(:)
!> the taper of the sphere of the error on rr(1:nx)
   real(dp), allocatable :: swin(:)
!> the points of the band-limited models: rr(1:nx), then ri (bohr)
   real(dp), allocatable :: xr(:)
!> the mesh of the integral below r_s, r_s exp((k-kin) al) (bohr)
   real(dp), allocatable :: ri(:)
!> the pseudo valence charge on the background, on ri
   real(dp), allocatable :: nvi(:)
!> work space on ri
   real(dp), allocatable :: wki(:)
!> the exchange-correlation energy per electron (Ha) on ri
   real(dp), allocatable :: epi(:)
!> exchange-correlation potential (Ha) on ri
   real(dp), allocatable :: vxi(:)
!> sin(dq xr), the rotation from one panel to the next
   real(dp), allocatable :: shr(:)
!> cos(dq xr)
   real(dp), allocatable :: chr(:)
!> sin(tk dq xr), the sines of the first panel (nt,ngk)
   real(dp), allocatable :: sr0(:,:)
!> cos(tk dq xr) (nt,ngk)
   real(dp), allocatable :: cr0(:,:)
!> the sines of the current panel, rotated from sr0 (nt,ngk)
   real(dp), allocatable :: sr(:,:)
!> the cosines of the current panel, rotated from cr0 (nt,ngk)
   real(dp), allocatable :: cr(:,:)
!> the model's transform at the nodes of the panels (npq,ngk)
   real(dp), allocatable :: ntk(:,:)
!> the band-limited models at xr at the grid cutoffs m*qstep (nt,0:nq)
   real(dp), allocatable :: rqs(:,:)
!> the band-limited model summed so far at xr; work space of the tables
   real(dp), allocatable :: rqp(:)
!> exchange-correlation potential of the background alone (Ha)
   real(dp) :: vbg
!> the needed density cutoff (bohr**-1)
   real(dp) :: qneed
!> the needed plane-wave cutoff, qneed**2/8 (Ha)
   real(dp) :: ecore
!> whether the search of icmod 6 was extended beyond q_b + quick
   logical :: full
!> icmod=7: the hardness at rhotps + rhoc (nprobe,nprobe)
   real(dp), allocatable :: d2ae(:,:)
!> the hardness at rhotps + the model (nprobe,nprobe)
   real(dp), allocatable :: d2mod(:,:)
!> zeros, the reference of the hardness at rhotps + rhoc (nprobe,nprobe)
   real(dp), allocatable :: d2zero(:,:)
!> the hardness mismatch (Ha)
   real(dp) :: mismatch
!> the radii of the basin (bohr)
   real(dp), allocatable :: crad(:)
!> their hardness mismatches (Ha)
   real(dp), allocatable :: cmis(:)
!> their fills, 0 for none
   integer, allocatable :: cfill(:)
!> the number of radii of the basin
   integer :: ncand
 end type budget_t

 contains

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
&                          weighted,errmax,slack,rcmin,mmax,rr,iexc,c_light, &
&                          zion,probes,lprobe,nprobe,iout,ifail)

!! @param icmod  6 for the smallest radius, 7 for the hardness-matched one
!! @param rhotps  pseudo valence charge density
!! @param rhoc  core charge density
!! @param rhomod  model core charge and its first 4 derivatives
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
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh (bohr)
!! @param iexc  exchange-correlation functional
!! @param c_light  speed of light (Hartree atomic units)
!! @param zion  ion charge (for der2exc)
!! @param probes  probe densities of the hardness (icmod=7), one electron each
!! @param lprobe  angular momenta of the probes
!! @param nprobe  number of probes (0 for icmod=6)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use modcore_m, only: core_pieces
 use pieces_m, only: piecewise_t

!Input variables
 integer, intent(in) :: icmod,mmax,iexc,nprobe
 real(dp), intent(in) :: c_light
 integer, intent(in) :: lprobe(nprobe)
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: probes(mmax,nprobe)
 real(dp), intent(in) :: fcfact,ecutw,errmax,slack,zion,rcmin
 logical, intent(in) :: weighted
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: rhomod(mmax,nrhomod)
 type(piecewise_t), intent(out) :: corepw(nrhomod)
 integer, intent(out) :: ifail

!Local variables
 type(budget_t) :: bc

 ifail=0
 call budget_init(bc,rhotps,rhoc,fcfact,ecutw,errmax,slack,rcmin,mmax,rr, &
&                 iexc,c_light,zion,probes,lprobe,nprobe,iout,ifail)
 if(ifail/=0) return
 call budget_tables(bc,iout,ifail)
 if(ifail/=0) return
 call budget_sphere(bc)
 call budget_search6(bc,iout,ifail)
 if(ifail/=0) return
 if(icmod==7) then
   call budget_hardness7(bc,iout,ifail)
   if(ifail/=0) return
 end if
 call budget_needed_cut(bc,iout,ifail)
 if(ifail/=0) return

! fill in model and derivatives, all-electron beyond r_s
 call core_pieces(rr,rhoc,mmax,bc%fill%rs,fill_form,fill_pack(bc%fill), &
&                  corepw,rhomod)

 call budget_report(bc,icmod,fcfact,weighted,ecutw,slack,rcmin,iout)

 return

 end subroutine modcore_budget

 subroutine budget_init(bc,rhotps,rhoc,fcfact,ecutw,errmax,slack,rcmin, &
&                       mmax,rr,iexc,c_light,zion,probes,lprobe,nprobe,iout, &
&                       ifail)

!! the cutoffs of the budget, the q grid of the band-limited
!! transforms, and the largest and smallest radii r_s; the
!! inputs are kept in bc
!!
!! @param bc  the state of the search
!! @param rhotps..nprobe  as in modcore_budget
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: mmax,iexc,nprobe
 integer, intent(in) :: lprobe(nprobe)
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: probes(mmax,nprobe)
 real(dp), intent(in) :: fcfact,ecutw,errmax,slack,zion,rcmin
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ii

 ifail=0
 bc%mmax=mmax
 bc%iexc=iexc
 bc%c_light=c_light
 bc%nprobe=nprobe
 bc%rr=rr
 bc%rhoc=rhoc
 bc%rhotps=rhotps
 bc%probes=probes
 bc%lprobe=lprobe
 bc%errmax=errmax
 bc%zion=zion

 if(fcfact>0.0d0) then
   bc%ecut=fcfact
 else
   bc%factor=1.0d0
   if(fcfact<0.0d0) bc%factor=-fcfact
   bc%ecut=bc%factor*ecutw
 end if
 bc%budget=slack*bc%ecut

 bc%al = 0.01d0 * log(bc%rr(101) / bc%rr(1))
 bc%dq = qstep/npanel
 bc%nq = nint(qmax/qstep)
 bc%npq = bc%nq*npanel

 bc%qlo=2.0d0*sqrt(2.0d0*bc%ecut)
 if(bc%qlo >= qmax) then
   write(iout,'(/a,f12.6,a)') 'modcore_budget: ERROR ecut beyond', &
&        qmax**2/8.0d0,' Ha'
   ifail=1
   return
 end if
 bc%qhi=min(3.0d0*bc%qlo,qmax)
! the first density cutoff that decides whether the model fits budget
 bc%qb=2.0d0*sqrt(2.0d0*bc%budget)

! the largest r_s: the smallest rc, or the last radius at which the core
! charge and the derivatives matched there are not zero
 bc%irtop=0
 do ii=1,bc%mmax-6
   if(.not. all(bc%rhoc(ii:ii+6)>0.0d0)) exit
   bc%irtop=ii
 end do
 if(bc%irtop<4) then
   write(iout,'(/a)') 'modcore_budget: ERROR no core charge'
   ifail=1
   return
 end if
 bc%rtop=min(rcmin,bc%rr(bc%irtop))
 bc%rlo=bc%rtop/50.0d0
 if(qmax*bc%rtop > 2.0d0*(nxg-50)) then
   write(iout,'(/a,f12.6,a)') 'modcore_budget: ERROR r_s up to', &
&        bc%rtop,' bohr needs more nodes of the transform'
   ifail=1
   return
 end if

 return
 end subroutine budget_init

 subroutine budget_tables(bc,iout,ifail)

!! the Gauss-Legendre nodes, the window of the minimum tail, the
!! tables of the core charge's transform beyond r_s and the
!! Legendre polynomials of the fill
!!
!! @param bc  the state of the search
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use logmesh_m, only: lm_interp
 use quad_m, only: quad_gauss_legendre

!Input variables
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ig,ii,ip,jj,kk
 real(dp) :: ca,sa,sb,tt

 ifail=0
 call quad_gauss_legendre(0.0d0,1.0d0,ngk,bc%tk,bc%wk)
 call quad_gauss_legendre(0.0d0,1.0d0,ngi,bc%tgi,bc%wgi)
 call quad_gauss_legendre(0.0d0,1.0d0,nxg,bc%xg,bc%wxg)

! the window of the minimum tail: nwin Gauss-Legendre panels on
! [qlo, qhi], the weights times q**2
 bc%nqn=ngq*nwin
 allocate(bc%qn(bc%nqn),bc%wn(bc%nqn))
 do ii=1,nwin
   call quad_gauss_legendre(bc%qlo+(ii-1)*(bc%qhi-bc%qlo)/nwin, &
&                      bc%qlo+ii*(bc%qhi-bc%qlo)/nwin, &
&                      ngq,bc%qn((ii-1)*ngq+1:ii*ngq), &
&                                bc%wn((ii-1)*ngq+1:ii*ngq))
 end do
 bc%wn(:)=bc%wn(:)*bc%qn(:)**2

! the core charge's transform beyond r_s: Gauss-Legendre in ln r on each
! mesh interval from below r_lo (and the steps of the central differences)
! to the last one where r**3 rhoc is above 1e-14 of its largest value,
! of the interpolant of rhoc; og the weights times r**2 rhoc, so that
! int r**2 rhoc j0(kr) dr = sum og sin(k ot)/k.  Their sums at the
! nodes of the band-limited transform and of the window, from each mesh
! point r_s can lie below on: tab(:,:,i), tabw(:,i) over the intervals
! from i on, the sines rotated from one panel to the next
 tt=maxval(bc%rr(:)**3*bc%rhoc(:))
 bc%iend=1
 do ii=1,bc%mmax-2
   if(bc%rr(ii)**3*bc%rhoc(ii)>1.0d-14*tt) bc%iend=ii+1
 end do
 bc%i0=max(1,floor(log(bc%rlo/(1.05d0*bc%rr(1)))/bc%al))
 if(bc%iend<=bc%i0+1) then
   write(iout,'(/a)') 'modcore_budget: ERROR no core charge'
   ifail=1
   return
 end if
 bc%nout=ngi*(bc%iend-bc%i0)
 allocate(bc%ot(bc%nout),bc%og(bc%nout),bc%sho(bc%nout),bc%cho(bc%nout))
 allocate(bc%so0(bc%nout,ngk),bc%co0(bc%nout,ngk))
 do ii=bc%i0,bc%iend-1
   do kk=1,ngi
     jj=(ii-bc%i0)*ngi+kk
     bc%ot(jj)=bc%rr(ii)*exp(bc%al*bc%tgi(kk))
     call lm_interp(bc%rr,bc%rhoc,bc%mmax,bc%ot(jj),tt)
     bc%og(jj)=bc%al*bc%wgi(kk)*bc%ot(jj)**2*tt
   end do
 end do
 bc%sho(:)=sin(bc%dq*bc%ot(:))
 bc%cho(:)=cos(bc%dq*bc%ot(:))
 do kk=1,ngk
   bc%so0(:,kk)=sin(bc%tk(kk)*bc%dq*bc%ot(:))
   bc%co0(:,kk)=cos(bc%tk(kk)*bc%dq*bc%ot(:))
 end do
 bc%irmx=min(bc%iend,floor(log(bc%rtop*exp(2.0d0*hmis)/bc%rr(1))/bc%al)+3)
 allocate(bc%tab(bc%npq,ngk,bc%i0:bc%irmx),bc%tabw(bc%nqn,bc%i0:bc%irmx))
 allocate(bc%rqp(bc%npq*ngk))
 bc%rqp(:)=0.0d0
 bc%tabw(:,bc%i0)=0.0d0
 do ii=bc%iend-1,bc%i0,-1
   if(ii+1<=bc%irmx) then
     bc%tab(:,:,ii+1)=reshape(bc%rqp,[bc%npq,ngk])
     bc%tabw(:,ii+1)=bc%tabw(:,bc%i0)
   end if
   do kk=1,ngi
     jj=(ii-bc%i0)*ngi+kk
     do ig=1,ngk
       sa=bc%so0(jj,ig)
       ca=bc%co0(jj,ig)
       do ip=1,bc%npq
         bc%rqp(ip+(ig-1)*bc%npq)=bc%rqp(ip+(ig-1)*bc%npq) + bc%og(jj)*sa
         sb=sa*bc%cho(jj) + ca*bc%sho(jj)
         ca=ca*bc%cho(jj) - sa*bc%sho(jj)
         sa=sb
       end do
     end do
     bc%tabw(:,bc%i0)=bc%tabw(:,bc%i0) + bc%og(jj)*sin(bc%qn(:)*bc%ot(jj))
   end do
 end do
 bc%tab(:,:,bc%i0)=reshape(bc%rqp,[bc%npq,ngk])
 deallocate(bc%rqp)
 allocate(bc%fx(nxg),bc%fw(nxg),bc%plg(nxg,0:npow))
! the Legendre polynomials of the fill at the nodes on [0, r_s]
 do ii=1,nxg
   call legendre_derivs(2.0d0*bc%xg(ii)**2-1.0d0,0,bc%plg(ii,:))
 end do

 return
 end subroutine budget_tables

 subroutine budget_sphere(bc)

!! the sphere of the error, the meshes of its integrals, and the
!! xc energy density of the valence charge on the core charge and
!! the potential of the background
!!
!! @param bc  the state of the search

 use cinf_taper_m, only: cinf_taper
 use logmesh_m, only: lm_fh_reach, lm_blend
 use xc_dispatch_m, only: exc_dispatch

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 integer :: ii,kk

! the sphere of the error, and the points its integrals' interpolant
! reads beyond it
 bc%nin=0
 do ii=1,bc%mmax
   if(bc%rr(ii)>radius) exit
   bc%nin=ii
 end do
 bc%nx=min(bc%mmax,bc%nin+lm_fh_reach()+lm_blend+10)
! the mesh of the integral below r_s (set_model): ri(k) = r_s
! exp((k-kin) al), from at least rr(1) to the points its interpolant reads
! beyond r_s (r_s from r_lo, less the steps of the central differences).
! The band-limited models are taken at xr: rr(1:nx), then ri
 bc%kin=floor((log(bc%rlo/bc%rr(1)) - 2.0d0*hmis)/bc%al)+1
 bc%ni=bc%kin+lm_fh_reach()+lm_blend+10
 bc%nt=bc%nx+bc%ni
 allocate(bc%core(bc%mmax),bc%nval(bc%nx),bc%work(bc%nx),bc%eps(bc%nx), &
&         bc%vxc(bc%nx),bc%rqp(bc%nt))
 allocate(bc%wf(bc%ni),bc%wc(bc%nx),bc%fext(bc%ni),bc%af(bc%ni), &
&         bc%ac(bc%nx),bc%swin(bc%nx))
 allocate(bc%xr(bc%nt),bc%ri(bc%ni),bc%nvi(bc%ni),bc%wki(bc%ni), &
&         bc%epi(bc%ni),bc%vxi(bc%ni))
 do ii=1,bc%nx
   bc%swin(ii)=cinf_taper(bc%rr(ii),taper*radius,radius)
 end do
 allocate(bc%shr(bc%nt),bc%chr(bc%nt),bc%sr0(bc%nt,ngk),bc%cr0(bc%nt,ngk), &
&         bc%sr(bc%nt,ngk),bc%cr(bc%nt,ngk))
 allocate(bc%ntk(bc%npq,ngk),bc%rqs(bc%nt,0:bc%nq))
 bc%xr(1:bc%nx)=bc%rr(1:bc%nx)
 bc%shr(1:bc%nx)=sin(bc%dq*bc%rr(1:bc%nx))
 bc%chr(1:bc%nx)=cos(bc%dq*bc%rr(1:bc%nx))
 do kk=1,ngk
   bc%sr0(1:bc%nx,kk)=sin(bc%tk(kk)*bc%dq*bc%rr(1:bc%nx))
   bc%cr0(1:bc%nx,kk)=cos(bc%tk(kk)*bc%dq*bc%rr(1:bc%nx))
 end do

 bc%nval(:)=pi4*background + bc%rhotps(1:bc%nx)

! the xc energy density of the valence charge on the core charge
 bc%work(:)=bc%nval(:) + bc%rhoc(1:bc%nx)
 call exc_dispatch(bc%iexc,bc%c_light,bc%al,bc%work,bc%vxc,bc%eps,bc%rr,bc%nx, &
&                  'modcore_budget')
 bc%wc(:)=bc%work(:)*bc%eps(:)

! xc potential of the background alone
 bc%work(:)=pi4*background
 call exc_dispatch(bc%iexc,bc%c_light,bc%al,bc%work,bc%vxc,bc%eps,bc%rr,bc%nx, &
&                  'modcore_budget')
 bc%vbg=bc%vxc(bc%nx/2)

 return
 end subroutine budget_sphere

 subroutine budget_search6(bc,iout,ifail)

!! the smallest radius r_s whose model fits the budget (icmod 6)
!!
!! @param bc  the state of the search
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Output variables
 integer, intent(out) :: ifail

 ifail=0
! smallest r_s, first with the cutoffs to q_b + quick; if a larger one
! has a larger error at that r_s, again with those up to it + quick: each
! r_s found is a root of the largest error over all cutoffs less errmax
! once it is confirmed, whatever the range it was found with
 bc%full=.false.
 bc%qend=min(qmax,bc%qb+quick)
 do
   bc%fill%rs=smallest_radius(bc,iout,ifail)
   if(ifail/=0) return
   call set_model(bc,bc%fill%rs)
   if(bc%qend>=qmax) exit
   if(err_max(bc,bc%qend,qmax,bc%errmax)<=bc%errmax) exit
   bc%full=.true.
   bc%qend=min(qmax,bc%qpk+quick)
 end do
 bc%rs6=bc%fill%rs
 bc%qext=bc%qend
 bc%qend=qmax

 return
 end subroutine budget_search6

 subroutine budget_hardness7(bc,iout,ifail)

!! the radius of icmod 7: the minimum of the hardness mismatch
!! in the basin of the least of the radii from the largest one
!! down to r_s of icmod 6
!!
!! @param bc  the state of the search
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use der2exc_m, only: der2exc

!Input variables
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ib,ii,kk
 real(dp) :: tt

 ifail=0
 kk=max(0,floor(log(bc%rtop/bc%rs6)/log(ratio)))
 allocate(bc%crad(kk+2),bc%cfill(kk+2),bc%cmis(kk+2))
 bc%ncand=0
 do ii=0,kk
   tt=bc%rtop/ratio**ii
   if(tt<=bc%rs6) exit
   bc%ncand=bc%ncand+1
   bc%crad(bc%ncand)=tt
 end do
 bc%ncand=bc%ncand+1
 bc%crad(bc%ncand)=bc%rs6

 allocate(bc%d2ae(bc%nprobe,bc%nprobe),bc%d2mod(bc%nprobe,bc%nprobe), &
&         bc%d2zero(bc%nprobe,bc%nprobe))
 bc%d2zero(:,:)=0.0d0
 call der2exc(bc%rhotps,bc%rhoc,bc%probes,bc%rr,bc%d2ae,bc%d2zero,bc%mismatch, &
&               bc%zion,bc%iexc,bc%c_light,0,bc%nprobe,bc%lprobe, &
&               bc%rr(bc%mmax),bc%mmax)
 ib=0
 do ii=1,bc%ncand
   bc%cmis(ii)=mismatch_at(bc,log(bc%crad(ii)))
   bc%cfill(ii)=bc%fill%ifill
   if(bc%fill%ifill==0) cycle
   if(ib==0) then
     ib=ii
   else if(bc%cmis(ii)<bc%cmis(ib)) then
     ib=ii
   end if
 end do
 if(ib==0) then
   call no_model(bc,iout,ifail)
   if(ifail/=0) return
 end if
 bc%fill%rs=minimum_radius(bc,ib,iout,ifail)
 if(ifail/=0) return
 call set_model(bc,bc%fill%rs)
 bc%mismatch=mismatch_at(bc,log(bc%fill%rs))
! r_s6 is a root of this error less errmax, so at r_s = r_s6 the two are
! equal to round-off; the margin is above that, far below a radius that
! does not fit
 if(err_max(bc,bc%qb,qmax,bc%errmax)>bc%errmax*(1.0d0+1.0d-7)) then
   write(iout,'(/a,f12.6,a)') 'modcore_budget: WARNING the radius',bc%fill%rs, &
&          ' of the least mismatch does not fit the budget'
 end if

 return
 end subroutine budget_hardness7

 subroutine budget_needed_cut(bc,iout,ifail)

!! the cutoff the model needs: q**2/8 for the last density
!! cutoff q below q_b whose error is errmax
!!
!! @param bc  the state of the search
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: kk
 real(dp) :: ee,tt

 ifail=0
 call set_model(bc,bc%fill%rs)
 call spectrum(bc,bc%qb)
 bc%qneed=0.0d0
 tt=bc%qb
 ee=err_at(bc,bc%qb)-bc%errmax
 do kk=ceiling(bc%qb/qstep)-1,1,-1
   if(ee>=0.0d0) then
     bc%qneed=tt
     exit
   end if
   if(err_grid(bc,kk)>bc%errmax) then
     bc%qneed=zeroin(excess_of,kk*qstep,tt,iout,ifail)
     if(ifail/=0) return
     exit
   end if
   tt=kk*qstep
 end do
 bc%ecore=bc%qneed**2/8.0d0

 return

 contains

 real(dp) function excess_of(xx)

!! excess of this search, for zeroin
!!
!! @param xx  the argument

!Input variables
 real(dp), intent(in) :: xx

 excess_of=excess(bc,xx)
 return
 end function excess_of

 end subroutine budget_needed_cut

 subroutine budget_report(bc,icmod,fcfact,weighted,ecutw,slack,rcmin,iout)

!! the report of the model core
!!
!! @param bc  the state of the search
!! @param icmod..rcmin  as in modcore_budget
!! @param iout  unit of the text output

!Input variables
 integer, intent(in) :: icmod
 real(dp), intent(in) :: fcfact,ecutw,slack,rcmin
 logical, intent(in) :: weighted
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 integer :: ii

 if(icmod==7) then
   write(iout,'(/a)') 'Model core correction: hardness-matched core (icmod 7)'
 else
   write(iout,'(/a)') 'Model core correction: budgeted core (icmod 6)'
 end if
 if(fcfact<=0.0d0 .and. weighted) then
   write(iout,'(a,1p,e10.2,a,0p,f12.4,a,f8.4)') &
&        '  pseudo-wave cutoff, weighted mean Ecut of the tables at', &
& bc%errmax, &
&        ' Ha/electron',ecutw,' Ha, times',bc%factor
 else if(fcfact<=0.0d0) then
   write(iout,'(a,1p,e10.2,a,0p,f12.4,a,f8.4)') &
&        '  pseudo-wave cutoff, largest Ecut of the tables at',bc%errmax, &
&        ' Ha/electron',ecutw,' Ha, times',bc%factor
 end if
 write(iout,'(a,f12.4,a,1p,e10.2,a)') '  cutoff budget',bc%ecut, &
&      ' Ha for an xc energy error of',bc%errmax,' Ha'
 write(iout,'(a,f8.4,a,f12.4,a)') '  slack',slack, &
&      ', the model core must fit',bc%budget,' Ha'
 write(iout,'(a,f12.6,a,f10.6)') '  model radius r_s',bc%fill%rs, &
&      '  smallest rc',rcmin
 if(bc%rr(bc%irtop)<rcmin) then
   write(iout,'(a,f12.6)') '  last radius with core charge',bc%rr(bc%irtop)
 end if
 if(bc%full) then
   write(iout,'(a,f10.4,a)') '  (the error peaks beyond the first &
&cutoffs: found with those up to',bc%qext,' bohr**-1)'
 end if
 if(bc%fill%ifill==1) then
   write(iout,'(a)') '  fill minimum tail'
   call legendre_powers(.true.,bc%fill%cleg,bc%fill%cpoly)
   write(iout,'(a,1p,9e16.8)') '  coefficients c0..c8',bc%fill%cpoly(:)
 else
   write(iout,'(a)') '  fill Martins'
   write(iout,'(a,1p,3e16.8)') '  a, b, c',bc%fill%aco,bc%fill%bco,bc%fill%cco
 end if
 write(iout,'(a,f12.4,a)') '  cutoff the model core needs',bc%ecore,' Ha'
 if(icmod==7) then
   write(iout,'(a,f12.6)') '  smallest fitting radius (icmod 6)',bc%rs6
   write(iout,'(a,1p,e12.4,a)') '  hardness mismatch',bc%mismatch,' Ha'
   write(iout,'(a)') '  radii of the basin: r_s, fill, hardness mismatch (Ha)'
   do ii=1,bc%ncand
     if(bc%cfill(ii)==0) then
       write(iout,'(f12.6,2x,a2,a12)') bc%crad(ii),fillname(0),'-'
     else
       write(iout,'(f12.6,2x,a2,1p,e12.4)') bc%crad(ii), &
&            fillname(bc%cfill(ii)), &
&            bc%cmis(ii)
     end if
   end do
   write(iout,'(a)') ''
   deallocate(bc%crad,bc%cfill,bc%cmis,bc%d2ae,bc%d2mod,bc%d2zero)
 end if

 return
 end subroutine budget_report

 real(dp) function smallest_radius(bc,iout,ifail)

!! the smallest r_s in [r_lo, r_top] whose model fits, with the density
!! cutoffs up to qend: bisection over the mesh points between them, then
!! Brent's method on the largest error less errmax in ln r_s (excess_rs)
!!
!! @param bc  the state of the search
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ia,ib,imid,ilo,ihi
 real(dp) :: ga
 ifail=0
 associate(rlo=>bc%rlo,rr=>bc%rr,rtop=>bc%rtop)
   if(excess_rs(bc,log(rtop))>0.0d0) then
     call no_model(bc,iout,ifail)
     if(ifail/=0) return
   end if
   smallest_radius=rlo
   ga=excess_rs(bc,log(rlo))
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
     if(excess_rs(bc,log(rr(imid)))>0.0d0) then
       ia=imid
     else
       ib=imid
     end if
   end do
   smallest_radius=exp(zeroin(excess_rs_of,log(radius_of(bc,ia)), &
&                             log(radius_of(bc,ib)),iout,ifail))
   if(ifail/=0) return
 end associate
 return

 contains

 real(dp) function excess_rs_of(xx)

!! excess_rs of this search, for zeroin
!!
!! @param xx  the argument

!Input variables
 real(dp), intent(in) :: xx

 excess_rs_of=excess_rs(bc,xx)
 return
 end function excess_rs_of

 end function smallest_radius

 real(dp) function radius_of(bc,ii)

!! rr(ii), r_lo below the points of the bisection and r_top above them
!!
!! @param bc  the state of the search
!! @param ii  index

!Input variables
 integer, intent(in) :: ii

!Input/output variables
 type(budget_t), intent(inout) :: bc
 associate(rlo=>bc%rlo,rr=>bc%rr,rtop=>bc%rtop)
   if(rr(ii)<=rlo) then
     radius_of=rlo
   else if(rr(ii)>=rtop) then
     radius_of=rtop
   else
     radius_of=rr(ii)
   end if
 end associate
 return
 end function radius_of

 real(dp) function excess_rs(bc,lnrs)

!! the largest error over the density cutoffs q_b..qend of the model of
!! radius exp(lnrs), less errmax; errmax if the radius has no model
!!
!! @param bc  the state of the search
!! @param lnrs  ln r_s

!Input variables
 real(dp), intent(in) :: lnrs

!Input/output variables
 type(budget_t), intent(inout) :: bc
 associate(errmax=>bc%errmax,ifill=>bc%fill%ifill,qb=>bc%qb,qend=>bc%qend)
   call set_model(bc,exp(lnrs))
   if(ifill==0) then
     excess_rs=errmax
     return
   end if
   excess_rs=err_max(bc,qb,qend,errmax) - errmax
 end associate
 return
 end function excess_rs

 real(dp) function excess(bc,qq)

!! the error at the density cutoff qq less errmax
!!
!! @param bc  the state of the search
!! @param qq  density cutoff (bohr**-1)

!Input variables
 real(dp), intent(in) :: qq

!Input/output variables
 type(budget_t), intent(inout) :: bc
 associate(errmax=>bc%errmax)
   excess=err_at(bc,qq) - errmax
 end associate
 return
 end function excess

 real(dp) function mismatch_at(bc,lnrs)

!! the hardness mismatch of the model of radius exp(lnrs); huge if the
!! radius has no model
!!
!! @param bc  the state of the search
!! @param lnrs  ln r_s

 use der2exc_m, only: der2exc

!Input variables
 real(dp), intent(in) :: lnrs

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp) :: mm
 associate(core=>bc%core,d2ae=>bc%d2ae,d2mod=>bc%d2mod,iexc=>bc%iexc, &
&          c_light=>bc%c_light,ifill=>bc%fill%ifill,lprobe=>bc%lprobe, &
&          mmax=>bc%mmax,nprobe=>bc%nprobe,probes=>bc%probes, &
&          rhotps=>bc%rhotps,rr=>bc%rr,zion=>bc%zion)
   call set_model(bc,exp(lnrs))
   mismatch_at=huge(1.0d0)
   if(ifill==0) return
   call der2exc(rhotps,core,probes,rr,d2mod,d2ae,mm, &
&               zion,iexc,c_light,0,nprobe,lprobe,rr(mmax),mmax)
   mismatch_at=mm
 end associate
 return
 end function mismatch_at

 real(dp) function slope_at(bc,lnrs)

!! the derivative of the hardness mismatch in ln r_s, by central
!! differences of step hmis
!!
!! @param bc  the state of the search
!! @param lnrs  ln r_s

!Input variables
 real(dp), intent(in) :: lnrs

!Input/output variables
 type(budget_t), intent(inout) :: bc
 slope_at=(mismatch_at(bc,lnrs+hmis) - mismatch_at(bc,lnrs-hmis))/(2.0d0*hmis)
 return
 end function slope_at

 real(dp) function minimum_radius(bc,ib,iout,ifail)

!! the radius of the least hardness mismatch in the basin of the radius
!! crad(ib): the root of its slope between the neighbours of crad(ib),
!! or an end of [r_s of icmod 6, r_top] where the slope points out of it
!!
!! @param bc  the state of the search
!! @param ib  index of the radius of the least mismatch of the list
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: ib
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Output variables
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: ua,ub,sa,sb
 integer :: ia,ic
 ifail=0
 associate(crad=>bc%crad,mismatch=>bc%mismatch,ncand=>bc%ncand)
   ia=min(ncand,ib+1)
   ic=max(1,ib-1)
   minimum_radius=crad(ib)
   if(ia==ib) then
     if(slope_at(bc,log(crad(ib)))>=0.0d0) return
   end if
   if(ic==ib) then
     if(slope_at(bc,log(crad(ib)))<=0.0d0) return
   end if
   ua=log(crad(ia))
   ub=log(crad(ic))
   sa=slope_at(bc,ua)
   sb=slope_at(bc,ub)
   if(sa>0.0d0 .or. sb<0.0d0) then
     write(iout,'(/a,f12.6)') 'modcore_budget: WARNING no minimum of the '// &
&         'mismatch bracketed about',crad(ib)
     return
   end if
   minimum_radius=exp(zeroin(slope_at_of,ua,ub,iout,ifail))
   if(ifail/=0) return
 end associate
 return

 contains

 real(dp) function slope_at_of(xx)

!! slope_at of this search, for zeroin
!!
!! @param xx  the argument

!Input variables
 real(dp), intent(in) :: xx

 slope_at_of=slope_at(bc,xx)
 return
 end function slope_at_of

 end function minimum_radius

 subroutine no_model(bc,iout,ifail)
!! reports that no model radius fits and stops
!!
!! @param bc  the state of the search
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: iout

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Output variables
 integer, intent(out) :: ifail
 ifail=0
 associate(budget=>bc%budget,errmax=>bc%errmax,rtop=>bc%rtop)

   write(iout,'(/a,f12.6,a,f12.6,a,es12.4,a)') &
&        'modcore_budget: ERROR no model core inside rc =',rtop, &
&        ' needs less than',budget,' Ha for an error of',errmax,' Ha'
   ifail=1
   return
 end associate
 end subroutine no_model

 subroutine set_model(bc,rad)

!! the model of radius rad: the fill (ifill 1 for the minimum tail in cpoly,
!! 2 for Martins' in aco, bco, cco, 0 for none), in core the model inside
!! rad and the core charge beyond, the nodes of its transform, the mesh
!! ri of the integral below rad, and the xc energy density of the valence
!! charge on the fill there
!!
!! @param bc  the state of the search
!! @param rad  model radius r_s (bohr)

 use logmesh_m, only: lm_derivs, lm_interp
 use xc_dispatch_m, only: exc_dispatch

!Input variables
 real(dp), value :: rad

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp) :: dm(3),l1,l2,xs,xe
 integer :: kk,is
 associate(aco=>bc%fill%aco,al=>bc%al,bco=>bc%fill%bco,cco=>bc%fill%cco, &
&          chr=>bc%chr,core=>bc%core,cr0=>bc%cr0,dq=>bc%dq,epi=>bc%epi, &
&          fext=>bc%fext,fw=>bc%fw,fx=>bc%fx,i0=>bc%i0,iexc=>bc%iexc, &
&          c_light=>bc%c_light,ifill=>bc%fill%ifill,irmx=>bc%irmx, &
&          irow=>bc%irow,js=>bc%js,kin=>bc%kin,mmax=>bc%mmax,msnap=>bc%msnap, &
&          ni=>bc%ni,npdone=>bc%npdone,nt=>bc%nt,nvi=>bc%nvi,nx=>bc%nx, &
&          pg=>bc%pg,pt=>bc%pt,rhoc=>bc%rhoc,rhotps=>bc%rhotps,ri=>bc%ri, &
&          rr=>bc%rr,rs=>bc%fill%rs,shr=>bc%shr,sr0=>bc%sr0,tgi=>bc%tgi, &
&          tk=>bc%tk,vxi=>bc%vxi,wf=>bc%wf,wgi=>bc%wgi,wki=>bc%wki, &
&          wxg=>bc%wxg,xg=>bc%xg,xr=>bc%xr)

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
     write(error_unit,'(/a)') &
&         'modcore_budget: ERROR r_s beyond the table'
     error stop 1
   end if
   xs=log(rs)
   xe=log(rr(is+1))
   do kk=1,ngi
     pt(kk)=exp(xs + (xe-xs)*tgi(kk))
     call lm_interp(rr,rhoc,mmax,pt(kk),l1)
     pg(kk)=(xe-xs)*wgi(kk)*pt(kk)**2*l1
   end do

   call min_tail(bc,dm)
   ifill=0
   if(tail_monotone(bc%fill)) then
     ifill=1
     do kk=1,mmax
       if(rr(kk)>=rs) exit
       core(kk)=poly(bc%fill,(rr(kk)/rs)**2)
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
     fw(kk)=rs*wxg(kk)*fx(kk)*fill_value(bc%fill,fx(kk))
   end do

! the mesh ri, its points of the band-limited models, the valence charge
! there and the xc energy density of the valence charge on the fill,
! extended beyond rs as far as the interpolant of the integral to rs reads
   do kk=1,ni
     ri(kk)=rs*exp((kk-kin)*al)
     call lm_interp(rr,rhotps,mmax,ri(kk),nvi(kk))
     fext(kk)=fill_value(bc%fill,ri(kk))
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
   call exc_dispatch(iexc,c_light,al,wki,vxi,epi,ri,ni,'modcore_budget')
   wf(:)=wki(:)*epi(:)
 end associate
 return
 end subroutine set_model

 real(dp) function fill_value(fill,rad)

!! the fill at rad <= rs
!!
!! @param fill  the fill
!! @param rad  radius (bohr)

!Input variables
 real(dp), intent(in) :: rad

!Input variables
 type(fill_t), intent(in) :: fill
 associate(aco=>fill%aco,bco=>fill%bco,cco=>fill%cco,ifill=>fill%ifill,rs=>fill%rs)
   if(ifill==1) then
     fill_value=poly(fill,(rad/rs)**2)
   else
     fill_value=exp(aco + bco*rad**2 + cco*rad**4)
   end if
 end associate
 return
 end function fill_value

 logical function tail_monotone(fill)

!! whether the minimum tail p, as a function of x = 2(r/rs)**2 - 1, does
!! not rise on [-1, 1]: p' <= 0 at -1, at 1, and at the roots of p'' in
!! between (companion matrix of p'' in powers of x, the roots polished by
!! Newton's method; a nearly real root counts by its real part, a point of
!! [-1, 1] all the same); p' from the Legendre series
!!
!! @param fill  the fill

 use lapack_m, only: dhseqr

!Input variables
 type(fill_t), intent(in) :: fill

!Local variables
 real(dp) :: cx(0:npow),d2(0:npow-2),hess(npow-2,npow-2)
 real(dp) :: wr(npow-2),wi(npow-2),zz(1,1),wkh(4*npow),tr,vv,vd
 integer :: kk,jj,nd,info,it
 associate(cleg=>fill%cleg)
   tail_monotone=.false.
   if(leg_value(fill,-1.0d0,1)>0.0d0) return
   if(leg_value(fill,1.0d0,1)>0.0d0) return
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
     write(error_unit,'(/a,i0)') 'modcore_budget: ERROR dhseqr info =',info
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
       if(leg_value(fill,tr,1)>0.0d0) then
         tail_monotone=.false.
         return
       end if
     end if
   end do
 end associate
 return
 end function tail_monotone

 real(dp) function leg_value(fill,xx,mm)

!! the mm-th derivative in x of the Legendre series cleg at xx
!!
!! @param fill  the fill
!! @param xx  argument, 2(r/rs)**2 - 1
!! @param mm  order of the derivative, 0..4

!Input variables
 real(dp), intent(in) :: xx
 integer, intent(in) :: mm

!Input variables
 type(fill_t), intent(in) :: fill

!Local variables
 real(dp) :: pl(0:npow,0:4)
 associate(cleg=>fill%cleg)
   call legendre_derivs(xx,mm,pl(:,0:mm))
   leg_value=sum(cleg(:)*pl(:,mm))
 end associate
 return
 end function leg_value

 subroutine legendre_derivs(xx,nder,pl)

!! the Legendre polynomials P_k(xx), k=0..npow, and their derivatives up
!! to order nder: P_(k+1) = ((2k+1) x P_k - k P_(k-1))/(k+1) and
!! P_(k+1)^(m) = P_(k-1)^(m) + (2k+1) P_k^(m-1)
!!
!! @param xx  argument
!! @param nder  highest derivative wanted, 0<=nder<=4
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

 real(dp) function poly(fill,tt)

!! the minimum tail at (r/r_s)**2 = tt
!!
!! @param fill  the fill
!! @param tt  (r/r_s)**2

!Input variables
 real(dp), intent(in) :: tt

!Input variables
 type(fill_t), intent(in) :: fill
 poly=leg_value(fill,2.0d0*tt-1.0d0,0)
 return
 end function poly

 subroutine min_tail(bc,dm)

!! the Legendre coefficients cleg of the minimum tail of radius rs: the least
!! squares of W (A c + b) with C c = d (LAPACK dgglse), where A c is the
!! transform of the polynomial at the window nodes (Gauss-Legendre on
!! [0, r_s]), b that of the core charge beyond r_s (tabw and the part of
!! the mesh interval at r_s), W
!! the square root of the weights, and C c = d the match of the core
!! charge and its first two derivatives at r_s
!!
!! @param bc  the state of the search
!! @param dm  the core charge and its first two derivatives at r_s

 use lapack_m, only: dgglse

!Input variables
 real(dp), intent(in) :: dm(3)

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp), allocatable :: amat(:,:),cvec(:),wk(:)
 real(dp) :: bmat(3,0:npow),dvec(3),tt,qx,jx,xx,wx
 integer :: ii,jj,kk,info,lwork
 associate(cleg=>bc%fill%cleg,irow=>bc%irow,nqn=>bc%nqn,pg=>bc%pg,plg=>bc%plg, &
&          pt=>bc%pt,qn=>bc%qn,rs=>bc%fill%rs,tabw=>bc%tabw,wn=>bc%wn, &
&          wxg=>bc%wxg,xg=>bc%xg)
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
     write(error_unit,'(/a,i0)') &
&         'modcore_budget: ERROR dgglse info =',info
     error stop 1
   end if
   deallocate(amat,cvec,wk)
 end associate
 return
 end subroutine min_tail

 real(dp) function outer_transform(bc,kk)

!! int_{r_s}^infinity r**2 rhoc j0(kk r) dr
!!
!! @param bc  the state of the search
!! @param kk  wave number (bohr**-1), > 0

!Input variables
 real(dp), intent(in) :: kk

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp) :: ss
 integer :: jj
 associate(js=>bc%js,nout=>bc%nout,og=>bc%og,ot=>bc%ot,pg=>bc%pg,pt=>bc%pt)
   ss=0.0d0
   do jj=1,ngi
     ss=ss + pg(jj)*sin(kk*pt(jj))
   end do
   do jj=js,nout
     ss=ss + og(jj)*sin(kk*ot(jj))
   end do
   outer_transform=ss/kk
 end associate
 return
 end function outer_transform

 real(dp) function model_transform(bc,kk)

!! int_0^infinity r**2 n j0(kk r) dr of the model n
!!
!! @param bc  the state of the search
!! @param kk  wave number (bohr**-1), > 0

!Input variables
 real(dp), intent(in) :: kk

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp) :: ss
 integer :: jj
 associate(fw=>bc%fw,fx=>bc%fx)
   ss=0.0d0
   do jj=1,nxg
     ss=ss + fw(jj)*sin(kk*fx(jj))
   end do
   model_transform=ss/kk + outer_transform(bc,kk)
 end associate
 return
 end function model_transform

 subroutine spectrum(bc,qe)

!! the model's transform ntk at the nodes of the panels of width dq up to
!! qe (and the next grid point), and the band-limited models rqs(:,m) at
!! the density cutoffs m*qstep: n_q(r) = (2/pi) int_0^q k**2 n(k) j0(kr)
!! dk.  The sines at the nodes of each family (the same node of every
!! panel) rotated from one panel to the next, those of the partial
!! interval at rs computed
!!
!! @param bc  the state of the search
!! @param qe  largest density cutoff needed (bohr**-1)

!Input variables
 real(dp), intent(in) :: qe

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp) :: kq,ss,tt
 integer :: mlast,ip,ig,jj,ll
 associate(chr=>bc%chr,cr=>bc%cr,cr0=>bc%cr0,dq=>bc%dq,fw=>bc%fw,fx=>bc%fx, &
&          irow=>bc%irow,msnap=>bc%msnap,npdone=>bc%npdone,nq=>bc%nq, &
&          nt=>bc%nt,ntk=>bc%ntk,pg=>bc%pg,pt=>bc%pt,rqp=>bc%rqp,rqs=>bc%rqs, &
&          shr=>bc%shr,sr=>bc%sr,sr0=>bc%sr0,tab=>bc%tab,tk=>bc%tk,wk=>bc%wk, &
&          xr=>bc%xr)
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
 end associate
 return
 end subroutine spectrum

 real(dp) function err_at(bc,qq)

!! the error of the density cutoff qq: the band-limited model from the
!! last grid point below qq, the whole panels after it and the part of
!! the last panel up to qq, at its own Gauss-Legendre nodes
!!
!! @param bc  the state of the search
!! @param qq  density cutoff (bohr**-1), 0 <= qq <= qmax

!Input variables
 real(dp), intent(in) :: qq

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp) :: bq(bc%nt),kq,tt,qa,ww
 integer :: m0,ip,ig,pf
 associate(dq=>bc%dq,msnap=>bc%msnap,nt=>bc%nt,ntk=>bc%ntk,rqs=>bc%rqs, &
&          tk=>bc%tk,wk=>bc%wk,xr=>bc%xr)
   call spectrum(bc,qq)
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
       tt=2.0d0/pi*wk(ig)*ww*kq*model_transform(bc,kq)
       bq(:)=bq(:) + tt*sin(kq*xr(:))/xr(:)
     end do
   end if
   err_at=xc_error(bc,bq)
 end associate
 return
 end function err_at

 real(dp) function err_grid(bc,mm)

!! the error of the density cutoff mm*qstep of the grid
!!
!! @param bc  the state of the search
!! @param mm  index of the cutoff

!Input variables
 integer, intent(in) :: mm

!Input/output variables
 type(budget_t), intent(inout) :: bc
 associate(rqs=>bc%rqs)
   call spectrum(bc,mm*qstep)
   err_grid=xc_error(bc,rqs(:,mm))
 end associate
 return
 end function err_grid

 real(dp) function err_max(bc,qa,qe,elim)

!! the largest error over the density cutoffs qa..qe: on the grid qa,
!! m*qstep, qe; if that is above elim, the grid's largest (only whether
!! the error exceeds elim is wanted); otherwise refined by Brent's method
!! between the neighbours of the grid's largest, and of each other local
!! maximum of the grid whose parabola peaks within 1% of the result.
!! qpk the grid point of the largest
!!
!! @param bc  the state of the search
!! @param qa  first density cutoff (bohr**-1)
!! @param qe  last density cutoff (bohr**-1)
!! @param elim  error limit (Ha)

 use ugrid_m, only: ug_fmin

!Input variables
 real(dp), intent(in) :: qa,qe,elim

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp), allocatable :: tq(:),eq(:)
 integer :: nn,mm,ii,ilo,ihi,imx
 associate(nq=>bc%nq,qpk=>bc%qpk)
   call spectrum(bc,qe)
   allocate(tq(nint((qe-qa)/qstep)+3),eq(nint((qe-qa)/qstep)+3))
   nn=1
   tq(1)=qa
   eq(1)=err_at(bc,qa)
   do mm=floor(qa/qstep)+1,nq
     if(mm*qstep>=qe) exit
     if(mm*qstep<=qa) cycle
     nn=nn+1
     tq(nn)=mm*qstep
     eq(nn)=err_grid(bc,mm)
   end do
   if(qe>qa) then
     nn=nn+1
     tq(nn)=qe
     eq(nn)=err_at(bc,qe)
   end if
   imx=maxloc(eq(1:nn),1)
   qpk=tq(imx)
   err_max=eq(imx)
   if(err_max>elim .or. nn==1) then
     deallocate(tq,eq)
     return
   end if
   err_max=max(err_max,-ug_fmin(neg_err_of,tq(max(1,imx-1)),tq(min(nn,imx+1))))
   do ii=1,nn
     if(ii==imx) cycle
     ilo=max(1,ii-1)
     ihi=min(nn,ii+1)
     if(eq(ii)<eq(ilo) .or. eq(ii)<eq(ihi)) cycle
     if(nn>=3) then
       if(parabola_peak(tq,eq,nn,ii)<0.99d0*err_max) cycle
     end if
     err_max=max(err_max,-ug_fmin(neg_err_of,tq(ilo),tq(ihi)))
   end do
   deallocate(tq,eq)
 end associate
 return

 contains

 real(dp) function neg_err_of(xx)

!! neg_err of this search, for ug_fmin
!!
!! @param xx  the argument

!Input variables
 real(dp), intent(in) :: xx

 neg_err_of=neg_err(bc,xx)
 return
 end function neg_err_of

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

 real(dp) function neg_err(bc,qq)

!! minus the error of the density cutoff qq
!!
!! @param bc  the state of the search
!! @param qq  density cutoff (bohr**-1)

!Input variables
 real(dp), intent(in) :: qq

!Input/output variables
 type(budget_t), intent(inout) :: bc
 neg_err=-err_at(bc,qq)
 return
 end function neg_err

 real(dp) function zeroin(fun,xa0,xb0,iout,ifail)

!! a root of fun between xa0 and xb0, where fun changes sign, by Brent's
!! method with no tolerance of its own (ug_brent): it stops when the
!! bracket is a few rounding units wide (as lm_root). Where fun jumps
!! across zero (the fill type changes, a discrete event of
!! doc/smoothness_plan.md) that is the root, and the jump is reported; the
!! iteration limit is an error
!!
!! @param fun  the function
!! @param xa0,xb0  the bracket
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use ugrid_m, only: ug_brent

!Input variables
 procedure(scalar_fn) :: fun
 real(dp), intent(in) :: xa0,xb0
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: fb,fc,fe(2)
 integer :: ierr
 logical :: conv

 ifail=0
 call ug_brent(fun,xa0,xb0,zeroin,ierr,fq=fb,fo=fc,fends=fe,converged=conv)
 if(ierr/=0) then
   write(iout,'(/a,2es12.4)') 'modcore_budget: ERROR no root between', &
&       xa0,xb0
!  a defined result, for the caller's expression before its test of ifail
   zeroin=xa0
   ifail=1
   return
 end if
 if(.not. conv) then
   write(iout,'(/a,es12.4,a,es12.4)') &
&        'modcore_budget: ERROR zeroin not converged in 200 steps, x=', &
&        zeroin,', f=',fb
   ifail=1
   return
 end if
! a continuous function is at its round-off on a bracket a few ulps wide
 if(fb/=0.0d0 .and. min(abs(fb),abs(fc))>1.0d-6*max(abs(fe(1)),abs(fe(2)))) &
& then
   write(iout,'(/2a,es12.4,a,es12.4,a,es12.4)') 'modcore_budget: WARNING', &
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

 subroutine tail_derivatives(fill,rad,dd)

!! the minimum tail and its first 4 derivatives at rad (beyond rs, its
!! continuation, for the pieces of the model): those of the
!! Legendre series in x = 2t - 1, t = (r/rs)**2, by the chain rule with
!! t' = 2r/rs**2, t'' = 2/rs**2 (d/dt = 2 d/dx)
!!
!! @param fill  the fill
!! @param rad  radius (bohr)
!! @param dd  the model and its first 4 derivatives at rad

!Input variables
 real(dp), intent(in) :: rad

!Input variables
 type(fill_t), intent(in) :: fill

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 real(dp) :: pl(0:npow,0:4),pt(0:4),t1,t2
 integer :: mm
 associate(cleg=>fill%cleg,rs=>fill%rs)
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
 end associate
 return
 end subroutine tail_derivatives

 function fill_pack(fill) result(par)

!! the fill as the parameters of fill_form: ifill, r_s, the Legendre
!! coefficients of the minimum tail and Martins' coefficients
!!
!! @param fill  the fill

!Input variables
 type(fill_t), intent(in) :: fill

!Output variables
 real(dp) :: par(npow+6)

 par(1)=dble(fill%ifill)
 par(2)=fill%rs
 par(3:npow+3)=fill%cleg(0:npow)
 par(npow+4:npow+6)=[fill%aco,fill%bco,fill%cco]
 return
 end function fill_pack

 subroutine fill_form(par,rad,dd)

!! the fill of fill_pack's parameters and its first 4 derivatives at rad
!! (fill_derivatives), the closed form of the model's pieces inside r_s
!!
!! @param par  the fill's parameters (fill_pack)
!! @param rad  radius (bohr)
!! @param dd  the fill and its first 4 derivatives at rad

!Input variables
 real(dp), intent(in) :: par(:),rad

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 type(fill_t) :: fill

 fill%ifill=nint(par(1))
 fill%rs=par(2)
 fill%cleg(0:npow)=par(3:npow+3)
 fill%aco=par(npow+4)
 fill%bco=par(npow+5)
 fill%cco=par(npow+6)
 call fill_derivatives(fill,rad,dd)
 return
 end subroutine fill_form

 subroutine fill_derivatives(fill,rad,dd)

!! the fill of the model and its first 4 derivatives at rad
!!
!! @param fill  the fill
!! @param rad  radius (bohr)
!! @param dd  the fill and its first 4 derivatives at rad

!Input variables
 real(dp), intent(in) :: rad

!Input variables
 type(fill_t), intent(in) :: fill

!Output variables
 real(dp), intent(out) :: dd(5)
 associate(ifill=>fill%ifill,rs=>fill%rs)

   if(ifill==1) then
     call tail_derivatives(fill,rad,dd)
   else
     call martins_derivatives(fill,rad,dd)
   end if
 end associate
 return
 end subroutine fill_derivatives

 subroutine martins_derivatives(fill,rad,dd)

!! the model and its first 4 derivatives at rad, from n' = p' n with
!! p = a + b r**2 + c r**4: n(k) = sum_j C(k-1,j) p(j+1) n(k-1-j)
!!
!! @param fill  the fill
!! @param rad  radius (bohr)
!! @param dd  the model and its first 4 derivatives at rad

!Input variables
 real(dp), intent(in) :: rad

!Input variables
 type(fill_t), intent(in) :: fill

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 real(dp) :: dp1(4)
 integer :: kk,jj
 integer, parameter :: binom(0:3,0:3)=reshape([1,1,1,1, 0,1,2,3, &
&                                              0,0,1,3, 0,0,0,1],[4,4])
 associate(aco=>fill%aco,bco=>fill%bco,cco=>fill%cco)
   dp1=[2.0d0*bco*rad + 4.0d0*cco*rad**3, 2.0d0*bco + 12.0d0*cco*rad**2, &
&       24.0d0*cco*rad, 24.0d0*cco]
   dd(1)=exp(aco + bco*rad**2 + cco*rad**4)
   do kk=1,4
     dd(kk+1)=0.0d0
     do jj=0,kk-1
       dd(kk+1)=dd(kk+1) + binom(kk-1,jj)*dp1(jj+1)*dd(kk-jj)
     end do
   end do
 end associate
 return
 end subroutine martins_derivatives

 real(dp) function xc_error(bc,bq)

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
!! @param bc  the state of the search
!! @param bq  band-limited model core charge at xr: rr(1:nx), then ri
!!             (size nt)

 use logmesh_m, only: lm_integral
 use xc_dispatch_m, only: exc_dispatch

!Input variables
 real(dp), intent(in) :: bq(:)

!Input/output variables
 type(budget_t), intent(inout) :: bc

!Local variables
 real(dp) :: si,so1
 associate(ac=>bc%ac,af=>bc%af,al=>bc%al,epi=>bc%epi,eps=>bc%eps, &
&          fext=>bc%fext,iexc=>bc%iexc,c_light=>bc%c_light,ni=>bc%ni, &
&          nt=>bc%nt,nval=>bc%nval,nvi=>bc%nvi,nx=>bc%nx,rhoc=>bc%rhoc, &
&          ri=>bc%ri,rr=>bc%rr,rs=>bc%fill%rs,swin=>bc%swin,vbg=>bc%vbg, &
&          vxc=>bc%vxc,vxi=>bc%vxi,wc=>bc%wc,wf=>bc%wf,wki=>bc%wki, &
&          work=>bc%work)
   work(:)=nval(:) + bq(1:nx)
   call exc_dispatch(iexc,c_light,al,work,vxc,eps,rr,nx,'modcore_budget')
   work(:)=work(:)*eps(:)
   ac(:)=swin(:)*rr(1:nx)**2*(work(:) - wc(:) - vbg*(bq(1:nx) - rhoc(1:nx)))
   wki(:)=nvi(:) + bq(nx+1:nt)
   call exc_dispatch(iexc,c_light,al,wki,vxi,epi,ri,ni,'modcore_budget')
   wki(:)=wki(:)*epi(:)
   af(:)=ri(:)**2*(wki(:) - wf(:) - vbg*(bq(nx+1:nt) - fext(:)))
   call lm_integral(ri,af,ni,2,0.0d0,rs,si)
   call lm_integral(rr,ac,nx,2,rs,radius,so1)
   xc_error=abs(si + so1)
 end associate
 return
 end function xc_error

 subroutine modcore_probes(lmax,lloc,nproj,epa,wellprj,vp,vlocpw,vkbpw, &
&                         evkb,rr,mmax,nkap,probes,lprobe,nprobe, &
&                         iout,ifail)

!! the probe densities of the hardness of icmod=7: for every projector
!! iprj of every l, the eigenstate with iprj-1 nodes of the screened
!! pseudo-Hamiltonian, as for run_diag, at unit occupation (the average
!! over J = l +/- 1/2 for nkap=2), if it is bound.  The occupied valence
!! states always are; an unoccupied one counts if bound, so that an atom
!! with one valence state still has a matrix to match
!!
!! @param lmax  maximum angular momentum
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param nproj  number of projectors for each l
!! @param epa  reference energies of the projectors
!! @param wellprj  whether each projector is a well state (its energy then
!!                 of either sign), else a bound state
!! @param vp  screened semi-local and local pseudopotentials
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param probes  probe densities, (u/r)**2 normalized to one electron
!! @param lprobe  angular momenta of the probes
!! @param nprobe  number of probes
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxl, mxlin, mxprj
 use kappa_m, only: kap_count
 use pieces_m, only: piecewise_t
 use lschvkbb_m, only: lschvkbb
 use pieces_m, only: piecewise_t

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,nkap
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap)
 real(dp), intent(in) :: epa(mxprj,mxlin+1,nkap),evkb(mxprj,mxl+1,nkap)
 logical, intent(in) :: wellprj(mxprj,mxlin+1)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: probes(mmax,mxprj*(lmax+1))
 integer, intent(out) :: lprobe(mxprj*(lmax+1))
 integer, intent(out) :: nprobe
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: uu(:),up(:),dens(:)
 real(dp) :: et,emin,emax,fj
 integer :: l1,ll,iprj,ikap,mkap,npr,ierr,mch
 logical :: bound

 ifail=0
 allocate(uu(mmax),up(mmax),dens(mmax))
 nprobe=0
 do l1=1,lmax+1
   ll=l1-1
   mkap=kap_count(nkap,l1-1)
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
&                    evkb(1,l1,ikap),uu,up,mmax,mch,iout,ifail)
       if(ifail/=0) return
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
&                                allwaves,errmax,ecut,iout,ifail)

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
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use logmesh_m, only: lm_fh_value
 use ugrid_m, only: ug_cross_down, ug_brent
 use ugrid_m, only: ug_brent

!Input variables
 integer, intent(in) :: ncurve,nprof
 real(dp), intent(in) :: ecurve(2,ncurve,nprof),wcurve(nprof),errmax
 logical, intent(in) :: weighted,allwaves
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: ecut
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: qerr,ecut1,qa,qb,hs,xa,xb,fa,fb
 real(dp) :: q1s(nprof),dqs(nprof)
 real(dp), allocatable :: lf(:,:)
 integer :: ip,kk,nk,ns,ierr
 integer :: nks(nprof)
 logical :: found,found1

 ifail=0
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
     write(iout,'(/a,es12.4,a)') &
&          'modcore_budget_table: ERROR no pseudo-wave reaches', &
&          errmax,' Ha per electron in its table'
     ifail=1
     return
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
     write(iout,'(/a,i0,a)') &
&          'modcore_budget_table: ERROR pseudo-wave ', &
&          ip,' has a weight and no convergence curve'
     ifail=1
     return
   end if
   qa=max(qa,q1s(ip))
   qb=max(qb,q1s(ip)+dble(nks(ip)-1)*dqs(ip))
   hs=min(hs,dqs(ip))
 end do
 if(.not. (hs<huge(1.0d0)) .or. qa>=qb) then
   write(iout,'(/a)') 'modcore_budget_table: ERROR no weighted curves'
   ifail=1
   return
 end if
 fb=excess_mean(qb)
 if(fb>0.0d0) then
   write(iout,'(/a,es12.4,a)') &
&        'modcore_budget_table: ERROR the weighted mean does not'// &
&        ' reach',errmax,' Ha per electron'
   ifail=1
   return
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
       write(iout,'(/a)') &
&            'modcore_budget_table: ERROR the bracket lost its sign change'
       ifail=1
       return
     end if
     ecut=0.5d0*qerr**2
     return
   end if
   xb=xa
 end do
 write(iout,'(/a,es12.4,a)') &
&      'modcore_budget_table: ERROR the weighted mean is below', &
&      errmax,' Ha per electron where the curves start'
 ifail=1
 return

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

 end module modcore_budget_m
