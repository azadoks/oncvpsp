!
! Copyright (c) 2026 by Austin Zadoks
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
!> The "phase shifts" of run_phsft and oncvpsp_ld.x, atan(r phi'/phi) of the
!> scattering solution phi = u/r at a set of radii, as continuous functions
!> of the energy, of the radius and of the potential
!> (doc/smoothness_plan.md, X-2). The radii are any radii: (u, al r u') are
!> interpolated to each from the mesh (lm_interp), all of them from the one
!> solution of each energy, which runs to the last. The angle is that of the
!> point (phi, r phi'), defined up to pi only, as the overall sign of a
!> solution is arbitrary (the outward solution of vkboutwf flips it as the
!> energy crosses a pole of its coefficients): its doubled angle is
!> sign-invariant. The phase is followed from one energy to the next by
!> continuation, taking the change of the doubled angle in (-pi, pi] and
!> halving the energy step until no step turns it by more than pi/2, at each
!> radius on its own (so a radius's phases are the same with any other
!> radii), so the result is the same continuous function of the energy
!> whatever steps were taken: no thresholds, no 2 pi or pi corrections, no
!> search for spurious jumps. A resonance narrower than the step turns the
!> doubled angle by 2 pi within it unseen (73_Ta's 4f at -0.92 Ha at 2.6
!> bohr, 0.05 Ha steps); the Pruefer angle, atan(r phi'/phi) less pi per
!> node of u inside the radius, sees it, as a node crosses the radius. It
!> is the same function up to a constant multiple of pi for a local
!> potential, so a step whose change differs from the Pruefer angle's by a
!> multiple of pi takes the Pruefer angle's (the all-electron potential).
!> With projectors a node can appear inside them, where the Pruefer angle
!> jumps by pi and the phase does not (Pb p with debl 4.7), so such a step
!> is halved instead, to maxdepth, which resolves a resonance and leaves
!> the jump to the continuation.

 module phsft_m
 use constants_m, only: dp, pi
 use outwf_m, only: outwf_t, outwf_setup, outwf_solve, outwf_free
 use adams_m, only: adams_t
 implicit none
 private
 public :: phsft_t, phsft_mch, phsft_ae_init, phsft_ps_init, phsft_free
 public :: phsft_pair, phsft_node_phase, phsft_node_phases

!> the deepest halving of an energy step: a phase that still turns by more
!> than pi/2 within depsh/2**maxdepth is an error
 integer, parameter :: maxdepth=40

 type phsft_t
!> 1 for the all-electron potential (scatterstate), 2 for the
!> pseudopotential (outwf_m)
   integer :: kind=0
!> 0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!> (kind 1)
   integer :: irel=0
!> angular momentum
   integer :: ll=0
!> Dirac kappa, l or -(l+1) for j = l -/+ 1/2 (kind 1)
   integer :: kap=0
!> size of the log radial mesh
   integer :: mmax=0
!> order and coefficients of the Adams-Moulton steps (kind 1; the
!> outward solver holds them for kind 2)
   type(adams_t) :: adams
!> the last point of the solutions, the largest of mchq (and of the
!> projectors' support, kind 2)
   integer :: mch=0
!> number of projectors (kind 2)
   integer :: nvkb=0
!> number of radii
   integer :: nq=0
!> atomic number (kind 1)
   real(dp) :: zz=0.0d0
!> log of the mesh ratio
   real(dp) :: al=0.0d0
!> the speed of light (irel>0)
   real(dp) :: c_light=0.0d0
!> the radii of the phases (bohr)
   real(dp), allocatable :: rq(:)
!> the last point of the solution each radius needs alone (mch is the
!> largest)
   integer, allocatable :: mchq(:)
!> log radial mesh (bohr)
   real(dp), allocatable :: rr(:)
!> the all-electron potential (Ha) (kind 1)
   real(dp), allocatable :: vv(:)
!> coefficients of the projectors (Ha) (kind 2)
   real(dp), allocatable :: evkb(:)
!> the outward solver of the pseudopotential (kind 2)
   type(outwf_t) :: ow
 end type phsft_t

 contains

 integer function phsft_mch(rr,mmax,rq)

!! the last mesh point a solution is needed to for its interpolant at rq:
!! the reach of the interpolant beyond the first point past rq; 0 if that
!! is beyond the mesh
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param rq  radius (bohr)

 use logmesh_m, only: lm_fh_reach

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rq

!Local variables
 integer :: ii

 phsft_mch=0
 do ii=1,mmax
   if(rr(ii)>rq) then
     if(ii+lm_fh_reach()<=mmax) phsft_mch=ii+lm_fh_reach()
     return
   end if
 end do
 return
 end function phsft_mch

 subroutine phsft_ae_init(irel,c_light,ll,kap,rr,vv,zz,mmax,adams,nq,rq,ps, &
&                         iout,ifail)

!! the all-electron solver, at the radii rq
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2
!!              fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2 (irel=2 only)
!! @param rr  log radial mesh (bohr)
!! @param vv  all-electron potential
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param nq  number of radii
!! @param rq  radii of the phases, each with phsft_mch(rr,mmax,rq)>0
!! @param ps  the solver
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax,nq
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: rr(mmax),vv(mmax),zz,rq(nq)
 integer, intent(in) :: iout

!Output variables
 type(phsft_t), intent(out) :: ps
 integer, intent(out) :: ifail

 ifail=0
 ps%kind=1
 ps%irel=irel
 ps%c_light=c_light
 ps%ll=ll
 ps%kap=kap
 ps%zz=zz
 ps%adams=adams
 call common_init(rr,mmax,nq,rq,ps,iout,ifail)
 if(ifail/=0) return
 ps%vv=vv
 return
 end subroutine phsft_ae_init

 subroutine phsft_ps_init(ll,nvkb,rr,vloc,vlocpw,vkbpw,evkb,mmax,adams,nq,rq, &
&      ps, &
&      iout,ifail)

!! the pseudopotential solver, local potential and nvkb projectors (none
!! for the semi-local potential that is the local one), at the radii rq
!!
!! @param ll  angular momentum
!! @param nvkb  number of projectors
!! @param rr  log radial mesh (bohr)
!! @param vloc  local potential on the mesh
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param nq  number of radii
!! @param rq  radii of the phases, each with phsft_mch(rr,mmax,rq)>0
!! @param ps  the solver
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use pieces_m, only: piecewise_t

!Input variables
 integer, intent(in) :: ll,nvkb,mmax,nq
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb),rq(nq)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)
 integer, intent(in) :: iout

!Output variables
 type(phsft_t), intent(out), target :: ps

 integer, intent(out) :: ifail
 ifail=0
 ps%kind=2
 ps%ll=ll
 ps%nvkb=nvkb
 call common_init(rr,mmax,nq,rq,ps,iout,ifail)
 if(ifail/=0) return
 ps%evkb=evkb
 call outwf_setup(ll,rr,mmax,adams,vloc,vlocpw,nvkb,vkbpw,ps%ow,iout,ifail)
 if(ifail/=0) return
!the projections need the solution over the projectors' support, also for
!radii inside it
 if(nvkb>0) then
   ps%mch=max(ps%mch,phsft_mch(rr,mmax,ps%ow%rsrc))
   ps%mchq(:)=max(ps%mchq(:),phsft_mch(rr,mmax,ps%ow%rsrc))
 end if
 return
 end subroutine phsft_ps_init

 subroutine common_init(rr,mmax,nq,rq,ps,iout,ifail)

!! the mesh, the radii and the last point of the solutions, that of the
!! outermost radius
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param nq  number of radii
!! @param rq  radii of the phases (bohr)
!! @param ps  the solver; its mesh, radii, al, mch and mchq set
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input
 use logmesh_m, only: lm_step

!Input variables
 integer, intent(in) :: mmax,nq
 real(dp), intent(in) :: rr(mmax),rq(nq)
 integer, intent(in) :: iout

!Input/output variables
 type(phsft_t), intent(inout) :: ps

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: iq,mq

 ifail=0
 ps%mmax=mmax
 ps%nq=nq
 ps%rq=rq
 ps%rr=rr
 ps%al=lm_step(rr,mmax)
 ps%mch=0
 allocate(ps%mchq(nq))
 do iq=1,nq
   mq=phsft_mch(rr,mmax,rq(iq))
   if(mq==0) then
     write(iout,'(/a/a,es12.4,a,es12.4)') &
&          'common_init: ERROR the radius is too close to the end '// &
&          'of the mesh', &
&          ' rq=',rq(iq),' rr(mmax)=',rr(mmax)
     ifail=fail_input
     return
   end if
   ps%mch=max(ps%mch,mq)
   ps%mchq(iq)=mq
 end do
 return
 end subroutine common_init

 subroutine phsft_free(ps)

!! releases the solver
!!
!! @param ps  the solver

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

 if(ps%kind==2) call outwf_free(ps%ow)
 if(allocated(ps%rr)) deallocate(ps%rr)
 if(allocated(ps%vv)) deallocate(ps%vv)
 if(allocated(ps%evkb)) deallocate(ps%evkb)
 if(allocated(ps%rq)) deallocate(ps%rq)
 if(allocated(ps%mchq)) deallocate(ps%mchq)
 ps%kind=0
 ps%nq=0
 return
 end subroutine phsft_free

 subroutine phsft_scan(ps,epsh2,depsh,npsh,eal,psh,pal,iout,ifail)

!! the phases at the energies epsh2, epsh2-depsh, ... (npsh of them), each
!! followed from the one before, the first in (-pi/2, pi/2]; and the
!! phases at the energy eal, followed from the nearest of them
!!
!! @param ps  the solver
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha),
!!               the first energy of the scan
!! @param depsh  energy step of the log derivative analysis (Ha), downwards
!! @param npsh  number of energies
!! @param eal  the energy of pal
!! @param psh  the phases, per energy and radius
!! @param pal  the phases at eal, per radius
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 integer, intent(in) :: npsh
 real(dp), intent(in) :: epsh2,depsh,eal
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: psh(npsh,ps%nq),pal(ps%nq)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: da(ps%nq),pr(ps%nq,npsh),pral(ps%nq)
 integer :: ii,kk,iq

 ifail=0
 call doubled_angle(ps,epsh2,da,pr(:,1),iout,ifail)
 psh(1,:)=0.5d0*da(:)
 if(ifail/=0) return
 do ii=2,npsh
   call doubled_angle(ps,energy(ii),da,pr(:,ii),iout,ifail)
   if(ifail/=0) return
   do iq=1,ps%nq
     call follow(ps,iq,psh(ii-1,iq),pr(iq,ii-1),energy(ii-1),energy(ii), &
&                da(iq),pr(iq,ii),0,psh(ii,iq),iout,ifail)
     if(ifail/=0) return
   end do
 end do
 kk=min(npsh,max(1,nint((epsh2-eal)/depsh)+1))
 call doubled_angle(ps,eal,da,pral,iout,ifail)
 if(ifail/=0) return
 do iq=1,ps%nq
   call follow(ps,iq,psh(kk,iq),pr(iq,kk),energy(kk),eal,da(iq),pral(iq), &
&              0,pal(iq),iout,ifail)
   if(ifail/=0) return
 end do
 return

 contains

 real(dp) function energy(ii)
 integer, intent(in) :: ii
 energy=epsh2-depsh*real(ii-1,dp)
 end function energy

 end subroutine phsft_scan

 subroutine phsft_pair(irel,c_light,ll,kap,nvkb,rr,vfull,zz,vloc,vlocpw,vkbpw, &
&                      evkb,mmax,adams,nq,rq,npsh,epsh2,depsh,eal,pshae,pshps, &
&                      iout,ifail)

!! the all-electron and the pseudopotential phases of one l and kappa at
!! the radii rq and the energies of phsft_scan, the pseudopotential's
!! shifted by the multiple of pi that makes them equal the all-electron
!! ones at the reference energy eal, where the two solutions agree beyond
!! the core radius
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2 (irel=2 only)
!! @param nvkb  number of projectors
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param zz  atomic number
!! @param vloc  local potential on the mesh
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param nq  number of radii
!! @param rq  radii of the phases
!! @param npsh  number of energies
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha),
!!               the first energy of the scan
!! @param depsh  energy step of the log derivative analysis (Ha), downwards
!! @param eal  the reference energy of the alignment
!! @param pshae  the all-electron phases, per energy and radius
!! @param pshps  the pseudopotential phases, likewise
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use pieces_m, only: piecewise_t

!Input variables
 integer, intent(in) :: irel,ll,kap,nvkb,mmax,nq,npsh
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: rr(mmax),vfull(mmax),zz,vloc(mmax),evkb(nvkb)
 real(dp), intent(in) :: rq(nq),epsh2,depsh,eal
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: pshae(npsh,nq),pshps(npsh,nq)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: palae(nq),palps(nq)
 integer :: iq
 type(phsft_t), target :: ps

 ifail=0
 call phsft_ae_init(irel,c_light,ll,kap,rr,vfull,zz,mmax,adams,nq,rq,ps,iout, &
&                   ifail)
 if(ifail/=0) return
 call phsft_scan(ps,epsh2,depsh,npsh,eal,pshae,palae,iout,ifail)
 if(ifail/=0) return
 call phsft_free(ps)
 call phsft_ps_init(ll,nvkb,rr,vloc,vlocpw,vkbpw,evkb,mmax,adams,nq,rq,ps, &
&                   iout,ifail)
 if(ifail/=0) return
 call phsft_scan(ps,epsh2,depsh,npsh,eal,pshps,palps,iout,ifail)
 if(ifail/=0) return
 call phsft_free(ps)
 do iq=1,nq
   pshps(:,iq)=pshps(:,iq)+pi*anint((palae(iq)-palps(iq))/pi)
 end do
 return
 end subroutine phsft_pair

 subroutine phsft_node_phase(ps,iq,ee,th,nodes,uq_out,iout,ifail)

!! the phase of radius iq at the energy ee with no continuation:
!! atan(r phi'/phi) - pi (the nodes of u inside the radius), the Pruefer
!! angle, continuous across a node's passing the radius. For a local
!! potential it falls as the energy rises, a monotonic function that a
!! root search brackets with no step to alias a turn by pi
!! (oncvpsp_confine.x). A pseudo-wave function can gain or lose a node
!! inside the projectors (Pb p with debl 4.7, at the same energy at every
!! radius), where the phase jumps by pi: with nodes false, only the angle
!! atan(r phi'/phi) in (-pi/2, pi/2], for the caller to count the poles of
!! r phi'/phi below ee instead (ws_eigenvalues), which r phi'/phi falls
!! between for any Hermitian potential
!!
!! @param ps  the solver
!! @param iq  the radius
!! @param ee  energy (Ha)
!! @param th  the phase
!! @param nodes  whether to subtract pi per node (default true)
!! @param uq_out  u at the radius (output)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 integer, intent(in) :: iq
 real(dp), intent(in) :: ee
 logical, intent(in), optional :: nodes
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: th
 real(dp), intent(out), optional :: uq_out
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp) :: uq
 logical :: count_nodes

 ifail=0
 th=0.0d0
 count_nodes=.true.
 if(present(nodes)) count_nodes=nodes
 allocate(uu(ps%mmax,2),up(ps%mmax,2))
 call solve(ps,ee,ps%mchq(iq),uu,up,iout,ifail)
 if(ifail/=0) return
 call angle_at(ps,iq,ps%mchq(iq),uu(:,1),up(:,1),count_nodes,th,uq)
 if(present(uq_out)) uq_out=uq
 deallocate(uu,up)
 return
 end subroutine phsft_node_phase

 subroutine phsft_node_phases(ps,ee,th,nodes,iout,ifail)

!! the phases of phsft_node_phase at every radius from one solution, out
!! to the outermost (the scan of oncvpsp_confine.x)
!!
!! @param ps  the solver
!! @param ee  energy (Ha)
!! @param th  the phase of each radius (output)
!! @param nodes  whether to subtract pi per node (default true)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 real(dp), intent(in) :: ee
 logical, intent(in), optional :: nodes
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: th(ps%nq)
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp) :: uq
 integer :: iq
 logical :: count_nodes

 ifail=0
 count_nodes=.true.
 if(present(nodes)) count_nodes=nodes
 allocate(uu(ps%mmax,2),up(ps%mmax,2))
 call solve(ps,ee,ps%mch,uu,up,iout,ifail)
 if(ifail/=0) return
 do iq=1,ps%nq
   call angle_at(ps,iq,ps%mch,uu(:,1),up(:,1),count_nodes,th(iq),uq)
 end do
 deallocate(uu,up)
 return
 end subroutine phsft_node_phases

 subroutine solve(ps,ee,mlast,uu,up,iout,ifail)
!! the solution at ee out to the point mlast
!!
!! @param ps  the solver
!! @param ee  the scattering energy (Ha)
!! @param mlast  the last point integrated to
!! @param uu  the radial wave function (times r), zero beyond mlast; the
!!            Dirac small component in uu(:,2) (kind 1)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use scatterstate_m, only: scatterstate

!Input variables
 real(dp), intent(in) :: ee
 integer, intent(in) :: mlast
 integer, intent(in) :: iout

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Output variables
 real(dp), intent(out) :: uu(ps%mmax,2),up(ps%mmax,2)
 integer, intent(out) :: ifail

!Local variables
 integer :: nn
 ifail=0
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 if(ps%kind==1) then
   nn=0
   call scatterstate(ps%irel,ps%c_light,nn,ps%ll,ps%kap,ee,ps%rr,ps%zz, &
&                    ps%vv,uu,up,ps%mmax,ps%adams,mlast)
 else
   call outwf_solve(ps%ow,ee,ps%evkb,mlast,mlast,uu(:,1),up(:,1), &
&        iout=iout,ifail=ifail)
   if(ifail/=0) return
 end if
 end subroutine solve

 subroutine angle_at(ps,iq,mlast,uu,up,count_nodes,th,uq)
!! the phase of radius iq of a solution out to the point mlast, and u there
!!
!! @param ps  the solver
!! @param iq  the index of the radius in ps%rq
!! @param mlast  the last point of the solution
!! @param uu  the radial wave function (times r)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param count_nodes  .true. to take pi off the phase for each node of
!!                     uu inside the radius
!! @param th  the phase
!! @param uq  uu interpolated to the radius
 use logmesh_m, only: lm_interp

!Input variables
 type(phsft_t), intent(in) :: ps
 integer, intent(in) :: iq,mlast
 real(dp), intent(in) :: uu(ps%mmax),up(ps%mmax)
 logical, intent(in) :: count_nodes

!Output variables
 real(dp), intent(out) :: th,uq

!Local variables
 real(dp) :: upq
 integer :: jj,nnode,jlast
 call lm_interp(ps%rr(1:mlast),uu(1:mlast),mlast,ps%rq(iq),uq)
 call lm_interp(ps%rr(1:mlast),up(1:mlast),mlast,ps%rq(iq),upq)
 nnode=0
 jlast=1
 do jj=2,mlast
   if(ps%rr(jj)>=ps%rq(iq)) exit
   if(uu(jj)*uu(jj-1)<0.0d0) nnode=nnode+1
   jlast=jj
 end do
 if(uq*uu(jlast)<0.0d0) nnode=nnode+1
! r phi'/phi = r u'/u - 1, and al r u' is the solvers' derivative
 th=atan2(upq/ps%al-uq,uq)
 if(th>0.5d0*pi) th=th-pi
 if(th<=-0.5d0*pi) th=th+pi
 if(.not. count_nodes) nnode=0
 th=th-pi*dble(nnode)
 end subroutine angle_at

 recursive subroutine follow(ps,iq,th0,pr0,e0,e1,da1,pr1,depth,th1,iout, &
&                             ifail)

!! the phase of radius iq at e1, continued from th0 at e0
!!
!! @param ps  the solver
!! @param iq  the radius
!! @param th0  the phase at e0
!! @param pr0  the Pruefer angle at e0
!! @param e0  the energy of th0
!! @param e1  the energy wanted
!! @param da1  the doubled angle at e1
!! @param pr1  the Pruefer angle at e1
!! @param depth  halvings so far
!! @param th1  the phase at e1
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_numerical

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 real(dp), intent(in) :: th0,pr0,e0,e1,da1,pr1
 integer, intent(in) :: iq,depth
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: th1
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: dd,em,thm,dam(ps%nq),prm(ps%nq)
 integer :: nturn

 ifail=0
 dd=da1-2.0d0*th0
 dd=dd-2.0d0*pi*anint(dd/(2.0d0*pi))
 if(abs(dd)>0.5d0*pi .and. depth>=maxdepth) then
   write(iout,'(/a,es12.4,a,es12.4)') &
&        'follow: ERROR the phase turns by more than pi/2 between',e0, &
&        ' and',e1
   th1=th0
   ifail=fail_numerical
   return
 end if
! the turns by pi that the step leaves out, by the Pruefer angle
 nturn=0
 if(abs(dd)<=0.5d0*pi) nturn=nint((pr1-pr0-0.5d0*dd)/pi)
 if(nturn/=0 .and. ps%kind==1) then
   th1=th0+(pr1-pr0)
 else if(abs(dd)>0.5d0*pi .or. (nturn/=0 .and. depth<maxdepth)) then
   em=0.5d0*(e0+e1)
   call doubled_angle(ps,em,dam,prm,iout,ifail)
   if(ifail/=0) return
   call follow(ps,iq,th0,pr0,e0,em,dam(iq),prm(iq),depth+1,thm,iout,ifail)
   if(ifail/=0) return
   call follow(ps,iq,thm,prm(iq),em,e1,da1,pr1,depth+1,th1,iout,ifail)
   if(ifail/=0) return
 else
   th1=th0+0.5d0*dd
 end if
 return
 end subroutine follow

 subroutine doubled_angle(ps,ee,da,pr,iout,ifail)

!! twice the angle of (phi, r phi') at each radius, in (-pi, pi]: that of
!! (u, u' - u/r) at the energy ee, the same for -u; and the Pruefer angle
!!
!! @param ps  the solver
!! @param ee  energy (Ha)
!! @param da  the doubled angle of each radius
!! @param pr  the Pruefer angle of each radius (phsft_node_phase)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use scatterstate_m, only: scatterstate
 use logmesh_m, only: lm_interp

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 real(dp), intent(in) :: ee
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: da(ps%nq),pr(ps%nq)
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp) :: uq,upq,xx,yy,hh,uq2
 integer :: nn,iq

 ifail=0
 allocate(uu(ps%mmax,2),up(ps%mmax,2))
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 if(ps%kind==1) then
   nn=0
   call scatterstate(ps%irel,ps%c_light,nn,ps%ll,ps%kap,ee,ps%rr,ps%zz, &
&                    ps%vv,uu,up,ps%mmax,ps%adams,ps%mch)
 else
   call outwf_solve(ps%ow,ee,ps%evkb,ps%mch,ps%mch,uu(:,1),up(:,1), &
&        iout=iout,ifail=ifail)
   if(ifail/=0) return
 end if
 do iq=1,ps%nq
   call lm_interp(ps%rr(1:ps%mch),uu(1:ps%mch,1),ps%mch,ps%rq(iq),uq)
   call lm_interp(ps%rr(1:ps%mch),up(1:ps%mch,1),ps%mch,ps%rq(iq),upq)
! r phi'/phi = r u'/u - 1, and al r u' is the solvers' derivative
   hh=hypot(uq,upq/ps%al-uq)
   xx=uq/hh
   yy=(upq/ps%al-uq)/hh
   da(iq)=atan2(2.0d0*xx*yy,(xx-yy)*(xx+yy))
   call angle_at(ps,iq,ps%mch,uu(:,1),up(:,1),.true.,pr(iq),uq2)
 end do
 deallocate(uu,up)
 return
 end subroutine doubled_angle

 end module phsft_m
