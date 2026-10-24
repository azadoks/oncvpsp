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
!> search for spurious jumps. Except across a resonance narrower than the
!> step, where the doubled angle turns by 2 pi within one step unseen, and
!> the phase below it is off by pi unless a halving lands in it (73_Ta's 4f
!> near -0.95 Ha at 2.6 bohr, both phases, the pseudopotential's through
!> the alignment).

 module phsft_m
 use constants_m, only: dp, pi
 use pieces_m, only: piecewise_t
 use outwf_m, only: outwf_t, outwf_setup, outwf_solve, outwf_free
 use logmesh_m, only: lm_interp, lm_fh_reach
 implicit none
 private
 public :: phsft_t, phsft_mch, phsft_ae_init, phsft_ps_init, phsft_free
 public :: phsft_scan, phsft_pair, phsft_node_phase

!> the deepest halving of an energy step: a phase that turns by more than
!> pi/2 within depsh/2**maxdepth is taken as it is found
 integer, parameter :: maxdepth=40

 type phsft_t
!> 1 for the all-electron potential (scatterstate), 2 for the
!> pseudopotential (outwf_m)
   integer :: kind=0
   integer :: irel=0,ll=0,kap=0,mmax=0,mch=0,nvkb=0,nq=0
   real(dp) :: zz=0.0d0,al=0.0d0
!> the radii of the phases
   real(dp), allocatable :: rq(:)
   real(dp), allocatable :: rr(:),vv(:),evkb(:)
   type(outwf_t) :: ow
 end type phsft_t

 contains

 integer function phsft_mch(rr,mmax,rq)

!! the last mesh point a solution is needed to for its interpolant at rq:
!! the reach of the interpolant beyond the first point past rq; 0 if that
!! is beyond the mesh
!!
!! @param rr  log radial mesh
!! @param mmax  size of the log radial mesh
!! @param rq  radius (bohr)

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

 subroutine phsft_ae_init(irel,ll,kap,rr,vv,zz,mmax,nq,rq,ps)

!! the all-electron solver, at the radii rq
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2
!!              fully-relativistic
!! @param ll  angular momentum
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param rr  log radial mesh
!! @param vv  all-electron potential
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param nq  number of radii
!! @param rq  radii of the phases, each with phsft_mch(rr,mmax,rq)>0
!! @param ps  the solver

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax,nq
 real(dp), intent(in) :: rr(mmax),vv(mmax),zz,rq(nq)

!Output variables
 type(phsft_t), intent(out) :: ps

 ps%kind=1
 ps%irel=irel
 ps%ll=ll
 ps%kap=kap
 ps%zz=zz
 call common_init(rr,mmax,nq,rq,ps)
 ps%vv=vv
 return
 end subroutine phsft_ae_init

 subroutine phsft_ps_init(ll,nvkb,rr,vloc,vlocpw,vkbpw,evkb,mmax,nq,rq,ps)

!! the pseudopotential solver, local potential and nvkb projectors (none
!! for the semi-local potential that is the local one), at the radii rq
!!
!! @param ll  angular momentum
!! @param nvkb  number of projectors
!! @param rr  log radial mesh
!! @param vloc  local potential on the mesh
!! @param vlocpw  its pieces
!! @param vkbpw  the projectors' pieces
!! @param evkb  their coefficients
!! @param mmax  size of the log radial mesh
!! @param nq  number of radii
!! @param rq  radii of the phases, each with phsft_mch(rr,mmax,rq)>0
!! @param ps  the solver

!Input variables
 integer, intent(in) :: ll,nvkb,mmax,nq
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb),rq(nq)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)

!Output variables
 type(phsft_t), intent(out), target :: ps

 ps%kind=2
 ps%ll=ll
 ps%nvkb=nvkb
 call common_init(rr,mmax,nq,rq,ps)
 ps%evkb=evkb
 call outwf_setup(ll,rr,mmax,vloc,vlocpw,nvkb,vkbpw,ps%ow)
!the projections need the solution over the projectors' support, also for
!radii inside it
 if(nvkb>0) ps%mch=max(ps%mch,phsft_mch(rr,mmax,ps%ow%rsrc))
 return
 end subroutine phsft_ps_init

 subroutine common_init(rr,mmax,nq,rq,ps)

!! the mesh, the radii and the last point of the solutions, that of the
!! outermost radius

!Input variables
 integer, intent(in) :: mmax,nq
 real(dp), intent(in) :: rr(mmax),rq(nq)

!Input/output variables
 type(phsft_t), intent(inout) :: ps

!Local variables
 integer :: iq,mq

 ps%mmax=mmax
 ps%nq=nq
 ps%rq=rq
 ps%rr=rr
 ps%al=0.01d0*log(rr(101)/rr(1))
 ps%mch=0
 do iq=1,nq
   mq=phsft_mch(rr,mmax,rq(iq))
   if(mq==0) then
     write(6,'(/a/a,1p,e12.4,a,e12.4)') &
&          'phsft_m: ERROR the radius is too close to the end of the mesh', &
&          ' rq=',rq(iq),' rr(mmax)=',rr(mmax)
     error stop 1
   end if
   ps%mch=max(ps%mch,mq)
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
 ps%kind=0
 ps%nq=0
 return
 end subroutine phsft_free

 subroutine phsft_scan(ps,epsh2,depsh,npsh,eal,psh,pal)

!! the phases at the energies epsh2, epsh2-depsh, ... (npsh of them), each
!! followed from the one before, the first in (-pi/2, pi/2]; and the
!! phases at the energy eal, followed from the nearest of them
!!
!! @param ps  the solver
!! @param epsh2  the first (highest) energy (Ha)
!! @param depsh  the energy step (downwards)
!! @param npsh  number of energies
!! @param eal  the energy of pal
!! @param psh  the phases, per energy and radius
!! @param pal  the phases at eal, per radius

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 integer, intent(in) :: npsh
 real(dp), intent(in) :: epsh2,depsh,eal

!Output variables
 real(dp), intent(out) :: psh(npsh,ps%nq),pal(ps%nq)

!Local variables
 real(dp) :: da(ps%nq)
 integer :: ii,kk,iq

 psh(1,:)=0.5d0*doubled_angle(ps,epsh2)
 do ii=2,npsh
   da(:)=doubled_angle(ps,energy(ii))
   do iq=1,ps%nq
     psh(ii,iq)=follow(ps,iq,psh(ii-1,iq),energy(ii-1),energy(ii),da(iq),0)
   end do
 end do
 kk=min(npsh,max(1,nint((epsh2-eal)/depsh)+1))
 da(:)=doubled_angle(ps,eal)
 do iq=1,ps%nq
   pal(iq)=follow(ps,iq,psh(kk,iq),energy(kk),eal,da(iq),0)
 end do
 return

 contains

 real(dp) function energy(ii)
 integer, intent(in) :: ii
 energy=epsh2-depsh*real(ii-1,dp)
 end function energy

 end subroutine phsft_scan

 subroutine phsft_pair(irel,ll,kap,nvkb,rr,vfull,zz,vloc,vlocpw,vkbpw, &
&                      evkb,mmax,nq,rq,npsh,epsh2,depsh,eal,pshae,pshps)

!! the all-electron and the pseudopotential phases of one l and kappa at
!! the radii rq and the energies of phsft_scan, the pseudopotential's
!! shifted by the multiple of pi that makes them equal the all-electron
!! ones at the reference energy eal, where the two solutions agree beyond
!! the core radius
!!
!! @param irel  0 non-, 1 scalar-, 2 fully-relativistic
!! @param ll  angular momentum
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param nvkb  number of projectors
!! @param rr  log radial mesh
!! @param vfull  all-electron potential
!! @param zz  atomic number
!! @param vloc  local potential on the mesh
!! @param vlocpw  its pieces
!! @param vkbpw  the projectors' pieces
!! @param evkb  their coefficients
!! @param mmax  size of the log radial mesh
!! @param nq  number of radii
!! @param rq  radii of the phases
!! @param npsh  number of energies
!! @param epsh2  the first (highest) energy (Ha)
!! @param depsh  the energy step (downwards)
!! @param eal  the reference energy of the alignment
!! @param pshae  the all-electron phases, per energy and radius
!! @param pshps  the pseudopotential phases, likewise

!Input variables
 integer, intent(in) :: irel,ll,kap,nvkb,mmax,nq,npsh
 real(dp), intent(in) :: rr(mmax),vfull(mmax),zz,vloc(mmax),evkb(nvkb)
 real(dp), intent(in) :: rq(nq),epsh2,depsh,eal
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)

!Output variables
 real(dp), intent(out) :: pshae(npsh,nq),pshps(npsh,nq)

!Local variables
 real(dp) :: palae(nq),palps(nq)
 integer :: iq
 type(phsft_t), target :: ps

 call phsft_ae_init(irel,ll,kap,rr,vfull,zz,mmax,nq,rq,ps)
 call phsft_scan(ps,epsh2,depsh,npsh,eal,pshae,palae)
 call phsft_free(ps)
 call phsft_ps_init(ll,nvkb,rr,vloc,vlocpw,vkbpw,evkb,mmax,nq,rq,ps)
 call phsft_scan(ps,epsh2,depsh,npsh,eal,pshps,palps)
 call phsft_free(ps)
 do iq=1,nq
   pshps(:,iq)=pshps(:,iq)+pi*anint((palae(iq)-palps(iq))/pi)
 end do
 return
 end subroutine phsft_pair

 real(dp) function phsft_node_phase(ps,iq,ee,nodes,uq_out)

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
!! @param nodes  whether to subtract pi per node (default true)
!! @param uq_out  u at the radius (output)

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 integer, intent(in) :: iq
 real(dp), intent(in) :: ee
 logical, intent(in), optional :: nodes

!Output variables
 real(dp), intent(out), optional :: uq_out

!Local variables
 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp) :: uq,upq
 integer :: nn,ierr,jj,nnode,jlast

 allocate(uu(ps%mmax,2),up(ps%mmax,2))
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 if(ps%kind==1) then
   nn=0
   call scatterstate(ps%irel,nn,ps%ll,ps%kap,ierr,ee,ps%rr,ps%zz,ps%vv, &
&                    uu,up,ps%mmax,ps%mch)
 else
   call outwf_solve(ps%ow,ee,ps%evkb,ps%mch,ps%mch,uu(:,1),up(:,1))
 end if
 call lm_interp(ps%rr(1:ps%mch),uu(1:ps%mch,1),ps%mch,ps%rq(iq),uq)
 call lm_interp(ps%rr(1:ps%mch),up(1:ps%mch,1),ps%mch,ps%rq(iq),upq)
 nnode=0
 jlast=1
 do jj=2,ps%mch
   if(ps%rr(jj)>=ps%rq(iq)) exit
   if(uu(jj,1)*uu(jj-1,1)<0.0d0) nnode=nnode+1
   jlast=jj
 end do
 if(uq*uu(jlast,1)<0.0d0) nnode=nnode+1
! r phi'/phi = r u'/u - 1, and al r u' is the solvers' derivative
 phsft_node_phase=atan2(upq/ps%al-uq,uq)
 if(phsft_node_phase>0.5d0*pi) phsft_node_phase=phsft_node_phase-pi
 if(phsft_node_phase<=-0.5d0*pi) phsft_node_phase=phsft_node_phase+pi
 if(present(uq_out)) uq_out=uq
 if(present(nodes)) then
   if(.not. nodes) nnode=0
 end if
 phsft_node_phase=phsft_node_phase-pi*dble(nnode)
 deallocate(uu,up)
 return
 end function phsft_node_phase

 recursive function follow(ps,iq,th0,e0,e1,da1,depth) result(th1)

!! the phase of radius iq at e1, continued from th0 at e0
!!
!! @param ps  the solver
!! @param iq  the radius
!! @param th0  the phase at e0
!! @param e0  the energy of th0
!! @param e1  the energy wanted
!! @param da1  the doubled angle at e1
!! @param depth  halvings so far

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 real(dp), intent(in) :: th0,e0,e1,da1
 integer, intent(in) :: iq,depth

!Output variables
 real(dp) :: th1

!Local variables
 real(dp) :: dd,em,dam(ps%nq)

 dd=da1-2.0d0*th0
 dd=dd-2.0d0*pi*anint(dd/(2.0d0*pi))
 if(abs(dd)>0.5d0*pi .and. depth<maxdepth) then
   em=0.5d0*(e0+e1)
   dam(:)=doubled_angle(ps,em)
   th1=follow(ps,iq,follow(ps,iq,th0,e0,em,dam(iq),depth+1),em,e1,da1, &
&             depth+1)
 else
   th1=th0+0.5d0*dd
 end if
 return
 end function follow

 function doubled_angle(ps,ee) result(da)

!! twice the angle of (phi, r phi') at each radius, in (-pi, pi]: that of
!! (u, u' - u/r) at the energy ee, the same for -u
!!
!! @param ps  the solver
!! @param ee  energy (Ha)

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Input variables
 real(dp), intent(in) :: ee

!Output variables
 real(dp) :: da(ps%nq)

!Local variables
 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp) :: uq,upq,xx,yy,hh
 integer :: nn,ierr,iq

 allocate(uu(ps%mmax,2),up(ps%mmax,2))
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 if(ps%kind==1) then
   nn=0
   call scatterstate(ps%irel,nn,ps%ll,ps%kap,ierr,ee,ps%rr,ps%zz,ps%vv, &
&                    uu,up,ps%mmax,ps%mch)
 else
   call outwf_solve(ps%ow,ee,ps%evkb,ps%mch,ps%mch,uu(:,1),up(:,1))
 end if
 do iq=1,ps%nq
   call lm_interp(ps%rr(1:ps%mch),uu(1:ps%mch,1),ps%mch,ps%rq(iq),uq)
   call lm_interp(ps%rr(1:ps%mch),up(1:ps%mch,1),ps%mch,ps%rq(iq),upq)
! r phi'/phi = r u'/u - 1, and al r u' is the solvers' derivative
   hh=hypot(uq,upq/ps%al-uq)
   xx=uq/hh
   yy=(upq/ps%al-uq)/hh
   da(iq)=atan2(2.0d0*xx*yy,(xx-yy)*(xx+yy))
 end do
 deallocate(uu,up)
 return
 end function doubled_angle

 end module phsft_m
