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
!> The "phase shifts" of run_phsft, atan(r phi'/phi) of the scattering
!> solution phi = u/r at a radius, as continuous functions of the energy,
!> of the radius and of the potential (doc/smoothness_plan.md, X-2). The
!> radius is any radius: (u, al r u') are interpolated to it from the mesh
!> (lm_interp). The angle is that of the point (phi, r phi'), defined up to
!> pi only, as the overall sign of a solution is arbitrary (the outward
!> solution of vkboutwf flips it as the energy crosses a pole of its
!> coefficients): its doubled angle is sign-invariant. The phase is
!> followed from one energy to the next by continuation, taking the change
!> of the doubled angle in (-pi, pi] and halving the energy step until no
!> step turns it by more than pi/2, so the result is the same continuous
!> function of the energy whatever steps were taken: no thresholds, no
!> 2 pi or pi corrections, no search for spurious jumps.

 module phsft_m
 use constants_m, only: dp, pi
 use pieces_m, only: piecewise_t
 use outwf_m, only: outwf_t, outwf_setup, outwf_solve, outwf_free
 use logmesh_m, only: lm_interp, lm_fh_reach
 implicit none
 private
 public :: phsft_t, phsft_mch, phsft_ae_init, phsft_ps_init, phsft_free
 public :: phsft_scan

!> the deepest halving of an energy step: a phase that turns by more than
!> pi/2 within depsh/2**maxdepth is taken as it is found
 integer, parameter :: maxdepth=40

 type phsft_t
!> 1 for the all-electron potential (scatterstate), 2 for the
!> pseudopotential (outwf_m)
   integer :: kind=0
   integer :: irel=0,ll=0,kap=0,mmax=0,mch=0,nvkb=0
   real(dp) :: zz=0.0d0,rq=0.0d0,al=0.0d0
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

 subroutine phsft_ae_init(irel,ll,kap,rr,vv,zz,mmax,rq,ps)

!! the all-electron solver, at the radius rq
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2
!!              fully-relativistic
!! @param ll  angular momentum
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param rr  log radial mesh
!! @param vv  all-electron potential
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param rq  radius of the phases, with phsft_mch(rr,mmax,rq)>0
!! @param ps  the solver

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax
 real(dp), intent(in) :: rr(mmax),vv(mmax),zz,rq

!Output variables
 type(phsft_t), intent(out) :: ps

 ps%kind=1
 ps%irel=irel
 ps%ll=ll
 ps%kap=kap
 ps%zz=zz
 call common_init(rr,mmax,rq,ps)
 ps%vv=vv
 return
 end subroutine phsft_ae_init

 subroutine phsft_ps_init(ll,nvkb,rr,vloc,vlocpw,vkbpw,evkb,mmax,rq,ps)

!! the pseudopotential solver, local potential and nvkb projectors (none
!! for the semi-local potential that is the local one), at the radius rq
!!
!! @param ll  angular momentum
!! @param nvkb  number of projectors
!! @param rr  log radial mesh
!! @param vloc  local potential on the mesh
!! @param vlocpw  its pieces
!! @param vkbpw  the projectors' pieces
!! @param evkb  their coefficients
!! @param mmax  size of the log radial mesh
!! @param rq  radius of the phases, with phsft_mch(rr,mmax,rq)>0
!! @param ps  the solver

!Input variables
 integer, intent(in) :: ll,nvkb,mmax
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb),rq
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)

!Output variables
 type(phsft_t), intent(out), target :: ps

 ps%kind=2
 ps%ll=ll
 ps%nvkb=nvkb
 call common_init(rr,mmax,rq,ps)
 ps%evkb=evkb
 call outwf_setup(ll,rr,mmax,vloc,vlocpw,nvkb,vkbpw,ps%ow)
 return
 end subroutine phsft_ps_init

 subroutine common_init(rr,mmax,rq,ps)

!! the mesh, the radius and the last point of the solutions

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rq

!Input/output variables
 type(phsft_t), intent(inout) :: ps

 ps%mmax=mmax
 ps%rq=rq
 ps%rr=rr
 ps%al=0.01d0*log(rr(101)/rr(1))
 ps%mch=phsft_mch(rr,mmax,rq)
 if(ps%mch==0) then
   write(6,'(/a/a,1p,e12.4,a,e12.4)') &
&        'phsft_m: ERROR the radius is too close to the end of the mesh', &
&        ' rq=',rq,' rr(mmax)=',rr(mmax)
   error stop 1
 end if
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
 ps%kind=0
 return
 end subroutine phsft_free

 subroutine phsft_scan(ps,epsh2,depsh,npsh,eal,psh,pal)

!! the phase at the energies epsh2, epsh2-depsh, ... (npsh of them), each
!! followed from the one before, the first in (-pi/2, pi/2]; and the phase
!! at the energy eal, followed from the nearest of them
!!
!! @param ps  the solver
!! @param epsh2  the first (highest) energy (Ha)
!! @param depsh  the energy step (downwards)
!! @param npsh  number of energies
!! @param eal  the energy of pal
!! @param psh  the phases
!! @param pal  the phase at eal

!Input variables
 integer, intent(in) :: npsh
 real(dp), intent(in) :: epsh2,depsh,eal

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Output variables
 real(dp), intent(out) :: psh(npsh),pal

!Local variables
 integer :: ii,kk

 psh(1)=0.5d0*doubled_angle(ps,epsh2)
 do ii=2,npsh
   psh(ii)=follow(ps,psh(ii-1),energy(ii-1),energy(ii),0)
 end do
 kk=min(npsh,max(1,nint((epsh2-eal)/depsh)+1))
 pal=follow(ps,psh(kk),energy(kk),eal,0)
 return

 contains

 real(dp) function energy(ii)
 integer, intent(in) :: ii
 energy=epsh2-depsh*real(ii-1,dp)
 end function energy

 end subroutine phsft_scan

 recursive function follow(ps,th0,e0,e1,depth) result(th1)

!! the phase at e1, continued from th0 at e0
!!
!! @param ps  the solver
!! @param th0  the phase at e0
!! @param e0  the energy of th0
!! @param e1  the energy wanted
!! @param depth  halvings so far

!Input variables
 real(dp), intent(in) :: th0,e0,e1
 integer, intent(in) :: depth

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Output variables
 real(dp) :: th1

!Local variables
 real(dp) :: dd,em

 dd=doubled_angle(ps,e1)-2.0d0*th0
 dd=dd-2.0d0*pi*anint(dd/(2.0d0*pi))
 if(abs(dd)>0.5d0*pi .and. depth<maxdepth) then
   em=0.5d0*(e0+e1)
   th1=follow(ps,follow(ps,th0,e0,em,depth+1),em,e1,depth+1)
 else
   th1=th0+0.5d0*dd
 end if
 return
 end function follow

 real(dp) function doubled_angle(ps,ee)

!! twice the angle of (phi, r phi') at the radius, in (-pi, pi]: that of
!! (u, u' - u/r) at the energy ee, the same for -u
!!
!! @param ps  the solver
!! @param ee  energy (Ha)

!Input variables
 real(dp), intent(in) :: ee

!Input/output variables
 type(phsft_t), intent(inout), target :: ps

!Local variables
 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp) :: uq,upq,xx,yy,hh
 integer :: nn,ierr,node

 allocate(uu(ps%mmax,2),up(ps%mmax,2))
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 if(ps%kind==1) then
   nn=0
   call scatterstate(ps%irel,nn,ps%ll,ps%kap,ierr,ee,ps%rr,ps%zz,ps%vv, &
&                    uu,up,ps%mmax,ps%mch)
 else
   call outwf_solve(ps%ow,ee,ps%evkb,ps%mch,ps%mch,uu(:,1),up(:,1),node)
 end if
 call lm_interp(ps%rr(1:ps%mch),uu(1:ps%mch,1),ps%mch,ps%rq,uq)
 call lm_interp(ps%rr(1:ps%mch),up(1:ps%mch,1),ps%mch,ps%rq,upq)
! r phi'/phi = r u'/u - 1, and al r u' is the solvers' derivative
 hh=hypot(uq,upq/ps%al-uq)
 xx=uq/hh
 yy=(upq/ps%al-uq)/hh
 doubled_angle=atan2(2.0d0*xx*yy,(xx-yy)*(xx+yy))
 deallocate(uu,up)
 return
 end function doubled_angle

 end module phsft_m
