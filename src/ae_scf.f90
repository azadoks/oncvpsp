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
! Convergence measure, energy terms and mixing of the all-electron SCF of
! sratom and relatom, which mix the screening (Hartree plus exchange-
! correlation) potential and add the nuclear -zz/r to it for the radial
! equations; psatom shares the convergence test and the one-step mixing

 module ae_scf_m
 implicit none
 private
 public :: ae_residual, ae_energies, ae_mix, ae_converged, ae_floor_warning
 public :: ae_anderson1, ae_outer_peak

 contains

 subroutine ae_residual(rho,vsi,vso,sf,rr,mmax,res)

!! root mean square over the electrons of the change of the screening
!! potential from input to output, the size of the first-order change of
!! the eigenvalues (the convergence measure of RadialAtomicSolvers.jl)
!!
!! @param rho  total charge density (r**2 rho integrates to sf)
!! @param vsi  input screening potential
!! @param vso  output screening potential
!! @param sf  number of electrons
!! @param rr  log radial mesh
!! @param mmax  dimension of log radial grid
!! @param res  residual (Ha)

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rho(mmax),vsi(mmax),vso(mmax),rr(mmax)
 real(dp), intent(in) :: sf

!Output variables
 real(dp), intent(out) :: res

!Local variables
 real(dp) :: al,sn
 integer :: ii

 al = 0.01d0 * log(rr(101) / rr(1))

 sn=0.0d0
 do ii=1,mmax
   sn=sn + rho(ii)*rr(ii)**3*(vso(ii)-vsi(ii))**2
 end do
 res=sqrt(al*sn/sf)

 return
 end subroutine ae_residual

 subroutine ae_energies(rho,vsi,vxc,zz,eeig,eeel,eexc,rr,mmax,etot,eterm, &
&                       vext)

!! total energy and its terms as RadialAtomicSolvers.jl and NIST define
!! them, the kinetic energy from the eigenvalue sum less the potential
!! energy in the input potential the orbitals are eigenstates of, so that
!! the total energy is stationary in that potential
!!
!! @param rho  total charge density of the orbitals of vsi
!! @param vsi  input screening potential
!! @param vxc  exchange-correlation potential of rho (vout)
!! @param zz  atomic number
!! @param eeig  sum of the occupied eigenvalues
!! @param eeel  electron-electron interaction of rho (vout)
!! @param eexc  exchange-correlation correction of rho (vout)
!! @param rr  log radial mesh
!! @param mmax  dimension of log radial grid
!! @param etot  total energy (Ha)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (Ha)
!! @param vext  external potential the orbitals also see (a confining wall,
!!              oncvpsp_confine.x; zero for the free atom), whose energy is
!!              in etot but in no term

 use constants_m, only: dp
 use logmesh_m, only: lm_head
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rho(mmax),vsi(mmax),vxc(mmax),rr(mmax),vext(mmax)
 real(dp), intent(in) :: zz,eeig,eeel,eexc

!Output variables
 real(dp), intent(out) :: etot,eterm(4)

!Local variables
 real(dp) :: al,ekin,enuc,ecoul,exc,escr,eext
 integer :: ii

 real(dp), allocatable :: ff(:)

 allocate(ff(mmax))

 al = 0.01d0 * log(rr(101) / rr(1))

! electron-nucleus energy; the integrand goes as rr**2 below rr(1)
 do ii=1,mmax
   ff(ii)=-zz*rho(ii)*rr(ii)**2
 end do
 enuc=radint(ff,2.0d0)

! potential energy in the input screening potential
 do ii=1,mmax
   ff(ii)=rho(ii)*vsi(ii)*rr(ii)**3
 end do
 escr=radint(ff,3.0d0)

! energy in the external potential
 eext=0.0d0
 if(any(vext(:)/=0.0d0)) then
   do ii=1,mmax
     ff(ii)=rho(ii)*vext(ii)*rr(ii)**3
   end do
   eext=radint(ff,3.0d0)
 end if

! exchange-correlation energy from vout's correction to the eigenvalue sum
 do ii=1,mmax
   ff(ii)=rho(ii)*vxc(ii)*rr(ii)**3
 end do
 exc=eexc + radint(ff,3.0d0)

 ecoul=0.5d0*eeel
 ekin=eeig - enuc - escr - eext

 etot=ekin + enuc + ecoul + exc + eext
 eterm(1)=ekin
 eterm(2)=enuc
 eterm(3)=ecoul
 eterm(4)=exc

 deallocate(ff)
 return

 contains

 function radint(ff,pw)

!! integral of ff/rr over rr on the log mesh, by vout's rule; the part
!! below rr(1) is that of ff going as rr**pw times the density (lm_head)

!Input variables
 real(dp), intent(in) :: ff(mmax)
 real(dp), intent(in) :: pw

!Output variables
 real(dp) :: radint

 radint=(9.0d0*ff(1) + 28.0d0*ff(2) + 23.0d0*ff(3))/24.0d0
 do ii=4,mmax
   radint=radint + ff(ii)
 end do
 radint=al*radint + lm_head(ff(1),rho(1),rho(2),al,pw)

 return
 end function radint

 end subroutine ae_energies

 subroutine ae_mix(it,vsi,vso,hv,hr,ndep,zz,rr,mmax,irel,c_light)

!! next input screening potential by Anderson mixing of the last ndep
!! iterations, with the residual filtered in ln(r) near the nucleus for the
!! scalar-relativistic equation (KerkerMixing of RadialAtomicSolvers.jl)
!!
!! Inside r* = (1-gamma)/zz ~ zz/(2c**2) the relativistic mass is set by
!! the potential, which with a GGA closes a loop potential -> orbital
!! exponent -> density gradient -> potential. Its modes there grow with the
!! number of mesh points per decade, and the one-step Anderson mixing of
!! sratom diverged for light atoms on fine meshes. The filter damps them;
!! it does not change the fixed point.
!!
!! @param it  iteration number, 1 for the first
!! @param vsi  input screening potential, the next one on output
!! @param vso  output screening potential
!! @param hv  history of input potentials, column 1 the newest
!! @param hr  history of residuals vso - vsi
!! @param ndep  depth of the history
!! @param zz  atomic number
!! @param rr  log radial mesh
!! @param mmax  dimension of log radial grid
!! @param irel  0 non-relativistic, 1 scalar-relativistic
!! @param c_light  speed of light (Hartree atomic units)

 use constants_m, only: dp
 use lapack_m, only: dgesv
 implicit none

!Input variables
 integer, intent(in) :: it,ndep,mmax,irel
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: vso(mmax),rr(mmax),zz

!Input/output variables
 real(dp), intent(inout) :: vsi(mmax),hv(mmax,ndep),hr(mmax,ndep)

!Local variables
! blend parameter of the residual
 real(dp), parameter :: bl=0.5d0
! filter strength bfil/zz**2 (in ln(r) squared) inside rfil*r*
 real(dp), parameter :: bfil=1.0d-4, rfil=10.0d0
 real(dp) :: aa(ndep,ndep),bb(ndep)
 real(dp) :: al,bk,rs
 integer :: ii,jj,kk,nh,info
 integer :: ipiv(ndep)

 real(dp), allocatable :: res(:),vb(:),ww(:)

 allocate(res(mmax),vb(mmax),ww(mmax))

 al = 0.01d0 * log(rr(101) / rr(1))

 do jj=ndep,2,-1
   hv(:,jj)=hv(:,jj-1)
   hr(:,jj)=hr(:,jj-1)
 end do
 hv(:,1)=vsi(:)
 hr(:,1)=vso(:)-vsi(:)
 nh=min(it,ndep)

! the combination of the history of least residual, in the r**2 weight
! of the one-step mixing this replaces
 do ii=1,mmax
   ww(ii)=rr(ii)**2
 end do
 res(:)=hr(:,1)
 vb(:)=hv(:,1)
 if(nh>1) then
   do jj=2,nh
     do kk=2,nh
       aa(jj-1,kk-1)=sum(ww(:)*(hr(:,1)-hr(:,jj))*(hr(:,1)-hr(:,kk)))
     end do
     bb(jj-1)=sum(ww(:)*hr(:,1)*(hr(:,1)-hr(:,jj)))
   end do
   do jj=1,nh-1
     aa(jj,jj)=aa(jj,jj)*(1.0d0+1.0d-10)
   end do
   call dgesv(nh-1,1,aa,ndep,ipiv,bb,ndep,info)
   if(info==0) then
     do jj=2,nh
       res(:)=res(:)-bb(jj-1)*(hr(:,1)-hr(:,jj))
       vb(:)=vb(:)-bb(jj-1)*(hv(:,1)-hv(:,jj))
     end do
   end if
 end if

 if(irel/=0) then
   rs=rfil*(1.0d0-sqrt(1.0d0-(zz/c_light)**2))/zz
   bk=bfil/zz**2
   call kerker_filter()
 end if

 do ii=1,mmax
   vsi(ii)=vb(ii)+bl*res(ii)
 end do

 deallocate(res,vb,ww)
 return

 contains

 subroutine kerker_filter()

!! res <- (1 + L)**(-1) res, L = -d/dx beta d/dx in x = ln(r), with
!! beta = bk/(1 + (r/rs)**2), by finite volumes with no flux at the ends

 real(dp), allocatable :: beta(:),dd(:)

 allocate(beta(mmax),dd(mmax))

! beta on the face between ii and ii+1, in mesh-index units
 do ii=1,mmax-1
   beta(ii)=bk/(1.0d0+rr(ii)*rr(ii+1)/rs**2)/al**2
 end do
 beta(mmax)=0.0d0

! tridiagonal solve, sub- and superdiagonal -beta
 dd(1)=1.0d0+beta(1)
 do ii=2,mmax
   dd(ii)=1.0d0+beta(ii)+beta(ii-1)
   dd(ii)=dd(ii)-beta(ii-1)**2/dd(ii-1)
   res(ii)=res(ii)+beta(ii-1)*res(ii-1)/dd(ii-1)
 end do
 res(mmax)=res(mmax)/dd(mmax)
 do ii=mmax-1,1,-1
   res(ii)=(res(ii)+beta(ii)*res(ii+1))/dd(ii)
 end do

 deallocate(beta,dd)
 return
 end subroutine kerker_filter

 end subroutine ae_mix

 subroutine ae_converged(it,res,tolres,resbest,nflat,convg,iout)

!! whether the SCF is converged: convg on entry says each eigenvalue has
!! settled; the residual must also be below tolres, or have stopped falling
!! at its round-off floor (not halved in nstall iterations, and below
!! resfloor). That floor depends on the mesh and the functional: 3-4e-12
!! for Nd and Bi on the default mesh, ~2e-11 for PBE N on amesh 1.004 out
!! to 45 bohr (rr1 7e-5), where a fixed 1e-11 was never reached; a
!! stop on the floor is reported as a warning
!!
!! @param it  iteration
!! @param res  residual of this iteration (ae_residual)
!! @param tolres  threshold of the residual
!! @param resbest  smallest residual so far by halvings (kept between calls)
!! @param nflat  iterations since it last halved (kept between calls)
!! @param convg  in: eigenvalues settled; out: converged
!! @param iout  unit of the text output

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: it
 real(dp), intent(in) :: res,tolres

!Input/output variables
 real(dp), intent(inout) :: resbest
 integer, intent(inout) :: nflat
 logical, intent(inout) :: convg
 integer, intent(in) :: iout

!Local variables
 integer, parameter :: nstall=5
 real(dp), parameter :: resfloor=1.0d-9

 if(it==1 .or. res<0.5d0*resbest) then
   resbest=res
   nflat=0
 else
   nflat=nflat+1
 end if
 if(it<2) convg=.false.
 if(res>tolres .and. (res>resfloor .or. nflat<nstall)) convg=.false.
 if(convg .and. res>tolres) call ae_floor_warning(res,tolres,iout)
 return
 end subroutine ae_converged

 subroutine ae_floor_warning(res,tolres,iout)

!! the warning of an SCF stopped on the round-off floor of its residual
!! (ae_converged; aeatom_cached repeats it for a solution it reuses)
!!
!! @param res  the residual at the stop
!! @param tolres  its threshold
!! @param iout  unit of the text output

 use constants_m, only: dp
 implicit none

!Input variables
 real(dp), intent(in) :: res,tolres
 integer, intent(in) :: iout

 write(iout,'(/2a,es12.4,a,es12.4)') 'ae_converged: WARNING scf stopped ', &
&      'on the round-off floor of the residual, res=',res,', tolres=',tolres
 return
 end subroutine ae_floor_warning

 subroutine ae_anderson1(mix,vi,vo,vi1,vo1,rr,mmax,bl)

!! the next input potential by D. G. Anderson's method with one previous
!! iteration (psatom, relatom): the input and output potentials vi, vo
!! and those before, vi1, vo1, blended with bl, and those saved for the
!! next; with mix false, plain linear mixing (thl = 0)
!!
!! @param mix  whether the previous iteration enters
!! @param vi  input potential, the next on output
!! @param vo  output potential
!! @param vi1,vo1  the previous input and output potentials, this
!!                 iteration's on output
!! @param rr  log radial mesh
!! @param mmax  size of the mesh
!! @param bl  the blend of the output potential

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax
 logical, intent(in) :: mix
 real(dp), intent(in) :: vo(mmax),rr(mmax),bl

!Input/output variables
 real(dp), intent(inout) :: vi(mmax),vi1(mmax),vo1(mmax)

!Local variables
 real(dp) :: thl,sn,sd,rl,rl1,dr,vn
 integer :: ii

 thl=0.0d0
 if(mix) then
   sn=0.0d0
   sd=0.0d0
   do ii=1,mmax
     rl=vo(ii)-vi(ii)
     rl1=vo1(ii)-vi1(ii)
     dr=rl-rl1
     sn=sn + rl*dr*rr(ii)**2
     sd=sd + dr*dr*rr(ii)**2
   end do
   thl=sn/sd
 end if

 do ii=1,mmax
   vn=(1.0d0-bl)*((1.0d0-thl)*vi(ii) + thl*vi1(ii)) &
&   + bl*((1.0d0-thl)*vo(ii) + thl*vo1(ii))
   vi1(ii)=vi(ii)
   vo1(ii)=vo(ii)
   vi(ii)=vn
 end do
 return
 end subroutine ae_anderson1

 subroutine ae_outer_peak(up,rr,mch,rpk)

!! the radius of the outermost extremum of a bound state inside mch, where
!! its derivative up changes sign; rpk unchanged if there is none
!!
!! @param up  the derivative of the state on the mesh
!! @param rr  log radial mesh
!! @param mch  the last point
!! @param rpk  the radius (output)

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mch
 real(dp), intent(in) :: up(mch),rr(mch)

!Input/output variables
 real(dp), intent(inout) :: rpk

!Local variables
 integer :: jj

 do jj=mch-1,1,-1
   if(up(jj)*up(jj+1)<0.0d0) then
     rpk=rr(jj)
     exit
   end if
 end do
 return
 end subroutine ae_outer_peak

 end module ae_scf_m
