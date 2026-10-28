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

 module sratom_m

 implicit none
 private
 public :: sratom

 contains

 subroutine sratom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho,rr,vi,zz,mmax,iexc,etot, &
&           eterm,res,ierr,irel,c_light,vext,warm,iout,ifail)

!! self-consistent scalar-relativistic all-electron atom
!! calculation using log mesh (non-relativistic when irel=0)
!!
!! @param na  principal quantum number array, dimension ncv
!! @param la  angular-momenta
!! @param ea  eigenvalues (output; input with warm)
!! @param fa  occupancies
!! @param rpk  radius of outermost peak of wave function
!! @param nc  number of core states
!! @param ncv  number of core+valence states
!! @param it  number of iterations (output)
!! @param rhoc  core charge (output)
!! @param rho  total core+valence charge (output)
!! @param rr  log radial mesh
!! @param vi  all-electron potential (output; input with warm), the one the
!!            output orbitals are eigenstates of
!! @param zz  atomic number
!! @param mmax  size of log grid
!! @param iexc  exchange-correlation function to be used
!! @param etot  all-electron total energy (output)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (output)
!! @param res  residual of the screening potential (output, ae_residual)
!! @param ierr  error flag
!! @param irel  0 non-relativistic, 1 scalar-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param vext  external potential the orbitals also see (aeatom)
!! @param warm  start from vi and ea (the solution of a nearby potential,
!!              oncvpsp_confine.x) instead of the Thomas-Fermi potential
!!              and hydrogenic energies
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, ae_tolres
 use ae_scf_m, only: ae_energies, ae_residual, ae_converged, ae_mix, &
&                    ae_outer_peak
 use tfapot_m, only: tfapot
 use vout_m, only: vout
 use lschfb_m, only: lschfb

!Input variables

 integer, intent(in) :: mmax,iexc,nc,ncv
 integer, intent(in) :: na(ncv),la(ncv)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(ncv),rr(mmax)
 integer, intent(in) :: irel
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: vext(mmax)

!Input/output variables
 integer, intent(inout) :: ierr
 real(dp), intent(inout) :: rpk(ncv)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot,eterm(4),res
 real(dp), intent(inout) :: ea(ncv)
 real(dp), intent(out) :: rho(mmax),rhoc(mmax)
 real(dp), intent(inout) :: vi(mmax)
 logical, intent(in) :: warm
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 integer :: nin,mch
 real(dp) :: amesh,al
 real(dp) :: eeel,eexc,et,sf,eeig
 real(dp) :: zion
 integer :: ii
 logical :: convg
 real(dp) :: resbest
 integer :: nflat

 real(dp), allocatable :: uu(:),up(:)
 real(dp), allocatable :: vsi(:),vso(:),vxc(:)
 real(dp), allocatable :: hv(:,:),hr(:,:)

! depth of the history of the Anderson mixing
 integer, parameter :: ndep=8
! convergence threshold of the residual and the eigenvalues (constants_m)
 real(dp), parameter ::  tolres=ae_tolres

 ifail=0
 allocate(uu(mmax),up(mmax))
 allocate(vsi(mmax),vso(mmax),vxc(mmax))
 allocate(hv(mmax,ndep),hr(mmax,ndep))

! why all this is necessary is unclear, but it seems to be
 uu(:)=0.0d0; up(:)=0.0d0; vso(:)=0.0d0
 vxc(:)=0.0d0; hv(:,:)=0.0d0; hr(:,:)=0.0d0
 eeel=0.0d0; eexc=0.0d0; et=0.0d0
 sf=0.0d0; eeig=0.0d0; zion=0.0d0
 nin=0; mch=0

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

 if(.not. warm) then
   do ii=1,mmax
     vi(ii)=tfapot(rr(ii),zz)
   end do
 end if

! the screening (Hartree plus exchange-correlation) potential is mixed,
! and the nuclear potential added to it for the radial equation
 do ii=1,mmax
   vsi(ii)=vi(ii) + zz/rr(ii)
 end do

! starting approximation for energies
 sf=0.0d0
 do ii=1,ncv
   sf=sf+fa(ii)
   zion=zz+1.0d0-sf
   if(warm) cycle
   ea(ii)=-0.5d0*(zion/na(ii))**2
   if(ea(ii)>vi(mmax)) ea(ii)=2.0d0*vi(mmax)
 end do

! big self  self-consietency loop

 do it=1,100
   convg=.true.

   rhoc(:) = 0.0d0
   rho(:)=0.0d0

! solve for bound states in turn
   eeig=0.0d0
   do ii=1,ncv

! skip unoccupied states
     if(fa(ii)==0.0d0) then
       ea(ii)=0.0d0
       cycle
     end if
     et=ea(ii)
     ierr = 0
     call lschfb(na(ii),la(ii),ierr,et, &
&                rr,vext+vi,uu,up,zz,mmax,mch,irel,c_light,iout,ifail)
     if(ifail/=0) return
! a state with no bound solution in this potential gets its box state;
! the SCF goes on but cannot converge with it
     if(ierr == 10) then
       convg=.false.
       ierr=0
     end if
     if(ierr /= 0) then
       write(iout,'(/a,3i4)') 'sratom123: lschfb convergence ERROR n,l,iter=', &
&       na(ii),la(ii),it
       ifail=1
       return
     end if

! each eigenvalue settled to tolres (see the convergence test below)
     if(abs(ea(ii)-et)>tolres*max(abs(et),1.0d0)) convg=.false.
     ea(ii)=et

! accumulate charge and eigenvalues
     eeig = eeig + fa(ii) * ea(ii)
     rho(:)=rho(:) + fa(ii)*(uu(:)/rr(:))**2
     if(ii<=nc) then
       rhoc(:)=rhoc(:) + fa(ii)*(uu(:)/rr(:))**2
     end if


! find outermost peak of wavefunction
     call ae_outer_peak(up,rr,mch,rpk(ii))

   end do

   if(ierr/=0) then
     exit
   end if


! output screening potential, and the energies of the orbitals of the
! input one
   call vout(0,rho,rhoc,vso,vxc,sf,eeel,eexc, &
&            rr,mmax,iexc,c_light)
   call ae_energies(rho,vsi,vxc,zz,eeig,eeel,eexc,rr,mmax,etot,eterm, &
&                   vext)

! converged when no eigenvalue has moved by more than tolres and the
! screening potential is self-consistent; the eigenvalues alone can stand
! still from a start far from self-consistency. Unchanged eigenvalues were
! once asked for, which a near-neutral mode of the scalar-relativistic GGA
! potential at the nucleus can deny for ever (1e-9 Ha jitter with E_tot
! converged to 1e-14 Ha, H at rr1 1e-5)
   call ae_residual(rho,vsi,vso,sf,rr,mmax,res)
   call ae_converged(it,res,tolres,resbest,nflat,convg,iout)
   if(convg) exit

! next input potential, by Anderson mixing of the last ndep iterations
   call ae_mix(it,vsi,vso,hv,hr,ndep,zz,rr,mmax,irel,c_light)
   do ii=1,mmax
     vi(ii)=vsi(ii) - zz/rr(ii)
   end do

   if(it==100 .and. .not. convg) then
     write(iout,'(/a)') 'sratom: WARNING failed to converge'
   end if

 end do !it

 if(.not. convg .and. ierr==0) then
   ierr=100
 end if

 deallocate(uu,up)
 deallocate(vsi,vso,vxc)
 deallocate(hv,hr)
 return

 end subroutine sratom

 end module sratom_m
