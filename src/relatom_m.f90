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

 module relatom_m

 implicit none
 private
 public :: relatom

 contains

 subroutine relatom(na,la,ea,fa,nc,ncv,it,rhoc,rho, &
&           rr,vi,zz,mmax,adams,iexc,c_light,etot,eterm,res,ierr,vext,warm, &
&           iout,ifail)

!! self-consistent fully-relativistic all-electron atom
!! calculation using log mesh
!! partially-occupied orbitals are weighted proportionally to 2j+1
!!
!! @param na  principal quantum numbers, dimension ncv
!! @param la  angular momenta
!! @param ea  eigenvalues (output; input with warm), ea(:,1) for j=l+1/2,
!!            ea(:,2) for j=l-1/2
!! @param fa  occupancies
!! @param nc  number of core states
!! @param ncv  number of core+valence states
!! @param it  number of iterations (output)
!! @param rhoc  core charge density (output)
!! @param rho  total core+valence charge (output)
!! @param rr  log radial mesh (bohr)
!! @param vi  all-electron potential (output; input with warm), the one the
!!            output orbitals are eigenstates of
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param iexc  exchange-correlation functional
!! @param c_light  speed of light (Hartree atomic units)
!! @param etot  all-electron total energy (output)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (Ha), output
!! @param res  residual of the screening potential (output, ae_residual)
!! @param ierr  error flag
!! @param vext  external potential the orbitals also see (aeatom)
!! @param warm  start from vi and ea (the solution of a nearby potential,
!!              oncvpsp_confine.x) instead of the Thomas-Fermi potential
!!              and hydrogenic energies
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_convergence
 use constants_m, only: dp, mxstate, mxkap, ae_tolres
 use kappa_m, only: kap_count, kap_of, kap_j_share
 use ae_scf_m, only: ae_energies, ae_residual, ae_converged, &
&                    ae_mix
 use tfapot_m, only: tfapot
 use vout_m, only: vout
 use ldiracfb_m, only: ldiracfb
 use logmesh_m, only: lm_step
 use adams_m, only: adams_t

!Input variables

 integer, intent(in) :: mmax,iexc,nc,ncv
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 integer, intent(in) :: na(ncv),la(ncv)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(ncv),rr(mmax),vext(mmax)
 integer, intent(in) :: iout
 logical, intent(in) :: warm

!Input/output variables
 integer, intent(inout) :: ierr
 real(dp), intent(inout) :: ea(mxstate,mxkap)
 real(dp), intent(inout) :: vi(mmax)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot,eterm(4),res
 real(dp), intent(out) :: rho(mmax),rhoc(mmax)
 integer, intent(out) :: ifail

!Local variables
 integer :: nin,mch
 real(dp) :: amesh,al
 real(dp) :: eeel,eexc,et,sf,eeig
 real(dp) :: zion,fj
 integer :: ii,ll,ikap,kap,mkap
 logical :: convg
 real(dp) :: resbest
 integer :: nflat

 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp), allocatable :: vsi(:),vso(:),vxc(:)
 real(dp), allocatable :: hv(:,:),hr(:,:)

! blend parameter for Anderson iterative potential mixing
!the depth of the mixing's history (ae_mix)
 integer, parameter :: ndep=8
! convergence threshold of the residual and the eigenvalues (constants_m)
 real(dp), parameter ::  tolres=ae_tolres

 ifail=0
 allocate(uu(mmax,2),up(mmax,2))
 allocate(vsi(mmax),vso(mmax),vxc(mmax))
 allocate(hv(mmax,ndep),hr(mmax,ndep))

! why all this is necessary is unclear, but it seems to be
 uu(:,:)=0.0d0; up(:,:)=0.0d0; vso(:)=0.0d0; hv(:,:)=0.0d0; hr(:,:)=0.0d0
 vxc(:)=0.0d0
 eeel=0.0d0; eexc=0.0d0; et=0.0d0
 sf=0.0d0; eeig=0.0d0; zion=0.0d0
 nin=0; mch=0

 al=lm_step(rr,mmax)
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
   do ikap=1,mxkap
     ea(ii,ikap)=-0.5d0*(zion/na(ii))**2
     if(ea(ii,ikap)>vi(mmax)) ea(ii,ikap)=2.0d0*vi(mmax)
   end do
   if(la(ii)==0) ea(ii,2)=0.0d0
 end do

! big self  self-consietency loop

 do it=1,100
   convg=.true.

   rhoc(:) = 0.0d0
   rho(:)=0.0d0

! solve for bound states in turn
   eeig=0.0d0
   do ii=1,ncv
     ll = la(ii)
     mkap=kap_count(2,ll)
! loop on J = ll +/- 1/2
     do ikap=1,mkap
       kap=kap_of(ll,ikap)
       fj=kap_j_share(ll,ikap,fa(ii))
       et=ea(ii,ikap)
       ierr = 0

       call ldiracfb(c_light,na(ii),ll,kap,ierr,et,rr,zz,vext+vi,uu,up,mmax, &
&                    adams, &
&                    mch,iout,ifail)
       if(ifail/=0) return
! a state with no bound solution in this potential gets its box state;
! the SCF goes on but cannot converge with it
       if(ierr == 10) then
         convg=.false.
         ierr=0
       end if

       if(ierr /= 0) then
         write(iout,'(/2a,5i4)') 'relatom: ldiracfb convergence ERROR', &
&             ' n,l,kap,iter,ierr=',na(ii),ll,kap,it,ierr
         ifail=fail_convergence
         return
       end if

! each eigenvalue settled to tolres (see the convergence test below)
       if(abs(ea(ii,ikap)-et)>tolres*max(abs(et),1.0d0)) convg=.false.
       ea(ii,ikap)=et

! accumulate charge and eigenvalues
       eeig = eeig + fj*ea(ii,ikap)
       rho(:)=rho(:) + fj*((uu(:,1)/rr(:))**2 + (uu(:,2)/rr(:))**2)
       if(ii<=nc) then
         rhoc(:)=rhoc(:) + fj*((uu(:,1)/rr(:))**2 + (uu(:,2)/rr(:))**2)
       end if

     end do
   end do

   if(ierr/=0) then
     exit
   end if


! output screening potential, and the energies of the orbitals of the
! input one
   call vout(.false.,rho,rhoc,vso,vxc,sf,eeel,eexc, &
&            rr,mmax,adams,iexc,c_light)
   call ae_energies(rho,vsi,vxc,zz,eeig,eeel,eexc,rr,mmax,etot,eterm, &
&                   vext)

! converged when no eigenvalue has moved by more than tolres and the
! screening potential is self-consistent; the eigenvalues alone can stand
! still from a start far from self-consistency. ldiracfb returns the
! energy of its last pass, never the guess unchanged, so unchanged
! eigenvalues cannot be asked for (as in sratom)
   call ae_residual(rho,vsi,vso,sf,rr,mmax,res)
   call ae_converged(it,res,tolres,resbest,nflat,convg,iout)
   if(convg) exit

! next input potential, by Anderson mixing of the last ndep iterations
   call ae_mix(it,vsi,vso,hv,hr,ndep,zz,rr,mmax,2,c_light)
   do ii=1,mmax
     vi(ii)=vsi(ii) - zz/rr(ii)
   end do

   if(it==100 .and. .not. convg) then
     write(iout,'(/a)') 'relatom: WARNING failed to converge'
   end if

 end do !it

 if(.not. convg .and. ierr==0) then
   ierr=100
 end if

 deallocate(uu,up)
 deallocate(vsi,vso,vxc)
 deallocate(hv,hr)
 return

 end subroutine relatom

 end module relatom_m
