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

 module ionic_tail_m

 implicit none
 private
 public :: ionic_screening,coulomb_tail

 contains

! The ionic tail of the pseudopotentials, after PseudoPotentialGenerators.jl
! (src/unscreening.jl: ionic_tail, ionic_switch (here cinf_taper, in
! cinf_taper.f90), ionic_tail_window, with_ionic_tail;
! src/verification.jl: coulomb_tail).
!
! Beyond the largest core radius rcmax the screened potentials vp are V_AE,
! screened by the all-electron valence density n_v.  Unscreening them with
! the pseudo valence density rho, which differs from n_v there (about 1e-4
! of it for a scalar- or fully relativistic reference, whose pseudo states
! miss the relativistic terms of the tail), leaves V_Hxc[n_v] - V_Hxc[rho]
! in vpuns: a long-ranged term that is not ionic, and that plane-wave codes
! see through the G=0 term alpha_z and wherever they cut the table.  Here
! the screened potentials are switched over a window (ra,rb) to
!
!   V_ion + V_Hxc[rho + rhomod],
!   V_ion = -Z/r + V_H[n_c] + V_xc[n_c + n_v] - V_xc[n_v + rhomod],
!
! the frozen-core ionic potential (-zion/r plus the core penetration and,
! without a model core, the core-valence xc nonlinearity), so that vpuns
! is vp - V_Hxc[rho + rhomod] up to ra and V_ion beyond rb.  The pseudo
! valence states see the switched local potential, so rho is iterated
! with the tail to a fixed point: the reference configuration stays a
! fixed point of the pseudo-atom SCF.  The projectors vanish beyond rc
! and are unchanged.  V_ion is built from the densities, not as
! V_AE - V_Hxc[n_v + rhomod]: V_AE is the input potential of the last
! all-electron SCF step.

 subroutine ionic_tail(zz,zion,rhoc,rhotae,rhomod,rr,mmax,adams,iexc,c_light, &
&                      vion,ifail)

!! the frozen-core ionic potential V_ion of the all-electron atom
!!
!! @param zz  atomic number
!! @param zion  ion charge (zz minus the core occupation)
!! @param rhoc  core charge density
!! @param rhotae  all-electron valence charge density
!! @param rhomod  model core charge
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param iexc  exchange-correlation functional
!! @param c_light  speed of light (Hartree atomic units)
!! @param vion  V_ion
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp
 use vout_m, only: vout
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: mmax,iexc
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: zz,zion
 real(dp), intent(in) :: rhoc(mmax),rhotae(mmax),rhomod(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: vion(mmax)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: eeel,eexc
 real(dp), allocatable :: vt(:),vxct(:)

 allocate(vt(mmax),vxct(mmax))
 vxct(:)=0.0d0

! V_H[n_c], its charge fixed at the end of the mesh to zz-zion so that
! -Z/r + V_H[n_c] is -zion/r there
 call vout(.false.,rhoc,rhoc,vt,vxct,zz-zion,eeel,eexc,rr,mmax,adams,iexc, &
&          c_light,ifail)
 if(ifail/=0) return
 vion(:)=-zz/rr(:)+vt(:)-vxct(:)

 call vout(.true.,rhotae,rhoc,vt,vxct,zion,eeel,eexc,rr,mmax,adams,iexc, &
&          c_light,ifail)
 if(ifail/=0) return
 vion(:)=vion(:)+vxct(:)
 call vout(.true.,rhotae,rhomod,vt,vxct,zion,eeel,eexc,rr,mmax,adams,iexc, &
&          c_light,ifail)
 if(ifail/=0) return
 vion(:)=vion(:)-vxct(:)

 deallocate(vt,vxct)
 return
 end subroutine ionic_tail

 subroutine ionic_window(rr,mmax,zion,vion,rcmax,tol,ra,rb,ifail)

!! the window (ra,rb) of the switch: ra is the largest core radius, beyond
!! which the screened potentials are V_AE, and rb is where V_ion has become
!! -zion/r to tol, but at least 2*ra (a switch over less than ra would be
!! harder than the pseudization)
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param zion  ion charge
!! @param vion  V_ion
!! @param rcmax  largest core radius (bohr)
!! @param tol  tolerance of V_ion + zion/r
!! @param ra  start of the window
!! @param rb  end of the window
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp
 use logmesh_m, only: lm_root

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),vion(mmax)
 real(dp), intent(in) :: zion,rcmax,tol

!Output variables
 real(dp), intent(out) :: ra,rb
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: gg(:)
 real(dp) :: rcore
 integer :: ii,ierr

! rcore: the outermost radius where (V_ion + zion/r)**2 = tol**2, the root
! of the interpolant between the mesh points where it crosses (R-13)
 ifail=0
 allocate(gg(mmax))
 do ii=1,mmax
   gg(ii)=(vion(ii)+zion/rr(ii))**2-tol**2
 end do
 rcore=0.0d0
 do ii=mmax-1,1,-1
   if(gg(ii)>=0.0d0) then
     call lm_root(rr,gg,mmax,rr(ii),rr(ii+1),rcore,ierr,ifail)
     if(ifail/=0) return
     if(ierr/=0) rcore=rr(ii+1)
     exit
   end if
 end do
 deallocate(gg)
 ra=rcmax
! a C-infinity maximum, within 0.05 bohr of the larger
 rb=0.5d0*(rcore+2.0d0*ra+sqrt((rcore-2.0d0*ra)**2+0.01d0))
 return
 end subroutine ionic_window

 subroutine ionic_valence(nkap,nc,nv,la,fa,ea,lloc,nproj,rr,vp,vlocpw, &
&                         vkbpw,evkb,mmax,adams,rho,rhops,uupsa,eeig,iout,ifail)

!! the pseudo valence states of the local potential vp(:,lloc+1,:) and the
!! VKB projectors, as oncvpsp_main solves them: rho, the state densities
!! rhops, the orbitals uupsa and the eigenvalue sum eeig
!!
!! @param nkap  number of kappas per l: 1, or 2 for irel=2 (both kappa)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular momenta
!! @param fa  occupancies
!! @param ea  all-electron eigenvalues (the starting guesses)
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param nproj  number of projectors for each l
!! @param rr  log radial mesh (bohr)
!! @param vp  screened semi-local and local pseudopotentials
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param rho  pseudo valence charge
!! @param rhops  state-by-state pseudo valence charge (one electron)
!! @param uupsa  pseudo valence wave functions
!! @param eeig  sum of the occupied eigenvalues
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_construction, fail_convergence, ierr_ghost, failed, &
&                    msg_ionic_valence_construction, msg_ionic_valence_conv
 use constants_m, only: dp, mxstate, mxl, mxlin, mxprj
 use kappa_m, only: kap_count, kap_j_share
 use pieces_m, only: piecewise_t
 use lschvkbb_m, only: lschvkbb
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: nkap,nc,nv,lloc,mmax
 type(adams_t), intent(in) :: adams
 integer, intent(in) :: la(mxstate),nproj(mxlin+1)
 real(dp), intent(in) :: fa(mxstate),ea(mxstate,nkap)
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: rho(mmax),rhops(mmax,nv),uupsa(mmax,nkap,nv)
 real(dp), intent(out) :: eeig
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: uu(:),up(:)
 real(dp) :: et,emin,emax,fj,de
 integer :: kk,ll,l1,ikap,mkap,mch,ierr
 integer :: nodes(mxl+1)

 ifail=0
 allocate(uu(mmax),up(mmax))
 if(nkap==1) then
   de=0.25d0
 else
   de=0.10d0
 end if
 rho(:)=0.0d0
 rhops(:,:)=0.0d0
 uupsa(:,:,:)=0.0d0
 eeig=0.0d0
 nodes(:)=0
 do kk=1,nv
   ll=la(nc+kk)
   l1=ll+1
   mkap=kap_count(nkap,ll)
   do ikap=1,mkap
     fj=1.0d0
     if(mkap==2) fj=kap_j_share(ll,ikap,1.0d0)
     et=ea(nc+kk,ikap)
     emax=(1.0d0-de)*et
     emin=(1.0d0+de)*et
     call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
&                  rr,vp(1,lloc+1,ikap),vlocpw,vkbpw(1:nproj(l1),l1,ikap), &
&                  evkb(1,l1,ikap), &
&                  uu,up,mmax,adams,mch,iout,ifail)
     if(ifail/=0) return
     if(ierr/=0) then
       write(iout,'(/a,3i0)') 'ionic_valence: ERROR lschvkbb n,l,ierr=', &
&            ll+nodes(l1)+1,ll,ierr
       if(ierr==ierr_ghost) write(iout,'(a,i2,a,f12.6,a/a)') &
&           'a state of l=',ll,' lies below',emin, &
&           ' Ha that the node count does not expect:', &
&           'a bound ghost of the projectors'
       if(ierr==ierr_ghost) then
         ifail=failed(fail_construction,msg_ionic_valence_construction,ierr)
       else
         ifail=failed(fail_convergence,msg_ionic_valence_conv,ierr)
       end if
       return
     end if
     uupsa(:,ikap,kk)=uu(:)
     rhops(:,kk)=rhops(:,kk)+fj*(uu(:)/rr(:))**2
     eeig=eeig+fa(nc+kk)*fj*et
   end do
   rho(:)=rho(:)+fa(nc+kk)*rhops(:,kk)
   nodes(l1)=nodes(l1)+1
 end do
 deallocate(uu,up)
 return
 end subroutine ionic_valence

 subroutine ionic_screening(nkap,nc,nv,la,fa,ea,lmax,lloc,nproj,rr,vp,vlocpw, &
&                           vkbpw,evkb,mmax,adams,zz,zion,zval,iexc,c_light, &
&                           rcmax,rhoc,rhotae,rhomod,rho,rhops,uupsa,eeig,vo, &
&                           vxc,eeel,eexc,vion,iout,ifail)

!! switches the screened potentials vp to V_ion + V_Hxc[rho + rhomod] over
!! the ionic window, iterating the pseudo valence states (rho, rhops,
!! uupsa, eeig) and their screening (vo, vxc, eeel, eexc) to a fixed point:
!! until the local potential, weighted by rho, changes by no more than
!! 1e-13 Ha.  The last states are those of the final vp
!!
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular momenta
!! @param fa  occupancies
!! @param ea  all-electron eigenvalues
!! @param lmax  maximum angular momentum
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param nproj  number of projectors for each l
!! @param rr  log radial mesh (bohr)
!! @param vp  screened semi-local and local pseudopotentials
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param zz  atomic number
!! @param zion  ion charge
!! @param zval  valence charge
!! @param iexc  exchange-correlation functional
!! @param c_light  speed of light (Hartree atomic units)
!! @param rcmax  largest core radius (bohr)
!! @param rhoc  core charge density
!! @param rhotae  all-electron valence charge density
!! @param rhomod  model core charge
!! @param rho  pseudo valence charge
!! @param rhops  state-by-state pseudo valence charge
!! @param uupsa  pseudo valence wave functions
!! @param eeig  sum of the occupied eigenvalues
!! @param vo  Hartree and exchange-correlation potential of rho + rhomod
!! @param vxc  exchange-correlation potential (Ha) of rho + rhomod
!! @param eeel  electron-electron interaction energy
!! @param eexc  exchange-correlation correction to the eigenvalue sum
!! @param vion  V_ion
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_convergence, failed, msg_ionic_screening_conv
 use constants_m, only: dp, mxstate, mxl, mxlin, mxprj
 use pieces_m, only: piecewise_t
 use cinf_taper_m, only: cinf_taper
 use vout_m, only: vout
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: nkap,nc,nv,lmax,lloc,mmax,iexc
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 integer, intent(in) :: la(mxstate),nproj(mxlin+1)
 real(dp), intent(in) :: fa(mxstate),ea(mxstate,nkap)
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 real(dp), intent(in) :: zz,zion,zval,rcmax
 real(dp), intent(in) :: rhoc(mmax),rhotae(mmax),rhomod(mmax)
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: vp(mmax,mxl+2,nkap)
 real(dp), intent(inout) :: rho(mmax),rhops(mmax,nv),uupsa(mmax,nkap,nv)
 real(dp), intent(inout) :: eeig,vo(mmax),vxc(mmax),eeel,eexc

!Output variables
 real(dp), intent(out) :: vion(mmax)
 integer, intent(out) :: ifail

!Local variables
 real(dp), parameter :: tol=1.0d-10
! the fixed point: it contracts 7-10 times a step to a floor of 1e-17 to
! 5e-15 Ha (the reference inputs), so 1e-13 leaves the result within
! 2e-14 of it, whatever the step count (doc/smoothness_plan.md, S-6)
 real(dp), parameter :: tolfix=1.0d-13
 integer, parameter :: maxiter=50
 real(dp), allocatable :: ww(:),vp0(:,:,:),vnew(:),vnext(:,:),vsw(:)
 real(dp) :: ra,rb,change,sn,sd
 integer :: ii,l1,ikap,iter

 ifail=0
 allocate(ww(mmax),vp0(mmax,mxl+2,nkap),vnew(mmax),vnext(mmax,nkap),vsw(mmax))

 call ionic_tail(zz,zion,rhoc,rhotae,rhomod,rr,mmax,adams,iexc,c_light,vion, &
&                ifail)
 if(ifail/=0) return
 call ionic_window(rr,mmax,zion,vion,rcmax,tol,ra,rb,ifail)
 if(ifail/=0) return
 do ii=1,mmax
   ww(ii)=cinf_taper(rr(ii),ra,rb)
 end do
 vp0(:,:,:)=vp(:,:,:)

 do ikap=1,nkap
   vp(:,lloc+1,ikap)=ww(:)*vp0(:,lloc+1,ikap)+(1.0d0-ww(:))*(vion(:)+vo(:))
 end do
 do iter=1,maxiter
   call ionic_valence(nkap,nc,nv,la,fa,ea,lloc,nproj,rr,vp,vlocpw,vkbpw, &
&                     evkb,mmax,adams,rho,rhops,uupsa,eeig,iout,ifail)
   if(ifail/=0) return
   call vout(.true.,rho,rhomod,vo,vxc,zval,eeel,eexc,rr,mmax,adams,iexc, &
&            c_light,ifail)
   if(ifail/=0) return
! change of the local potential where the valence electrons are; once it
! is below tolfix, vp is kept as the potential whose states rho, eeig and
! (through rho) vo are, so that the pseudo-atom SCF of the reference
! configuration starts at its own fixed point
   sn=0.0d0
   sd=0.0d0
   do ikap=1,nkap
     vnew(:)=ww(:)*vp0(:,lloc+1,ikap)+(1.0d0-ww(:))*(vion(:)+vo(:))
     sn=sn+sum(rho(:)*rr(:)**3*(vnew(:)-vp(:,lloc+1,ikap))**2)
     sd=sd+sum(rho(:)*rr(:)**3)
     vnext(:,ikap)=vnew(:)
   end do
   change=sqrt(sn/sd)
   if(change<=tolfix) exit
   if(iter==maxiter) then
     write(iout,'(/a,i0,a,es12.4)') &
&          'ionic_screening: ERROR the ionic tail did not converge in', &
&          maxiter,' steps, change =',change
     ifail=failed(fail_convergence,msg_ionic_screening_conv)
     return
   end if
   vp(:,lloc+1,:)=vnext(:,:)
 end do

! the semi-local potentials, switched with the screening their local
! potential was switched with
 do ikap=1,nkap
   vsw(:)=vp(:,lloc+1,ikap)-ww(:)*vp0(:,lloc+1,ikap)
   do l1=1,max(lmax+1,lloc+1)
     if(l1==lloc+1) cycle
     vp(:,l1,ikap)=ww(:)*vp0(:,l1,ikap)+vsw(:)
   end do
 end do

 write(iout,'(/a,f8.4,a,f8.4,a,i3,a)') 'Ionic tail switched over r =',ra, &
&      ' to',rb,' a.u., converged in',iter,' steps'

 deallocate(ww,vp0,vnew,vnext,vsw)
 return
 end subroutine ionic_screening

 subroutine coulomb_tail(vloc,vion,rr,mmax,zion,rcmax,iout,ifail)

!! reports how closely the unscreened local potential vloc is -zion/r
!! beyond rcmax, as plane-wave codes see it (coulomb_tail of
!! PseudoPotentialGenerators.jl): alpha_z = 4 pi int r**2 (vloc + zion/r)
!! dr over the whole table (the G=0 term of the energy, Ha bohr**3), its
!! change when the table is cut at 6 and 10 bohr, the radius beyond which
!! |vloc + zion/r| < 1e-10, and the largest |vloc - V_ion| beyond rcmax
!! (the noise floor of the unscreening, below which that radius means
!! nothing). The integrals are lm_integral's, cut at the radii exactly,
!! and the radius is the root of the interpolant of (vloc + zion/r)**2 -
!! 1e-10**2 past the last mesh point above it (lm_root), so all are
!! smooth functions of the inputs (doc/smoothness_plan.md, X-6)
!!
!! @param vloc  unscreened local potential
!! @param vion  V_ion
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param zion  ion charge
!! @param rcmax  largest core radius (bohr)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp, pi
 use logmesh_m, only: lm_integral, lm_root

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: vloc(mmax),vion(mmax),rr(mmax)
 real(dp), intent(in) :: zion,rcmax
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 real(dp), parameter :: tol=1.0d-10
 real(dp), parameter :: radii(2)=[6.0d0,10.0d0]
 real(dp) :: alpha,cut,rcoul,dvmax,dv
 real(dp), allocatable :: ff(:)
 integer :: ii,jj,ierr

 ifail=0
! r**2 (vloc + zion/r), which goes as r at the origin
 allocate(ff(mmax))
 ff(:)=rr(:)**2*vloc(:)+zion*rr(:)
 call lm_integral(rr,ff,mmax,1,0.0d0,rr(mmax),alpha)
 alpha=4.0d0*pi*alpha

 rcoul=rr(1)
 do ii=mmax,1,-1
   if(abs(vloc(ii)+zion/rr(ii))>=tol) then
     rcoul=rr(min(ii+1,mmax))
     if(ii<mmax) then
       ff(:)=(vloc(:)+zion/rr(:))**2-tol**2
       call lm_root(rr,ff,mmax,rr(ii),rr(ii+1),rcoul,ierr,ifail)
       if(ifail/=0) return
       if(ierr/=0) rcoul=rr(ii+1)
     end if
     exit
   end if
 end do

 dvmax=0.0d0
 do ii=1,mmax
   if(rr(ii)<rcmax) cycle
   dv=abs(vloc(ii)-vion(ii))
   dvmax=max(dvmax,dv)
 end do

 write(iout,'(/a)') 'Coulomb tail of the local potential'
 write(iout,'(a,1p,e18.10,a)') &
&      '  alpha_z = 4 pi int r**2 (V_loc + zion/r) dr', &
&      alpha,' Ha bohr**3'
 ff(:)=rr(:)**2*vloc(:)+zion*rr(:)
 do jj=1,2
   call lm_integral(rr,ff,mmax,1,min(radii(jj),rr(mmax)),rr(mmax),cut)
   write(iout,'(a,f6.2,a,1p,e12.4)') '  alpha_z change, table cut at r =', &
&        radii(jj),':',-4.0d0*pi*cut
 end do
 write(iout,'(a,f10.4)') '  |V_loc + zion/r| < 1e-10 beyond r =',rcoul
 write(iout,'(a,1p,e12.4)') '  max |V_loc - V_ion| beyond rcmax',dvmax
 deallocate(ff)

 return
 end subroutine coulomb_tail

 end module ionic_tail_m
