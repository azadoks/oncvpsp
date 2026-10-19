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
! The ionic tail of the pseudopotentials, after PseudoPotentialGenerators.jl
! (src/unscreening.jl: ionic_tail, ionic_switch, ionic_tail_window,
! with_ionic_tail; src/verification.jl: coulomb_tail).
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

 subroutine ionic_tail(zz,zion,rhoc,rhotae,rhomod,rr,mmax,iexc,vion)

!! the frozen-core ionic potential V_ion of the all-electron atom
!!
!! @param zz  nuclear charge
!! @param zion  ionic charge (zz minus the core occupation)
!! @param rhoc  all-electron core charge
!! @param rhotae  all-electron valence charge
!! @param rhomod  model core charge
!! @param rr  log radial grid
!! @param mmax  dimension of log grid
!! @param iexc  exchange-correlation function to be used
!! @param vion  V_ion

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax,iexc
 real(dp), intent(in) :: zz,zion
 real(dp), intent(in) :: rhoc(mmax),rhotae(mmax),rhomod(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: vion(mmax)

!Local variables
 real(dp) :: eeel,eexc
 real(dp), allocatable :: vt(:),vxct(:)

 allocate(vt(mmax),vxct(mmax))
 vxct(:)=0.0d0

! V_H[n_c], its charge fixed at the end of the mesh to zz-zion so that
! -Z/r + V_H[n_c] is -zion/r there
 call vout(0,rhoc,rhoc,vt,vxct,zz-zion,eeel,eexc,rr,mmax,iexc)
 vion(:)=-zz/rr(:)+vt(:)-vxct(:)

 call vout(1,rhotae,rhoc,vt,vxct,zion,eeel,eexc,rr,mmax,iexc)
 vion(:)=vion(:)+vxct(:)
 call vout(1,rhotae,rhomod,vt,vxct,zion,eeel,eexc,rr,mmax,iexc)
 vion(:)=vion(:)-vxct(:)

 deallocate(vt,vxct)
 return
 end subroutine ionic_tail

 subroutine ionic_window(rr,mmax,zion,vion,rcmax,tol,ra,rb)

!! the window (ra,rb) of the switch: ra is the largest core radius, beyond
!! which the screened potentials are V_AE, and rb is where V_ion has become
!! -zion/r to tol, but at least 2*ra (a switch over less than ra would be
!! harder than the pseudization)
!!
!! @param rr  log radial grid
!! @param mmax  dimension of log grid
!! @param zion  ionic charge
!! @param vion  V_ion
!! @param rcmax  largest core radius
!! @param tol  tolerance of V_ion + zion/r
!! @param ra  start of the window
!! @param rb  end of the window

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),vion(mmax)
 real(dp), intent(in) :: zion,rcmax,tol

!Output variables
 real(dp), intent(out) :: ra,rb

!Local variables
 real(dp) :: rcore
 integer :: ii

 rcore=0.0d0
 do ii=mmax,1,-1
   if(abs(vion(ii)+zion/rr(ii))>=tol) then
     rcore=rr(min(ii+1,mmax))
     exit
   end if
 end do
 ra=rcmax
 rb=max(rcore,2.0d0*ra)
 return
 end subroutine ionic_window

 real(dp) function ionic_switch(rad,ra,rb)

!! C-infinity switch from 1 (r <= ra) to 0 (r >= rb): f(1-x)/(f(1-x)+f(x)),
!! f(t) = exp(-1/t), x = (r-ra)/(rb-ra); its transform decays faster than
!! any power of q, so it adds no plane-wave tail of its own
!!
!! @param rad  radius
!! @param ra  start of the window
!! @param rb  end of the window
!! @return the switch at rad

 use constants_m, only: dp
 implicit none

!Input variables
 real(dp), intent(in) :: rad,ra,rb

!Local variables
 real(dp) :: xx,f0,f1

 if(rad<=ra) then
   ionic_switch=1.0d0
 else if(rad>=rb) then
   ionic_switch=0.0d0
 else
   xx=(rad-ra)/(rb-ra)
   f0=exp(-1.0d0/(1.0d0-xx))
   f1=exp(-1.0d0/xx)
   ionic_switch=f0/(f0+f1)
 end if
 return
 end function ionic_switch

 subroutine ionic_valence(nkap,nc,nv,la,fa,ea,lloc,nproj,rr,vp,vkb,evkb, &
&                         mmax,mxprj,rho,rhops,uupsa,eeig)

!! the pseudo valence states of the local potential vp(:,lloc+1,:) and the
!! VKB projectors, as oncvpsp_main solves them: rho, the state densities
!! rhops, the orbitals uupsa and the eigenvalue sum eeig
!!
!! @param nkap  1 for scalar-relativistic or non-relativistic, 2 for
!!              relativistic (both kappa)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular momenta
!! @param fa  occupancies
!! @param ea  all-electron eigenvalues (the starting guesses)
!! @param lloc  index of local potential
!! @param nproj  number of projectors for each l
!! @param rr  log radial grid
!! @param vp  screened semi-local and local pseudopotentials
!! @param vkb  Vanderbilt-Kleinman-Bylander projectors
!! @param evkb  coefficients of the projectors
!! @param mmax  dimension of log grid
!! @param mxprj  dimension of number of projector
!! @param rho  pseudo valence charge
!! @param rhops  state-by-state pseudo valence charge (one electron)
!! @param uupsa  pseudo valence wave functions
!! @param eeig  sum of the occupied eigenvalues

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: nkap,nc,nv,lloc,mmax,mxprj
 integer, intent(in) :: la(30),nproj(6)
 real(dp), intent(in) :: fa(30),ea(30,2)
 real(dp), intent(in) :: rr(mmax),vp(mmax,5,nkap)
 real(dp), intent(in) :: vkb(mmax,mxprj,4,nkap),evkb(mxprj,4,nkap)

!Output variables
 real(dp), intent(out) :: rho(mmax),rhops(mmax,nv),uupsa(mmax,nkap,nv)
 real(dp), intent(out) :: eeig

!Local variables
 real(dp), allocatable :: uu(:),up(:)
 real(dp) :: et,emin,emax,fj,de
 integer :: kk,ll,l1,ikap,mkap,mch,ierr
 integer :: nodes(4)

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
   mkap=nkap
   if(ll==0) mkap=1
   do ikap=1,mkap
     if(mkap==1) then
       fj=1.0d0
     else if(ikap==1) then
       fj=dble(ll+1)/dble(2*ll+1)
     else
       fj=dble(ll)/dble(2*ll+1)
     end if
     et=ea(nc+kk,ikap)
     emax=(1.0d0-de)*et
     emin=(1.0d0+de)*et
     call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
&                  rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                  uu,up,mmax,mch)
     if(ierr/=0) then
       write(6,'(/a,3i4)') 'ionic_valence: ERROR lschvkbb n,l,ierr=', &
&            ll+nodes(l1)+1,ll,ierr
       error stop 1
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

 subroutine ionic_screening(nkap,nc,nv,la,fa,ea,lmax,lloc,nproj,rr,vp,vkb, &
&                           evkb,mmax,mxprj,zz,zion,zval,iexc,rcmax,rhoc, &
&                           rhotae,rhomod,rho,rhops,uupsa,eeig,vo,vxc, &
&                           eeel,eexc,vion)

!! switches the screened potentials vp to V_ion + V_Hxc[rho + rhomod] over
!! the ionic window, iterating the pseudo valence states (rho, rhops,
!! uupsa, eeig) and their screening (vo, vxc, eeel, eexc) to a fixed point:
!! until the local potential, weighted by rho, changes by no more than
!! 1e-10.  The last states are those of the final vp
!!
!! @param nkap  1 for scalar-relativistic or non-relativistic, 2 for
!!              relativistic
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular momenta
!! @param fa  occupancies
!! @param ea  all-electron eigenvalues
!! @param lmax  maximum angular momentum
!! @param lloc  index of local potential
!! @param nproj  number of projectors for each l
!! @param rr  log radial grid
!! @param vp  screened semi-local and local pseudopotentials
!! @param vkb  Vanderbilt-Kleinman-Bylander projectors
!! @param evkb  coefficients of the projectors
!! @param mmax  dimension of log grid
!! @param mxprj  dimension of number of projector
!! @param zz  nuclear charge
!! @param zion  ionic charge
!! @param zval  valence charge
!! @param iexc  exchange-correlation function to be used
!! @param rcmax  largest core radius
!! @param rhoc  all-electron core charge
!! @param rhotae  all-electron valence charge
!! @param rhomod  model core charge
!! @param rho  pseudo valence charge
!! @param rhops  state-by-state pseudo valence charge
!! @param uupsa  pseudo valence wave functions
!! @param eeig  sum of the occupied eigenvalues
!! @param vo  Hartree and exchange-correlation potential of rho + rhomod
!! @param vxc  exchange-correlation potential of rho + rhomod
!! @param eeel  electron-electron interaction energy
!! @param eexc  exchange-correlation correction to the eigenvalue sum
!! @param vion  V_ion

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: nkap,nc,nv,lmax,lloc,mmax,mxprj,iexc
 integer, intent(in) :: la(30),nproj(6)
 real(dp), intent(in) :: fa(30),ea(30,2)
 real(dp), intent(in) :: rr(mmax),vkb(mmax,mxprj,4,nkap),evkb(mxprj,4,nkap)
 real(dp), intent(in) :: zz,zion,zval,rcmax
 real(dp), intent(in) :: rhoc(mmax),rhotae(mmax),rhomod(mmax)

!Input/output variables
 real(dp), intent(inout) :: vp(mmax,5,nkap)
 real(dp), intent(inout) :: rho(mmax),rhops(mmax,nv),uupsa(mmax,nkap,nv)
 real(dp), intent(inout) :: eeig,vo(mmax),vxc(mmax),eeel,eexc

!Output variables
 real(dp), intent(out) :: vion(mmax)

!Local variables
 real(dp), parameter :: tol=1.0d-10
 integer, parameter :: maxiter=50
 real(dp), allocatable :: ww(:),vp0(:,:,:),vnew(:),vnext(:,:),vsw(:)
 real(dp) :: ra,rb,change,sn,sd
 real(dp), external :: ionic_switch
 integer :: ii,l1,ik,iter

 allocate(ww(mmax),vp0(mmax,5,nkap),vnew(mmax),vnext(mmax,nkap),vsw(mmax))

 call ionic_tail(zz,zion,rhoc,rhotae,rhomod,rr,mmax,iexc,vion)
 call ionic_window(rr,mmax,zion,vion,rcmax,tol,ra,rb)
 do ii=1,mmax
   ww(ii)=ionic_switch(rr(ii),ra,rb)
 end do
 vp0(:,:,:)=vp(:,:,:)

 do ik=1,nkap
   vp(:,lloc+1,ik)=ww(:)*vp0(:,lloc+1,ik)+(1.0d0-ww(:))*(vion(:)+vo(:))
 end do
 do iter=1,maxiter
   call ionic_valence(nkap,nc,nv,la,fa,ea,lloc,nproj,rr,vp,vkb,evkb, &
&                     mmax,mxprj,rho,rhops,uupsa,eeig)
   call vout(1,rho,rhomod,vo,vxc,zval,eeel,eexc,rr,mmax,iexc)
! change of the local potential where the valence electrons are; once it
! is below tol, vp is kept as the potential whose states rho, eeig and
! (through rho) vo are, so that the pseudo-atom SCF of the reference
! configuration starts at its own fixed point
   sn=0.0d0
   sd=0.0d0
   do ik=1,nkap
     vnew(:)=ww(:)*vp0(:,lloc+1,ik)+(1.0d0-ww(:))*(vion(:)+vo(:))
     sn=sn+sum(rho(:)*rr(:)**3*(vnew(:)-vp(:,lloc+1,ik))**2)
     sd=sd+sum(rho(:)*rr(:)**3)
     vnext(:,ik)=vnew(:)
   end do
   change=sqrt(sn/sd)
   if(change<=tol) exit
   if(iter==maxiter) then
     write(6,'(/a,i4,a,1p,e10.2)') &
&          'ionic_screening: ERROR the ionic tail did not converge in', &
&          maxiter,' steps, change =',change
     error stop 1
   end if
   vp(:,lloc+1,:)=vnext(:,:)
 end do

! the semi-local potentials, switched with the screening their local
! potential was switched with
 do ik=1,nkap
   vsw(:)=vp(:,lloc+1,ik)-ww(:)*vp0(:,lloc+1,ik)
   do l1=1,max(lmax+1,lloc+1)
     if(l1==lloc+1) cycle
     vp(:,l1,ik)=ww(:)*vp0(:,l1,ik)+vsw(:)
   end do
 end do

 write(6,'(/a,f8.4,a,f8.4,a,i3,a)') 'Ionic tail switched over r =',ra, &
&      ' to',rb,' a.u., converged in',iter,' steps'

 deallocate(ww,vp0,vnew,vnext,vsw)
 return
 end subroutine ionic_screening

 subroutine coulomb_tail(vloc,vion,rr,mmax,zion,rcmax)

!! reports how closely the unscreened local potential vloc is -zion/r
!! beyond rcmax, as plane-wave codes see it (coulomb_tail of
!! PseudoPotentialGenerators.jl): alpha_z = 4 pi int r**2 (vloc + zion/r)
!! dr over the whole table (the G=0 term of the energy, Ha bohr**3), its
!! change when the table is cut at 6 and 10 bohr, the radius beyond which
!! |vloc + zion/r| < 1e-10, and the largest |vloc - V_ion| beyond rcmax
!!
!! @param vloc  unscreened local potential
!! @param vion  V_ion
!! @param rr  log radial grid
!! @param mmax  dimension of log grid
!! @param zion  ionic charge
!! @param rcmax  largest core radius

 use constants_m, only: dp, pi
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: vloc(mmax),vion(mmax),rr(mmax)
 real(dp), intent(in) :: zion,rcmax

!Local variables
 real(dp), parameter :: tol=1.0d-10
 real(dp), parameter :: radii(2)=[6.0d0,10.0d0]
 real(dp) :: al,alpha,cut,rcoul,dvmax,dv
 integer :: ii,jj

 al = 0.01d0 * log(rr(101) / rr(1))

 alpha=0.0d0
 do ii=1,mmax
   alpha=alpha+rr(ii)**3*(vloc(ii)+zion/rr(ii))
 end do
 alpha=4.0d0*pi*al*alpha

 rcoul=rr(1)
 do ii=mmax,1,-1
   if(abs(vloc(ii)+zion/rr(ii))>=tol) then
     rcoul=rr(min(ii+1,mmax))
     exit
   end if
 end do

 dvmax=0.0d0
 do ii=1,mmax
   if(rr(ii)<rcmax) cycle
   dv=abs(vloc(ii)-vion(ii))
   dvmax=max(dvmax,dv)
 end do

 write(6,'(/a)') 'Coulomb tail of the local potential'
 write(6,'(a,1p,e18.10,a)') '  alpha_z = 4 pi int r**2 (V_loc + zion/r) dr', &
&      alpha,' Ha bohr**3'
 do jj=1,2
   cut=0.0d0
   do ii=1,mmax
     if(rr(ii)>radii(jj)) exit
     cut=cut+rr(ii)**3*(vloc(ii)+zion/rr(ii))
   end do
   cut=4.0d0*pi*al*cut
   write(6,'(a,f6.2,a,1p,e12.4)') '  alpha_z change, table cut at r =', &
&        radii(jj),':',cut-alpha
 end do
 write(6,'(a,f10.4)') '  |V_loc + zion/r| < 1e-10 beyond r =',rcoul
 write(6,'(a,1p,e12.4)') '  max |V_loc - V_ion| beyond rcmax',dvmax

 return
 end subroutine coulomb_tail
