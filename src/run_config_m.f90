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

 module run_config_m

 implicit none
 private
 public :: run_config_ps

 contains

 subroutine run_config_ps(jj,aec,nc,rhov,rhomod,rr,rcmax,mmax,mxprj,iexc, &
&                         etot,epstot,nproj,vpuns,lloc,vlocpw,vkbpw,evkb, &
&                         irel,c_light,nkap,ntc,itc,etc,dtc,iout,ifail)

!! self-consistent pseudoatom calculation for test configuration jj, compared
!! with the all-electron results for the reference and test configurations
!!
!! @param jj  index of current configuration
!! @param aec  its all-electron half (run_config_ae)
!! @param nc  number of core states
!! @param rhov  valence pseudo-charge of reference configuration
!! @param rhomod  model core charge
!! @param rr  log radial mesh
!! @param rcmax  maximum core radius for psp
!! @param mmax  size of log grid
!! @param mxprj  dimension of number of projectors
!! @param iexc  exchange-correlation function to be used
!! @param etot  reference configuration total energy
!! @param epstot  reference configuration pseudoatom total energy
!! @param nproj  number of VKB projectors to use for each l
!! @param vpuns  unscreened semi-local pseudopotentials (plus different vloc if
!!               lloc==4)
!! @param lloc  index-1 of local potential
!! @param vlocpw  the screened local potential's pieces
!! @param vkbpw  Vanderbilt-Kleinman-Bylander projectors as functions of
!!               pieces, last index 1 for j=l+1/2, 2 for j=l-1/2 when irel=2
!! @param evkb  VKB projector coefficients, last index as for vkbpw
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!! @param ntc  rows of the table of eigenvalues, -1 if the configuration
!!             failed (output)
!! @param itc  n, l and ikap of each row, ikap 1 for j=l+1/2 (or an
!!             l-channel), 2 for j=l-1/2 (output)
!! @param etc  occupation, all-electron and pseudo eigenvalue (Ha) of each
!!             row (output)
!! @param dtc  the all-electron and pseudo total-energy differences from
!!             the reference configuration and their difference, the
!!             excitation error (Ha) (output)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use config_ae_m, only: config_ae_warning, ae_config_t
 use constants_m, only: dp, mxstate, mxl, mxlin, nrhomod, nlev
 use pieces_m, only: piecewise_t
 use kappa_m, only: kap_count, kap_of, j_share
 use psatom_m, only: psatom
 implicit none

!Input variables

 type(ae_config_t), intent(in) :: aec
 integer, intent(in) :: jj,mmax,mxprj,iexc,nc,lloc,irel,nkap
 real(dp), intent(in) :: c_light
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: etot,epstot,rcmax
 real(dp), intent(in) :: rhov(mmax),rr(mmax)
 real(dp), intent(in) :: vpuns(mmax,mxl+2)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 real(dp), intent(in) :: rhomod(mmax,nrhomod)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)

!Output variables
 integer, intent(out) :: ntc,itc(3,2*mxstate)
 real(dp), intent(out) :: etc(3,2*mxstate),dtc(3)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,it,kk,ll,l1,ierr,nvt
 integer :: ikap,kap,mkap
 integer :: nat(mxstate),lat(mxstate),natp(mxstate),latp(mxstate),nav(mxl+1)
 real(dp) :: eaetst,etsttot,fj
 real(dp) :: eat(mxstate,nkap,nlev),fat(mxstate,nlev),eatp(mxstate,nkap)
 real(dp) :: fatp(mxstate,nlev)

 real(dp),allocatable :: rho(:),rhocps(:),vfree(:)

 ifail=0
! atom tests comparing reference state and excited configuration for
! all-electron and pseudo atoms.  Total excitation energies and excited-
! -configuration eigenvalues are compared

! For multi-projector potentials, the bound-state solver lschvkbb has poor
! stability, and requires a trial energy sufficiently close to the correct
! result for a given local potential. This appears to be due to the fact
! that intermediate stages towards the solution of the radial Schroedinger
! equation with a non-local potential may not increase their node count
! monotonically with increasing energy (the same phenonenon which gives
! rise to ghost states).

! To aviod these problems, the pseudoatom calculation is started with the
! screened local potential and valence eigenvalues of the reference
! configuration.  This is adiabatically changed in two steps.  In the first,
! the occupation numbers are changed by 2% within the potential iteration
! loop for each of the first 50 steps to those of a maximally-ionized
! intermediate configuration.

! The second step is started with the local potential of this maximally-
! ionized psuedo-atom, and eigenvalues from an all-electron
! calculation of this configuration including unoccupied states that
! are to be adiabatically filled to produce the final configuration.
! The same 2% strategy is used filling pseudo-atom states in this step.

! The all-electron atoms for the maximally-ionized and excited
! configurations are computed by run_config_ae.

 ntc=-1
 itc(:,:)=0
 etc(:,:)=0.0d0
 dtc(:)=0.0d0
 if(aec%ist/=0) then
   call config_ae_warning(jj,irel,aec,iout)
   return
 end if

 nvt=aec%nvt
 nat(:)=aec%nat(:)
 lat(:)=aec%lat(:)
 fat(:,:)=aec%fat(:,:)
 eat(:,:,:)=aec%eat(:,:,:)
 eaetst=aec%etst

 allocate(rho(mmax),rhocps(mmax),vfree(mmax))
 vfree(:)=0.0d0

!first pseudoatom run from reference to maximally-ionized configuration
 do l1=1,4
   nav(l1)=l1-1
 end do

 do kk=1,nvt
   latp(kk)=lat(nc+kk)
   l1=latp(kk)+1
   nav(l1)=nav(l1)+1
   natp(kk)=nav(l1)
   eatp(kk,:)=eat(nc+kk,:,1)
   fatp(kk,1)=fat(nc+kk,1)
   fatp(kk,2)=fat(nc+kk,2)
 end do

 rho(:)=rhov(:)

 rhocps(:)=rhomod(:,1)

 call psatom(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
&           rr,rcmax,mmax,mxprj,iexc,c_light,etsttot,nproj,vpuns,lloc, &
&           vlocpw,vkbpw,evkb,ierr,nkap,vfree,iout,ifail)
 if(ifail/=0) return

 if(ierr/=0) then
   if(irel==2) then
     write(iout,'(a/a,i2)') &
&          'run_config_r: WARNING for fully non-local PS atom,', &
&          ' WARNING no output for configuration',jj
   else
     write(iout,'(a,a/a,i2)') &
&          'run_config: WARNING for fully non-local PS atom,', &
&          ' stg. 1', &
&          ' WARNING no output for configuration',jj
   end if
   deallocate(rho,rhocps,vfree)
   return
 end if

!second pseudoatom run from maximally-ionized to excited configuration

 do kk=1,nvt
   eatp(kk,:)=eat(nc+kk,:,2)
   fatp(kk,1)=fat(nc+kk,2)
   fatp(kk,2)=fat(nc+kk,3)
 end do

 call psatom(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
&           rr,rcmax,mmax,mxprj,iexc,c_light,etsttot,nproj,vpuns,lloc, &
&           vlocpw,vkbpw,evkb,ierr,nkap,vfree,iout,ifail)
 if(ifail/=0) return

 if(ierr/=0) then
   if(irel==2) then
     write(iout,'(a/a,i2)') &
&          'run_config_r: WARNING for fully non-local PS atom,', &
&          ' WARNING no output for configuration',jj
   else
     write(iout,'(a,a/a,i2)') &
&          'run_config: WARNING for fully non-local PS atom,', &
&          ' stg. 2', &
&          ' WARNING no output for configuration',jj
   end if
   deallocate(rho,rhocps,vfree)
   return
 end if



 if(irel==2) then
   call write_results_r
 else
   call write_results
 end if

 deallocate(rho,rhocps,vfree)

 return

 contains

 subroutine write_results
!! the eigenvalues and energies of the all-electron and pseudo atoms
 write(iout,'(/a)') '   n   l     f        eae           eps        diff'
 ntc=0
 do ii=1,nc
   write(iout,'(2i4,f8.4,f14.8)') nat(ii),lat(ii),fat(ii,3),eat(ii,1,3)
 end do
 do ii=1,nvt
   if(fat(ii+nc,3)==0.0d0) cycle
   write(iout,'(2i4,f8.4,2f14.8,1p,d12.2)') nat(ii+nc),lat(ii+nc), &
&        fat(ii+nc,3), &
&   eat(ii+nc,1,3),eatp(ii,1),eatp(ii,1)-eat(ii+nc,1,3)
   ntc=ntc+1
   itc(:,ntc)=[nat(ii+nc),lat(ii+nc),1]
   etc(:,ntc)=[fat(ii+nc,3),eat(ii+nc,1,3),eatp(ii,1)]
 end do
 dtc(:)=[eaetst-etot,etsttot-epstot,eaetst-etot-etsttot+epstot]

 write(iout,'(/a)') '    Total energies and differences'
 write(iout,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      AE_ref=',etot, &
& '  AE_tst=',eaetst,'  dif=',eaetst-etot
 write(iout,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      PS_ref=',epstot, &
& '  PS_tst=',etsttot,'  dif=',etsttot-epstot
 write(iout,'(a,1p,d10.2)') '      PSP excitation error=', &
& eaetst-etot-etsttot+epstot
 return
 end subroutine write_results

 subroutine write_results_r
!! write_results for irel=2, a line for each kappa with its share of the
!! occupation
 write(iout,'(/a)') '   n   l kap    f        eae           eps        diff'
 ntc=0
 do ii=1,nc
   ll = lat(ii)
   mkap=kap_count(2,ll)
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     kap=kap_of(ll,ikap)
     fj=j_share(ll,ikap,fat(ii,3))
     write(iout,'(3i4,f8.4,f14.8)') nat(ii),lat(ii),kap,fj,eat(ii,ikap,3)
   end do
 end do

 do ii=1,nvt
   ll = lat(ii+nc)
   mkap=kap_count(2,ll)
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     kap=kap_of(ll,ikap)
     fj=j_share(ll,ikap,fat(ii+nc,3))
     write(iout,'(3i4,f8.4,2f14.8,1p,e12.2)') nat(ii+nc),lat(ii+nc),kap,fj, &
&     eat(ii+nc,ikap,3),eatp(ii,ikap),eatp(ii,ikap)-eat(ii+nc,ikap,3)
     ntc=ntc+1
     itc(:,ntc)=[nat(ii+nc),lat(ii+nc),ikap]
     etc(:,ntc)=[fj,eat(ii+nc,ikap,3),eatp(ii,ikap)]
   end do
 end do
 dtc(:)=[eaetst-etot,etsttot-epstot,eaetst-etot-etsttot+epstot]

 write(iout,'(/a)') '    Total energies and differences'
 write(iout,'(a,1p,e16.8,a,e16.8,a,e10.2)') '      AE_ref=',etot, &
& '  AE_tst=',eaetst,'  dif=',eaetst-etot
 write(iout,'(a,1p,e16.8,a,e16.8,a,e10.2)') '      PS_ref=',epstot, &
& '  PS_tst=',etsttot,'  dif=',etsttot-epstot
 write(iout,'(a,1p,e10.2)') '      PSP excitation error=', &
& eaetst-etot-etsttot+epstot
 return
 end subroutine write_results_r

 end subroutine run_config_ps

 end module run_config_m
