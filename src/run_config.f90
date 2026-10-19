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
 subroutine run_config_ps(jj,istcnf,ierrcnf,itcnf,kerrcnf,nvtcnf, &
&                         natcnf,latcnf,fatcnf,eatcnf,etstcnf,nc,rhov,rhomod, &
&                         rr,rcmax,mmax,mxprj,iexc,etot,epstot,nproj,vpuns, &
&                         lloc,vkb,evkb,irel,nkap)

!! self-consistent pseudoatom calculation for test configuration jj, compared
!! with the all-electron results for the reference and test configurations
!!
!! @param jj  index of current configuration
!! @param istcnf..etstcnf  all-electron results, all configurations, from
!!                         run_config_ae
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
!! @param vkb  Vanderbilt-Kleinman-Bylander projectors, last index 1 for
!!             j=l+1/2, 2 for j=l-1/2 when irel=2
!! @param evkb  VKB projector coefficients, last index as for vkb
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappa per l: 1, or 2 for irel=2

 use config_ae_m, only: config_ae_warning
 use constants_m, only: dp
 implicit none

!Input variables

 integer, intent(in) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(in) :: nvtcnf(5),natcnf(30,5),latcnf(30,5)
 real(dp), intent(in) :: fatcnf(30,3,5),eatcnf(30,nkap,3,5),etstcnf(5)
 integer, intent(in) :: jj,mmax,mxprj,iexc,nc,lloc,irel,nkap
 integer, intent(in) :: nproj(5)
 real(dp), intent(in) :: etot,epstot,rcmax
 real(dp), intent(in) :: rhov(mmax),rr(mmax)
 real(dp), intent(in) :: vpuns(mmax,5),vkb(mmax,mxprj,4,nkap)
 real(dp), intent(in) :: evkb(mxprj,4,nkap)
 real(dp), intent(in) :: rhomod(mmax,5)

!Local variables
 integer :: ii,it,kk,ll,l1,ierr,nvt
 integer :: ikap,kap,mkap
 integer :: nat(30),lat(30),natp(30),latp(30),nav(4)
 real(dp) :: eaetst,etsttot,fj
 real(dp) :: eat(30,nkap,3),fat(30,3),eatp(30,nkap),fatp(30,3)

 real(dp),allocatable :: rho(:),rhocps(:)

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

 if(istcnf(jj)/=0) then
   call config_ae_warning(jj,irel,istcnf,ierrcnf,itcnf,kerrcnf,natcnf,latcnf)
   return
 end if

 nvt=nvtcnf(jj)
 nat(:)=natcnf(:,jj)
 lat(:)=latcnf(:,jj)
 fat(:,:)=fatcnf(:,:,jj)
 eat(:,:,:)=eatcnf(:,:,:,jj)
 eaetst=etstcnf(jj)

 allocate(rho(mmax),rhocps(mmax))

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
&           rr,rcmax,mmax,mxprj,iexc,etsttot,nproj,vpuns,lloc, &
&           vkb,evkb,ierr,nkap)

 if(ierr/=0) then
   if(irel==2) then
     write(6,'(a/a,i2)') &
&          'run_config_r: WARNING for fully non-local PS atom,', &
&          ' WARNING no output for configuration',jj
   else
     write(6,'(a,a/a,i2)') &
&          'run_config: WARNING for fully non-local PS atom,', &
&          ' stg. 1', &
&          ' WARNING no output for configuration',jj
   end if
   deallocate(rho,rhocps)
   return
 end if

!second pseudoatom run from maximally-ionized to excited configuration

 do kk=1,nvt
   eatp(kk,:)=eat(nc+kk,:,2)
   fatp(kk,1)=fat(nc+kk,2)
   fatp(kk,2)=fat(nc+kk,3)
 end do

 call psatom(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
&           rr,rcmax,mmax,mxprj,iexc,etsttot,nproj,vpuns,lloc, &
&           vkb,evkb,ierr,nkap)

 if(ierr/=0) then
   if(irel==2) then
     write(6,'(a/a,i2)') &
&          'run_config_r: WARNING for fully non-local PS atom,', &
&          ' WARNING no output for configuration',jj
   else
     write(6,'(a,a/a,i2)') &
&          'run_config: WARNING for fully non-local PS atom,', &
&          ' stg. 2', &
&          ' WARNING no output for configuration',jj
   end if
   deallocate(rho,rhocps)
   return
 end if



 if(irel==2) then
   call write_results_r
 else
   call write_results
 end if

 deallocate(rho,rhocps)

 return

 contains

 subroutine write_results
!! the eigenvalues and energies of the all-electron and pseudo atoms
 write(6,'(/a)') '   n   l     f        eae           eps        diff'
 do ii=1,nc
   write(6,'(2i4,f8.4,f14.8)') nat(ii),lat(ii),fat(ii,3),eat(ii,1,3)
 end do
 do ii=1,nvt
   if(fat(ii+nc,3)==0.0d0) cycle
   write(6,'(2i4,f8.4,2f14.8,1p,d12.2)') nat(ii+nc),lat(ii+nc),fat(ii+nc,3), &
&   eat(ii+nc,1,3),eatp(ii,1),eatp(ii,1)-eat(ii+nc,1,3)
 end do

 write(6,'(/a)') '    Total energies and differences'
 write(6,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      AE_ref=',etot, &
& '  AE_tst=',eaetst,'  dif=',eaetst-etot
 write(6,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      PS_ref=',epstot, &
& '  PS_tst=',etsttot,'  dif=',etsttot-epstot
 write(6,'(a,1p,d10.2)') '      PSP excitation error=', &
& eaetst-etot-etsttot+epstot
 return
 end subroutine write_results

 subroutine write_results_r
!! write_results for irel=2, a line for each kappa with its share of the
!! occupation
 write(6,'(/a)') '   n   l kap    f        eae           eps        diff'
 do ii=1,nc
   ll = lat(ii)
   if(ll==0) then
     mkap=1
   else
     mkap=2
   end if
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     if(ikap==1) then
       kap=-(ll+1)
       fj=dble(ll+1)*fat(ii,3)/dble(2*ll+1)
     else
       kap=  ll
       fj=dble(ll)*fat(ii,3)/dble(2*ll+1)
     end if
     write(6,'(3i4,f8.4,f14.8)') nat(ii),lat(ii),kap,fj,eat(ii,ikap,3)
   end do
 end do

 do ii=1,nvt
   ll = lat(ii+nc)
   if(ll==0) then
     mkap=1
   else
     mkap=2
   end if
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     if(ikap==1) then
       kap=-(ll+1)
       fj=dble(ll+1)*fat(ii+nc,3)/dble(2*ll+1)
     else
       kap=  ll
       fj=dble(ll)*fat(ii+nc,3)/dble(2*ll+1)
     end if
     write(6,'(3i4,f8.4,2f14.8,1p,e12.2)') nat(ii+nc),lat(ii+nc),kap,fj, &
&     eat(ii+nc,ikap,3),eatp(ii,ikap),eatp(ii,ikap)-eat(ii+nc,ikap,3)
   end do
 end do

 write(6,'(/a)') '    Total energies and differences'
 write(6,'(a,1p,e16.8,a,e16.8,a,e10.2)') '      AE_ref=',etot, &
& '  AE_tst=',eaetst,'  dif=',eaetst-etot
 write(6,'(a,1p,e16.8,a,e16.8,a,e10.2)') '      PS_ref=',epstot, &
& '  PS_tst=',etsttot,'  dif=',etsttot-epstot
 write(6,'(a,1p,e10.2)') '      PSP excitation error=', &
& eaetst-etot-etsttot+epstot
 return
 end subroutine write_results_r

 end subroutine run_config_ps
