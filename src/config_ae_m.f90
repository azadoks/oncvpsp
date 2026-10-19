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
! all-electron half of the test-configuration calculations in run_config:
! the all-electron atom for the maximally-ionized intermediate and the
! excited configuration.  These depend only on the configurations, not
! on the pseudopotential, so they can be computed once and reused
! (see ae_file_m), for each of the non-, scalar- and fully-relativistic
! codes.

 module config_ae_m
 use constants_m, only: dp
 implicit none
 private
 public :: run_config_ae, config_ae_warning

 contains

 subroutine run_config_ae(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                         mmax,iexc,ea,irel,nkap,istcnf,ierrcnf,itcnf, &
&                         kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)

! solves the all-electron atom in the maximally-ionized and excited states
! of test configuration jj, the all-electron half of the test-configuration
! calculations.  These depend only on the configurations, not on the
! pseudopotential, so oncvpsp_ae.x can compute them once for oncvpsp_ps.x
! (see ae_file_m).  run_config_ps does the pseudo-atom half.

!jj  index of current configufation
!nacnf  principal quantum number array, all configurations
!lacnf  angular-momenta array, all config.
!facnf  occupancies, all config.
!nc  number of core states
!nvcnf  number of valence states, all config.
!rr  log radial mesh
!zz  atomic number
!mmax  size of log grid
!iexc  exchange-correlation function to be used
!ea  reference configurattion eigenvalues, ea(:,1) for j=l+1/2 and
!    ea(:,2) for j=l-1/2 when irel=2
!irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!nkap  number of kappa per level: 1, or 2 for irel=2
!istcnf  0 success, 1 all-electron atom failed, 2 the bound-state solver
!        failed for an empty level (with j=l+1/2 when irel=2), 3 the
!        same with j=l-1/2, all config.
!ierrcnf  aeatom error code when istcnf==1, boundstate error code when
!         istcnf>=2, all config.
!itcnf  aeatom iteration count when istcnf/=0, all config.
!kerrcnf  index in natcnf/latcnf of the level whose bound state failed when
!         istcnf>=2, all config.
!nvtcnf  number of merged valence levels, all config.
!natcnf  merged level principal quantum numbers, reference levels first,
!        all config.
!latcnf  merged level angular momenta, all config.
!fatcnf  merged level occupancies (:,1:3,jj): reference,
!        maximally-ionized, excited, all config.
!eatcnf  merged level all-electron eigenvalues (:,ikap,1:3,jj): reference,
!        maximally-ionized, excited, all config.
!etstcnf  all-electron total energy of the excited configuration, all
!         config.
!Only the column jj of istcnf through etstcnf is set.

!Input variables
 integer, intent(in) :: jj,mmax,iexc,nc,irel,nkap
 integer, intent(in) :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: facnf(30,5),ea(30,nkap),rr(mmax)

!Input/output variables
 integer, intent(inout) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(inout) :: nvtcnf(5),natcnf(30,5),latcnf(30,5)
 real(dp), intent(inout) :: fatcnf(30,3,5),eatcnf(30,nkap,3,5),etstcnf(5)

!Local variables
 integer :: ii,it,kk,ll,ierr,mch,nvt
 integer :: ikap,kap,mkap
 integer :: nat(30),lat(30)
 integer ::  indxr(30),indxe(30)
 real(dp) :: et,eaetst
 real(dp) :: eat(30,nkap,3),fat(30,3),rpk(30,nkap)

 real(dp),allocatable :: rho(:),rhoc(:),vfull(:)
 real(dp),allocatable :: uu(:,:),up(:,:)

 allocate(rho(mmax),rhoc(mmax),vfull(mmax))
 allocate(uu(mmax,2),up(mmax,2))

 istcnf(jj)=0
 ierrcnf(jj)=0
 itcnf(jj)=0
 kerrcnf(jj)=0

 it=0

 fat(:,:)=0.0d0
 nat(:)=0
 lat(:)=0

! set up cross indces of reference (r) and excited (e) cofigurations
 indxr(:)=0
 indxe(:)=0
 do ii=1,nc+nvcnf(1)
   do kk=1,nc+nvcnf(jj)
     if(nacnf(kk,jj)==nacnf(ii,1) .and. &
&       lacnf(kk,jj)==lacnf(ii,1)) then
       indxr(ii)=kk
       indxe(kk)=ii
     end if
   end do
 end do

! set up merged list with reference levels first
! levels not present in one confguration or the other get zero occupancy
 eat(:,:,:)=0.0d0
 do ii=1,nc+nvcnf(1)
   nat(ii)=nacnf(ii,1)
   lat(ii)=lacnf(ii,1)
   fat(ii,1)=facnf(ii,1)
   eat(ii,:,1)=ea(ii,:)
   if(indxr(ii)>0) then
     fat(ii,2)=facnf(indxr(ii),jj)
   else
     fat(ii,2)=0.0d0
   end if
   fat(ii,3)=fat(ii,2)
 end do

 nvt=nvcnf(1)
 do kk=1,nc+nvcnf(jj)
   if(indxe(kk)==0) then
     nvt=nvt+1
     nat(nc+nvt)=nacnf(kk,jj)
     lat(nc+nvt)=lacnf(kk,jj)
     fat(nc+nvt,1)=0.0d0
     fat(nc+nvt,3)=facnf(kk,jj)
   end if
 end do

!set up array for maximally ionized intermediate state
 do kk=1,nc+nvt
   fat(kk,2)=min(fat(kk,1),fat(kk,3))
 end do



 eat(:,:,2)=0.d0; rpk(:,:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
 eaetst=0.d0

!all-electron atom solution for maximally-ionized state

 call aeatom(irel,nkap,nat,lat,eat(1,1,2),fat(1,2),rpk,nc,nc+nvt,it, &
&            rhoc,rho,rr,vfull,zz,mmax,iexc,eaetst,ierr)
 if(ierr/=0) then
   call store_results
   istcnf(jj)=1
   ierrcnf(jj)=ierr
   itcnf(jj)=it
   deallocate(rho,rhoc,vfull)
   deallocate(uu,up)
   return
 end if

!fill in energies for empty levels of maximally ionized state, for each
!kappa when irel=2 (the j=l-1/2 search starting from the j=l+1/2 energy)

 do kk=1,nc+nvt
   if(eat(kk,1,2)==0.0d0) then
     ll=lat(kk)
     if(nkap==1 .or. ll==0) then
       mkap=1
     else
       mkap=2
     end if
     et=0.0d0
     do ikap=1,mkap
       if(ikap==1) then
         kap=-(ll+1)
       else
         kap=ll
       end if
       call boundstate(irel,nat(kk),ll,kap,ierr,et, &
&                      rr,zz,vfull,uu,up,mmax,mch)
!a negative ierr (no classical turning point) is accepted for j=l-1/2
       if(ierr>0 .or. (ierr/=0 .and. ikap==1)) then
         call store_results
         istcnf(jj)=1+ikap
         ierrcnf(jj)=ierr
         kerrcnf(jj)=kk
         itcnf(jj)=it
         deallocate(rho,rhoc,vfull)
         deallocate(uu,up)
         return
       end if
       eat(kk,ikap,2)=et
     end do
     if(nkap==2 .and. ll==0) eat(kk,2,2)=eat(kk,1,2)
   end if
 end do

 eat(:,:,3)=0.d0; rpk(:,:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
 eaetst=0.d0

!all-electron atom solution for excited state

 call aeatom(irel,nkap,nat,lat,eat(1,1,3),fat(1,3),rpk,nc,nc+nvt,it, &
&            rhoc,rho,rr,vfull,zz,mmax,iexc,eaetst,ierr)
 if(ierr/=0) then
   call store_results
   istcnf(jj)=1
   ierrcnf(jj)=ierr
   itcnf(jj)=it
   deallocate(rho,rhoc,vfull)
   deallocate(uu,up)
   return
 end if

 call store_results
 deallocate(rho,rhoc,vfull)
 deallocate(uu,up)
 return

 contains

 subroutine store_results
! copies the results so far into column jj of the output arrays
 nvtcnf(jj)=nvt
 natcnf(:,jj)=nat(:)
 latcnf(:,jj)=lat(:)
 fatcnf(:,:,jj)=fat(:,:)
 eatcnf(:,:,:,jj)=eat(:,:,:)
 etstcnf(jj)=eaetst
 return
 end subroutine store_results

 end subroutine run_config_ae

 subroutine config_ae_warning(jj,irel,istcnf,ierrcnf,itcnf,kerrcnf, &
&                             natcnf,latcnf)
! prints the warning run_config_ps gives when the all-electron atom for test
! configuration jj failed

!jj  index of the configuration (1 for the reference)
!irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!istcnf, ierrcnf, itcnf, kerrcnf, natcnf, latcnf  run_config_ae results, all
!  configurations

!Input variables
 integer, intent(in) :: jj,irel
 integer, intent(in) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(in) :: natcnf(30,5),latcnf(30,5)

!Local variables
 integer :: ll

 if(istcnf(jj)==1) then
   if(irel==2) then
     write(6,'(a/a,i2)') 'run_config_r: WARNING  for AE atom,', &
&          ' WARNING no output for configuration',jj
   else
     write(6,'(a/a,i2)') 'run_config: WARNING  for AE atom,', &
&          ' WARNING no output for configuration',jj
   end if
   if(ierrcnf(jj)==-1 .and. irel==2) then
     write(6,'(a,i4)') &
&          'run_config_r: WARNING no classical turning point ' // &
&          'error, iteration', &
&          itcnf(jj)
   else if(ierrcnf(jj)==-1) then
     write(6,'(a,i4)') 'no classical turning point error, iteration',itcnf(jj)
   else if(irel==2) then
     write(6,'(a)') 'run_config_r: WARNING self-consistency failed to converge'
   else
     write(6,'(a)') 'run_config: WARNING self-consistency failed to converge'
   end if
 else if(istcnf(jj)>=2) then
   ll=latcnf(kerrcnf(jj),jj)
   if(irel/=2) then
     write(6,'(/a,3i4)') &
&         'runconfig: WARNING lschfb convergence ERROR n,l,iter=', &
&         natcnf(kerrcnf(jj),jj),ll,itcnf(jj)
   else if(istcnf(jj)==2) then
     write(6,'(/a,4i4)') &
&         'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,ierr=', &
&         natcnf(kerrcnf(jj),jj),ll,-(ll+1),ierrcnf(jj)
   else
     write(6,'(/a,4i4)') &
&         'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=', &
&         natcnf(kerrcnf(jj),jj),ll,ll,itcnf(jj)
   end if
 end if
 return
 end subroutine config_ae_warning

 end module config_ae_m
