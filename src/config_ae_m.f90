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
! (see ae_file_m).  run_config_ae_r is the fully-relativistic version,
! the all-electron half of run_config_r.

 module config_ae_m
 use constants_m, only: dp
 implicit none
 private
 public :: run_config_ae, config_ae_warning
 public :: run_config_ae_r, config_ae_r_warning

 contains

 subroutine run_config_ae(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                         mmax,iexc,ea,srel,istcnf,ierrcnf,itcnf,kerrcnf, &
&                         nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)

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
!ea  reference configurattion eigenvalues
!srel .true. for scalar-relativistic, .false. for non-relativistic
!istcnf  0 success, 1 all-electron atom failed, 2 lschfb failed for an
!        empty level, all config.
!ierrcnf  sratom error code when istcnf==1, all config.
!itcnf  sratom iteration count when istcnf/=0, all config.
!kerrcnf  index in natcnf/latcnf of the level whose lschfb failed when
!         istcnf==2, all config.
!nvtcnf  number of merged valence levels, all config.
!natcnf  merged level principal quantum numbers, reference levels first,
!        all config.
!latcnf  merged level angular momenta, all config.
!fatcnf  merged level occupancies (:,1:3,jj): reference,
!        maximally-ionized, excited, all config.
!eatcnf  merged level all-electron eigenvalues (:,1:3,jj): reference,
!        maximally-ionized, excited, all config.
!etstcnf  all-electron total energy of the excited configuration, all
!         config.
!Only the column jj of istcnf through etstcnf is set.

!Input variables
 integer, intent(in) :: jj,mmax,iexc,nc
 integer, intent(in) :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: facnf(30,5),ea(30),rr(mmax)
 logical, intent(in) :: srel

!Input/output variables
 integer, intent(inout) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(inout) :: nvtcnf(5),natcnf(30,5),latcnf(30,5)
 real(dp), intent(inout) :: fatcnf(30,3,5),eatcnf(30,3,5),etstcnf(5)

!Local variables
 integer :: ii,it,kk,ierr,mch,nvt
 integer :: nat(30),lat(30)
 integer ::  indxr(30),indxe(30)
 real(dp) :: et,eaetst
 real(dp) :: eat(30,3),fat(30,3),rpk(30)

 real(dp),allocatable :: rho(:),rhoc(:),vfull(:)
 real(dp),allocatable :: uu(:),up(:)

 allocate(rho(mmax),rhoc(mmax),vfull(mmax))
 allocate(uu(mmax),up(mmax))

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
 eat(:,1)=0.0d0
 do ii=1,nc+nvcnf(1)
   nat(ii)=nacnf(ii,1)
   lat(ii)=lacnf(ii,1)
   fat(ii,1)=facnf(ii,1)
   eat(ii,1)=ea(ii)
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



 eat(:,2)=0.d0; rpk(:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
 eaetst=0.d0

!all-electron atom solution for maximally-ionized state

 call sratom(nat,lat,eat(1,2),fat(1,2),rpk,nc,nc+nvt,it,rhoc,rho, &
&            rr,vfull,zz,mmax,iexc,eaetst,ierr,srel)
 if(ierr/=0) then
   call store_results
   istcnf(jj)=1
   ierrcnf(jj)=ierr
   itcnf(jj)=it
   deallocate(rho,rhoc,vfull)
   deallocate(uu,up)
   return
 end if

!fill in energies for empty levels of maximally ionized state

 do kk=1,nc+nvt
   if(eat(kk,2)==0.0d0) then
     et=0.0d0
     call lschfb(nat(kk),lat(kk),ierr,et, &
&                rr,vfull,uu,up,zz,mmax,mch,srel)
     if(ierr /= 0) then
       call store_results
       istcnf(jj)=2
       kerrcnf(jj)=kk
       itcnf(jj)=it
       deallocate(rho,rhoc,vfull)
       deallocate(uu,up)
       return
     end if
     eat(kk,2)=et
   end if
 end do

 eat(:,3)=0.d0; rpk(:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
 eaetst=0.d0

!all-electron atom solution for excited state

 call sratom(nat,lat,eat(1,3),fat(1,3),rpk,nc,nc+nvt,it,rhoc,rho, &
&            rr,vfull,zz,mmax,iexc,eaetst,ierr,srel)
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
 eatcnf(:,:,jj)=eat(:,:)
 etstcnf(jj)=eaetst
 return
 end subroutine store_results

 end subroutine run_config_ae

 subroutine config_ae_warning(jj,istcnf,ierrcnf,itcnf,kerrcnf,natcnf,latcnf)
! prints the warning run_config gives when the all-electron atom for test
! configuration jj failed

!jj  index of the configuration (1 for the reference)
!istcnf, ierrcnf, itcnf, kerrcnf, natcnf, latcnf  run_config_ae results, all
!  configurations

!Input variables
 integer, intent(in) :: jj
 integer, intent(in) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(in) :: natcnf(30,5),latcnf(30,5)

 select case(istcnf(jj))
 case(1)
   write(6,'(a/a,i2)') 'run_config: WARNING  for AE atom,', &
&        ' WARNING no output for configuration',jj
   if(ierrcnf(jj)==-1) then
     write(6,'(a,i4)') 'no classical turning point error, iteration',itcnf(jj)
   else
     write(6,'(a)') 'run_config: WARNING self-consistency failed to converge'
   end if
 case(2)
   write(6,'(/a,3i4)') &
&       'runconfig: WARNING lschfb convergence ERROR n,l,iter=', &
&       natcnf(kerrcnf(jj),jj),latcnf(kerrcnf(jj),jj),itcnf(jj)
 case default
! the atom converged: nothing to report
 end select
 return
 end subroutine config_ae_warning

 subroutine run_config_ae_r(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                           mmax,iexc,ea,istcnf,ierrcnf,itcnf,kerrcnf, &
&                           nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)

! fully-relativistic run_config_ae, the all-electron half of the
! test-configuration calculations of run_config_ps_r

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
!ea  reference configurattion eigenvalues
!istcnf  0 success, 1 all-electron atom failed, 2 ldiracfb failed for an
!        empty level with kap=l+1, 3 with kap=l, all config.
!ierrcnf  relatom error code when istcnf==1, ldiracfb error code when
!         istcnf>=2, all config.
!itcnf  relatom iteration count when istcnf/=0, all config.
!kerrcnf  index in natcnf/latcnf of the level whose ldiracfb failed when
!         istcnf>=2, all config.
!nvtcnf, natcnf, latcnf, fatcnf, etstcnf  as in run_config_ae
!eatcnf  merged level all-electron eigenvalues (:,ikap,1:3,jj) for
!        kap=-(l+1) (ikap=1) and kap=l (ikap=2): reference,
!        maximally-ionized, excited, all config.
!Only the column jj of istcnf through etstcnf is set.

!Input variables
 integer, intent(in) :: jj,mmax,iexc,nc
 integer, intent(in) :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: facnf(30,5),ea(30,2),rr(mmax)

!Input/output variables
 integer, intent(inout) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(inout) :: nvtcnf(5),natcnf(30,5),latcnf(30,5)
 real(dp), intent(inout) :: fatcnf(30,3,5),eatcnf(30,2,3,5),etstcnf(5)

!Local variables
 integer :: ii,it,kk,ll,ierr,mch,nvt
 integer :: kap
 integer :: nat(30),lat(30)
 integer ::  indxr(30),indxe(30)
 real(dp) :: et,eaetst
 real(dp) :: eat(30,2,3),fat(30,3),rpk(30,2)

 real(dp),allocatable :: rho(:),rhoc(:),vfull(:)
 real(dp),allocatable :: uu(:,:),up(:,:)

 allocate(rho(mmax),rhoc(mmax),vfull(mmax))
 allocate(uu(mmax,2),up(mmax,2))

 istcnf(jj)=0
 ierrcnf(jj)=0
 itcnf(jj)=0
 kerrcnf(jj)=0

 it=0

! set up pseudo atom arrays
!arrange data arrays for adiabatic configufation switching of pseudo atom

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

 call relatom(nat,lat,eat(1,1,2),fat(1,2),rpk,nc,nc+nvt,it,rhoc,rho, &
&            rr,vfull,zz,mmax,iexc,eaetst,ierr)
 if(ierr/=0) then
   call store_results
   istcnf(jj)=1
   ierrcnf(jj)=ierr
   itcnf(jj)=it
   deallocate(rho,rhoc,vfull)
   deallocate(uu,up)
   return
 end if

!fill in energies for empty levels of maximally ionized state

 do kk=1,nc+nvt
   if(eat(kk,1,2)==0.0d0) then
     et=0.0d0
     ll=lat(kk)
     kap=ll+1
     call ldiracfb(nat(kk),ll,kap,ierr,et, &
&                rr,zz,vfull,uu,up,mmax,mch)
     if(ierr /= 0) then
       call store_results
       istcnf(jj)=2
       ierrcnf(jj)=ierr
       kerrcnf(jj)=kk
       itcnf(jj)=it
       deallocate(rho,rhoc,vfull)
       deallocate(uu,up)
       return
     end if
     eat(kk,1,2)=et
     if(ll==0) then
       eat(kk,2,2)=et
     else
       kap=ll
       call ldiracfb(nat(kk),ll,kap,ierr,et, &
&                  rr,zz,vfull,uu,up,mmax,mch)
       if(ierr>0) then
         call store_results
         istcnf(jj)=3
         ierrcnf(jj)=ierr
         kerrcnf(jj)=kk
         itcnf(jj)=it
         deallocate(rho,rhoc,vfull)
         deallocate(uu,up)
         return
       end if
       eat(kk,2,2)=et
     end if
   end if
 end do

 eat(:,:,3)=0.d0; rpk(:,:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
 eaetst=0.d0

!all-electron atom solution for excited state

 call relatom(nat,lat,eat(1,1,3),fat(1,3),rpk,nc,nc+nvt,it,rhoc,rho, &
&            rr,vfull,zz,mmax,iexc,eaetst,ierr)
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

 end subroutine run_config_ae_r

 subroutine config_ae_r_warning(jj,istcnf,ierrcnf,itcnf,kerrcnf,natcnf,latcnf)
! prints the warning run_config_r gives when the all-electron atom for
! test configuration jj failed

!jj  index of the configuration (1 for the reference)
!istcnf, ierrcnf, itcnf, kerrcnf, natcnf, latcnf  run_config_ae_r results, all
!  configurations

!Input variables
 integer, intent(in) :: jj
 integer, intent(in) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(in) :: natcnf(30,5),latcnf(30,5)

!Local variables
 integer :: ll

 select case(istcnf(jj))
 case(1)
   write(6,'(a/a,i2)') 'run_config_r: WARNING  for AE atom,', &
&        ' WARNING no output for configuration',jj
   if(ierrcnf(jj)==-1) then
     write(6,'(a,i4)') &
&          'run_config_r: WARNING no classical turning point ' // &
&          'error, iteration', &
&          itcnf(jj)
   else
     write(6,'(a)') 'run_config_r: WARNING self-consistency failed to converge'
   end if
 case(2)
   ll=latcnf(kerrcnf(jj),jj)
   write(6,'(/a,4i4)') &
&       'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,ierr=', &
&       natcnf(kerrcnf(jj),jj),ll,ll+1,ierrcnf(jj)
 case(3)
   ll=latcnf(kerrcnf(jj),jj)
   write(6,'(/a,4i4)') &
&       'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=', &
&       natcnf(kerrcnf(jj),jj),ll,ll,itcnf(jj)
 case default
! the atom converged: nothing to report
 end select
 return
 end subroutine config_ae_r_warning

 end module config_ae_m
