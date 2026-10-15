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
 implicit none
 private
 public :: config_ae_t, run_config_ae, config_ae_warning
 public :: config_ae_r_t, run_config_ae_r, config_ae_r_warning

 integer, parameter :: dp=kind(1.0d0)

!all-electron results for one test configuration, merged with the
!reference configuration (reference levels first)
 type :: config_ae_t
!stat  0 success, 1 all-electron atom failed, 2 lschfb failed for an
!      empty level
!ierr  sratom error code when stat==1
!it  sratom iteration count
!kerr  index in nat/lat of the level whose lschfb failed when stat==2
!nvt  number of merged valence levels
!nat, lat  merged level quantum numbers
!fat(:,1:3)  occupations: reference, maximally-ionized, excited
!eat(:,1:3)  all-electron eigenvalues: reference, maximally-ionized,
!            excited
!eaetst  all-electron total energy of the excited configuration
  integer :: stat=0,ierr=0,it=0,kerr=0,nvt=0
  integer :: nat(30)=0,lat(30)=0
  real(dp) :: fat(30,3)=0.0d0,eat(30,3)=0.0d0
  real(dp) :: eaetst=0.0d0
 end type config_ae_t

!fully-relativistic counterpart, from run_config_ae_r
 type :: config_ae_r_t
!stat  0 success, 1 all-electron atom failed, 2 ldiracfb failed for an
!      empty level with kap=l+1, 3 with kap=l
!ierr  relatom error code when stat==1, ldiracfb error code when stat==2
!it  relatom iteration count
!kerr  index in nat/lat of the level whose ldiracfb failed when stat>=2
!nvt, nat, lat, fat, eaetst  as in config_ae_t
!eat(:,ikap,1:3)  all-electron eigenvalues for kap=-(l+1) (ikap=1) and
!                 kap=l (ikap=2): reference, maximally-ionized, excited
  integer :: stat=0,ierr=0,it=0,kerr=0,nvt=0
  integer :: nat(30)=0,lat(30)=0
  real(dp) :: fat(30,3)=0.0d0,eat(30,2,3)=0.0d0
  real(dp) :: eaetst=0.0d0
 end type config_ae_r_t

 contains

 subroutine run_config_ae(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                         mmax,iexc,ea,srel,cae)

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
!cae  all-electron results for this configuration

!Input variables
 integer :: jj,mmax,iexc,nc
 integer :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 real(dp) :: zz
 real(dp) :: facnf(30,5),ea(30),rr(mmax)
 logical :: srel

!Output variables
 type(config_ae_t), intent(out) :: cae

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
&         lacnf(kk,jj)==lacnf(ii,1)) then
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
     fat(ii,2)=dmin1(fat(ii,1),fat(ii,3))
   end do



 eat(:,2)=0.d0; rpk(:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
 eaetst=0.d0

!all-electron atom solution for maximally-ionized state

 call sratom(nat,lat,eat(1,2),fat(1,2),rpk,nc,nc+nvt,it,rhoc,rho, &
&            rr,vfull,zz,mmax,iexc,eaetst,ierr,srel)
 if(ierr/=0) then
  call store_results
  cae%stat=1
  cae%ierr=ierr
  cae%it=it
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
     if(ierr .ne. 0) then
       call store_results
       cae%stat=2
       cae%kerr=kk
       cae%it=it
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
  cae%stat=1
  cae%ierr=ierr
  cae%it=it
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
  cae%nvt=nvt
  cae%nat(:)=nat(:)
  cae%lat(:)=lat(:)
  cae%fat(:,:)=fat(:,:)
  cae%eat(:,:)=eat(:,:)
  cae%eaetst=eaetst
 end subroutine store_results

 end subroutine run_config_ae

 subroutine config_ae_warning(jj,cae)
! prints the warning run_config gives when the all-electron atom for test
! configuration jj failed
 integer :: jj
 type(config_ae_t), intent(in) :: cae

 select case(cae%stat)
 case(1)
  write(6,'(a/a,i2)') 'run_config: WARNING  for AE atom,', &
&       ' WARNING no output for configuration',jj
  if(cae%ierr==-1) then
    write(6,'(a,i4)') 'no classical turning point error, iteration',cae%it
  else
    write(6,'(a)') 'run_config: WARNING self-consistency failed to converge'
  end if
 case(2)
  write(6,'(/a,3i4)') &
&      'runconfig: WARNING lschfb convergence ERROR n,l,iter=', &
&      cae%nat(cae%kerr),cae%lat(cae%kerr),cae%it
 end select
 end subroutine config_ae_warning

 subroutine run_config_ae_r(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                           mmax,iexc,ea,cae)

! fully-relativistic run_config_ae, the all-electron half of run_config_r

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
!cae  all-electron results for this configuration

!Input variables
 integer :: jj,mmax,iexc,nc
 integer :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 real(dp) :: zz
 real(dp) :: facnf(30,5),ea(30,2),rr(mmax)

!Output variables
 type(config_ae_r_t), intent(out) :: cae

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
&         lacnf(kk,jj)==lacnf(ii,1)) then
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
     fat(ii,2)=dmin1(fat(ii,1),fat(ii,3))
   end do



 eat(:,:,2)=0.d0; rpk(:,:)=0.d0; rhoc(:)=0.d0; rho(:)=0.d0; vfull(:)=0.d0
 eaetst=0.d0

!all-electron atom solution for maximally-ionized state

 call relatom(nat,lat,eat(1,1,2),fat(1,2),rpk,nc,nc+nvt,it,rhoc,rho, &
&            rr,vfull,zz,mmax,iexc,eaetst,ierr)
 if(ierr/=0) then
  call store_results
  cae%stat=1
  cae%ierr=ierr
  cae%it=it
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
     if(ierr .ne. 0) then
       call store_results
       cae%stat=2
       cae%ierr=ierr
       cae%kerr=kk
       cae%it=it
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
         cae%stat=3
         cae%ierr=ierr
         cae%kerr=kk
         cae%it=it
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
  cae%stat=1
  cae%ierr=ierr
  cae%it=it
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
  cae%nvt=nvt
  cae%nat(:)=nat(:)
  cae%lat(:)=lat(:)
  cae%fat(:,:)=fat(:,:)
  cae%eat(:,:,:)=eat(:,:,:)
  cae%eaetst=eaetst
 end subroutine store_results

 end subroutine run_config_ae_r

 subroutine config_ae_r_warning(jj,cae)
! prints the warning run_config_r gives when the all-electron atom for
! test configuration jj failed
 integer :: jj
 type(config_ae_r_t), intent(in) :: cae

 integer :: ll

 select case(cae%stat)
 case(1)
  write(6,'(a/a,i2)') 'run_config_r: WARNING  for AE atom,', &
&       ' WARNING no output for configuration',jj
  if(cae%ierr==-1) then
    write(6,'(a,i4)') &
&         'run_config_r: WARNING no classical turning point error, iteration', &
&         cae%it
  else
    write(6,'(a)') 'run_config_r: WARNING self-consistency failed to converge'
  end if
 case(2)
  ll=cae%lat(cae%kerr)
  write(6,'(/a,4i4)') &
&      'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,ierr=', &
&      cae%nat(cae%kerr),ll,ll+1,cae%ierr
 case(3)
  ll=cae%lat(cae%kerr)
  write(6,'(/a,4i4)') &
&      'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=', &
&      cae%nat(cae%kerr),ll,ll,cae%it
 end select
 end subroutine config_ae_r_warning

 end module config_ae_m
