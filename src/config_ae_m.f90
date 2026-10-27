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
 use constants_m, only: dp, mxstate, mxcnf, nlev, ae_tolres
 use kappa_m, only: kap_count, kap_of
 use ae_scf_m, only: ae_floor_warning
 implicit none
 private
 public :: run_config_ae, config_ae_warning, aeatom_cached

!> a successful aeatom solution, with the inputs that determine it: sratom
!> and relatom overwrite the eigenvalues' guess, so the atom is a function
!> of these alone
 type :: ae_entry_t
   integer :: irel,nkap,nc,ncv,mmax,iexc,it
   real(dp) :: zz,etot,eterm(4),res
   integer, allocatable :: na(:),la(:)
   real(dp), allocatable :: fa(:),rr(:),ea(:,:),rpk(:,:),rhoc(:),rho(:),vi(:)
 end type ae_entry_t

!> the solutions of this run (aeatom_cached)
 type(ae_entry_t), save :: ae_cache(2*mxcnf+1)
 integer, save :: nae_cache=0

 contains

 subroutine aeatom_cached(irel,nkap,na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&                         rr,vi,zz,mmax,iexc,etot,eterm,res,ierr,iout,ifail)

!! aeatom, or the copy of its earlier solution with the same inputs. The
!! reference configuration recurs as test configuration 0, at both levels
!! of run_config_ae, and the excited level of a configuration that only
!! takes charge away is its intermediate one; each was solved again, to
!! the same bits. A failed solution is not kept, so that its warnings
!! recur, and a reused one stopped on the residual's round-off floor
!! repeats that warning. rpk comes back as that solution left it (rows
!! 1..ncv). The arguments are those of aeatom

 use aeatom_m, only: aeatom

!Input variables
 integer, intent(in) :: irel,nkap,mmax,iexc,nc,ncv
 integer, intent(in) :: na(ncv),la(ncv)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(ncv),rr(mmax)

!Input/output variables
 integer, intent(inout) :: ierr
 real(dp), intent(inout) :: ea(mxstate,nkap),rpk(mxstate,nkap)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot,eterm(4),res
 real(dp), intent(out) :: rho(mmax),rhoc(mmax),vi(mmax)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 integer :: ii
 real(dp), allocatable :: vfree(:)

 ifail=0
 do ii=1,nae_cache
   associate(ce=>ae_cache(ii))
     if(ce%irel/=irel .or. ce%nkap/=nkap .or. ce%nc/=nc .or. &
&       ce%ncv/=ncv .or. ce%mmax/=mmax .or. ce%iexc/=iexc .or. &
&       ce%zz/=zz) cycle
     if(any(ce%na/=na) .or. any(ce%la/=la) .or. any(ce%fa/=fa) .or. &
&       any(ce%rr/=rr)) cycle
     ea(1:ncv,:)=ce%ea
     rpk(1:ncv,:)=ce%rpk
     it=ce%it
     rhoc(:)=ce%rhoc
     rho(:)=ce%rho
     vi(:)=ce%vi
     etot=ce%etot
     eterm(:)=ce%eterm
     res=ce%res
     ierr=0
     if(res>ae_tolres) call ae_floor_warning(res,ae_tolres,iout)
     return
   end associate
 end do

 allocate(vfree(mmax))
 vfree(:)=0.0d0
 call aeatom(irel,nkap,na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&            rr,vi,zz,mmax,iexc,etot,eterm,res,ierr,vfree,.false.,iout,ifail)
 if(ifail/=0) return
 deallocate(vfree)
 if(ierr/=0 .or. nae_cache==size(ae_cache)) return

 nae_cache=nae_cache+1
 associate(ce=>ae_cache(nae_cache))
   ce%irel=irel
   ce%nkap=nkap
   ce%nc=nc
   ce%ncv=ncv
   ce%mmax=mmax
   ce%iexc=iexc
   ce%zz=zz
   ce%na=na
   ce%la=la
   ce%fa=fa
   ce%rr=rr
   ce%ea=ea(1:ncv,:)
   ce%rpk=rpk(1:ncv,:)
   ce%it=it
   ce%rhoc=rhoc
   ce%rho=rho
   ce%vi=vi
   ce%etot=etot
   ce%eterm=eterm
   ce%res=res
 end associate
 return
 end subroutine aeatom_cached

 subroutine run_config_ae(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                         mmax,iexc,ea,irel,nkap,istcnf,ierrcnf,itcnf, &
&                         kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf, &
&                         iout,ifail)

!! solves the all-electron atom in the maximally-ionized and excited states
!! of test configuration jj, the all-electron half of the test-configuration
!! calculations.  These depend only on the configurations, not on the
!! pseudopotential, so oncvpsp_ae.x can compute them once for oncvpsp_ps.x
!! (see ae_file_m).  run_config_ps does the pseudo-atom half.
!!
!! @param jj  index of current configufation
!! @param nacnf  principal quantum number array, all configurations
!! @param lacnf  angular-momenta array, all config.
!! @param facnf  occupancies, all config.
!! @param nc  number of core states
!! @param nvcnf  number of valence states, all config.
!! @param rr  log radial mesh
!! @param zz  atomic number
!! @param mmax  size of log grid
!! @param iexc  exchange-correlation function to be used
!! @param ea  reference configurattion eigenvalues, ea(:,1) for j=l+1/2 and
!!            ea(:,2) for j=l-1/2 when irel=2
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappa per level: 1, or 2 for irel=2
!! @param istcnf  0 success, 1 all-electron atom failed, 2 the bound-state
!!                solver failed for an empty level (with j=l+1/2 when irel=2), 3
!!                the same with j=l-1/2, all config.
!! @param ierrcnf  aeatom error code when istcnf==1, boundstate error code when
!!                 istcnf>=2, all config.
!! @param itcnf  aeatom iteration count when istcnf/=0, all config.
!! @param kerrcnf  index in natcnf/latcnf of the level whose bound state failed
!!                 when istcnf>=2, all config.
!! @param nvtcnf  number of merged valence levels, all config.
!! @param natcnf  merged level principal quantum numbers, reference levels
!!                first, all config.
!! @param latcnf  merged level angular momenta, all config.
!! @param fatcnf  merged level occupancies (:,1:3,jj): reference,
!!                maximally-ionized, excited, all config.
!! @param eatcnf  merged level all-electron eigenvalues (:,ikap,1:3,jj):
!!                reference, maximally-ionized, excited, all config.
!! @param etstcnf  all-electron total energy of the excited configuration, all
!!                 config.
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written
!!
!! Only the column jj of istcnf through etstcnf is set.

!Input variables
 use boundstate_m, only: boundstate
 integer, intent(in) :: jj,mmax,iexc,nc,irel,nkap
 integer, intent(in) :: nacnf(mxstate,mxcnf),lacnf(mxstate,mxcnf),nvcnf(mxcnf)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: facnf(mxstate,mxcnf),ea(mxstate,nkap),rr(mmax)

!Input/output variables
 integer, intent(inout) :: istcnf(mxcnf),ierrcnf(mxcnf),itcnf(mxcnf)
 integer, intent(inout) :: kerrcnf(mxcnf)
 integer, intent(inout) :: nvtcnf(mxcnf),natcnf(mxstate,mxcnf)
 integer, intent(inout) :: latcnf(mxstate,mxcnf)
 real(dp), intent(inout) :: fatcnf(mxstate,nlev,mxcnf)
 real(dp), intent(inout) :: eatcnf(mxstate,nkap,nlev,mxcnf)
 real(dp), intent(inout) :: etstcnf(mxcnf)

!Output variables
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,it,kk,ll,ierr,mch,nvt
 integer :: ikap,kap,mkap
 integer :: nat(mxstate),lat(mxstate)
 integer ::  indxr(mxstate),indxe(mxstate)
 real(dp) :: et,eaetst,eterm(4),res
 real(dp) :: eat(mxstate,nkap,nlev),fat(mxstate,nlev),rpk(mxstate,nkap)

 real(dp),allocatable :: rho(:),rhoc(:),vfull(:)
 real(dp),allocatable :: uu(:,:),up(:,:)

 ifail=0
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

 call aeatom_cached(irel,nkap,nat,lat,eat(1,1,2),fat(1,2),rpk,nc, &
&                   nc+nvt,it,rhoc,rho,rr,vfull,zz,mmax,iexc,eaetst,eterm, &
&                   res,ierr,iout,ifail)
 if(ifail/=0) return
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
     mkap=kap_count(nkap,ll)
     et=0.0d0
     do ikap=1,mkap
       kap=kap_of(ll,ikap)
       call boundstate(irel,nat(kk),ll,kap,ierr,et, &
&                      rr,zz,vfull,uu,up,mmax,mch,iout,ifail)
       if(ifail/=0) return
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

 call aeatom_cached(irel,nkap,nat,lat,eat(1,1,3),fat(1,3),rpk,nc, &
&                   nc+nvt,it,rhoc,rho,rr,vfull,zz,mmax,iexc,eaetst,eterm, &
&                   res,ierr,iout,ifail)
 if(ifail/=0) return
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
!! copies the results so far into column jj of the output arrays
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
&                             natcnf,latcnf,iout)
!! prints the warning run_config_ps gives when the all-electron atom for test
!! configuration jj failed
!!
!! @param jj  index of the configuration (1 for the reference)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param istcnf,ierrcnf,itcnf,kerrcnf,natcnf,latcnf  run_config_ae results, all
!!        configurations
!! @param iout  unit of the text output

!Input variables
 integer, intent(in) :: jj,irel
 integer, intent(in) :: istcnf(mxcnf),ierrcnf(mxcnf),itcnf(mxcnf),kerrcnf(mxcnf)
 integer, intent(in) :: natcnf(mxstate,mxcnf),latcnf(mxstate,mxcnf)
 integer, intent(in) :: iout

!Local variables
 integer :: ll

 if(istcnf(jj)==1) then
   if(irel==2) then
     write(iout,'(a/a,i2)') 'run_config_r: WARNING  for AE atom,', &
&          ' WARNING no output for configuration',jj
   else
     write(iout,'(a/a,i2)') 'run_config: WARNING  for AE atom,', &
&          ' WARNING no output for configuration',jj
   end if
   if(ierrcnf(jj)==-1 .and. irel==2) then
     write(iout,'(a,i4)') &
&          'run_config_r: WARNING no classical turning point ' // &
&          'error, iteration', &
&          itcnf(jj)
   else if(ierrcnf(jj)==-1) then
     write(iout,'(a,i4)') 'no classical turning point error, iteration', &
&          itcnf(jj)
   else if(irel==2) then
     write(iout,'(a)') &
&          'run_config_r: WARNING self-consistency failed to converge'
   else
     write(iout,'(a)') 'run_config: WARNING self-consistency failed to converge'
   end if
 else if(istcnf(jj)>=2) then
   ll=latcnf(kerrcnf(jj),jj)
   if(irel/=2) then
     write(iout,'(/a,3i4)') &
&         'runconfig: WARNING lschfb convergence ERROR n,l,iter=', &
&         natcnf(kerrcnf(jj),jj),ll,itcnf(jj)
   else if(istcnf(jj)==2) then
     write(iout,'(/a,4i4)') &
&         'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,ierr=', &
&         natcnf(kerrcnf(jj),jj),ll,-(ll+1),ierrcnf(jj)
   else
     write(iout,'(/a,4i4)') &
&         'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=', &
&         natcnf(kerrcnf(jj),jj),ll,ll,itcnf(jj)
   end if
 end if
 return
 end subroutine config_ae_warning

 end module config_ae_m
