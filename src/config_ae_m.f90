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
 implicit none
 private
 public :: run_config_ae, config_ae_warning, aeatom_cached, ae_config_t, &
&           ae_cache_t

!> the all-electron half of a test configuration (run_config_ae): its
!> levels merged with the reference's, the reference's first
 type :: ae_config_t
!> 0 success, 1 the all-electron atom failed, 2 the bound-state solver
!> failed for an empty level (with j=l+1/2 when irel=2), 3 the same with
!> j=l-1/2
   integer :: ist=0
!> aeatom's error code when ist=1, boundstate's when ist>=2
   integer :: ierr=0
!> aeatom's iteration count when ist/=0
   integer :: it=0
!> the index in nat, lat of the level whose bound state failed when ist>=2
   integer :: kerr=0
!> the number of merged valence levels
   integer :: nvt=0
   integer :: nat(mxstate)=0  !< their principal quantum numbers
   integer :: lat(mxstate)=0  !< their angular momenta
!> their occupancies (:,1:3): reference, maximally-ionized, excited
   real(dp) :: fat(mxstate,nlev)=0.0d0
!> their all-electron eigenvalues (mxstate,nkap,nlev), the same
   real(dp), allocatable :: eat(:,:,:)
!> the all-electron total energy of the excited configuration
   real(dp) :: etst=0.0d0
 end type ae_config_t

!> a successful aeatom solution, with the inputs that determine it: sratom
!> and relatom overwrite the eigenvalues' guess, so the atom is a function
!> of these alone
 type :: ae_entry_t
!> 0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
   integer :: irel=0
   integer :: nkap=0  !< number of kappas per l: 1, or 2 for irel=2
   integer :: nc=0  !< number of core states
   integer :: ncv=0  !< number of core and valence states
   integer :: mmax=0  !< size of the log radial mesh
   integer :: iexc=0  !< exchange-correlation functional
   integer :: it=0  !< aeatom's iteration count
   real(dp) :: zz=0.0d0  !< atomic number
   real(dp) :: etot=0.0d0  !< total energy (Ha)
!> kinetic, electron-nucleus, Hartree and exchange-correlation energies (Ha)
   real(dp) :: eterm(4)=0.0d0
   real(dp) :: res=0.0d0  !< residual of the screening potential (ae_residual)
   integer, allocatable :: na(:)  !< principal quantum numbers (ncv)
   integer, allocatable :: la(:)  !< angular momenta (ncv)
   real(dp), allocatable :: fa(:)  !< occupancies (ncv)
   real(dp), allocatable :: rr(:)  !< log radial mesh (bohr)
   real(dp), allocatable :: ea(:,:)  !< eigenvalues (ncv,nkap) (Ha)
!> radius of the outermost peak of each wave function (ncv,nkap) (bohr)
   real(dp), allocatable :: rpk(:,:)
   real(dp), allocatable :: rhoc(:)  !< core charge density
   real(dp), allocatable :: rho(:)  !< total charge density
   real(dp), allocatable :: vi(:)  !< all-electron potential (Ha)
 end type ae_entry_t

!> the all-electron atoms of a run, the reference and two per test
!> configuration (aeatom_cached), which a run owns and passes down
 type :: ae_cache_t
   integer :: nentry=0  !< number of entries in use
   type(ae_entry_t) :: entry(2*mxcnf+1)  !< the solutions, first nentry
 end type ae_cache_t

 contains

 subroutine aeatom_cached(cache,irel,c_light,nkap,na,la,ea,fa,rpk,nc,ncv,it, &
&                         rhoc,rho,rr,vi,zz,mmax,iexc,etot,eterm,res,ierr, &
&                         iout,ifail)

!! aeatom, or the copy of its earlier solution with the same inputs. The
!! reference configuration recurs as test configuration 0, at both levels
!! of run_config_ae, and the excited level of a configuration that only
!! takes charge away is its intermediate one; each was solved again, to
!! the same bits. A failed solution is not kept, so that its warnings
!! recur, and a reused one stopped on the residual's round-off floor
!! repeats that warning. rpk comes back as that solution left it (rows
!! 1..ncv). cache holds the solutions of the run; the other arguments
!! are those of aeatom
!!
!! @param cache  the all-electron atoms of the run; a new solution added
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param na  principal quantum numbers
!! @param la  angular momenta
!! @param ea  eigenvalues (Ha), the guess on input and the solution's on
!!            output; ea(:,1) for j=l+1/2, ea(:,2) for j=l-1/2 when irel=2
!! @param fa  occupancies
!! @param rpk  radius of the outermost peak of each wave function (bohr),
!!             indexed as ea
!! @param nc  number of core states
!! @param ncv  number of core and valence states
!! @param it  number of iterations
!! @param rhoc  core charge density
!! @param rho  total charge density
!! @param rr  log radial mesh (bohr)
!! @param vi  all-electron potential (Ha), output
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param iexc  exchange-correlation functional
!! @param etot  total energy (Ha)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (Ha)
!! @param res  residual of the screening potential (ae_residual)
!! @param ierr  0, or non-zero after a failure
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use aeatom_m, only: aeatom
 use ae_scf_m, only: ae_floor_warning

!Input variables
 integer, intent(in) :: irel,nkap,mmax,iexc,nc,ncv
 real(dp), intent(in) :: c_light
 integer, intent(in) :: na(ncv),la(ncv)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(ncv),rr(mmax)
 integer, intent(in) :: iout

!Input/output variables
 type(ae_cache_t), intent(inout) :: cache
 integer, intent(inout) :: ierr
 real(dp), intent(inout) :: ea(mxstate,nkap),rpk(mxstate,nkap)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot,eterm(4),res
 real(dp), intent(out) :: rho(mmax),rhoc(mmax),vi(mmax)
 integer, intent(out) :: ifail

!Local variables
 integer :: ii
 real(dp), allocatable :: vfree(:)

 ifail=0
 do ii=1,cache%nentry
   associate(ce=>cache%entry(ii))
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
 call aeatom(irel,c_light,nkap,na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&            rr,vi,zz,mmax,iexc,etot,eterm,res,ierr,vfree,.false.,iout,ifail)
 if(ifail/=0) return
 deallocate(vfree)
 if(ierr/=0 .or. cache%nentry==size(cache%entry)) return

 cache%nentry=cache%nentry+1
 associate(ce=>cache%entry(cache%nentry))
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

 subroutine run_config_ae(cache,jj,cnf,nc,rr,zz,mmax,iexc,ea,irel,c_light, &
&                         nkap,aec,iout,ifail)

!! solves the all-electron atom in the maximally-ionized and excited states
!! of test configuration jj, the all-electron half of the test-configuration
!! calculations.  These depend only on the configurations, not on the
!! pseudopotential, so oncvpsp_ae.x can compute them once for oncvpsp_ps.x
!! (see ae_file_m).  run_config_ps does the pseudo-atom half.
!!
!! @param cache  the all-electron atoms of the run
!! @param jj  index of current configufation
!! @param cnf  the reference and test configurations (input_m)
!! @param nc  number of core states
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param iexc  exchange-correlation functional
!! @param ea  reference configurattion eigenvalues, ea(:,1) for j=l+1/2 and
!!            ea(:,2) for j=l-1/2 when irel=2
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param aec  the all-electron half of configuration jj
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use boundstate_m, only: boundstate
 use input_m, only: config_t
 use kappa_m, only: kap_count, kap_of

!Input variables
 integer, intent(in) :: jj,mmax,iexc,nc,irel,nkap
 real(dp), intent(in) :: c_light
 type(config_t), intent(in) :: cnf(mxcnf)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: ea(mxstate,nkap),rr(mmax)
 integer, intent(in) :: iout

!Input/output variables
 type(ae_cache_t), intent(inout) :: cache

!Output variables
 type(ae_config_t), intent(out) :: aec
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

 aec%ist=0
 aec%ierr=0
 aec%it=0
 aec%kerr=0

 it=0

 fat(:,:)=0.0d0
 nat(:)=0
 lat(:)=0

! set up cross indces of reference (r) and excited (e) cofigurations
 indxr(:)=0
 indxe(:)=0
 do ii=1,nc+cnf(1)%nv
   do kk=1,nc+cnf(jj)%nv
     if(cnf(jj)%na(kk)==cnf(1)%na(ii) .and. &
&       cnf(jj)%la(kk)==cnf(1)%la(ii)) then
       indxr(ii)=kk
       indxe(kk)=ii
     end if
   end do
 end do

! set up merged list with reference levels first
! levels not present in one confguration or the other get zero occupancy
 eat(:,:,:)=0.0d0
 do ii=1,nc+cnf(1)%nv
   nat(ii)=cnf(1)%na(ii)
   lat(ii)=cnf(1)%la(ii)
   fat(ii,1)=cnf(1)%fa(ii)
   eat(ii,:,1)=ea(ii,:)
   if(indxr(ii)>0) then
     fat(ii,2)=cnf(jj)%fa(indxr(ii))
   else
     fat(ii,2)=0.0d0
   end if
   fat(ii,3)=fat(ii,2)
 end do

 nvt=cnf(1)%nv
 do kk=1,nc+cnf(jj)%nv
   if(indxe(kk)==0) then
     nvt=nvt+1
     nat(nc+nvt)=cnf(jj)%na(kk)
     lat(nc+nvt)=cnf(jj)%la(kk)
     fat(nc+nvt,1)=0.0d0
     fat(nc+nvt,3)=cnf(jj)%fa(kk)
   end if
 end do

!set up array for maximally ionized intermediate state
 do kk=1,nc+nvt
   fat(kk,2)=min(fat(kk,1),fat(kk,3))
 end do



 eat(:,:,2)=0.0d0; rpk(:,:)=0.0d0; rhoc(:)=0.0d0; rho(:)=0.0d0; vfull(:)=0.0d0
 eaetst=0.0d0

!all-electron atom solution for maximally-ionized state

 call aeatom_cached(cache,irel,c_light,nkap,nat,lat,eat(1,1,2),fat(1,2),rpk, &
&                   nc,nc+nvt,it,rhoc,rho,rr,vfull,zz,mmax,iexc,eaetst,eterm, &
&                   res,ierr,iout,ifail)
 if(ifail/=0) return
 if(ierr/=0) then
   call store_results
   aec%ist=1
   aec%ierr=ierr
   aec%it=it
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
       call boundstate(irel,c_light,nat(kk),ll,kap,ierr,et, &
&                      rr,zz,vfull,uu,up,mmax,mch,iout,ifail)
       if(ifail/=0) return
!a negative ierr (no classical turning point) is accepted for j=l-1/2
       if(ierr>0 .or. (ierr/=0 .and. ikap==1)) then
         call store_results
         aec%ist=1+ikap
         aec%ierr=ierr
         aec%kerr=kk
         aec%it=it
         deallocate(rho,rhoc,vfull)
         deallocate(uu,up)
         return
       end if
       eat(kk,ikap,2)=et
     end do
     if(nkap==2 .and. ll==0) eat(kk,2,2)=eat(kk,1,2)
   end if
 end do

 eat(:,:,3)=0.0d0; rpk(:,:)=0.0d0; rhoc(:)=0.0d0; rho(:)=0.0d0; vfull(:)=0.0d0
 eaetst=0.0d0

!all-electron atom solution for excited state

 call aeatom_cached(cache,irel,c_light,nkap,nat,lat,eat(1,1,3),fat(1,3),rpk, &
&                   nc,nc+nvt,it,rhoc,rho,rr,vfull,zz,mmax,iexc,eaetst,eterm, &
&                   res,ierr,iout,ifail)
 if(ifail/=0) return
 if(ierr/=0) then
   call store_results
   aec%ist=1
   aec%ierr=ierr
   aec%it=it
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
!! copies the results so far into aec
 aec%nvt=nvt
 aec%nat(:)=nat(:)
 aec%lat(:)=lat(:)
 aec%fat(:,:)=fat(:,:)
 aec%eat=eat
 aec%etst=eaetst
 return
 end subroutine store_results

 end subroutine run_config_ae

 subroutine config_ae_warning(jj,irel,aec,iout)
!! prints the warning run_config_ps gives when the all-electron atom for test
!! configuration jj failed
!!
!! @param jj  index of the configuration (1 for the reference)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param aec  its all-electron half (run_config_ae)
!! @param iout  unit of the text output

!Input variables
 integer, intent(in) :: jj,irel
 type(ae_config_t), intent(in) :: aec
 integer, intent(in) :: iout

!Local variables
 integer :: ll

 if(aec%ist==1) then
   if(irel==2) then
     write(iout,'(a/a,i2)') 'run_config_r: WARNING  for AE atom,', &
&          ' WARNING no output for configuration',jj
   else
     write(iout,'(a/a,i2)') 'run_config: WARNING  for AE atom,', &
&          ' WARNING no output for configuration',jj
   end if
   if(aec%ierr==-1 .and. irel==2) then
     write(iout,'(a,i4)') &
&          'run_config_r: WARNING no classical turning point ' // &
&          'error, iteration', &
&          aec%it
   else if(aec%ierr==-1) then
     write(iout,'(a,i4)') 'no classical turning point error, iteration', &
&          aec%it
   else if(irel==2) then
     write(iout,'(a)') &
&          'run_config_r: WARNING self-consistency failed to converge'
   else
     write(iout,'(a)') 'run_config: WARNING self-consistency failed to converge'
   end if
 else if(aec%ist>=2) then
   ll=aec%lat(aec%kerr)
   if(irel/=2) then
     write(iout,'(/a,3i4)') &
&         'runconfig: WARNING lschfb convergence ERROR n,l,iter=', &
&         aec%nat(aec%kerr),ll,aec%it
   else if(aec%ist==2) then
     write(iout,'(/a,4i4)') &
&         'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,ierr=', &
&         aec%nat(aec%kerr),ll,-(ll+1),aec%ierr
   else
     write(iout,'(/a,4i4)') &
&         'runconfig: WARNING ldiracfb convergence ERROR n,l,kap,iter=', &
&         aec%nat(aec%kerr),ll,ll,aec%it
   end if
 end if
 return
 end subroutine config_ae_warning

 end module config_ae_m
