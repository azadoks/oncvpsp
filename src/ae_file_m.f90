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
! all-electron data file passed from oncvpsp_ae.x to oncvpsp_ps.x
! (oncvpspnr_ae.x to oncvpspnr_ps.x, oncvpspr_ae.x to oncvpspr_ps.x)
!
! oncvpsp_ae.x solves the reference all-electron atom and the all-electron
! halves of the test configurations, which do not depend on any
! pseudopotential input, and writes them here.  oncvpsp_ps.x reads them
! instead of solving the all-electron atom again, so a loop that only
! changes pseudopotential inputs solves the all-electron atom once.
!
! The file is unformatted stream I/O in the native byte order:
!   magic, format version, irel                (0 non-, 1 scalar-,
!                                               2 fully-relativistic)
!   atsym, zz, nc, nv, iexc, na, la, fa         (checked against the input)
!   mmax, rr                                   (log radial mesh, checked
!                                               against the input mesh)
!   it, etot, ea, rhoc, rho, vfull             (reference atom, sratom or
!                                               relatom; ea(:,1:2) if irel=2)
!   ncfg, then for each test configuration:
!     nvcnf, nacnf, lacnf, facnf               (key, matched to the input)
!     istcnf through etstcnf for it            (run_config_ae[_r] results)

 module ae_file_m
 use, intrinsic :: iso_fortran_env, only: int8
 use constants_m, only: dp
 implicit none
 private
 public :: write_ae_file, read_ae_file
 public :: write_ae_file_r, read_ae_file_r

 character(len=8), parameter :: magic='ONCVAE  '
 integer, parameter :: version=1
 integer, parameter :: irel_fr=2
!relative tolerance for comparing real input (zz, occupations) and the
!radial mesh with the file: both programs parse the same input text, but
!may be built by different compilers or with different flags
 real(dp), parameter :: rtol=1.0d-12

 contains

 subroutine write_ae_file(fname,srel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                         mmax,rr,it,etot,ea,rhoc,rho,vfull, &
&                         ncnf,nvcnf,nacnf,lacnf,facnf,istcnf,ierrcnf,itcnf, &
&                         kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)

! writes the all-electron data file: the sratom results for the reference
! atom and the run_config_ae results for all configurations

!fname  file to write
!srel  .true. for scalar-relativistic, .false. for non-relativistic
!atsym  atomic symbol
!zz  atomic number
!nc  number of core states
!nv  number of valence states
!iexc  exchange-correlation function to be used
!na  principal quantum numbers, reference configuration
!la  angular momenta, reference configuration
!fa  occupancies, reference configuration
!mmax  size of log grid
!rr  log radial mesh
!it  sratom iteration count
!etot  all-electron total energy (Ha)
!ea  all-electron eigenvalues (Ha)
!rhoc  all-electron core charge density
!rho  all-electron valence charge density
!vfull  all-electron potential
!ncnf  number of test configurations (not counting the reference)
!nvcnf  number of valence states, all configurations
!nacnf  principal quantum numbers, all configurations
!lacnf  angular momenta, all configurations
!facnf  occupancies, all configurations
!istcnf through etstcnf  run_config_ae results, all configurations

!Input variables
 character(len=*), intent(in) :: fname
 logical, intent(in) :: srel
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax,it,ncnf
 integer, intent(in) :: na(30),la(30),nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp), intent(in) :: zz,etot
 real(dp), intent(in) :: fa(30),ea(30),facnf(30,5)
 real(dp), intent(in) :: rr(mmax),rhoc(mmax),rho(mmax),vfull(mmax)
 integer, intent(in) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(in) :: nvtcnf(5),natcnf(30,5),latcnf(30,5)
 real(dp), intent(in) :: fatcnf(30,3,5),eatcnf(30,3,5),etstcnf(5)

!Local variables
 integer :: unit,jj

 call write_head(fname,unit,merge(1,0,srel),atsym,zz,nc,nv,iexc,na,la,fa, &
&                mmax,rr)
 write(unit) it,etot,ea(1:nc+nv)
 write(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)
 write(unit) ncnf+1
 do jj=1,ncnf+1
   call write_config_key(unit,jj,nc,nvcnf,nacnf,lacnf,facnf)
   write(unit) istcnf(jj),ierrcnf(jj),itcnf(jj),kerrcnf(jj), &
&              nvtcnf(jj),natcnf(:,jj),latcnf(:,jj),fatcnf(:,:,jj), &
&              eatcnf(:,:,jj),etstcnf(jj)
 end do

 close(unit)
 return
 end subroutine write_ae_file

 subroutine write_ae_file_r(fname,atsym,zz,nc,nv,iexc,na,la,fa, &
&                           mmax,rr,it,etot,ea,rhoc,rho,vfull, &
&                           ncnf,nvcnf,nacnf,lacnf,facnf,istcnf,ierrcnf,itcnf, &
&                           kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)

! fully-relativistic write_ae_file: ea holds both kap, and the results
! come from run_config_ae

!fname  file to write
!atsym  atomic symbol
!zz  atomic number
!nc  number of core states
!nv  number of valence states
!iexc  exchange-correlation function to be used
!na  principal quantum numbers, reference configuration
!la  angular momenta, reference configuration
!fa  occupancies, reference configuration
!mmax  size of log grid
!rr  log radial mesh
!it  relatom iteration count
!etot  all-electron total energy (Ha)
!ea  all-electron eigenvalues (Ha), kap=-(l+1) and kap=l
!rhoc  all-electron core charge density
!rho  all-electron valence charge density
!vfull  all-electron potential
!ncnf  number of test configurations (not counting the reference)
!nvcnf  number of valence states, all configurations
!nacnf  principal quantum numbers, all configurations
!lacnf  angular momenta, all configurations
!facnf  occupancies, all configurations
!istcnf through etstcnf  run_config_ae results, all configurations

!Input variables
 character(len=*), intent(in) :: fname
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax,it,ncnf
 integer, intent(in) :: na(30),la(30),nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp), intent(in) :: zz,etot
 real(dp), intent(in) :: fa(30),ea(30,2),facnf(30,5)
 real(dp), intent(in) :: rr(mmax),rhoc(mmax),rho(mmax),vfull(mmax)
 integer, intent(in) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(in) :: nvtcnf(5),natcnf(30,5),latcnf(30,5)
 real(dp), intent(in) :: fatcnf(30,3,5),eatcnf(30,2,3,5),etstcnf(5)

!Local variables
 integer :: unit,jj

 call write_head(fname,unit,irel_fr,atsym,zz,nc,nv,iexc,na,la,fa, &
&                mmax,rr)
 write(unit) it,etot,ea(1:nc+nv,1:2)
 write(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)
 write(unit) ncnf+1
 do jj=1,ncnf+1
   call write_config_key(unit,jj,nc,nvcnf,nacnf,lacnf,facnf)
   write(unit) istcnf(jj),ierrcnf(jj),itcnf(jj),kerrcnf(jj), &
&              nvtcnf(jj),natcnf(:,jj),latcnf(:,jj),fatcnf(:,:,jj), &
&              eatcnf(:,:,:,jj),etstcnf(jj)
 end do

 close(unit)
 return
 end subroutine write_ae_file_r

 subroutine read_ae_file(fname,srel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                        mmax,rr,ncnf,nvcnf,nacnf,lacnf,facnf, &
&                        it,etot,ea,rhoc,rho,vfull,istcnf,ierrcnf,itcnf, &
&                        kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)

! reads the all-electron data file: the reference atom and the
! run_config_ae results for each configuration of this run, found in the
! file by matching the configuration itself, so the file may hold more
! test configurations than the input, in any order.  Stops if the file
! was written for a different atom or relativistic treatment, or lacks a
! configuration.  srel through facnf are this run's input, checked
! against the file; it through etstcnf are read from it.

!fname  file to read
!srel  .true. for scalar-relativistic, .false. for non-relativistic
!atsym  atomic symbol
!zz  atomic number
!nc  number of core states
!nv  number of valence states
!iexc  exchange-correlation function to be used
!na  principal quantum numbers, reference configuration
!la  angular momenta, reference configuration
!fa  occupancies, reference configuration
!mmax  size of log grid
!rr  log radial mesh
!ncnf  number of test configurations (not counting the reference)
!nvcnf  number of valence states, all configurations
!nacnf  principal quantum numbers, all configurations
!lacnf  angular momenta, all configurations
!facnf  occupancies, all configurations
!it  sratom iteration count
!etot  all-electron total energy (Ha)
!ea  all-electron eigenvalues (Ha)
!rhoc  all-electron core charge density
!rho  all-electron valence charge density
!vfull  all-electron potential
!istcnf through etstcnf  run_config_ae results, all configurations

!Input variables
 character(len=*), intent(in) :: fname
 logical, intent(in) :: srel
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax,ncnf
 integer, intent(in) :: na(30),la(30),nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(30),rr(mmax),facnf(30,5)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot
 real(dp), intent(out) :: ea(30),rhoc(mmax),rho(mmax),vfull(mmax)
 integer, intent(out) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(out) :: nvtcnf(5),natcnf(30,5),latcnf(30,5)
 real(dp), intent(out) :: fatcnf(30,3,5),eatcnf(30,3,5),etstcnf(5)

!Local variables
 integer :: unit,jj,kk,ncfg
 integer :: ist,ierr,itc,kerr,nvt,nat(30),lat(30)
 real(dp) :: eaetst,fat(30,3),eat(30,3)
 logical :: found(5),match(5)

 call read_head(fname,unit,merge(1,0,srel),atsym,zz,nc,nv,iexc,na,la,fa, &
&               mmax,rr)
 read(unit) it,etot,ea(1:nc+nv)
 read(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)

! copy the results of each configuration in the file to the configurations
! of this run that it matches
 found(:)=.false.
 read(unit) ncfg
 do kk=1,ncfg
   call match_config(unit,nc,ncnf,nvcnf,nacnf,lacnf,facnf,found,match)
   read(unit) ist,ierr,itc,kerr,nvt,nat,lat,fat,eat,eaetst
   do jj=1,ncnf+1
     if(.not. match(jj)) cycle
     istcnf(jj)=ist
     ierrcnf(jj)=ierr
     itcnf(jj)=itc
     kerrcnf(jj)=kerr
     nvtcnf(jj)=nvt
     natcnf(:,jj)=nat(:)
     latcnf(:,jj)=lat(:)
     fatcnf(:,:,jj)=fat(:,:)
     eatcnf(:,:,jj)=eat(:,:)
     etstcnf(jj)=eaetst
     found(jj)=.true.
   end do
 end do
 call check_found(ncnf,found)

 close(unit)
 return
 end subroutine read_ae_file

 subroutine read_ae_file_r(fname,atsym,zz,nc,nv,iexc,na,la,fa, &
&                          mmax,rr,ncnf,nvcnf,nacnf,lacnf,facnf, &
&                          it,etot,ea,rhoc,rho,vfull,istcnf,ierrcnf,itcnf, &
&                          kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)

! fully-relativistic read_ae_file.  atsym through facnf are this run's
! input, checked against the file; it through etstcnf are read from it.

!fname  file to read
!atsym  atomic symbol
!zz  atomic number
!nc  number of core states
!nv  number of valence states
!iexc  exchange-correlation function to be used
!na  principal quantum numbers, reference configuration
!la  angular momenta, reference configuration
!fa  occupancies, reference configuration
!mmax  size of log grid
!rr  log radial mesh
!ncnf  number of test configurations (not counting the reference)
!nvcnf  number of valence states, all configurations
!nacnf  principal quantum numbers, all configurations
!lacnf  angular momenta, all configurations
!facnf  occupancies, all configurations
!it  relatom iteration count
!etot  all-electron total energy (Ha)
!ea  all-electron eigenvalues (Ha), kap=-(l+1) and kap=l
!rhoc  all-electron core charge density
!rho  all-electron valence charge density
!vfull  all-electron potential
!istcnf through etstcnf  run_config_ae results, all configurations

!Input variables
 character(len=*), intent(in) :: fname
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax,ncnf
 integer, intent(in) :: na(30),la(30),nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(30),rr(mmax),facnf(30,5)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot
 real(dp), intent(out) :: ea(30,2),rhoc(mmax),rho(mmax),vfull(mmax)
 integer, intent(out) :: istcnf(5),ierrcnf(5),itcnf(5),kerrcnf(5)
 integer, intent(out) :: nvtcnf(5),natcnf(30,5),latcnf(30,5)
 real(dp), intent(out) :: fatcnf(30,3,5),eatcnf(30,2,3,5),etstcnf(5)

!Local variables
 integer :: unit,jj,kk,ncfg
 integer :: ist,ierr,itc,kerr,nvt,nat(30),lat(30)
 real(dp) :: eaetst,fat(30,3),eat(30,2,3)
 logical :: found(5),match(5)

 call read_head(fname,unit,irel_fr,atsym,zz,nc,nv,iexc,na,la,fa, &
&               mmax,rr)
 read(unit) it,etot,ea(1:nc+nv,1:2)
 read(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)

! copy the results of each configuration in the file to the configurations
! of this run that it matches
 found(:)=.false.
 read(unit) ncfg
 do kk=1,ncfg
   call match_config(unit,nc,ncnf,nvcnf,nacnf,lacnf,facnf,found,match)
   read(unit) ist,ierr,itc,kerr,nvt,nat,lat,fat,eat,eaetst
   do jj=1,ncnf+1
     if(.not. match(jj)) cycle
     istcnf(jj)=ist
     ierrcnf(jj)=ierr
     itcnf(jj)=itc
     kerrcnf(jj)=kerr
     nvtcnf(jj)=nvt
     natcnf(:,jj)=nat(:)
     latcnf(:,jj)=lat(:)
     fatcnf(:,:,jj)=fat(:,:)
     eatcnf(:,:,:,jj)=eat(:,:,:)
     etstcnf(jj)=eaetst
     found(jj)=.true.
   end do
 end do
 call check_found(ncnf,found)

 close(unit)
 return
 end subroutine read_ae_file_r

 subroutine write_head(fname,unit,irel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                      mmax,rr)
! opens fname and writes everything before the reference atom results

!fname  file name
!unit  unit number of the opened file
!irel  0 non-, 1 scalar-, 2 fully-relativistic
!atsym  atomic symbol
!zz  atomic number
!nc  number of core states
!nv  number of valence states
!iexc  exchange-correlation function to be used
!na  principal quantum numbers, reference configuration
!la  angular momenta, reference configuration
!fa  occupancies, reference configuration
!mmax  size of log grid
!rr  log radial mesh

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: irel
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax
 integer, intent(in) :: na(30),la(30)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(30),rr(mmax)

!Output variables
 integer, intent(out) :: unit

!Local variables
 integer :: ios

 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='replace',action='write',iostat=ios)
 if(ios/=0) then
   write(6,'(2a)') 'write_ae_file: ERROR cannot open ',trim(fname)
   error stop 1
 end if

 write(unit) magic,version,irel
 write(unit) atsym,zz,nc,nv,iexc
 write(unit) na(1:nc+nv),la(1:nc+nv),fa(1:nc+nv)
 write(unit) mmax,rr(1:mmax)
 return
 end subroutine write_head

 subroutine read_head(fname,unit,irel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                     mmax,rr)
! opens fname, checks it against this run's input and reads rr; stops if
! the file was written for a different atom or relativistic treatment

!fname  file name
!unit  unit number of the opened file
!irel  expected treatment: 0 non-, 1 scalar-, 2 fully-relativistic
!atsym  atomic symbol
!zz  atomic number
!nc  number of core states
!nv  number of valence states
!iexc  exchange-correlation function to be used
!na  principal quantum numbers, reference configuration
!la  angular momenta, reference configuration
!fa  occupancies, reference configuration
!mmax  size of log grid
!rr  log radial mesh

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: irel
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax
 integer, intent(in) :: na(30),la(30)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(30),rr(mmax)

!Output variables
 integer, intent(out) :: unit

!Local variables
 character(len=8) :: magic_f
 character(len=2) :: atsym_f
 integer :: ios,version_f,irel_f,nc_f,nv_f,iexc_f,mmax_f
 integer :: na_f(30),la_f(30)
 real(dp) :: zz_f,fa_f(30)
 real(dp), allocatable :: rr_f(:)

 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='old',action='read',iostat=ios)
 if(ios/=0) then
   write(6,'(2a)') 'read_ae_file: ERROR cannot open ',trim(fname)
   error stop 1
 end if

 read(unit,iostat=ios) magic_f,version_f,irel_f
 if(ios/=0 .or. magic_f/=magic) then
   write(6,'(2a)') 'read_ae_file: ERROR not an all-electron data file: ', &
&        trim(fname)
   error stop 1
 end if
!the magic string reads the same in either byte order, so a file from a
!machine of the other byte order shows up as a byte-swapped version
 if(version_f/=version .and. byteswap(version_f)==version) then
   write(6,'(2a)') 'read_ae_file: ERROR file was written on a machine', &
&        ' with the other byte order'
   error stop 1
 end if
 if(version_f/=version) then
   write(6,'(a,i4,a,i4)') 'read_ae_file: ERROR file format version', &
&        version_f,' but this program reads version',version
   error stop 1
 end if
 if(irel_f/=irel) then
   write(6,'(2a)') 'read_ae_file: ERROR file was written with a different', &
&        ' relativistic treatment'
   error stop 1
 end if

 read(unit) atsym_f,zz_f,nc_f,nv_f,iexc_f
 if(atsym_f/=atsym .or. differ(zz_f,zz) .or. nc_f/=nc .or. nv_f/=nv &
&   .or. iexc_f/=iexc) then
   write(6,'(2a)') 'read_ae_file: ERROR atsym, z, nc, nv or iexc differ', &
&        ' from the input'
   error stop 1
 end if
 read(unit) na_f(1:nc+nv),la_f(1:nc+nv),fa_f(1:nc+nv)
 if(any(na_f(1:nc+nv)/=na(1:nc+nv)) .or. any(la_f(1:nc+nv)/=la(1:nc+nv)) &
&   .or. any(differ(fa_f(1:nc+nv),fa(1:nc+nv)))) then
   write(6,'(2a)') 'read_ae_file: ERROR reference configuration differs', &
&        ' from the input'
   error stop 1
 end if

 read(unit) mmax_f
 if(mmax_f/=mmax) then
   write(6,'(a,2i8)') 'read_ae_file: ERROR radial mesh size differs', &
&        mmax_f,mmax
   error stop 1
 end if
 allocate(rr_f(mmax))
 read(unit) rr_f(1:mmax)
 if(any(abs(rr_f(1:mmax)-rr(1:mmax))>rtol*rr(1:mmax))) then
   write(6,'(2a)') 'read_ae_file: ERROR radial mesh differs from the', &
&        ' one this program builds'
   error stop 1
 end if
 return
 end subroutine read_head

 subroutine write_config_key(unit,jj,nc,nvcnf,nacnf,lacnf,facnf)
! writes configuration jj, the key by which match_config matches it

!unit  unit number of the open file
!jj  index of the configuration in this run (1 for the reference)
!nc  number of core states
!nvcnf  number of valence states, all configurations
!nacnf  principal quantum numbers, all configurations
!lacnf  angular momenta, all configurations
!facnf  occupancies, all configurations

!Input variables
 integer, intent(in) :: unit,jj,nc
 integer, intent(in) :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp), intent(in) :: facnf(30,5)

!Local variables
 integer :: nt

 nt=nc+nvcnf(jj)
 write(unit) nvcnf(jj),nacnf(1:nt,jj),lacnf(1:nt,jj),facnf(1:nt,jj)
 return
 end subroutine write_config_key

 subroutine match_config(unit,nc,ncnf,nvcnf,nacnf,lacnf,facnf,found,match)
! reads the key of the next configuration in the file and flags the
! configurations of this run, not found before, that it matches

!unit  unit number of the open file
!nc  number of core states
!ncnf  number of test configurations (not counting the reference)
!nvcnf  number of valence states, all configurations
!nacnf  principal quantum numbers, all configurations
!lacnf  angular momenta, all configurations
!facnf  occupancies, all configurations
!found  .true. for the configurations already found in the file
!match  .true. for the configurations this key matches

!Input variables
 integer, intent(in) :: unit,nc,ncnf
 integer, intent(in) :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp), intent(in) :: facnf(30,5)
 logical, intent(in) :: found(5)

!Output variables
 logical, intent(out) :: match(5)

!Local variables
 integer :: jj,nt,nvcnf_f
 integer :: nacnf_f(30),lacnf_f(30)
 real(dp) :: facnf_f(30)

 read(unit) nvcnf_f
 nt=nc+nvcnf_f
 read(unit) nacnf_f(1:nt),lacnf_f(1:nt),facnf_f(1:nt)

 match(:)=.false.
 do jj=1,ncnf+1
   if(found(jj)) cycle
   if(nvcnf_f/=nvcnf(jj)) cycle
   if(any(nacnf_f(1:nt)/=nacnf(1:nt,jj))) cycle
   if(any(lacnf_f(1:nt)/=lacnf(1:nt,jj))) cycle
   if(any(differ(facnf_f(1:nt),facnf(1:nt,jj)))) cycle
   match(jj)=.true.
 end do
 return
 end subroutine match_config

 subroutine check_found(ncnf,found)
! stops if a configuration of this run is not in the file

!ncnf  number of test configurations (not counting the reference)
!found  .true. for the configurations found in the file

!Input variables
 integer, intent(in) :: ncnf
 logical, intent(in) :: found(5)

!Local variables
 integer :: jj

 do jj=1,ncnf+1
   if(.not. found(jj)) then
     write(6,'(a,i2,a)') 'read_ae_file: ERROR test configuration',jj-1, &
&         ' is not in the all-electron data file'
     error stop 1
   end if
 end do
 return
 end subroutine check_found

 elemental logical function differ(aa,bb)
! .true. if aa and bb differ by more than rtol relative to the larger,
! or absolutely for values below 1 (occupations of 0); not for the mesh,
! whose points near the origin are far below 1

!aa  first value
!bb  second value

!Input variables
 real(dp), intent(in) :: aa,bb
 differ=abs(aa-bb)>rtol*max(1.0d0,abs(aa),abs(bb))
 return
 end function differ

 integer function byteswap(ii)
! ii with its bytes in reverse order

!ii  integer to byte-swap

!Input variables
 integer, intent(in) :: ii

!Local variables
 integer(int8) :: bytes(storage_size(ii)/8)
 bytes=transfer(ii,bytes)
 byteswap=transfer(bytes(size(bytes):1:-1),ii)
 return
 end function byteswap

 end module ae_file_m
