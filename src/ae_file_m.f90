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
!   mmax, rr                                   (log radial mesh)
!   it, etot, ea, rhoc, rho, vfull             (reference atom, sratom or
!                                               relatom; ea(:,1:2) if irel=2)
!   ncfg, then for each test configuration:
!     nvcnf, nacnf, lacnf, facnf               (key, matched to the input)
!     config_ae_t or config_ae_r_t fields      (run_config_ae[_r] results)

 module ae_file_m
 use config_ae_m, only: config_ae_t, config_ae_r_t
 implicit none
 private
 public :: write_ae_file, read_ae_file, ae_file_config
 public :: write_ae_file_r, read_ae_file_r, ae_file_config_r

 integer, parameter :: dp=kind(1.0d0)
 character(len=8), parameter :: magic='ONCVAE  '
 integer, parameter :: version=1
 integer, parameter :: irel_fr=2

!test configurations read from the file
 integer :: ncfg_file=0
 integer :: nvcnf_file(5)
 integer :: nacnf_file(30,5),lacnf_file(30,5)
 real(dp) :: facnf_file(30,5)
 type(config_ae_t) :: cae_file(5)
 type(config_ae_r_t) :: cae_r_file(5)

 contains

 subroutine write_ae_file(fname,srel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                         mmax,rr,it,etot,ea,rhoc,rho,vfull, &
&                         ncnf,nvcnf,nacnf,lacnf,facnf,cae)

!fname  file to write
!srel .true. for scalar-relativistic, .false. for non-relativistic
!atsym, zz, nc, nv, iexc, na, la, fa  reference atom input
!mmax  size of log grid
!rr  log radial mesh
!it, etot, ea, rhoc, rho, vfull  sratom results for the reference atom
!ncnf  number of test configurations (not counting the reference)
!nvcnf, nacnf, lacnf, facnf  all configurations
!cae  run_config_ae results for all configurations

!Input variables
 character(len=*) :: fname
 logical :: srel
 character*2 :: atsym
 integer :: nc,nv,iexc,mmax,it,ncnf
 integer :: na(30),la(30),nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp) :: zz,etot
 real(dp) :: fa(30),ea(30),facnf(30,5)
 real(dp) :: rr(mmax),rhoc(mmax),rho(mmax),vfull(mmax)
 type(config_ae_t) :: cae(5)

!Local variables
 integer :: unit,jj

 call write_head(fname,unit,merge(1,0,srel),atsym,zz,nc,nv,iexc,na,la,fa, &
&                mmax,rr)
 write(unit) it,etot,ea(1:nc+nv)
 write(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)
 write(unit) ncnf+1
 do jj=1,ncnf+1
  call write_config_key(unit,jj,nc,nvcnf,nacnf,lacnf,facnf)
  write(unit) cae(jj)%stat,cae(jj)%ierr,cae(jj)%it,cae(jj)%kerr, &
&             cae(jj)%nvt,cae(jj)%nat,cae(jj)%lat,cae(jj)%fat, &
&             cae(jj)%eat,cae(jj)%eaetst
 end do

 close(unit)
 return
 end subroutine write_ae_file

 subroutine write_ae_file_r(fname,atsym,zz,nc,nv,iexc,na,la,fa, &
&                           mmax,rr,it,etot,ea,rhoc,rho,vfull, &
&                           ncnf,nvcnf,nacnf,lacnf,facnf,cae)

! fully-relativistic write_ae_file: ea holds both kap, and cae comes from
! run_config_ae_r

!Input variables
 character(len=*) :: fname
 character*2 :: atsym
 integer :: nc,nv,iexc,mmax,it,ncnf
 integer :: na(30),la(30),nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp) :: zz,etot
 real(dp) :: fa(30),ea(30,2),facnf(30,5)
 real(dp) :: rr(mmax),rhoc(mmax),rho(mmax),vfull(mmax)
 type(config_ae_r_t) :: cae(5)

!Local variables
 integer :: unit,jj

 call write_head(fname,unit,irel_fr,atsym,zz,nc,nv,iexc,na,la,fa, &
&                mmax,rr)
 write(unit) it,etot,ea(1:nc+nv,1:2)
 write(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)
 write(unit) ncnf+1
 do jj=1,ncnf+1
  call write_config_key(unit,jj,nc,nvcnf,nacnf,lacnf,facnf)
  write(unit) cae(jj)%stat,cae(jj)%ierr,cae(jj)%it,cae(jj)%kerr, &
&             cae(jj)%nvt,cae(jj)%nat,cae(jj)%lat,cae(jj)%fat, &
&             cae(jj)%eat,cae(jj)%eaetst
 end do

 close(unit)
 return
 end subroutine write_ae_file_r

 subroutine read_ae_file(fname,srel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                        mmax,rr,it,etot,ea,rhoc,rho,vfull)

! reads the reference atom and keeps the test configurations for
! ae_file_config; stops if the file was written for a different atom
! or relativistic treatment

!fname  file to read
!srel, atsym, zz, nc, nv, iexc, na, la, fa  this run's input, checked
!                                            against the file
!mmax  size of log grid, checked against the file

!Output variables
!rr  log radial mesh
!it, etot, ea, rhoc, rho, vfull  sratom results for the reference atom

 character(len=*) :: fname
 logical :: srel
 character*2 :: atsym
 integer :: nc,nv,iexc,mmax,it
 integer :: na(30),la(30)
 real(dp) :: zz,etot
 real(dp) :: fa(30),ea(30)
 real(dp) :: rr(mmax),rhoc(mmax),rho(mmax),vfull(mmax)

!Local variables
 integer :: unit,jj

 call read_head(fname,unit,merge(1,0,srel),atsym,zz,nc,nv,iexc,na,la,fa, &
&               mmax,rr)
 read(unit) it,etot,ea(1:nc+nv)
 read(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)

 read(unit) ncfg_file
 do jj=1,ncfg_file
  call read_config_key(unit,jj,nc)
  read(unit) cae_file(jj)%stat,cae_file(jj)%ierr,cae_file(jj)%it, &
&            cae_file(jj)%kerr,cae_file(jj)%nvt,cae_file(jj)%nat, &
&            cae_file(jj)%lat,cae_file(jj)%fat,cae_file(jj)%eat, &
&            cae_file(jj)%eaetst
 end do

 close(unit)
 return
 end subroutine read_ae_file

 subroutine read_ae_file_r(fname,atsym,zz,nc,nv,iexc,na,la,fa, &
&                          mmax,rr,it,etot,ea,rhoc,rho,vfull)

! fully-relativistic read_ae_file

 character(len=*) :: fname
 character*2 :: atsym
 integer :: nc,nv,iexc,mmax,it
 integer :: na(30),la(30)
 real(dp) :: zz,etot
 real(dp) :: fa(30),ea(30,2)
 real(dp) :: rr(mmax),rhoc(mmax),rho(mmax),vfull(mmax)

!Local variables
 integer :: unit,jj

 call read_head(fname,unit,irel_fr,atsym,zz,nc,nv,iexc,na,la,fa, &
&               mmax,rr)
 read(unit) it,etot,ea(1:nc+nv,1:2)
 read(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)

 read(unit) ncfg_file
 do jj=1,ncfg_file
  call read_config_key(unit,jj,nc)
  read(unit) cae_r_file(jj)%stat,cae_r_file(jj)%ierr,cae_r_file(jj)%it, &
&            cae_r_file(jj)%kerr,cae_r_file(jj)%nvt,cae_r_file(jj)%nat, &
&            cae_r_file(jj)%lat,cae_r_file(jj)%fat,cae_r_file(jj)%eat, &
&            cae_r_file(jj)%eaetst
 end do

 close(unit)
 return
 end subroutine read_ae_file_r

 subroutine ae_file_config(jj,nc,nvcnf,nacnf,lacnf,facnf,cae)

! run_config_ae results for test configuration jj of this run, found in
! the file by matching the configuration itself, so the file may hold
! more test configurations than the input, in any order; stops if the
! configuration is not in the file

!Input variables
 integer :: jj,nc
 integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp) :: facnf(30,5)

!Output variables
 type(config_ae_t), intent(out) :: cae

 cae=cae_file(find_config(jj,nc,nvcnf,nacnf,lacnf,facnf))
 end subroutine ae_file_config

 subroutine ae_file_config_r(jj,nc,nvcnf,nacnf,lacnf,facnf,cae)

! fully-relativistic ae_file_config

!Input variables
 integer :: jj,nc
 integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp) :: facnf(30,5)

!Output variables
 type(config_ae_r_t), intent(out) :: cae

 cae=cae_r_file(find_config(jj,nc,nvcnf,nacnf,lacnf,facnf))
 end subroutine ae_file_config_r

 subroutine write_head(fname,unit,irel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                      mmax,rr)
! opens fname and writes everything before the reference atom results

 character(len=*) :: fname
 integer :: unit,irel
 character*2 :: atsym
 integer :: nc,nv,iexc,mmax
 integer :: na(30),la(30)
 real(dp) :: zz
 real(dp) :: fa(30),rr(mmax)

 integer :: ios

 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='replace',action='write',iostat=ios)
 if(ios/=0) then
  write(6,'(2a)') 'write_ae_file: ERROR cannot open ',trim(fname)
  stop 1
 end if

 write(unit) magic,version,irel
 write(unit) atsym,zz,nc,nv,iexc
 write(unit) na(1:nc+nv),la(1:nc+nv),fa(1:nc+nv)
 write(unit) mmax,rr(1:mmax)
 end subroutine write_head

 subroutine read_head(fname,unit,irel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                     mmax,rr)
! opens fname, checks it against this run's input and reads rr; stops if
! the file was written for a different atom or relativistic treatment

 character(len=*) :: fname
 integer :: unit,irel
 character*2 :: atsym
 integer :: nc,nv,iexc,mmax
 integer :: na(30),la(30)
 real(dp) :: zz
 real(dp) :: fa(30),rr(mmax)

 character(len=8) :: magic_f
 character*2 :: atsym_f
 integer :: ios,version_f,irel_f,nc_f,nv_f,iexc_f,mmax_f
 integer :: na_f(30),la_f(30)
 real(dp) :: zz_f,fa_f(30)

 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='old',action='read',iostat=ios)
 if(ios/=0) then
  write(6,'(2a)') 'read_ae_file: ERROR cannot open ',trim(fname)
  stop 1
 end if

 read(unit,iostat=ios) magic_f,version_f,irel_f
 if(ios/=0 .or. magic_f/=magic) then
  write(6,'(2a)') 'read_ae_file: ERROR not an all-electron data file: ', &
&       trim(fname)
  stop 1
 end if
 if(version_f/=version) then
  write(6,'(a,i4,a,i4)') 'read_ae_file: ERROR file format version', &
&       version_f,' but this program reads version',version
  stop 1
 end if
 if(irel_f/=irel) then
  write(6,'(2a)') 'read_ae_file: ERROR file was written with a different', &
&       ' relativistic treatment'
  stop 1
 end if

 read(unit) atsym_f,zz_f,nc_f,nv_f,iexc_f
 if(atsym_f/=atsym .or. zz_f/=zz .or. nc_f/=nc .or. nv_f/=nv &
&   .or. iexc_f/=iexc) then
  write(6,'(2a)') 'read_ae_file: ERROR atsym, z, nc, nv or iexc differ', &
&       ' from the input'
  stop 1
 end if
 read(unit) na_f(1:nc+nv),la_f(1:nc+nv),fa_f(1:nc+nv)
 if(any(na_f(1:nc+nv)/=na(1:nc+nv)) .or. any(la_f(1:nc+nv)/=la(1:nc+nv)) &
&   .or. any(fa_f(1:nc+nv)/=fa(1:nc+nv))) then
  write(6,'(2a)') 'read_ae_file: ERROR reference configuration differs', &
&       ' from the input'
  stop 1
 end if

 read(unit) mmax_f
 if(mmax_f/=mmax) then
  write(6,'(a,2i8)') 'read_ae_file: ERROR radial mesh size differs', &
&       mmax_f,mmax
  stop 1
 end if
 read(unit) rr(1:mmax)
 end subroutine read_head

 subroutine write_config_key(unit,jj,nc,nvcnf,nacnf,lacnf,facnf)
 integer :: unit,jj,nc
 integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp) :: facnf(30,5)

 integer :: nt

 nt=nc+nvcnf(jj)
 write(unit) nvcnf(jj),nacnf(1:nt,jj),lacnf(1:nt,jj),facnf(1:nt,jj)
 end subroutine write_config_key

 subroutine read_config_key(unit,jj,nc)
 integer :: unit,jj,nc

 integer :: nt

 read(unit) nvcnf_file(jj)
 nt=nc+nvcnf_file(jj)
 read(unit) nacnf_file(1:nt,jj),lacnf_file(1:nt,jj),facnf_file(1:nt,jj)
 end subroutine read_config_key

 integer function find_config(jj,nc,nvcnf,nacnf,lacnf,facnf)
! index in the file of test configuration jj of this run; stops if it is
! not there
 integer :: jj,nc
 integer :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp) :: facnf(30,5)

 integer :: kk,nt

 nt=nc+nvcnf(jj)
 do kk=1,ncfg_file
  if(nvcnf_file(kk)/=nvcnf(jj)) cycle
  if(any(nacnf_file(1:nt,kk)/=nacnf(1:nt,jj))) cycle
  if(any(lacnf_file(1:nt,kk)/=lacnf(1:nt,jj))) cycle
  if(any(facnf_file(1:nt,kk)/=facnf(1:nt,jj))) cycle
  find_config=kk
  return
 end do

 write(6,'(a,i2,a)') 'ae_file_config: ERROR test configuration',jj-1, &
&     ' is not in the all-electron data file'
 stop 1
 end function find_config

 end module ae_file_m
