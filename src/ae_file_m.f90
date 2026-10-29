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
!   c_light                                    (the run's, checked)
!   it, etot, eterm, res, ea, rhoc, rho, vfull (reference atom, aeatom;
!                                               ea(:,1:2) if irel=2)
!   ncfg, then for each test configuration:
!     nvcnf, nacnf, lacnf, facnf               (key, matched to the input)
!     istcnf through etstcnf for it            (run_config_ae results)

 module ae_file_m
 use, intrinsic :: iso_fortran_env, only: int8
 use constants_m, only: dp, mxstate, mxcnf, nlev
 implicit none
 private
 public :: write_ae_file, read_ae_file, read_ae_potential, byteswap, differ

 character(len=8), parameter :: magic='ONCVAE  ' !< first record of the file
 integer, parameter :: version=3 !< format version
!> relative tolerance for comparing real input (zz, occupations) and the
!> radial mesh with the file: both programs parse the same input text, but
!> may be built by different compilers or with different flags
 real(dp), parameter :: rtol=1.0d-12

 contains

 subroutine write_ae_file(fname,irel,c_light,nkap,atsym,zz,nc,nv,iexc,na,la, &
&                         fa,mmax,rr,it,etot,eterm,res,ea,rhoc,rho,vfull,ncnf, &
&                         cnf,aec,iout,ifail)

!! writes the all-electron data file: the aeatom results for the reference
!! atom and the run_config_ae results for all configurations
!!
!! @param fname  file to write
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional
!! @param na  principal quantum numbers, reference configuration
!! @param la  angular momenta, reference configuration
!! @param fa  occupancies, reference configuration
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh (bohr)
!! @param it  aeatom iteration count
!! @param etot  all-electron total energy (Ha)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (Ha)
!! @param res  aeatom residual of the screening potential (Ha)
!! @param ea  all-electron eigenvalues (Ha), kap=-(l+1) and kap=l if irel=2
!! @param rhoc  core charge density
!! @param rho  all-electron valence charge density
!! @param vfull  all-electron potential (Ha)
!! @param ncnf  number of test configurations (not counting the reference)
!! @param cnf  the reference and test configurations (input_m)
!! @param aec  the run_config_ae results, all configurations
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use config_ae_m, only: ae_config_t
 use input_m, only: config_t

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: irel,nkap
 real(dp), intent(in) :: c_light
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax,it,ncnf
 type(config_t), intent(in) :: cnf(mxcnf)
 integer, intent(in) :: na(mxstate),la(mxstate)
 real(dp), intent(in) :: zz,etot,eterm(4),res
 real(dp), intent(in) :: fa(mxstate),ea(mxstate,nkap)
 real(dp), intent(in) :: rr(mmax),rhoc(mmax),rho(mmax),vfull(mmax)
 type(ae_config_t), intent(in) :: aec(mxcnf)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: unit,jj

 ifail=0
 call write_head(fname,unit,irel,c_light,atsym,zz,nc,nv,iexc,na,la,fa, &
&                mmax,rr,iout,ifail)
 if(ifail/=0) return
 write(unit) it,etot,eterm,res,ea(1:nc+nv,1:nkap)
 write(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)
 write(unit) ncnf+1
 do jj=1,ncnf+1
   call write_config_key(unit,jj,nc,cnf)
   write(unit) aec(jj)%ist,aec(jj)%ierr,aec(jj)%it,aec(jj)%kerr, &
&              aec(jj)%nvt,aec(jj)%nat(:),aec(jj)%lat(:),aec(jj)%fat(:,:), &
&              aec(jj)%eat(:,:,:),aec(jj)%etst
 end do

 close(unit)
 return
 end subroutine write_ae_file

 subroutine read_ae_file(fname,irel,c_light,nkap,atsym,zz,nc,nv,iexc,na,la,fa, &
&                        mmax,rr,ncnf,cnf, &
&                        it,etot,eterm,res,ea,rhoc,rho,vfull,aec,iout,ifail)

!! reads the all-electron data file: the reference atom and the
!! run_config_ae results for each configuration of this run, found in the
!! file by matching the configuration itself, so the file may hold more
!! test configurations than the input, in any order.  Stops if the file
!! was written for a different atom or relativistic treatment, or lacks a
!! configuration.  irel through cnf are this run's input, checked
!! against the file; it through aec are read from it.
!!
!! @param fname  file to read
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional
!! @param na  principal quantum numbers, reference configuration
!! @param la  angular momenta, reference configuration
!! @param fa  occupancies, reference configuration
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh (bohr)
!! @param ncnf  number of test configurations (not counting the reference)
!! @param cnf  the reference and test configurations (input_m)
!! @param it  aeatom iteration count
!! @param etot  all-electron total energy (Ha)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (Ha)
!! @param res  aeatom residual of the screening potential (Ha)
!! @param ea  all-electron eigenvalues (Ha), kap=-(l+1) and kap=l if irel=2
!! @param rhoc  core charge density
!! @param rho  all-electron valence charge density
!! @param vfull  all-electron potential (Ha)
!! @param aec  the run_config_ae results, all configurations
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use config_ae_m, only: ae_config_t
 use input_m, only: config_t

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: irel,nkap
 real(dp), intent(in) :: c_light
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax,ncnf
 type(config_t), intent(in) :: cnf(mxcnf)
 integer, intent(in) :: na(mxstate),la(mxstate)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(mxstate),rr(mmax)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot,eterm(4),res
 real(dp), intent(out) :: ea(mxstate,nkap),rhoc(mmax),rho(mmax),vfull(mmax)
 type(ae_config_t), intent(out) :: aec(mxcnf)
 integer, intent(out) :: ifail

!Local variables
 integer :: unit,jj,kk,ncfg
 integer :: ist,ierr,itc,kerr,nvt,nat(mxstate),lat(mxstate)
 real(dp) :: eaetst,fat(mxstate,nlev),eat(mxstate,nkap,nlev)
 logical :: found(5),match(5)

 ifail=0
 call read_head(fname,unit,irel,c_light,atsym,zz,nc,nv,iexc,na,la,fa, &
&               mmax,rr,iout,ifail)
 if(ifail/=0) return
 read(unit) it,etot,eterm,res,ea(1:nc+nv,1:nkap)
 read(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)

! copy the results of each configuration in the file to the configurations
! of this run that it matches
 found(:)=.false.
 read(unit) ncfg
 do kk=1,ncfg
   call match_config(unit,nc,ncnf,cnf,found,match)
   read(unit) ist,ierr,itc,kerr,nvt,nat,lat,fat,eat,eaetst
   do jj=1,ncnf+1
     if(.not. match(jj)) cycle
     aec(jj)%ist=ist
     aec(jj)%ierr=ierr
     aec(jj)%it=itc
     aec(jj)%kerr=kerr
     aec(jj)%nvt=nvt
     aec(jj)%nat(:)=nat(:)
     aec(jj)%lat(:)=lat(:)
     aec(jj)%fat(:,:)=fat(:,:)
     aec(jj)%eat=eat
     aec(jj)%etst=eaetst
     found(jj)=.true.
   end do
 end do
 call check_found(ncnf,found,iout,ifail)
 if(ifail/=0) then
   close(unit)
   return
 end if

 close(unit)
 return
 end subroutine read_ae_file

 subroutine write_head(fname,unit,irel,c_light,atsym,zz,nc,nv,iexc,na,la,fa, &
&                      mmax,rr,iout,ifail)
!! opens fname and writes everything before the reference atom results
!!
!! @param fname  file name
!! @param unit  unit number of the opened file
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional
!! @param na  principal quantum numbers, reference configuration
!! @param la  angular momenta, reference configuration
!! @param fa  occupancies, reference configuration
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh (bohr)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_io, failed, msg_ae_file_open_write

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: irel
 real(dp), intent(in) :: c_light
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax
 integer, intent(in) :: na(mxstate),la(mxstate)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(mxstate),rr(mmax)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: unit
 integer, intent(out) :: ifail

!Local variables
 integer :: ios

 ifail=0
 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='replace',action='write',iostat=ios)
 if(ios/=0) then
   write(iout,'(/2a)') 'write_ae_file: ERROR cannot open ',trim(fname)
   ifail=failed(fail_io,msg_ae_file_open_write)
   return
 end if

 write(unit) magic,version,irel
 write(unit) atsym,zz,nc,nv,iexc
 write(unit) na(1:nc+nv),la(1:nc+nv),fa(1:nc+nv)
 write(unit) mmax,rr(1:mmax)
 write(unit) c_light
 return
 end subroutine write_head

 subroutine read_head(fname,unit,irel,c_light,atsym,zz,nc,nv,iexc,na,la,fa, &
&                     mmax,rr,iout,ifail)
!! opens fname, checks it against this run's input and reads rr; fails if
!! the file was written for a different atom, relativistic treatment or
!! speed of light
!!
!! @param fname  file name
!! @param unit  unit number of the opened file
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic:
!!              the treatment expected
!! @param c_light  speed of light (Hartree atomic units)
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional
!! @param na  principal quantum numbers, reference configuration
!! @param la  angular momenta, reference configuration
!! @param fa  occupancies, reference configuration
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh (bohr)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_io, failed, msg_ae_file_atom, msg_ae_file_c_light, &
&                    msg_ae_file_irel, msg_ae_file_mesh, msg_ae_file_mmax, &
&                    msg_ae_file_reference

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: irel
 real(dp), intent(in) :: c_light
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nc,nv,iexc,mmax
 integer, intent(in) :: na(mxstate),la(mxstate)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(mxstate),rr(mmax)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: unit
 integer, intent(out) :: ifail

!Local variables
 character(len=2) :: atsym_f
 integer :: irel_f,nc_f,nv_f,iexc_f,mmax_f
 integer :: na_f(mxstate),la_f(mxstate)
 real(dp) :: zz_f,fa_f(mxstate),c_light_f
 real(dp), allocatable :: rr_f(:)

 ifail=0
 call open_ae(fname,unit,irel_f,iout,ifail)
 if(ifail/=0) return
 if(irel_f/=irel) then
   write(iout,'(/2a)') &
&        'read_ae_file: ERROR file was written with a different', &
&        ' relativistic treatment'
   ifail=failed(fail_io,msg_ae_file_irel)
   close(unit)
   return
 end if

 read(unit) atsym_f,zz_f,nc_f,nv_f,iexc_f
 if(atsym_f/=atsym .or. differ(zz_f,zz) .or. nc_f/=nc .or. nv_f/=nv &
&   .or. iexc_f/=iexc) then
   write(iout,'(/2a)') 'read_ae_file: ERROR atsym, z, nc, nv or iexc differ', &
&        ' from the input'
   ifail=failed(fail_io,msg_ae_file_atom)
   close(unit)
   return
 end if
 read(unit) na_f(1:nc+nv),la_f(1:nc+nv),fa_f(1:nc+nv)
 if(any(na_f(1:nc+nv)/=na(1:nc+nv)) .or. any(la_f(1:nc+nv)/=la(1:nc+nv)) &
&   .or. any(differ(fa_f(1:nc+nv),fa(1:nc+nv)))) then
   write(iout,'(/2a)') 'read_ae_file: ERROR reference configuration differs', &
&        ' from the input'
   ifail=failed(fail_io,msg_ae_file_reference)
   close(unit)
   return
 end if

 read(unit) mmax_f
 if(mmax_f/=mmax) then
   write(iout,'(/a,2i0)') 'read_ae_file: ERROR radial mesh size differs', &
&        mmax_f,mmax
   ifail=failed(fail_io,msg_ae_file_mmax)
   close(unit)
   return
 end if
 allocate(rr_f(mmax))
 read(unit) rr_f(1:mmax)
 if(any(abs(rr_f(1:mmax)-rr(1:mmax))>rtol*rr(1:mmax))) then
   write(iout,'(/2a)') 'read_ae_file: ERROR radial mesh differs from the', &
&        ' one this program builds'
   ifail=failed(fail_io,msg_ae_file_mesh)
   close(unit)
   return
 end if
 read(unit) c_light_f
 if(differ(c_light_f,c_light)) then
   write(iout,'(/2a)') 'read_ae_file: ERROR c_light differs from the input'
   ifail=failed(fail_io,msg_ae_file_c_light)
   close(unit)
   return
 end if
 return
 end subroutine read_head

 subroutine read_ae_potential(fname,irel,atsym,zz,mmax,rr,c_light_f,vfull, &
&                             nc,nv,iexc,na,la,fa,ea,etot,iout,ifail)

!! reads the atom, the mesh and the all-electron potential of the
!! reference atom from the all-electron data file, without an input to
!! check them against (oncvpsp_ld.x), and optionally its configuration,
!! eigenvalues and total energy (oncvpsp_confine.x)
!!
!! @param fname  file to read
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh (bohr) (allocated here)
!! @param c_light_f  the speed of light of the file
!! @param vfull  all-electron potential (Ha) (allocated here)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional
!! @param na  principal quantum numbers, reference configuration
!! @param la  angular momenta
!! @param fa  occupancies
!! @param ea  eigenvalues, ea(:,2) those of j=l-1/2 if irel=2
!! @param etot  total energy (Ha)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: mxkap

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: irel,mmax
 character(len=2), intent(out) :: atsym
 real(dp), intent(out) :: zz,c_light_f
 real(dp), allocatable, intent(out) :: rr(:),vfull(:)
 integer, intent(out), optional :: nc,nv,iexc,na(mxstate),la(mxstate)
 real(dp), intent(out), optional :: fa(mxstate),ea(mxstate,mxkap),etot
 integer, intent(out) :: ifail

!Local variables
 integer :: unit,nc_f,nv_f,iexc_f,nt,nkap,it
 integer :: na_f(mxstate),la_f(mxstate)
 real(dp) :: fa_f(mxstate),etot_f,eterm(4),res,ea_f(mxstate,mxkap)
 real(dp), allocatable :: rhoc(:),rho(:)

 ifail=0
 call open_ae(fname,unit,irel,iout,ifail)
 if(ifail/=0) return
 read(unit) atsym,zz,nc_f,nv_f,iexc_f
 nt=nc_f+nv_f
 na_f(:)=0
 la_f(:)=0
 fa_f(:)=0.0d0
 read(unit) na_f(1:nt),la_f(1:nt),fa_f(1:nt)
 read(unit) mmax
 allocate(rr(mmax),vfull(mmax),rhoc(mmax),rho(mmax))
 read(unit) rr(1:mmax)
 read(unit) c_light_f
 nkap=1
 if(irel==2) nkap=2
 ea_f(:,:)=0.0d0
 read(unit) it,etot_f,eterm,res,ea_f(1:nt,1:nkap)
 read(unit) rhoc(1:mmax),rho(1:mmax),vfull(1:mmax)
 close(unit)
 deallocate(rhoc,rho)
 if(present(nc)) nc=nc_f
 if(present(nv)) nv=nv_f
 if(present(iexc)) iexc=iexc_f
 if(present(na)) na(:)=na_f(:)
 if(present(la)) la(:)=la_f(:)
 if(present(fa)) fa(:)=fa_f(:)
 if(present(ea)) ea(:,:)=ea_f(:,:)
 if(present(etot)) etot=etot_f
 return
 end subroutine read_ae_potential

 subroutine open_ae(fname,unit,irel_f,iout,ifail)
!! opens fname and reads its first record, failing if it is not an
!! all-electron data file of this format version and byte order
!!
!! @param fname  file name
!! @param unit  unit number of the opened file
!! @param irel_f  the relativistic treatment of the file
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_io, failed, msg_ae_file_byte_order, &
&                    msg_ae_file_not_ae, msg_ae_file_open_read, &
&                    msg_ae_file_version

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: unit,irel_f
 integer, intent(out) :: ifail

!Local variables
 character(len=8) :: magic_f
 integer :: ios,version_f

 ifail=0
 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='old',action='read',iostat=ios)
 if(ios/=0) then
   write(iout,'(/2a)') 'read_ae_file: ERROR cannot open ',trim(fname)
   ifail=failed(fail_io,msg_ae_file_open_read)
   return
 end if

 read(unit,iostat=ios) magic_f,version_f,irel_f
 if(ios/=0 .or. magic_f/=magic) then
   write(iout,'(/2a)') 'read_ae_file: ERROR not an all-electron data file: ', &
&        trim(fname)
   ifail=failed(fail_io,msg_ae_file_not_ae)
   close(unit)
   return
 end if
!the magic string reads the same in either byte order, so a file from a
!machine of the other byte order shows up as a byte-swapped version
 if(version_f/=version .and. byteswap(version_f)==version) then
   write(iout,'(/2a)') 'read_ae_file: ERROR file was written on a machine', &
&        ' with the other byte order'
   ifail=failed(fail_io,msg_ae_file_byte_order)
   close(unit)
   return
 end if
 if(version_f/=version) then
   write(iout,'(/a,i0,a,i0)') 'read_ae_file: ERROR file format version ', &
&        version_f,' but this program reads version ',version
   ifail=failed(fail_io,msg_ae_file_version)
   close(unit)
   return
 end if
 return
 end subroutine open_ae

 subroutine write_config_key(unit,jj,nc,cnf)
!! writes configuration jj, the key by which match_config matches it
!!
!! @param unit  unit number of the open file
!! @param jj  index of the configuration in this run (1 for the reference)
!! @param nc  number of core states
!! @param cnf  the reference and test configurations (input_m)

 use input_m, only: config_t

!Input variables
 integer, intent(in) :: unit,jj,nc
 type(config_t), intent(in) :: cnf(mxcnf)

!Local variables
 integer :: nt

 nt=nc+cnf(jj)%nv
 write(unit) cnf(jj)%nv,cnf(jj)%na(1:nt),cnf(jj)%la(1:nt),cnf(jj)%fa(1:nt)
 return
 end subroutine write_config_key

 subroutine match_config(unit,nc,ncnf,cnf,found,match)
!! reads the key of the next configuration in the file and flags the
!! configurations of this run, not found before, that it matches
!!
!! @param unit  unit number of the open file
!! @param nc  number of core states
!! @param ncnf  number of test configurations (not counting the reference)
!! @param cnf  the reference and test configurations (input_m)
!! @param found  .true. for the configurations already found in the file
!! @param match  .true. for the configurations this key matches

 use input_m, only: config_t

!Input variables
 integer, intent(in) :: unit,nc,ncnf
 type(config_t), intent(in) :: cnf(mxcnf)
 logical, intent(in) :: found(5)

!Output variables
 logical, intent(out) :: match(5)

!Local variables
 integer :: jj,nt,nv_f
 integer :: na_f(mxstate),la_f(mxstate)
 real(dp) :: fa_f(mxstate)

 read(unit) nv_f
 nt=nc+nv_f
 read(unit) na_f(1:nt),la_f(1:nt),fa_f(1:nt)

 match(:)=.false.
 do jj=1,ncnf+1
   if(found(jj)) cycle
   if(nv_f/=cnf(jj)%nv) cycle
   if(any(na_f(1:nt)/=cnf(jj)%na(1:nt))) cycle
   if(any(la_f(1:nt)/=cnf(jj)%la(1:nt))) cycle
   if(any(differ(fa_f(1:nt),cnf(jj)%fa(1:nt)))) cycle
   match(jj)=.true.
 end do
 return
 end subroutine match_config

 subroutine check_found(ncnf,found,iout,ifail)
!! fails if a configuration of this run is not in the file
!!
!! @param ncnf  number of test configurations (not counting the reference)
!! @param found  .true. for the configurations found in the file
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_io, failed, msg_ae_file_test_config

!Input variables
 integer, intent(in) :: ncnf
 logical, intent(in) :: found(5)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: jj

 ifail=0
 do jj=1,ncnf+1
   if(.not. found(jj)) then
     write(iout,'(/a,i0,a)') 'read_ae_file: ERROR test configuration ',jj-1, &
&         ' is not in the all-electron data file'
     ifail=failed(fail_io,msg_ae_file_test_config)
     return
   end if
 end do
 return
 end subroutine check_found

 elemental logical function differ(aa,bb)
!! .true. if aa and bb differ by more than rtol relative to the larger,
!! or absolutely for values below 1 (occupations of 0); not for the mesh,
!! whose points near the origin are far below 1
!!
!! @param aa  first value
!! @param bb  second value

!Input variables
 real(dp), intent(in) :: aa,bb
 differ=abs(aa-bb)>rtol*max(1.0d0,abs(aa),abs(bb))
 return
 end function differ

 integer function byteswap(ii)
!! ii with its bytes in reverse order
!!
!! @param ii  integer to byte-swap

!Input variables
 integer, intent(in) :: ii

!Local variables
 integer(int8) :: bytes(storage_size(ii)/8)
 bytes=transfer(ii,bytes)
 byteswap=transfer(bytes(size(bytes):1:-1),ii)
 return
 end function byteswap

 end module ae_file_m
