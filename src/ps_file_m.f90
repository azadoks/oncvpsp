!
! Copyright (c) 2026 by Austin Zadoks
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
! pseudopotential data file passed from oncvpsp.x or oncvpsp_ps.x (and the
! nr and fr ones, --ps-file) to oncvpsp_ld.x and oncvpsp_confine.x
!
! It holds what the log derivatives of run_phsft need of the
! pseudopotential: the screened local potential and its pieces, the
! projectors' pieces and coefficients, the reference energies at which the
! pseudopotential's phases are aligned with the all-electron ones, and the
! input's energies and radius of the log derivatives, oncvpsp_ld.x's
! defaults; and what the self-consistent pseudo-atom of psatom needs
! besides (oncvpsp_confine.x's atom in a wall): the unscreened local
! potential, the reference valence pseudo-density, the model core charge
! and the exchange-correlation functional.  The all-electron potential and
! the valence configuration are those of the all-electron data file
! (ae_file_m), which both programs read with it.
!
! The file is unformatted stream I/O in the native byte order:
!   magic, format version, irel                (0 non-, 1 scalar-,
!                                               2 fully-relativistic)
!   atsym, zz, nkap
!   mmax, rr                                   (log radial mesh)
!   c_light                                    (relativity_m)
!   lmax, lloc, nproj(1:4), rc(1:mxlin+1)
!   epsh1, epsh2, depsh, rxpsh
!   epa(1:4,1:nkap)                            (the first projector's
!                                               reference energies)
!   vloc(1:mmax,1:nkap), vlocpw                (vp(:,lloc+1,:), screened)
!   for l=0..3 and each kappa of l (one for l=0 or irel<2):
!     nvkb, evkb(1:nvkb), vkbpw(1:nvkb)        (nvkb 0 for l=lloc)
!   iexc, rcmax
!   vpuns(1:mmax), rhov(1:mmax), rhomod(1:mmax) (the unscreened local
!                                               potential, vpuns(:,lloc+1),
!                                               the valence pseudo-density
!                                               and the model core charge)
! a piecewise function (pieces_m) as nseg, bound(0:nseg), vals(:,1:nseg)

 module ps_file_m
 use constants_m, only: dp, mxprj, mxl, mxlin
 use pieces_m, only: piecewise_t, pw_maxseg, pw_alloc, pw_free
 use relativity_m, only: c_light
 use ae_file_m, only: byteswap
 implicit none
 private
 public :: ps_file_t, write_ps_file, read_ps_file, ps_file_free

 character(len=8), parameter :: magic='ONCVPS  ' !< first record of the file
 integer, parameter :: version=2 !< format version

!> the contents of the file
 type ps_file_t
   integer :: irel=0,nkap=0,mmax=0,lmax=0,lloc=0,iexc=0
   integer :: nproj(4)=0
!> nvkb(l+1,ikap): the projectors of l and kappa in the file
   integer :: nvkb(4,2)=0
   character(len=2) :: atsym='  '
   real(dp) :: zz=0.0d0,c_light=0.0d0
   real(dp) :: rc(mxlin+1)=0.0d0
   real(dp) :: epsh1=0.0d0,epsh2=0.0d0,depsh=0.0d0,rxpsh=0.0d0,rcmax=0.0d0
   real(dp) :: epa(4,2)=0.0d0
   real(dp) :: evkb(mxprj,4,2)=0.0d0
   real(dp), allocatable :: rr(:),vloc(:,:)
!> the unscreened local potential, the valence pseudo-density and the model
!> core charge
   real(dp), allocatable :: vpuns(:),rhov(:),rhomod(:)
   type(piecewise_t) :: vlocpw
   type(piecewise_t) :: vkbpw(mxprj,4,2)
 end type ps_file_t

 contains

 subroutine write_ps_file(fname,irel,nkap,atsym,zz,mmax,rr,lmax,lloc, &
&                         nproj,rc,epsh1,epsh2,depsh,rxpsh,epa,vp,vlocpw, &
&                         evkb,vkbpw,iexc,rcmax,vpuns,rhov,rhomod,ifail)

!! writes the pseudopotential data file
!!
!! @param fname  file to write
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh
!! @param lmax  maximum angular momentum of the projectors
!! @param lloc  l of the local potential (4 for the polynomial)
!! @param nproj  number of projectors of each l
!! @param rc  core radii, rc(lloc+1) that of the local potential
!! @param epsh1  low energy of the log derivatives (Ha)
!! @param epsh2  high energy of the log derivatives (Ha)
!! @param depsh  their energy step
!! @param rxpsh  their radius, negative for the default of run_phsft
!! @param epa  reference energies of the projectors
!! @param vp  semi-local pseudopotentials, vp(:,lloc+1,:) the local one
!! @param vlocpw  the local potential's pieces
!! @param evkb  coefficients of the projectors
!! @param vkbpw  the projectors' pieces
!! @param iexc  exchange-correlation function
!! @param rcmax  maximum core radius
!! @param vpuns  unscreened local potential
!! @param rhov  valence pseudo-density of the reference configuration
!! @param rhomod  model core charge
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 character(len=*), intent(in) :: fname
 integer, intent(in) :: irel,nkap,mmax,lmax,lloc,iexc
 character(len=2), intent(in) :: atsym
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: zz,rr(mmax),rc(mxlin+1)
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh,rcmax
 real(dp), intent(in) :: vpuns(mmax),rhov(mmax),rhomod(mmax)
 real(dp), intent(in) :: epa(mxprj,mxlin+1,nkap),vp(mmax,mxl+2,nkap)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: unit,ios,l1,ikap,jj,nvkb

 ifail=0
 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='replace',action='write',iostat=ios)
 if(ios/=0) then
   write(6,'(2a)') 'write_ps_file: ERROR cannot open ',trim(fname)
   ifail=1
   return
 end if

 write(unit) magic,version,irel
 write(unit) atsym,zz,nkap
 write(unit) mmax,rr(1:mmax)
 write(unit) c_light
 write(unit) lmax,lloc,nproj(1:4),rc(1:mxlin+1)
 write(unit) epsh1,epsh2,depsh,rxpsh
 write(unit) epa(1,1:4,1:nkap)
 write(unit) vp(1:mmax,lloc+1,1:nkap)
 call write_pw(unit,vlocpw)
 do l1=1,4
   do ikap=1,nkap_of(nkap,l1)
     nvkb=nvkb_of(l1,lloc,nproj)
     write(unit) nvkb
     if(nvkb==0) cycle
     write(unit) evkb(1:nvkb,l1,ikap)
     do jj=1,nvkb
       call write_pw(unit,vkbpw(jj,l1,ikap))
     end do
   end do
 end do
 write(unit) iexc,rcmax
 write(unit) vpuns(1:mmax),rhov(1:mmax),rhomod(1:mmax)
 close(unit)
 return
 end subroutine write_ps_file

 subroutine read_ps_file(fname,pf,ifail)

!! reads the pseudopotential data file
!!
!! @param fname  file to read
!! @param pf  its contents
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 character(len=*), intent(in) :: fname

!Output variables
 type(ps_file_t), intent(out) :: pf
 integer, intent(out) :: ifail

!Local variables
 character(len=8) :: magic_f
 integer :: unit,ios,version_f,l1,ikap,jj,nvkb

 ifail=0
 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='old',action='read',iostat=ios)
 if(ios/=0) then
   write(6,'(2a)') 'read_ps_file: ERROR cannot open ',trim(fname)
   ifail=1
   return
 end if
 read(unit,iostat=ios) magic_f,version_f,pf%irel
 if(ios/=0 .or. magic_f/=magic) then
   write(6,'(2a)') 'read_ps_file: ERROR not a pseudopotential data file: ', &
&        trim(fname)
   ifail=1
   close(unit)
   return
 end if
 if(version_f/=version .and. byteswap(version_f)==version) then
   write(6,'(2a)') 'read_ps_file: ERROR file was written on a machine', &
&        ' with the other byte order'
   ifail=1
   close(unit)
   return
 end if
 if(version_f/=version) then
   write(6,'(a,i4,a,i4)') 'read_ps_file: ERROR file format version', &
&        version_f,' but this program reads version',version
   ifail=1
   close(unit)
   return
 end if

 read(unit) pf%atsym,pf%zz,pf%nkap
 read(unit) pf%mmax
 allocate(pf%rr(pf%mmax),pf%vloc(pf%mmax,pf%nkap))
 read(unit) pf%rr(:)
 read(unit) pf%c_light
 read(unit) pf%lmax,pf%lloc,pf%nproj(1:4),pf%rc(1:mxlin+1)
 read(unit) pf%epsh1,pf%epsh2,pf%depsh,pf%rxpsh
 read(unit) pf%epa(1:4,1:pf%nkap)
 read(unit) pf%vloc(:,:)
 call read_pw(unit,pf%mmax,pf%vlocpw,ifail)
 if(ifail/=0) then
   close(unit)
   return
 end if
 do l1=1,4
   do ikap=1,nkap_of(pf%nkap,l1)
     read(unit) nvkb
     pf%nvkb(l1,ikap)=nvkb
     if(nvkb==0) cycle
     read(unit) pf%evkb(1:nvkb,l1,ikap)
     do jj=1,nvkb
       call read_pw(unit,pf%mmax,pf%vkbpw(jj,l1,ikap),ifail)
       if(ifail/=0) then
         close(unit)
         return
       end if
     end do
   end do
 end do
 read(unit) pf%iexc,pf%rcmax
 allocate(pf%vpuns(pf%mmax),pf%rhov(pf%mmax),pf%rhomod(pf%mmax))
 read(unit) pf%vpuns(:),pf%rhov(:),pf%rhomod(:)
 close(unit)
 return
 end subroutine read_ps_file

 subroutine ps_file_free(pf)

!! releases the contents
!!
!! @param pf  the contents of a file

!Input/output variables
 type(ps_file_t), intent(inout) :: pf

!Local variables
 integer :: l1,ikap,jj

 if(allocated(pf%rr)) deallocate(pf%rr)
 if(allocated(pf%vloc)) deallocate(pf%vloc)
 if(allocated(pf%vpuns)) deallocate(pf%vpuns,pf%rhov,pf%rhomod)
 call pw_free(pf%vlocpw)
 do ikap=1,2
   do l1=1,4
     do jj=1,mxprj
       call pw_free(pf%vkbpw(jj,l1,ikap))
     end do
   end do
 end do
 return
 end subroutine ps_file_free

 integer function nkap_of(nkap,l1)
!! the kappas of l=l1-1: one for l=0
 integer, intent(in) :: nkap,l1
 nkap_of=nkap
 if(l1==1) nkap_of=1
 end function nkap_of

 integer function nvkb_of(l1,lloc,nproj)
!! the projectors of l=l1-1: none for the local potential's l
 integer, intent(in) :: l1,lloc,nproj(mxlin+1)
 nvkb_of=nproj(l1)
 if(l1-1==lloc) nvkb_of=0
 end function nvkb_of

 subroutine write_pw(unit,pw)
!! writes a piecewise function
 integer, intent(in) :: unit
 type(piecewise_t), intent(in) :: pw
 write(unit) pw%nseg,pw%bound(0:pw%nseg)
 write(unit) pw%vals(:,1:pw%nseg)
 end subroutine write_pw

 subroutine read_pw(unit,mmax,pw,ifail)
!! reads a piecewise function written by write_pw
 integer, intent(in) :: unit,mmax
 type(piecewise_t), intent(out) :: pw
 integer :: nseg
 real(dp) :: bound(0:pw_maxseg)
 integer, intent(out) :: ifail
 ifail=0
 read(unit) nseg
 if(nseg<1 .or. nseg>pw_maxseg) then
   write(6,'(a,i0)') 'read_ps_file: ERROR number of segments ',nseg
   ifail=1
   return
 end if
 read(unit) bound(0:nseg)
 call pw_alloc(pw,mmax,nseg,bound(0:nseg))
 read(unit) pw%vals(:,:)
 end subroutine read_pw

 end module ps_file_m
