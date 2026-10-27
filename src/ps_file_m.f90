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
!   c_light                                    (the run's)
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
 use constants_m, only: dp, mxprj, mxl, mxlin, nrhomod
 use pieces_m, only: piecewise_t, pw_maxseg, pw_alloc
 use pseudo_m, only: pseudo_t
 use ae_file_m, only: byteswap
 implicit none
 private
 public :: ps_file_t, write_ps_file, read_ps_file, ps_file_free

 character(len=8), parameter :: magic='ONCVPS  ' !< first record of the file
 integer, parameter :: version=2 !< format version

!> the contents of the file
 type ps_file_t
!> the pseudopotential: of its potentials only the local one,
!> vpuns(:,lloc+1), and of the model core only the charge, rhomod(:,1);
!> neither vkb nor the fully-relativistic split nor corepw
   type(pseudo_t) :: pot
!> nvkb(l+1,ikap): the projectors of l and kappa in the file
   integer :: nvkb(4,2)=0
!> the core radii, and the energies (Ha) and radius of the log
!> derivatives of the input
   real(dp) :: rc(mxlin+1)=0.0d0
   real(dp) :: epsh1=0.0d0,epsh2=0.0d0,depsh=0.0d0,rxpsh=0.0d0
!> the first projector's reference energy of each l and kappa
   real(dp) :: epa(4,2)=0.0d0
!> the screened local potential of each kappa (mmax,nkap)
   real(dp), allocatable :: vloc(:,:)
 end type ps_file_t

 contains

 subroutine write_ps_file(fname,pot,rc,epsh1,epsh2,depsh,rxpsh,epa,vp, &
&                         iout,ifail)

!! writes the pseudopotential data file
!!
!! @param fname  file to write
!! @param pot  the pseudopotential (pseudo_m)
!! @param rc  core radii, rc(lloc+1) that of the local potential
!! @param epsh1  low energy of the log derivatives (Ha)
!! @param epsh2  high energy of the log derivatives (Ha)
!! @param depsh  their energy step
!! @param rxpsh  their radius, negative for the default of run_phsft
!! @param epa  reference energies of the projectors
!! @param vp  semi-local pseudopotentials, vp(:,lloc+1,:) the local one
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 character(len=*), intent(in) :: fname
 type(pseudo_t), intent(in) :: pot
 real(dp), intent(in) :: rc(mxlin+1)
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh
 real(dp), intent(in) :: epa(mxprj,mxlin+1,pot%nkap)
 real(dp), intent(in) :: vp(pot%mmax,mxl+2,pot%nkap)

!Output variables
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 integer :: unit,ios,l1,ikap,jj,nvkb

 ifail=0
 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='replace',action='write',iostat=ios)
 if(ios/=0) then
   write(iout,'(2a)') 'write_ps_file: ERROR cannot open ',trim(fname)
   ifail=1
   return
 end if

 write(unit) magic,version,pot%irel
 write(unit) pot%atsym,pot%zz,pot%nkap
 write(unit) pot%mmax,pot%rr(1:pot%mmax)
 write(unit) pot%c_light
 write(unit) pot%lmax,pot%lloc,pot%nproj(1:4),rc(1:mxlin+1)
 write(unit) epsh1,epsh2,depsh,rxpsh
 write(unit) epa(1,1:4,1:pot%nkap)
 write(unit) vp(1:pot%mmax,pot%lloc+1,1:pot%nkap)
 call write_pw(unit,pot%vlocpw)
 do l1=1,4
   do ikap=1,nkap_of(pot%nkap,l1)
     nvkb=nvkb_of(l1,pot%lloc,pot%nproj)
     write(unit) nvkb
     if(nvkb==0) cycle
     write(unit) pot%evkb(1:nvkb,l1,ikap)
     do jj=1,nvkb
       call write_pw(unit,pot%vkbpw(jj,l1,ikap))
     end do
   end do
 end do
 write(unit) pot%iexc,pot%rcmax
 write(unit) pot%vpuns(1:pot%mmax,pot%lloc+1),pot%rho(1:pot%mmax), &
&            pot%rhomod(1:pot%mmax,1)
 close(unit)
 return
 end subroutine write_ps_file

 subroutine read_ps_file(fname,pf,iout,ifail)

!! reads the pseudopotential data file
!!
!! @param fname  file to read
!! @param pf  its contents
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 character(len=*), intent(in) :: fname

!Output variables
 type(ps_file_t), intent(out) :: pf
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 character(len=8) :: magic_f
 integer :: unit,ios,version_f,irel,l1,ikap,jj,nvkb

 ifail=0
 open(newunit=unit,file=fname,access='stream',form='unformatted', &
&     status='old',action='read',iostat=ios)
 if(ios/=0) then
   write(iout,'(2a)') 'read_ps_file: ERROR cannot open ',trim(fname)
   ifail=1
   return
 end if
 read(unit,iostat=ios) magic_f,version_f,irel
 if(ios/=0 .or. magic_f/=magic) then
   write(iout,'(2a)') 'read_ps_file: ERROR not a pseudopotential data file: ', &
&        trim(fname)
   ifail=1
   close(unit)
   return
 end if
 if(version_f/=version .and. byteswap(version_f)==version) then
   write(iout,'(2a)') 'read_ps_file: ERROR file was written on a machine', &
&        ' with the other byte order'
   ifail=1
   close(unit)
   return
 end if
 if(version_f/=version) then
   write(iout,'(a,i4,a,i4)') 'read_ps_file: ERROR file format version', &
&        version_f,' but this program reads version',version
   ifail=1
   close(unit)
   return
 end if

 pf%pot%irel=irel
 read(unit) pf%pot%atsym,pf%pot%zz,pf%pot%nkap
 read(unit) pf%pot%mmax
 associate(pot=>pf%pot,mmax=>pf%pot%mmax,nkap=>pf%pot%nkap)
   allocate(pot%rr(mmax),pf%vloc(mmax,nkap))
   read(unit) pot%rr(:)
   read(unit) pot%c_light
   read(unit) pot%lmax,pot%lloc,pot%nproj(1:4),pf%rc(1:mxlin+1)
   read(unit) pf%epsh1,pf%epsh2,pf%depsh,pf%rxpsh
   read(unit) pf%epa(1:4,1:nkap)
   read(unit) pf%vloc(:,:)
   call read_pw(unit,mmax,pot%vlocpw,iout,ifail)
   if(ifail/=0) then
     close(unit)
     return
   end if
   allocate(pot%evkb(mxprj,mxl+1,nkap),pot%vkbpw(mxprj,mxl+1,nkap))
   pot%evkb(:,:,:)=0.0d0
   do l1=1,4
     do ikap=1,nkap_of(nkap,l1)
       read(unit) nvkb
       pf%nvkb(l1,ikap)=nvkb
       if(nvkb==0) cycle
       read(unit) pot%evkb(1:nvkb,l1,ikap)
       do jj=1,nvkb
         call read_pw(unit,mmax,pot%vkbpw(jj,l1,ikap),iout,ifail)
         if(ifail/=0) then
           close(unit)
           return
         end if
       end do
     end do
   end do
   read(unit) pot%iexc,pot%rcmax
   allocate(pot%vpuns(mmax,mxl+2),pot%rho(mmax),pot%rhomod(mmax,nrhomod))
   pot%vpuns(:,:)=0.0d0
   pot%rhomod(:,:)=0.0d0
   read(unit) pot%vpuns(:,pot%lloc+1),pot%rho(:),pot%rhomod(:,1)
 end associate
 close(unit)
 return
 end subroutine read_ps_file

 subroutine ps_file_free(pf)

!! releases the contents
!!
!! @param pf  the contents of a file, deallocated on entry

!Output variables
 type(ps_file_t), intent(out) :: pf

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

 subroutine read_pw(unit,mmax,pw,iout,ifail)
!! reads a piecewise function written by write_pw
 integer, intent(in) :: unit,mmax
 type(piecewise_t), intent(out) :: pw
 integer :: nseg
 real(dp) :: bound(0:pw_maxseg)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail
 ifail=0
 read(unit) nseg
 if(nseg<1 .or. nseg>pw_maxseg) then
   write(iout,'(a,i0)') 'read_ps_file: ERROR number of segments ',nseg
   ifail=1
   return
 end if
 read(unit) bound(0:nseg)
 call pw_alloc(pw,mmax,nseg,bound(0:nseg))
 read(unit) pw%vals(:,:)
 end subroutine read_pw

 end module ps_file_m
