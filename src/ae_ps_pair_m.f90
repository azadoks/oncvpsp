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
!> The all-electron data file of an _ae.x (ae_file_m) and the
!> pseudopotential data file of a run with --ps-file (ps_file_m) of
!> oncvpsp_ld.x and oncvpsp_confine.x: both read, checked against each
!> other, and the build written for their relativistic treatment.
 module ae_ps_pair_m
 use constants_m, only: dp, mxstate
 implicit none
 private
 public :: load_pair

!> relative tolerance of the comparison of the two files' meshes, as in
!> ae_file_m
 real(dp), parameter :: rtol=1.0d-12

 contains

 subroutine load_pair(prog,ae_filename,ps_filename,irel,c_light,zz,mmax,rr, &
&                     vfull,pf,nc,nv,iexc,na,la,fa,ea,iout,ifail)
!! the two files of program prog, which must be of the same relativistic
!! treatment, atom, mesh and c_light (and functional, with iexc), and the
!! all-electron atom's configuration if asked; writes the build line
!!
!! @param prog  the program name, for the messages
!! @param ae_filename  the all-electron data file (ae_file_m)
!! @param ps_filename  the pseudopotential data file (ps_file_m)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units), that of the files
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh (bohr); allocated here
!! @param vfull  all-electron potential (Ha); allocated here
!! @param pf  the contents of the pseudopotential data file
!! @param nc  number of core states (optional)
!! @param nv  number of valence states (optional)
!! @param iexc  exchange-correlation functional (optional); when present,
!!              the two files must have the same
!! @param na  principal quantum numbers (optional)
!! @param la  angular momenta (optional)
!! @param fa  occupancies (optional)
!! @param ea  the all-electron eigenvalues (Ha) of each state and kappa
!!            (optional)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written
 use ae_file_m, only: read_ae_potential, differ
 use build_info_m, only: write_build_info
 use ps_file_m, only: ps_file_t, read_ps_file
 use constants_m, only: mxkap

!Input variables
 character(len=*), intent(in) :: prog,ae_filename,ps_filename
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: irel,mmax
 real(dp), intent(out) :: zz
!> the speed of light of the files
 real(dp), intent(out) :: c_light
 real(dp), allocatable, intent(out) :: rr(:),vfull(:)
 type(ps_file_t), intent(out) :: pf
 integer, intent(out), optional :: nc,nv,iexc,na(mxstate),la(mxstate)
 real(dp), intent(out), optional :: fa(mxstate),ea(mxstate,mxkap)
 integer, intent(out) :: ifail

!Local variables
 character(len=2) :: atsym
 character(len=19) :: treatment
 real(dp) :: c_light_f
 logical :: bad

 ifail=0
 call read_ae_potential(ae_filename,irel,atsym,zz,mmax,rr,c_light_f, &
&                       vfull,nc=nc,nv=nv,iexc=iexc,na=na,la=la,fa=fa,ea=ea, &
&                       iout=iout,ifail=ifail)
 if(ifail/=0) return
 call read_ps_file(ps_filename,pf,iout,ifail)
 if(ifail/=0) return
 bad=pf%pot%irel/=irel .or. pf%pot%atsym/=atsym .or. differ(pf%pot%zz,zz) &
&   .or. pf%pot%mmax/=mmax .or. differ(pf%pot%c_light,c_light_f)
 if(present(iexc)) then
   if(bad .or. pf%pot%iexc/=iexc) then
     write(iout,'(3a)') prog,': ERROR the files differ in ', &
&         'relativistic treatment, atom, mesh size, c_light or functional'
     ifail=1
     return
   end if
 else if(bad) then
   write(iout,'(3a)') prog,': ERROR the files differ in relativistic ', &
&       'treatment, atom, mesh size or c_light'
   ifail=1
   return
 end if
 if(any(abs(pf%pot%rr(:)-rr(:))>rtol*rr(:))) then
   write(iout,'(2a)') prog,': ERROR the files'' radial meshes differ'
   ifail=1
   return
 end if
 c_light=c_light_f

 if(irel==0) then
   treatment='non-relativistic'
 else if(irel==1) then
   treatment='scalar-relativistic'
 else
   treatment='fully-relativistic'
 end if
 call write_build_info(6,prog//'.x',trim(treatment))
 end subroutine load_pair

 end module ae_ps_pair_m
