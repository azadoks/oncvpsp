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
 use ae_file_m, only: read_ae_potential, differ
 use ps_file_m, only: ps_file_t, read_ps_file
 use relativity_m, only: set_relativity
 use build_info_m, only: write_build_info
 implicit none
 private
 public :: load_pair

!> relative tolerance of the comparison of the two files' meshes, as in
!> ae_file_m
 real(dp), parameter :: rtol=1.0d-12

 contains

 subroutine load_pair(prog,ae_filename,ps_filename,irel,zz,mmax,rr,vfull, &
&                     pf,nc,nv,iexc,na,la,fa,ea)
!! the two files of program prog, which must be of the same relativistic
!! treatment, atom, mesh and c_light (and functional, with iexc), and the
!! all-electron atom's configuration if asked; sets c_light and writes
!! the build line
 character(len=*), intent(in) :: prog,ae_filename,ps_filename
 integer, intent(out) :: irel,mmax
 real(dp), intent(out) :: zz
 real(dp), allocatable, intent(out) :: rr(:),vfull(:)
 type(ps_file_t), intent(out) :: pf
 integer, intent(out), optional :: nc,nv,iexc,na(mxstate),la(mxstate)
 real(dp), intent(out), optional :: fa(mxstate),ea(mxstate,2)
 character(len=2) :: atsym
 character(len=19) :: treatment
 real(dp) :: c_light_f
 logical :: bad

 call read_ae_potential(ae_filename,irel,atsym,zz,mmax,rr,c_light_f, &
&                       vfull,nc=nc,nv=nv,iexc=iexc,na=na,la=la,fa=fa,ea=ea)
 call read_ps_file(ps_filename,pf)
 bad=pf%irel/=irel .or. pf%atsym/=atsym .or. differ(pf%zz,zz) &
&   .or. pf%mmax/=mmax .or. differ(pf%c_light,c_light_f)
 if(present(iexc)) then
   if(bad .or. pf%iexc/=iexc) then
     write(6,'(3a)') prog,': ERROR the files differ in ', &
&         'relativistic treatment, atom, mesh size, c_light or functional'
     error stop 1
   end if
 else if(bad) then
   write(6,'(3a)') prog,': ERROR the files differ in relativistic ', &
&       'treatment, atom, mesh size or c_light'
   error stop 1
 end if
 if(any(abs(pf%rr(:)-rr(:))>rtol*rr(:))) then
   write(6,'(2a)') prog,': ERROR the files'' radial meshes differ'
   error stop 1
 end if
 call set_relativity(c_light_f)

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
