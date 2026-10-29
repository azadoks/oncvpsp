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
!> The C interface of the library (include/oncvpsp.h): a run handle that
!> holds one calculation, its input given as TOML text, its results and its
!> failure. The handle is the only state: two handles are two independent
!> calculations, in one thread or in several. Every function returns 0 or
!> the ifail of the failure (failure_m); a binding raises from the kind and
!> the message of the record (oncv_failure).

 module oncvpsp_c_m
 use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_null_ptr, c_char, &
&                 c_size_t, c_associated, c_loc, c_f_pointer
 use input_m, only: input_t
 use results_m, only: results_t
 use failure_m, only: failure_t
 use constants_m, only: mxlin
 implicit none
 private
 public :: oncv_api_version, oncv_new, oncv_free, oncv_set_output, &
&          oncv_run_toml, oncv_failure, oncv_dims, dims_t

!> the version of the interface: a change to a function or a struct
!> raises it
 integer(c_int), parameter :: api_version=1

!> the sizes of a run's results (oncv_dims), for the caller's arrays
 type, bind(c) :: dims_t
!> the version of the interface that filled it (api_version)
   integer(c_int) :: version
!> 0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
   integer(c_int) :: irel
!> number of kappas per l: 1, or 2 for irel=2
   integer(c_int) :: nkap
!> the functional (iexc)
   integer(c_int) :: iexc
!> size of the log radial mesh
   integer(c_int) :: mmax
!> the largest l of the projectors
   integer(c_int) :: lmax
!> the l of the local potential, 4 for the polynomial one
   integer(c_int) :: lloc
!> the projectors of each l, l=0..mxlin
   integer(c_int) :: nproj(mxlin+1)
!> points of the linear output mesh
   integer(c_int) :: nrl
 end type dims_t

!> a calculation: its settings, input, results and failure
 type :: handle_t
!> the file of the text output
   character(len=:), allocatable :: out_path
!> the input values, completed by the run
   type(input_t) :: inp
!> the results
   type(results_t) :: psp
!> where the run was when it failed
   type(failure_t) :: at=failure_t(0,0,0,0,-1,0,0)
!> the run's ifail, 0 before a run
   integer :: ifail=0
!> the run's relativistic treatment
   integer :: irel=0
!> whether a run finished, so that the results are there
   logical :: done=.false.
 end type handle_t

 contains

 function oncv_api_version() result(version) bind(c,name='oncv_api_version')

!! the version of the interface (api_version)

!Output variables
 integer(c_int) :: version

 version=api_version
 return
 end function oncv_api_version

 function oncv_new() result(handle) bind(c,name='oncv_new')

!! a new run handle, its text output going nowhere (/dev/null); free it
!! with oncv_free

!Output variables
 type(c_ptr) :: handle

!Local variables
 type(handle_t), pointer :: hh

 allocate(hh)
 hh%out_path='/dev/null'
 handle=c_loc(hh)
 return
 end function oncv_new

 subroutine oncv_free(handle) bind(c,name='oncv_free')

!! frees a run handle and all it holds; a null handle is ignored
!!
!! @param handle  the run handle

!Input variables
 type(c_ptr), value :: handle

!Local variables
 type(handle_t), pointer :: hh

 if(.not. c_associated(handle)) return
 call c_f_pointer(handle,hh)
 deallocate(hh)
 return
 end subroutine oncv_free

 function oncv_set_output(handle,path,npath) result(ifail) &
&         bind(c,name='oncv_set_output')

!! the file the runs of the handle write their text output to, replaced
!! at each run (/dev/null for none)
!!
!! @param handle  the run handle
!! @param path  the file name, npath characters, not terminated
!! @param npath  its length

!Input variables
 type(c_ptr), value :: handle
 character(kind=c_char), intent(in) :: path(npath)
 integer(c_size_t), value :: npath

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh

 ifail=0
 call c_f_pointer(handle,hh)
 hh%out_path=from_c(path,npath)
 return
 end function oncv_set_output

 function oncv_run_toml(handle,irel,text,ntext,skip_ld) result(ifail) &
&         bind(c,name='oncv_run_toml')

!! the pseudopotential and its tests of the TOML input text, with the
!! relativistic treatment irel, into the handle, whose earlier results it
!! replaces (oncvpsp_run)
!!
!! @param handle  the run handle
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param text  the TOML input, ntext characters, not terminated
!! @param ntext  its length
!! @param skip_ld  non-zero to leave out the log derivatives

 use oncvpsp_run_m, only: oncvpsp_run
#ifdef WITH_TOML
 use input_toml_m, only: read_input_toml_string
#endif
 use failure_m, only: fail_input, fail_io, failed, msg_c_output, &
&                    msg_c_no_toml, stage_input
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 type(c_ptr), value :: handle
 integer(c_int), value :: irel,skip_ld
 character(kind=c_char), intent(in) :: text(ntext)
 integer(c_size_t), value :: ntext

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 integer :: iout,ios,jfail

 call c_f_pointer(handle,hh)
 hh%done=.false.
 hh%irel=irel
 hh%at=failure_t(0,0,0,stage_input,-1,0,0)
 open(newunit=iout,file=hh%out_path,action='write',status='replace', &
&     iostat=ios)
 if(ios/=0) then
   write(error_unit,'(/2a)') &
&       'oncv_run_toml: ERROR cannot open the text output ',hh%out_path
   hh%ifail=failed(fail_io,msg_c_output)
   ifail=hh%ifail
   return
 end if
#ifdef WITH_TOML
 call read_input_toml_string(from_c(text,ntext),hh%inp,iout,jfail)
 if(jfail==0) then
   call oncvpsp_run(irel,0,hh%inp,'','',skip_ld/=0,hh%psp,hh%at,iout,jfail)
 end if
#else
 write(iout,'(/a)') 'oncv_run_toml: ERROR TOML input needs a build with TOML'
 jfail=failed(fail_input,msg_c_no_toml)
#endif
 close(iout)
 hh%ifail=jfail
 hh%done=jfail==0
 ifail=jfail
 return
 end function oncv_run_toml

 function oncv_failure(handle,failure) result(ifail) &
&         bind(c,name='oncv_failure')

!! the failure of the handle's last run (fail_context): its kind, message
!! and status, and where the run was; kind 0 if it did not fail
!!
!! @param handle  the run handle
!! @param failure  the record

 use failure_m, only: fail_context

!Input variables
 type(c_ptr), value :: handle

!Output variables
 type(failure_t), intent(out) :: failure
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh

 call c_f_pointer(handle,hh)
 failure=fail_context(hh%ifail,hh%at)
 ifail=0
 return
 end function oncv_failure

 function oncv_dims(handle,dims) result(ifail) bind(c,name='oncv_dims')

!! the sizes of the handle's results, for the caller's arrays; all 0 but
!! the version without a finished run
!!
!! @param handle  the run handle
!! @param dims  the sizes

!Input variables
 type(c_ptr), value :: handle

!Output variables
 type(dims_t), intent(out) :: dims
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh

 call c_f_pointer(handle,hh)
 dims=dims_t(api_version,0,0,0,0,0,0,0,0)
 ifail=0
 if(.not. hh%done) return
 associate(pot=>hh%psp%pot)
   dims%irel=pot%irel
   dims%nkap=pot%nkap
   dims%iexc=pot%iexc
   dims%mmax=pot%mmax
   dims%lmax=pot%lmax
   dims%lloc=pot%lloc
   dims%nproj(:)=pot%nproj(:)
 end associate
 dims%nrl=hh%psp%nrl
 return
 end function oncv_dims

 function from_c(chars,nchars) result(string)

!! a Fortran string of nchars C characters
!!
!! @param chars  the characters
!! @param nchars  their number

!Input variables
 character(kind=c_char), intent(in) :: chars(nchars)
 integer(c_size_t), intent(in) :: nchars

!Output variables
 character(len=:), allocatable :: string

!Local variables
 integer(c_size_t) :: ii

 allocate(character(len=nchars) :: string)
 do ii=1,nchars
   string(ii:ii)=chars(ii)
 end do
 return
 end function from_c

 end module oncvpsp_c_m
