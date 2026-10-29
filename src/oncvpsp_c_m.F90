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
 use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char, c_size_t, &
&                 c_double, c_associated, c_loc, c_f_pointer
 use input_m, only: input_t
 use results_m, only: results_t
 use failure_m, only: failure_t
 use constants_m, only: dp, mxlin
 implicit none
 private
 public :: oncv_api_version, oncv_new, oncv_free, oncv_set_output, &
&          oncv_run_toml, oncv_failure, oncv_dims, dims_t, oncv_field_name, &
&          oncv_shape, oncv_get_double, oncv_get_int

!> the version of the interface: a change to a function or a struct
!> raises it
 integer(c_int), parameter :: api_version=1

!> 1 if c_double is the kind of the reals, else a division by zero that
!> stops the compilation
 integer, parameter :: same_double=1/merge(1,0,c_double==dp)

!> the largest rank of a field
 integer, parameter :: mxrank=5

!> the fields of the results that oncv_shape and oncv_get_* copy out, each
!> its kind (r real, i integer, a logical as 0 or 1) and its name, the
!> component of results_t (pot.*, tests.*: pseudo_m, test_results_m)
 character(len=32), parameter :: fields(*)=[character(len=32) :: &
&  'r pot.zz', 'r pot.zion', 'r pot.c_light', 'r pot.rcmax', 'r pot.rr', &
&  'i pot.nproj', 'r pot.evkb', 'r pot.vkb', 'r pot.esr', 'r pot.eso', &
&  'r pot.vpuns', 'r pot.rho', 'r pot.rhomod', 'r epstot', 'r ea', 'r rhotae', &
&  'r rhoc', 'r rcore', 'r vpsml', 'r uupsa', 'i tests.logder.ne', &
&  'r tests.logder.ee', 'r tests.logder.rl', 'r tests.logder.ae', &
&  'r tests.logder.ps', 'i tests.ghosts.nrow', 'i tests.ghosts.nn', &
&  'r tests.ghosts.ee', 'r tests.ghosts.least', 'i tests.configs.ncfg', &
&  'i tests.configs.nrow', 'i tests.configs.state', 'r tests.configs.ee', &
&  'r tests.configs.de', 'r tests.convergence.cross', &
&  'r tests.convergence.curve', 'r tests.convergence.tail', &
&  'r tests.convergence.opt', 'r tests.convergence.wshare', &
&  'i tests.weighted_cut.used', 'r tests.weighted_cut.errw', &
&  'r tests.weighted_cut.ecutw', 'r tests.weighted_cut.weight', &
&  'r tests.vloc_cut.cut', 'r tests.vloc_cut.check', 'r tests.bbeig', &
&  'i tests.ghosts.nminus', 'r tests.ghosts.minus', 'i tests.ghosts.nplus', &
&  'r tests.ghosts.plus', 'r tests.ghosts.ecut', 'i tests.ghosts.basis', &
&  'r core.teter', 'r core.budget.ecut', 'r core.budget.error', &
&  'r core.budget.slack', 'r core.budget.ecutw', 'r core.budget.factor', &
&  'r core.budget.r_s', 'r core.budget.rc_min', 'r core.budget.r_core', &
&  'i core.budget.fill', 'r core.budget.cpoly', 'r core.budget.martins', &
&  'r core.budget.ecut_needed', 'r core.budget.r_s6', &
&  'r core.budget.mismatch', 'r core.budget.crad', 'i core.budget.cfill', &
&  'r core.budget.cmis']

!> a field copied out of the results: its values, flat, and its shape
 type :: field_t
!> whether the name is a field
   logical :: found=.false.
!> whether its values are integers (iv), else reals (rv)
   logical :: isint=.false.
!> its rank, 0 for a scalar
   integer :: rank=0
!> its shape, shp(1:rank)
   integer :: shp(mxrank)=0
!> its real values, in Fortran order
   real(dp), allocatable :: rv(:)
!> its integer values, in Fortran order
   integer, allocatable :: iv(:)
 end type field_t

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

!> a field's values and shape from an array or scalar of the results
 interface put
   module procedure put_r0, put_r1, put_r2, put_r3, put_r4, put_r5, &
&                   put_i0, put_i1, put_i2, put_i3, put_l0
 end interface put

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

 function oncv_write_files(handle,format,nformat,path,npath) result(ifail) &
&         bind(c,name='oncv_write_files')

!! writes the pseudopotential of the handle's run as the file format
!! ("psp8", "upf" or "psml") to path: the lines the executables append to
!! the text output between its Begin line and END_PSP (write_psp_file),
!! without the blank line before END_PSP
!!
!! @param handle  the run handle
!! @param format  the format, nformat characters, not terminated
!! @param nformat  its length
!! @param path  the file name, npath characters, not terminated
!! @param npath  its length

 use psp_files_m, only: write_psp_file
 use failure_m, only: fail_input, fail_io, failed, msg_c_no_results, &
&                    msg_c_format, msg_c_file
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 type(c_ptr), value :: handle
 character(kind=c_char), intent(in) :: format(nformat),path(npath)
 integer(c_size_t), value :: nformat,npath

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 character(len=:), allocatable :: kind,file,line,held
 integer :: iunit,iscr,ios,jfail
 logical :: inside,hold

 call c_f_pointer(handle,hh)
 if(.not. hh%done) then
   write(error_unit,'(/a)') 'oncv_write_files: ERROR the handle has no results'
   ifail=failed(fail_input,msg_c_no_results)
   return
 end if
 kind=from_c(format,nformat)
 if(kind/='psp8' .and. kind/='upf' .and. kind/='psml') then
   write(error_unit,'(/2a)') &
&       'oncv_write_files: ERROR the format is psp8, upf or psml, not ',kind
   ifail=failed(fail_input,msg_c_format)
   return
 end if
 file=from_c(path,npath)
 open(newunit=iunit,file=file,action='write',status='replace',iostat=ios)
 if(ios/=0) then
   write(error_unit,'(/2a)') 'oncv_write_files: ERROR cannot open ',file
   ifail=failed(fail_io,msg_c_file)
   return
 end if
 open(newunit=iscr,status='scratch',action='readwrite')
 call write_psp_file(kind,hh%inp,hh%psp,iscr,jfail)
 if(jfail==0) then
   rewind(iscr)
   inside=.false.
   hold=.false.
   do
     call read_line(iscr,line,ios)
     if(ios/=0) exit
     if(.not. inside) then
       inside=index(line,'Begin ')==1
       cycle
     end if
     if(line=='END_PSP') exit
     if(hold) write(iunit,'(a)') held
     held=line
     hold=.true.
   end do
   if(hold .and. len(held)>0) write(iunit,'(a)') held
 end if
 close(iscr)
 close(iunit)
 ifail=jfail
 return
 end function oncv_write_files

 subroutine read_line(iunit,line,ios)

!! the next line of iunit, whole, its trailing blanks kept
!!
!! @param iunit  the unit
!! @param line  the line
!! @param ios  0, or the iostat at the end of the file

!Input variables
 integer, intent(in) :: iunit

!Output variables
 character(len=:), allocatable, intent(out) :: line
 integer, intent(out) :: ios

!Local variables
 character(len=256) :: buf
 integer :: nread

 line=''
 do
   read(iunit,'(a)',advance='no',size=nread,iostat=ios) buf
   line=line//buf(1:nread)
   if(is_iostat_eor(ios)) then
     ios=0
     return
   end if
   if(ios/=0) return
 end do
 end subroutine read_line

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

 function oncv_field_name(ifield,name,nname) result(nout) &
&         bind(c,name='oncv_field_name')

!! the field ifield (1, 2, ...) of the results that oncv_shape and
!! oncv_get_* copy out: its kind (r real, i integer), a blank and its name,
!! into name, nname bytes at most, not terminated; the number of bytes
!! written, 0 past the last field
!!
!! @param ifield  the field
!! @param name  its kind and name
!! @param nname  the size of name

!Input variables
 integer(c_int), value :: ifield
 integer(c_size_t), value :: nname

!Input/output variables
 character(kind=c_char), intent(inout) :: name(nname)

!Output variables
 integer(c_int) :: nout

!Local variables
 integer :: ii

 nout=0
 if(ifield<1 .or. ifield>size(fields)) return
 nout=int(min(int(len_trim(fields(ifield)),c_size_t),nname),c_int)
 do ii=1,nout
   name(ii)=fields(ifield)(ii:ii)
 end do
 return
 end function oncv_field_name

 function oncv_shape(handle,name,nname,rank,shp) result(ifail) &
&         bind(c,name='oncv_shape')

!! the rank and shape of a field of the handle's results (oncv_field_name),
!! in Fortran order (the first index fastest); a field the run did not
!! compute has the size 0
!!
!! @param handle  the run handle
!! @param name  the field's name, nname characters, not terminated
!! @param nname  its length
!! @param rank  its rank, 0 for a scalar
!! @param shp  its shape, shp(1:rank)

!Input variables
 type(c_ptr), value :: handle
 integer(c_size_t), value :: nname
 character(kind=c_char), intent(in) :: name(nname)

!Output variables
 integer(c_int), intent(out) :: rank
 integer(c_int), intent(out) :: shp(mxrank)
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 type(field_t) :: ff

 rank=0
 shp(:)=0
 call c_f_pointer(handle,hh)
 call lookup(hh,from_c(name,nname),ff,ifail)
 if(ifail/=0) return
 rank=ff%rank
 shp(:)=ff%shp(:)
 return
 end function oncv_shape

 function oncv_get_double(handle,name,nname,values,nvalues) result(ifail) &
&         bind(c,name='oncv_get_double')

!! the values of a real field of the handle's results, in Fortran order,
!! into values, whose size nvalues must be the field's (oncv_shape)
!!
!! @param handle  the run handle
!! @param name  the field's name, nname characters, not terminated
!! @param nname  its length
!! @param values  the values
!! @param nvalues  their number

!Input variables
 type(c_ptr), value :: handle
 integer(c_size_t), value :: nname,nvalues
 character(kind=c_char), intent(in) :: name(nname)

!Input/output variables
 real(c_double), intent(inout) :: values(nvalues)

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 type(field_t) :: ff

 call c_f_pointer(handle,hh)
 call lookup(hh,from_c(name,nname),ff,ifail)
 if(ifail/=0) return
 call check_size(from_c(name,nname),ff,.false.,nvalues,ifail)
 if(ifail/=0) return
 values(:)=ff%rv(:)
 return
 end function oncv_get_double

 function oncv_get_int(handle,name,nname,values,nvalues) result(ifail) &
&         bind(c,name='oncv_get_int')

!! the values of an integer field of the handle's results (a logical as 0
!! or 1), in Fortran order, into values, whose size nvalues must be the
!! field's (oncv_shape)
!!
!! @param handle  the run handle
!! @param name  the field's name, nname characters, not terminated
!! @param nname  its length
!! @param values  the values
!! @param nvalues  their number

!Input variables
 type(c_ptr), value :: handle
 integer(c_size_t), value :: nname,nvalues
 character(kind=c_char), intent(in) :: name(nname)

!Input/output variables
 integer(c_int), intent(inout) :: values(nvalues)

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 type(field_t) :: ff

 call c_f_pointer(handle,hh)
 call lookup(hh,from_c(name,nname),ff,ifail)
 if(ifail/=0) return
 call check_size(from_c(name,nname),ff,.true.,nvalues,ifail)
 if(ifail/=0) return
 values(:)=ff%iv(:)
 return
 end function oncv_get_int

 subroutine check_size(name,ff,isint,nvalues,ifail)

!! whether a caller's buffer of nvalues values of the kind isint fits the
!! field ff
!!
!! @param name  the field's name
!! @param ff  the field
!! @param isint  whether the buffer holds integers
!! @param nvalues  its size
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_c_kind, msg_c_size
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 character(len=*), intent(in) :: name
 type(field_t), intent(in) :: ff
 logical, intent(in) :: isint
 integer(c_size_t), intent(in) :: nvalues

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: nff

 ifail=0
 if(ff%isint.neqv.isint) then
   write(error_unit,'(/3a)') 'check_size: ERROR the field ',trim(name), &
&       ' is of the other kind (oncv_get_double, oncv_get_int)'
   ifail=failed(fail_input,msg_c_kind)
   return
 end if
 nff=product(ff%shp(1:ff%rank))
 if(nvalues/=int(nff,c_size_t)) then
   write(error_unit,'(/3a,i0,a,i0)') 'check_size: ERROR the field ', &
&       trim(name),' has ',nff,' values, the buffer ',nvalues
   ifail=failed(fail_input,msg_c_size)
   return
 end if
 return
 end subroutine check_size

 subroutine lookup(hh,name,ff,ifail)

!! the field name of the handle's results, copied (fields)
!!
!! @param hh  the handle
!! @param name  the field's name
!! @param ff  the field
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_c_field, msg_c_no_results
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 type(handle_t), intent(in) :: hh
 character(len=*), intent(in) :: name

!Output variables
 type(field_t), intent(out) :: ff
 integer, intent(out) :: ifail

 ifail=0
 if(.not. hh%done) then
   write(error_unit,'(/a)') 'lookup: ERROR the handle has no results'
   ifail=failed(fail_input,msg_c_no_results)
   return
 end if
 associate(pot=>hh%psp%pot,tt=>hh%psp%tests,bb=>hh%psp%core%budget)
   select case(name)
   case('pot.zz')
     call put(ff,pot%zz)
   case('pot.zion')
     call put(ff,pot%zion)
   case('pot.c_light')
     call put(ff,pot%c_light)
   case('pot.rcmax')
     call put(ff,pot%rcmax)
   case('pot.rr')
     if(allocated(pot%rr)) then
       call put(ff,pot%rr)
     else
       call put_none(ff,1)
     end if
   case('pot.nproj')
     call put(ff,pot%nproj)
   case('pot.evkb')
     if(allocated(pot%evkb)) then
       call put(ff,pot%evkb)
     else
       call put_none(ff,3)
     end if
   case('pot.vkb')
     if(allocated(pot%vkb)) then
       call put(ff,pot%vkb)
     else
       call put_none(ff,4)
     end if
   case('pot.esr')
     if(allocated(pot%esr)) then
       call put(ff,pot%esr)
     else
       call put_none(ff,2)
     end if
   case('pot.eso')
     if(allocated(pot%eso)) then
       call put(ff,pot%eso)
     else
       call put_none(ff,2)
     end if
   case('pot.vpuns')
     if(allocated(pot%vpuns)) then
       call put(ff,pot%vpuns)
     else
       call put_none(ff,2)
     end if
   case('pot.rho')
     if(allocated(pot%rho)) then
       call put(ff,pot%rho)
     else
       call put_none(ff,1)
     end if
   case('pot.rhomod')
     if(allocated(pot%rhomod)) then
       call put(ff,pot%rhomod)
     else
       call put_none(ff,2)
     end if
   case('epstot')
     call put(ff,hh%psp%epstot)
   case('ea')
     call put(ff,hh%psp%ea)
   case('rhotae')
     if(allocated(hh%psp%rhotae)) then
       call put(ff,hh%psp%rhotae)
     else
       call put_none(ff,1)
     end if
   case('rhoc')
     if(allocated(hh%psp%rhoc)) then
       call put(ff,hh%psp%rhoc)
     else
       call put_none(ff,1)
     end if
   case('rcore')
     call put(ff,hh%psp%rcore)
   case('vpsml')
     if(allocated(hh%psp%vpsml)) then
       call put(ff,hh%psp%vpsml)
     else
       call put_none(ff,3)
     end if
   case('uupsa')
     if(allocated(hh%psp%uupsa)) then
       call put(ff,hh%psp%uupsa)
     else
       call put_none(ff,3)
     end if
   case('tests.logder.ne')
     call put(ff,tt%logder%ne)
   case('tests.logder.ee')
     if(allocated(tt%logder%ee)) then
       call put(ff,tt%logder%ee)
     else
       call put_none(ff,1)
     end if
   case('tests.logder.rl')
     call put(ff,tt%logder%rl)
   case('tests.logder.ae')
     if(allocated(tt%logder%ae)) then
       call put(ff,tt%logder%ae)
     else
       call put_none(ff,3)
     end if
   case('tests.logder.ps')
     if(allocated(tt%logder%ps)) then
       call put(ff,tt%logder%ps)
     else
       call put_none(ff,3)
     end if
   case('tests.ghosts.nrow')
     call put(ff,tt%ghosts%nrow)
   case('tests.ghosts.nn')
     call put(ff,tt%ghosts%nn)
   case('tests.ghosts.ee')
     call put(ff,tt%ghosts%ee)
   case('tests.ghosts.least')
     call put(ff,tt%ghosts%least)
   case('tests.configs.ncfg')
     call put(ff,tt%configs%ncfg)
   case('tests.configs.nrow')
     call put(ff,tt%configs%nrow)
   case('tests.configs.state')
     call put(ff,tt%configs%state)
   case('tests.configs.ee')
     call put(ff,tt%configs%ee)
   case('tests.configs.de')
     call put(ff,tt%configs%de)
   case('tests.convergence.cross')
     if(allocated(tt%convergence%cross)) then
       call put(ff,tt%convergence%cross)
     else
       call put_none(ff,5)
     end if
   case('tests.convergence.curve')
     if(allocated(tt%convergence%curve)) then
       call put(ff,tt%convergence%curve)
     else
       call put_none(ff,5)
     end if
   case('tests.convergence.tail')
     if(allocated(tt%convergence%tail)) then
       call put(ff,tt%convergence%tail)
     else
       call put_none(ff,4)
     end if
   case('tests.convergence.opt')
     if(allocated(tt%convergence%opt)) then
       call put(ff,tt%convergence%opt)
     else
       call put_none(ff,4)
     end if
   case('tests.convergence.wshare')
     if(allocated(tt%convergence%wshare)) then
       call put(ff,tt%convergence%wshare)
     else
       call put_none(ff,4)
     end if
   case('tests.weighted_cut.used')
     call put(ff,tt%weighted_cut%used)
   case('tests.weighted_cut.errw')
     call put(ff,tt%weighted_cut%errw)
   case('tests.weighted_cut.ecutw')
     call put(ff,tt%weighted_cut%ecutw)
   case('tests.weighted_cut.weight')
     if(allocated(tt%weighted_cut%weight)) then
       call put(ff,tt%weighted_cut%weight)
     else
       call put_none(ff,3)
     end if
   case('tests.vloc_cut.cut')
     call put(ff,tt%vloc_cut%cut)
   case('tests.vloc_cut.check')
     call put(ff,tt%vloc_cut%check)
   case('tests.bbeig')
     if(allocated(tt%bbeig)) then
       call put(ff,tt%bbeig)
     else
       call put_none(ff,3)
     end if
   case('tests.ghosts.nminus')
     call put(ff,tt%ghosts%nminus)
   case('tests.ghosts.minus')
     call put(ff,tt%ghosts%minus)
   case('tests.ghosts.nplus')
     call put(ff,tt%ghosts%nplus)
   case('tests.ghosts.plus')
     call put(ff,tt%ghosts%plus)
   case('tests.ghosts.ecut')
     call put(ff,tt%ghosts%ecut)
   case('tests.ghosts.basis')
     call put(ff,tt%ghosts%basis)
   case('core.teter')
     call put(ff,hh%psp%core%teter)
   case('core.budget.ecut')
     call put(ff,bb%ecut)
   case('core.budget.error')
     call put(ff,bb%error)
   case('core.budget.slack')
     call put(ff,bb%slack)
   case('core.budget.ecutw')
     call put(ff,bb%ecutw)
   case('core.budget.factor')
     call put(ff,bb%factor)
   case('core.budget.r_s')
     call put(ff,bb%r_s)
   case('core.budget.rc_min')
     call put(ff,bb%rc_min)
   case('core.budget.r_core')
     call put(ff,bb%r_core)
   case('core.budget.fill')
     call put(ff,bb%fill)
   case('core.budget.cpoly')
     call put(ff,bb%cpoly)
   case('core.budget.martins')
     call put(ff,bb%martins)
   case('core.budget.ecut_needed')
     call put(ff,bb%ecut_needed)
   case('core.budget.r_s6')
     call put(ff,bb%r_s6)
   case('core.budget.mismatch')
     call put(ff,bb%mismatch)
   case('core.budget.crad')
     if(allocated(bb%crad)) then
       call put(ff,bb%crad)
     else
       call put_none(ff,1)
     end if
   case('core.budget.cfill')
     if(allocated(bb%cfill)) then
       call put(ff,bb%cfill)
     else
       call put(ff,[integer ::])
     end if
   case('core.budget.cmis')
     if(allocated(bb%cmis)) then
       call put(ff,bb%cmis)
     else
       call put_none(ff,1)
     end if
   case default
     write(error_unit,'(/3a)') 'lookup: ERROR no field ',name, &
&         ' (oncv_field_name)'
     ifail=failed(fail_input,msg_c_field)
     return
   end select
 end associate
 return
 end subroutine lookup

 subroutine put_none(ff,rank)

!! a field the run did not compute: rank rank, size 0
!!
!! @param ff  the field
!! @param rank  its rank

!Input variables
 integer, intent(in) :: rank

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=rank
 allocate(ff%rv(0))
 return
 end subroutine put_none

 subroutine put_r0(ff,xx)

!! a real field of rank 0
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=0
 ff%rv=[xx]
 return
 end subroutine put_r0

 subroutine put_r1(ff,xx)

!! a real field of rank 1
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=1
 ff%shp(1:1)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r1

 subroutine put_r2(ff,xx)

!! a real field of rank 2
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=2
 ff%shp(1:2)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r2

 subroutine put_r3(ff,xx)

!! a real field of rank 3
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:,:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=3
 ff%shp(1:3)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r3

 subroutine put_r4(ff,xx)

!! a real field of rank 4
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:,:,:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=4
 ff%shp(1:4)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r4

 subroutine put_r5(ff,xx)

!! a real field of rank 5
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:,:,:,:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=5
 ff%shp(1:5)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r5

 subroutine put_i0(ff,xx)

!! an integer field of rank 0
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 integer, intent(in) :: xx

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=0
 ff%iv=[xx]
 return
 end subroutine put_i0

 subroutine put_i1(ff,xx)

!! an integer field of rank 1
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 integer, intent(in) :: xx(:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=1
 ff%shp(1:1)=shape(xx)
 ff%iv=reshape(xx,[size(xx)])
 return
 end subroutine put_i1

 subroutine put_i2(ff,xx)

!! an integer field of rank 2
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 integer, intent(in) :: xx(:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=2
 ff%shp(1:2)=shape(xx)
 ff%iv=reshape(xx,[size(xx)])
 return
 end subroutine put_i2

 subroutine put_i3(ff,xx)

!! an integer field of rank 3
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 integer, intent(in) :: xx(:,:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=3
 ff%shp(1:3)=shape(xx)
 ff%iv=reshape(xx,[size(xx)])
 return
 end subroutine put_i3

 subroutine put_l0(ff,xx)

!! a logical field, as the integer 0 or 1
!!
!! @param ff  the field
!! @param xx  its value

!Input variables
 logical, intent(in) :: xx

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%iv=[merge(1,0,xx)]
 return
 end subroutine put_l0

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
