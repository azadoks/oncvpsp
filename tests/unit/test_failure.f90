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
!> Unit test of the failure record (failure_m): failed's ifail split back
!> into its kind, message and status, and two runs that fail, at the input
!> and in a channel, with fail_context saying which message and where.
!>
!> usage: test_failure INPUT.dat (08_O: lloc 4, an occupied l=1)

 program test_failure
 use unit_check_m, only: check, check_done
 use input_m, only: input_t
 use input_text_m, only: read_input_text
 use oncvpsp_run_m, only: oncvpsp_run
 use results_m, only: results_t
 use failure_m, only: failure_t, failed, fail_context, fail_kind, &
&                     fail_msg, fail_ierr, fail_input, fail_convergence, &
&                     ierr_not_converged, stage_input, stage_channel, &
&                     msg_sratom_conv, msg_run_well_norm_free, &
&                     msg_ecut_weight_unoccupied
 implicit none

 type(results_t) :: psp
 type(input_t) :: inp
 type(failure_t) :: at,fl
 integer :: iout,ifail
 character(len=1024) :: input_filename

 call get_command_argument(1,input_filename)
 open(newunit=iout,file='/dev/null',action='write')

! the parts of an ifail
 ifail=failed(fail_convergence,msg_sratom_conv,ierr_not_converged)
 call check('kind',fail_kind(ifail)==fail_convergence)
 call check('message',fail_msg(ifail)==msg_sratom_conv)
 call check('status',fail_ierr(ifail)==ierr_not_converged)
 ifail=failed(fail_input,msg_run_well_norm_free)
 call check('no status',fail_ierr(ifail)==0 .and. &
&              fail_msg(ifail)==msg_run_well_norm_free)
 ifail=failed(fail_input,msg_run_well_norm_free,-3)
 call check('negative status dropped',fail_ierr(ifail)==0 .and. &
&              fail_kind(ifail)==fail_input)

! a failure of the input: the free tail's norm over all space
 call read_input(inp)
 inp%well_tail=1
 inp%well_norm=1
 call oncvpsp_run(0,0,inp,'','',.false.,psp,at,iout,ifail)
 fl=fail_context(ifail,at)
 call check('input: fails',ifail/=0)
 call check('input: kind',fl%kind==fail_input)
 call check('input: message',fl%msg==msg_run_well_norm_free)
 call check('input: stage',fl%stage==stage_input .and. fl%ll==-1)

! a failure in channel l=1: an ecut_weight for an occupied channel
 call read_input(inp)
 inp%ecut_weight(2)=1.0d0
 call oncvpsp_run(0,0,inp,'','',.false.,psp,at,iout,ifail)
 fl=fail_context(ifail,at)
 call check('channel: fails',ifail/=0)
 call check('channel: kind',fl%kind==fail_input)
 call check('channel: message',fl%msg==msg_ecut_weight_unoccupied)
 call check('channel: stage',fl%stage==stage_channel .and. fl%ll==1)
 close(iout)

 call check_done('failure')

 contains

 subroutine read_input(inp)
!! the input of the test
!!
!! @param inp  the input values (input_m)
 type(input_t), intent(out) :: inp
 integer :: unit,inline
 open(newunit=unit,file=trim(input_filename),status='old',action='read')
 call read_input_text(unit,inline,inp,iout,ifail)
 close(unit)
 if(ifail/=0) error stop 'test_failure: read_input_text failed'
 end subroutine read_input

 end program test_failure
