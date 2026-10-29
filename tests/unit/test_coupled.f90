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
!> Unit test of run_coupled called again on one run's stages
!> (run_stages): with the input's model core, then with none (icmod 0),
!> then with the input's again, each bit for bit what a whole run with
!> that model core gives, so the stages' results are not changed and the
!> model core can be redone alone (the budget's probe and its run).
!>
!> usage: test_coupled IREL INPUT.dat

 program test_coupled
 use constants_m, only: mxcnf
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check_done, nfail
 use input_m, only: input_t, config_t
 use input_text_m, only: read_input_text
 use oncvpsp_run_m, only: oncvpsp_run, run_stages, run_coupled, stages_t
 use results_m, only: results_t
 use failure_m, only: failure_t
 implicit none

 type(results_t) :: psp0,psp1,psp2,psp3,psp4
 type(input_t) :: inp,inp1,inp2
 type(config_t) :: cnf(mxcnf)
 type(stages_t) :: st
 type(failure_t) :: at
 integer :: irel,iout,ifail
 character(len=1024) :: arg,input_filename

 call get_command_argument(1,arg)
 read(arg,*) irel
 call get_command_argument(2,input_filename)
 open(newunit=iout,file='/dev/null',action='write')

! whole runs, with the input's model core and with none
 call read_input(inp)
 call oncvpsp_run(irel,0,inp,'','',.false.,psp0,at,iout,ifail)
 if(ifail/=0) error stop 'test_coupled: oncvpsp_run failed'
 call read_input(inp)
 inp%icmod=0
 call oncvpsp_run(irel,0,inp,'','',.false.,psp3,at,iout,ifail)
 if(ifail/=0) error stop 'test_coupled: oncvpsp_run (icmod 0) failed'

! the stages once, and the model core and the tests from them three times
 call read_input(inp)
 call run_stages(irel,0,inp,'',cnf,st,at,iout,ifail)
 if(ifail/=0) error stop 'test_coupled: run_stages failed'
 inp1=inp
 call run_coupled(irel,0,inp1,cnf,st,'',.false.,psp1,at,iout,ifail)
 if(ifail/=0) error stop 'test_coupled: run_coupled failed'
 inp2=inp
 inp2%icmod=0
 inp2%fcfact=0.0d0
 call run_coupled(irel,0,inp2,cnf,st,'',.false.,psp2,at,iout,ifail)
 if(ifail/=0) error stop 'test_coupled: run_coupled (icmod 0) failed'
 inp1=inp
 call run_coupled(irel,0,inp1,cnf,st,'',.false.,psp4,at,iout,ifail)
 if(ifail/=0) error stop 'test_coupled: run_coupled (again) failed'
 close(iout)

 call check_same(psp0,psp1,'the input''s model core')
 call check_same(psp3,psp2,'no model core')
 call check_same(psp0,psp4,'the input''s model core again')

 call check_done('coupled')

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
 if(ifail/=0) error stop 'test_coupled: read_input_text failed'
 end subroutine read_input

 subroutine check_same(pa,pb,what)
!! the results of two runs, bit for bit
!!
!! @param pa  those of the whole run
!! @param pb  those of run_coupled
!! @param what  the case, for the messages
 type(results_t), intent(in) :: pa,pb
 character(len=*), intent(in) :: what
 call differ(any(pa%pot%evkb/=pb%pot%evkb),what,'evkb')
 call differ(any(pa%pot%vkb/=pb%pot%vkb),what,'vkb')
 call differ(any(pa%pot%vpuns/=pb%pot%vpuns),what,'vpuns')
 call differ(any(pa%pot%rho/=pb%pot%rho),what,'rho')
 call differ(any(pa%pot%rhomod/=pb%pot%rhomod),what,'rhomod')
 call differ(pa%epstot/=pb%epstot,what,'pseudo-atom energy')
 call differ(any(pa%tests%logder%ae/=pb%tests%logder%ae),what, &
&            'log derivatives ae')
 call differ(any(pa%tests%logder%ps/=pb%tests%logder%ps),what, &
&            'log derivatives ps')
 call differ(any(pa%tests%ghosts%ee/=pb%tests%ghosts%ee),what,'ghosts')
 call differ(any(pa%tests%configs%ee/=pb%tests%configs%ee),what, &
&            'test eigenvalues')
 call differ(any(pa%tests%configs%de/=pb%tests%configs%de),what, &
&            'test energies')
 call differ(any(pa%tests%convergence%cross/= &
&                pb%tests%convergence%cross),what,'cutoffs')
 call differ(any(pa%tests%vloc_cut%cut/=pb%tests%vloc_cut%cut),what, &
&            'vloc cutoffs')
 call differ(pa%tests%weighted_cut%ecutw/=pb%tests%weighted_cut%ecutw,what, &
&            'weighted cutoff')
 end subroutine check_same

 subroutine differ(dd,what,res)
!! fails if the two runs differ in the result res
!!
!! @param dd  whether they differ
!! @param what  the case
!! @param res  the result
 logical, intent(in) :: dd
 character(len=*), intent(in) :: what,res
 if(dd) then
   nfail=nfail+1
   write(stderr,'(4a)') 'FAIL: ',what,': ',res
 end if
 end subroutine differ

 end program test_coupled
