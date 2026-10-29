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
!> Unit test of check_ae: the all-electron atom solved once (solve_ae)
!> is the input's, and each change of the input check_ae looks at is
!> refused with its message: another irel, nuclear charge, configuration
!> (n, l, f, nv), functional, mesh (amesh, rr1, rmax, adams_order) and
!> speed of light.
!>
!> usage: test_check_ae IREL INPUT.dat

 program test_check_ae
 use unit_check_m, only: check, check_done
 use input_m, only: input_t, config_t
 use input_text_m, only: read_input_text
 use oncvpsp_run_m, only: complete_input
 use ae_state_m, only: ae_state_t, solve_ae, check_ae
 use failure_m, only: fail_input, fail_kind, fail_msg, msg_ae_given_irel, &
&                    msg_ae_given_atom, msg_ae_given_mesh, msg_ae_given_c_light
 use constants_m, only: mxcnf
 implicit none

 type(input_t) :: inp,other
 type(config_t) :: cnf(mxcnf)
 type(ae_state_t) :: ae
 integer :: irel,nkap,iout,ifail,ierr,unit,inline
 character(len=1024) :: arg,input_filename

 call get_command_argument(1,arg)
 read(arg,*) irel
 call get_command_argument(2,input_filename)
 nkap=1
 if(irel==2) nkap=2

 open(newunit=iout,file='/dev/null',action='write')
 open(newunit=unit,file=trim(input_filename),status='old',action='read')
 call read_input_text(unit,inline,inp,iout,ifail)
 close(unit)
 if(ifail/=0) error stop 'test_check_ae: read_input_text failed'
 call complete_input(inp,cnf,iout,ifail)
 if(ifail/=0) error stop 'test_check_ae: complete_input failed'
 call solve_ae(irel,nkap,0,inp,cnf,'',ae,ierr,iout,ifail)
 if(ifail/=0 .or. ierr/=0) error stop 'test_check_ae: solve_ae failed'

 call check_ae(ae,irel,inp,iout,ifail)
 call check('the input',ifail==0)
 call check_ae(ae,mod(irel+1,3),inp,iout,ifail)
 call refused('irel',msg_ae_given_irel)

 other=inp
 other%zz=inp%zz+1.0d0
 call check_input(other,'zz',msg_ae_given_atom)
 other=inp
 other%na(inp%nc+1)=inp%na(inp%nc+1)+1
 call check_input(other,'n',msg_ae_given_atom)
 other=inp
 other%la(inp%nc+inp%nv)=inp%la(inp%nc+inp%nv)+1
 call check_input(other,'l',msg_ae_given_atom)
 other=inp
 other%fa(inp%nc+inp%nv)=inp%fa(inp%nc+inp%nv)*0.5d0
 call check_input(other,'f',msg_ae_given_atom)
 other=inp
 other%nv=inp%nv-1
 call check_input(other,'nv',msg_ae_given_atom)
 other=inp
 other%iexc=7-inp%iexc
 call check_input(other,'iexc',msg_ae_given_atom)

 other=inp
 other%amesh=1.02d0
 call check_input(other,'amesh',msg_ae_given_mesh)
 other=inp
 other%rr1=2.0d-4
 call check_input(other,'rr1',msg_ae_given_mesh)
 other=inp
 other%rmax=0.8d0*ae%rr(ae%mmax)
 call check_input(other,'rmax',msg_ae_given_mesh)
 other=inp
 other%adams_order=6
 call check_input(other,'adams_order',msg_ae_given_mesh)

 other=inp
 other%c_light=100.0d0
 call check_input(other,'c_light',msg_ae_given_c_light)
 close(iout)

 call check_done('check_ae')

 contains

 subroutine check_input(other,what,msg)
!! check_ae of the input other refused with the message msg
 type(input_t), intent(in) :: other
 character(len=*), intent(in) :: what
 integer, intent(in) :: msg
 call check_ae(ae,irel,other,iout,ifail)
 call refused(what,msg)
 end subroutine check_input

 subroutine refused(what,msg)
!! ifail an input failure with the message msg
 character(len=*), intent(in) :: what
 integer, intent(in) :: msg
 call check('another '//what,fail_kind(ifail)==fail_input .and. &
&              fail_msg(ifail)==msg)
 end subroutine refused

 end program test_check_ae
