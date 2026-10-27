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

 module oncvpsp_main_m

 implicit none
 private
 public :: oncvpsp_main

 contains

 subroutine oncvpsp_main(irel,istage)
!! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
!! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
!! and references therein.
!!
!!   D. R. Hamann
!!   Mat-Sim Research LLC
!!   P.O. Box 742
!!   Murray Hill, NJ 07974
!!   USA
!!
!!   Developed from original "gncpp" code of March 8,1987
!!
!!   Output format for ABINIT pspcod=8 and upf format for quantumespresso
!!
!!   The whole of each executable: oncvpsp.x, oncvpspnr.x and oncvpspr.x and
!!   their _ae.x and _ps.x stages are short programs calling this with their
!!   irel and istage.  With istage=1 (the _ae.x), it solves only the
!!   all-electron atom and the all-electron halves of the test
!!   configurations, and writes them to the file given by --ae-file.  With
!!   istage=2 (the _ps.x), it reads that file instead of solving the
!!   all-electron atom and otherwise runs as with istage=0, with the same
!!   output.  Both read the same input.  See ae_file_m.  Except in the
!!   all-electron stage, --ps-file writes the pseudopotential for
!!   oncvpsp_ld.x, which takes the log derivatives from it and the
!!   all-electron file at any radii and energies (ps_file_m), and
!!   --no-log-derivatives leaves them out of this run.
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param istage  0 the whole calculation, 1 the all-electron stage only, 2 the
!!                pseudopotential stage only
 use, intrinsic :: iso_fortran_env, only: stdin => input_unit
 use input_text_m, only: read_input_text
#if (defined WITH_TOML)
 use input_toml_m, only: read_input_toml
#endif
 use constants_m, only: well_onset_default
 use main_setup_m, only: main_names, parse_command_line, write_banner, &
&                        INPUT_STDIN, INPUT_TEXT, INPUT_TOML
 use input_m, only: input_t
 use oncvpsp_run_m, only: oncvpsp_run
 use results_m, only: results_t
 use psp_files_m, only: write_psp_files
 implicit none
 integer, parameter :: iout=6

 type(input_t) :: inp
 type(results_t) :: psp

!Input variables
 integer, intent(in) :: irel,istage

!Local variables
 integer :: ifail
 integer :: ios,inline,nkap

!names in the messages: the equation solved and the program
 character(len=19) :: treatment
 character(len=9) :: pname
 character(len=16) :: exe

 integer :: input_mode
 integer :: unit
 character(len=1024) :: input_filename
 character(len=1024) :: ae_filename
 character(len=1024) :: ps_filename
 logical :: skip_ld

 call main_names(irel,istage,nkap,treatment,pname,exe)
 call parse_command_line(istage,treatment,pname,exe,input_mode, &
&                        input_filename,ae_filename,ps_filename,skip_ld)
 call write_banner(irel,treatment,exe,iout)

 if(istage/=0 .and. len_trim(ae_filename) == 0) then
   write(iout,'(/4a)') trim(pname),': ERROR ',trim(exe),' requires --ae-file'
   error stop 1
 end if

 ios = 0
 inp%ep_as_given=.false.
 inp%well_onset=well_onset_default
 inp%well_tail=0
 inp%well_n(:)=0
 inp%well_norm=0
 inp%ecut_wt(:)=-1.0d0
 inp%wreg=0.0d0
 select case(input_mode)
 case(INPUT_STDIN)
   unit=stdin
 case(INPUT_TEXT,INPUT_TOML)
   open(newunit=unit, file=input_filename, status='old', action='read', &
&       iostat=ios)
   if(ios /= 0) then
     write(iout, '(/3a)') trim(pname),': ERROR cannot open input file, file=', &
&       trim(input_filename)
     error stop 1
   end if
 case default
   write(iout,'(/2a,i0)') trim(pname), &
&        ': ERROR unknown input mode, input_mode=', &
&        input_mode
   error stop 1
 end select
 if(input_mode==INPUT_TOML) then
#if (defined WITH_TOML)
   call read_input_toml(unit, inp, iout, ifail)
   if(ifail/=0) error stop 1
#endif
 else
   call read_input_text(unit, inline, inp, iout, ifail)
   if(ifail/=0) error stop 1
 end if
 if(input_mode/=INPUT_STDIN) close(unit)

 call oncvpsp_run(irel,istage,inp,ae_filename,ps_filename,skip_ld,psp, &
&                 iout,ifail)
 if(ifail/=0) error stop 1
 if(istage==1) stop

 call write_psp_files(inp,psp,iout,ifail)
 if(ifail/=0) error stop 1

 stop

 end subroutine oncvpsp_main

 end module oncvpsp_main_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
