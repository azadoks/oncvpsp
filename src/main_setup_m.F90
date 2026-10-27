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
!> The setup of oncvpsp_main before the pseudopotential: the names of the
!> executable, its command line and help, the banner, the log mesh, the
!> reference atom's results and the mesh points of the core radii.
 module main_setup_m
 use constants_m, only: dp
 use build_info_m, only: write_build_info
 implicit none
 private
 public :: main_names, parse_command_line, write_banner, make_log_mesh
 public :: write_reference_results, core_radius_points

!> where the input is read from
 integer, parameter, public :: INPUT_STDIN=1, INPUT_TEXT=2, INPUT_TOML=3

 contains

 subroutine main_names(irel,istage,nkap,treatment,pname,exe)
!! the number of kappa of each l (see doc/fr_merge_plan.md: 1 for j=l+1/2
!! and 2 for j=l-1/2 when fully relativistic, only 1 otherwise), and the
!! names in the messages: the equation solved, the program and the
!! executable of irel and istage
 integer, intent(in) :: irel,istage
 integer, intent(out) :: nkap
 character(len=19), intent(out) :: treatment
 character(len=9), intent(out) :: pname
 character(len=16), intent(out) :: exe
 if(irel==2) then
   nkap=2
   treatment='fully-relativistic'
   pname='oncvpsp_r'
 else
   nkap=1
   if(irel==0) then
     treatment='non-relativistic'
   else
     treatment='scalar-relativistic'
   end if
   pname='oncvpsp'
 end if
 if(irel==0) then
   exe='oncvpspnr'
 else if(irel==1) then
   exe='oncvpsp'
 else
   exe='oncvpspr'
 end if
 if(istage==1) then
   exe=trim(exe)//'_ae'
 else if(istage==2) then
   exe=trim(exe)//'_ps'
 end if
 exe=trim(exe)//'.x'
 end subroutine main_names

 subroutine parse_command_line(istage,treatment,pname,exe,input_mode, &
&                              input_filename,ae_filename,ps_filename, &
&                              skip_ld)
!! the options of the executable exe of stage istage: where the input is
!! read from, the data files and whether to leave out the log derivatives
 integer, intent(in) :: istage
 character(len=*), intent(in) :: treatment,pname,exe
 integer, intent(out) :: input_mode
 character(len=1024), intent(out) :: input_filename,ae_filename,ps_filename
 logical, intent(out) :: skip_ld
 integer :: iarg
 character(len=256) :: arg
 input_mode = INPUT_STDIN
 ae_filename=''
 ps_filename=''
 skip_ld=.false.
 do iarg = 1, command_argument_count()
   call get_command_argument(iarg, arg)
   select case(arg)
   case('-h', '--help')
     call write_help(istage,treatment,exe)
     stop
   case('-v', '--version')
     write(6,'(3a)') 'ONCVPSP (',trim(treatment),') version 4.0.1'
     call write_build_info(6,trim(exe),trim(treatment))
     stop
   case('-i', '--input')
     if(iarg + 1 > command_argument_count()) then
       write(6,'(/3a)') trim(pname),': ERROR --input requires a file name'
       error stop 1
     end if
     call get_command_argument(iarg + 1, input_filename)
     input_mode = INPUT_TEXT
#if (defined WITH_TOML)
   case('-t', '--toml-input')
     if(iarg + 1 > command_argument_count()) then
       write(6,'(/3a)') trim(pname), &
&            ': ERROR --toml-input requires a file name'
       error stop 1
     end if
     call get_command_argument(iarg + 1, input_filename)
     input_mode = INPUT_TOML
#else
   case('-t', '--toml-input')
     write(6,'(/3a)') trim(pname),': ERROR TOML input is not enabled', &
&          ' in this build'
     error stop 1
#endif
   case('-a', '--ae-file')
!only the stages take an all-electron data file; the whole program ignores it
     if(istage/=0) then
       if(iarg + 1 > command_argument_count()) then
         write(6,'(/3a)') trim(pname), &
&              ': ERROR --ae-file requires a file name'
         error stop 1
       end if
       call get_command_argument(iarg + 1, ae_filename)
     end if
   case('-p', '--ps-file')
!the all-electron stage has no pseudopotential
     if(istage/=1) then
       if(iarg + 1 > command_argument_count()) then
         write(6,'(/3a)') trim(pname), &
&              ': ERROR --ps-file requires a file name'
         error stop 1
       end if
       call get_command_argument(iarg + 1, ps_filename)
     end if
   case('--no-log-derivatives')
     skip_ld=(istage/=1)
   case default
   end select
 end do
 end subroutine parse_command_line

 subroutine write_banner(irel,treatment,exe,iout)
!! the version, the citation and the build
 integer, intent(in) :: irel
 character(len=*), intent(in) :: treatment,exe
 integer, intent(in) :: iout
 if(irel==2) then
   write(iout,'(a/a//)') &
&        'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&        'relativistic version 4.0.1 03/01/2019'
 else
   write(iout,'(a/a//)') &
&        'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&        trim(treatment)//' version 4.0.1 03/01/2019'
 end if

 write(iout,'(a/a/a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&      'in any publication utilizing these pseudopotentials.'

 call write_build_info(6,trim(exe),trim(treatment))
 write(iout,'(/)')
 end subroutine write_banner

 subroutine make_log_mesh(irel,zz,amesh_in,rr1_in,rmax_in,al,mmax,rr, &
&                         log_mesh)
!! the log mesh rr of mmax points, ratio exp(al): the input values, or the
!! built-in ones for those 0 or less; and the mesh for the echoes of the
!! input, only when the input sets it (else 0), so that the output of an
!! input without it is unchanged
 integer, intent(in) :: irel
 real(dp), intent(in) :: zz,amesh_in,rr1_in,rmax_in
 real(dp), intent(out) :: al,log_mesh(3)
 integer, intent(out) :: mmax
 real(dp), allocatable, intent(out) :: rr(:)
 real(dp) :: amesh,rr1,rmax
 integer :: ii
 if(amesh_in>0.0d0) then
   amesh=amesh_in
 else
!Alternative built-in meshes: amesh=1.012d0 (coarser), or 1.003d0 and
!1.0015d0 (finer)
   amesh=1.006d0
 end if

 al=log(amesh)
 if(rr1_in>0.0d0) then
   rr1=rr1_in
 else
   rr1=0.0005d0/zz
   if(irel/=2) rr1=min(rr1,0.0005d0/10)
 end if
 if(rmax_in>0.0d0) then
   rmax=rmax_in
 else
   rmax=45.0d0
 end if
 mmax=log(rmax/rr1)/al
 log_mesh(:)=0.0d0
 if(amesh_in>0.0d0 .or. rr1_in>0.0d0 .or. rmax_in>0.0d0) then
   log_mesh(:)=[amesh,rr1,rmax]
 end if

 allocate(rr(mmax))
 do ii=1,mmax
   rr(ii)=rr1*exp(al*(ii-1))
 end do
 end subroutine make_log_mesh

 subroutine write_reference_results(pname,it,etot,res,eterm,iout)
!! the all-electron reference atom's iterations, total energy, residual
!! and energy terms; an error if it did not converge
 character(len=*), intent(in) :: pname
 integer, intent(in) :: it
 real(dp), intent(in) :: etot,res,eterm(4)
 integer, intent(in) :: iout
 write(iout,'(//a)') 'Reference configufation results'
 write(iout,'(a,i6)') '  iterations',it
 if(it >= 100) then
   write(iout,'(2a)') trim(pname), &
&        ': ERROR all-electron reference atom not converged'
   error stop 1
 end if
 write(iout,'(a,1p,d18.8)') '  all-electron total energy (Ha)',etot
 write(iout,'(a,1p,e10.2)') '  residual of the screening potential (Ha)',res
 write(iout,'(a)') '  energy terms (Ha)'
 write(iout,'(a,f22.10)') '    total               ',etot
 write(iout,'(a,f22.10)') '    kinetic             ',eterm(1)
 write(iout,'(a,f22.10)') '    electron-nucleus    ',eterm(2)
 write(iout,'(a,f22.10)') '    Hartree             ',eterm(3)
 write(iout,'(a,f22.10)') '    exchange-correlation',eterm(4)
 end subroutine write_reference_results

 subroutine core_radius_points(lmax,lloc,rc,rr,mmax,irc,rcmax)
!! the last log mesh point at or inside each rc, for the bookkeeping of the
!! routines that index the mesh (R-2); the radii themselves stay as given;
!! and the largest
 integer, intent(in) :: lmax,lloc,mmax
 real(dp), intent(in) :: rc(:),rr(mmax)
 integer, intent(out) :: irc(:)
 real(dp), intent(out) :: rcmax
 integer :: ii,l1
 rcmax=0.0d0
 irc(:)=0
 do l1=1,max(lmax+1,lloc+1)
   irc(l1)=0
   do ii=1,mmax
     if(rr(ii)<=rc(l1)) irc(l1)=ii
   end do
   if(rc(l1)>0.0d0 .and. irc(l1)<2) then
     write(6,'(a,i2,a,es12.4)') &
&         'oncvpsp: ERROR rc inside the mesh start for l=',l1-1,': rc=',rc(l1)
     error stop 1
   end if
   rcmax=max(rcmax,rc(l1))
 end do
 end subroutine core_radius_points

 subroutine write_help(istage,treatment,exe)
!! the usage of this executable: what it does and its options
 integer, intent(in) :: istage
 character(len=*), intent(in) :: treatment,exe
 character(len=16) :: full

 full=exe(1:index(exe,'_')-1)
 if(istage==0) full=exe
 write(6,'(3a)') 'Usage: ',trim(exe),' [options] < INPUT.dat'
 write(6,'(a)') ''
 if(istage==0) then
   write(6,'(3a)') 'Generates a ',trim(treatment),' optimized norm-conserving'
   write(6,'(a)') 'Vanderbilt pseudopotential from the input (D. R. Hamann, PRB'
   write(6,'(a)') '88, 085117 (2013)) and writes its tests and the'
   write(6,'(a)') 'pseudopotential file(s) of the input''s psfile to stdout.'
 else if(istage==1) then
   write(6,'(3a)') 'The all-electron stage of ',trim(full),'.x: solves the'
   write(6,'(3a)') trim(treatment),' all-electron atom and the all-electron'
   write(6,'(a)') 'halves of the test configurations, which depend on no'
   write(6,'(a)') 'pseudopotential input, and writes them to --ae-file for'
   write(6,'(3a)') trim(full),'_ps.x. Its stdout is the all-electron part of'
   write(6,'(3a)') 'the output of ',trim(full),'.x.'
 else
   write(6,'(3a)') 'The pseudopotential stage of ',trim(full),'.x: the same'
   write(6,'(a)') 'input and output, but the all-electron atom is read from'
   write(6,'(3a)') 'the --ae-file of ',trim(full),'_ae.x instead of solved,'
   write(6,'(a)') 'so a loop over pseudopotential inputs solves it once.'
 end if
 write(6,'(a)') ''
 write(6,'(a)') 'Options:'
 write(6,'(a)') '  -h, --help'
 write(6,'(a)') '      Show this help and exit'
 write(6,'(a)') '  -v, --version'
 write(6,'(a)') '      Show the version and the build, and exit'
 write(6,'(a)') '  -i, --input FILE'
 write(6,'(a)') '      Read the input, in the .dat text format, from FILE'
 write(6,'(a)') '      instead of stdin'
 write(6,'(a)') '  -t, --toml-input FILE'
 write(6,'(a)') '      Read the input from the TOML file FILE'
#if !(defined WITH_TOML)
 write(6,'(a)') '      (not available: this build has no TOML)'
#endif
 if(istage==1) then
   write(6,'(a)') '  -a, --ae-file FILE'
   write(6,'(a)') '      Write the all-electron data file FILE (required)'
 else if(istage==2) then
   write(6,'(a)') '  -a, --ae-file FILE'
   write(6,'(3a)') '      Read the all-electron data file FILE, written by ', &
&                  trim(full),'_ae.x'
   write(6,'(a)') '      for the same atom, reference configuration and mesh'
   write(6,'(a)') '      (required); the input''s test configurations must'
   write(6,'(a)') '      be among those the file was written for'
 end if
 if(istage/=1) then
   write(6,'(a)') '  -p, --ps-file FILE'
   write(6,'(a)') '      Also write the pseudopotential data file FILE, from'
   write(6,'(a)') '      which oncvpsp_ld.x computes the log derivatives at'
   write(6,'(a)') '      any radii and energies'
   write(6,'(a)') '  --no-log-derivatives'
   write(6,'(a)') '      Leave out the log derivatives (their plot data and'
   write(6,'(a)') '      the DIAGNOSTICS ld lines); oncvpsp_ld.x gives them'
   write(6,'(a)') '      from --ps-file'
 end if
 write(6,'(a)') ''
 write(6,'(a)') 'Errors are reported on stdout, with a nonzero exit status.'
 end subroutine write_help

 end module main_setup_m
