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
! sections of the input echo whose form depends on the input
! optional inputs are written only when they are set, so inputs without
! them echo as they always have. The *_lines routines return the lines,
! for the PSML echo; write_icmod and write_epsh write them to unit 6.
! write_input_echo writes the whole echo of the psp8 and UPF files.

 module input_echo_m

 use constants_m, only: dp
 implicit none
 private
 public :: icmod_lines, epsh_lines, write_icmod, write_epsh
 public :: write_input_echo

!> length of the lines returned by the *_lines routines
 integer, parameter, public :: echo_len=100
!> maximum number of lines returned by the *_lines routines
 integer, parameter, public :: echo_max=5

 contains

 subroutine icmod_lines(icmod,fcfact,rcfact,slack,fcgrid,rcgrid,text,ntext)

!! MODEL CORE CHARGE section
!! the icmod=4 search grid is written only when it differs from the default
!!
!! @param icmod  model core charge option
!! @param fcfact  prefactor for model amplitude
!! @param rcfact  prefactor for model scale
!! @param slack  share of the cutoff the model core must fit (icmod=6, 7)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param text  the lines of the section
!! @param ntext  number of lines

!Input variables
 integer, intent(in) :: icmod
 real(dp), intent(in) :: fcfact,rcfact
 real(dp), intent(in) :: slack
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Output variables
 character(len=echo_len), intent(out) :: text(echo_max)
 integer, intent(out) :: ntext

!Local variables
!defaults, as set by the input readers
 real(dp), parameter :: fcgrid_def(3)=[1.5d0,6.0d0,0.5d0]
 real(dp), parameter :: rcgrid_def(3)=[1.0d0,1.9d0,0.1d0]

 text(1)='#'
 text(2)='# MODEL CORE CHARGE'
 if(icmod==4 .and. (any(fcgrid/=fcgrid_def) &
&   .or. any(rcgrid/=rcgrid_def))) then
   text(3)='# icmod, fcfact, rcfact,'
   text(4)='#   fcfact_min, fcfact_max, fcfact_step, '// &
&   'rcfact_min, rcfact_max, rcfact_step'
   write(text(5),'(i5,8f10.5)') icmod,fcfact,rcfact,fcgrid(:),rcgrid(:)
   ntext=5
 else if((icmod==6 .or. icmod==7) .and. slack/=1.0d0) then
!icmod=6, 7 rcfact is an error tolerance, slack is written if not 1
   text(3)='# icmod, fcfact, rcfact, slack'
   write(text(4),'(i5,f10.5,es12.4,f10.5)') icmod,fcfact,rcfact,slack
   ntext=4
 else
   text(3)='# icmod, fcfact, rcfact'
   if(icmod==6 .or. icmod==7) then
     write(text(4),'(i5,f10.5,es12.4)') icmod,fcfact,rcfact
   else
     write(text(4),'(i5,2f10.5)') icmod,fcfact,rcfact
   end if
   ntext=4
 end if

 return
 end subroutine icmod_lines

 subroutine epsh_lines(epsh1,epsh2,depsh,rxpsh,text,ntext)

!! LOG DERIVATIVE ANALYSIS section
!! rxpsh is written only when the input sets it (the readers set -1 otherwise)
!!
!! @param epsh1  lower limit of the log derivative analysis energy range (Ha)
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha)
!! @param depsh  energy step (Ha)
!! @param rxpsh  radius of the analysis, negative for the default
!! @param text  the lines of the section
!! @param ntext  number of lines

!Input variables
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh

!Output variables
 character(len=echo_len), intent(out) :: text(echo_max)
 integer, intent(out) :: ntext

 text(1)='#'
 text(2)='# LOG DERIVATIVE ANALYSIS'
 if(rxpsh>=0.0d0) then
   text(3)='# epsh1, epsh2, depsh, rxpsh'
   write(text(4),'(3f8.2,f10.5)') epsh1,epsh2,depsh,rxpsh
 else
   text(3)='# epsh1, epsh2, depsh'
   write(text(4),'(3f8.2)') epsh1,epsh2,depsh
 end if
 ntext=4

 return
 end subroutine epsh_lines

 subroutine write_icmod(icmod,fcfact,rcfact,slack,fcgrid,rcgrid)

!! writes the MODEL CORE CHARGE section to unit 6
!!
!! @param icmod  model core charge option
!! @param fcfact  prefactor for model amplitude
!! @param rcfact  prefactor for model scale
!! @param slack  share of the cutoff the model core must fit (icmod=6, 7)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)

!Input variables
 integer, intent(in) :: icmod
 real(dp), intent(in) :: fcfact,rcfact
 real(dp), intent(in) :: slack
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Local variables
 character(len=echo_len) :: text(echo_max)
 integer :: ii,ntext

 call icmod_lines(icmod,fcfact,rcfact,slack,fcgrid,rcgrid,text,ntext)
 write(6,'(a)') (trim(text(ii)),ii=1,ntext)

 return
 end subroutine write_icmod

 subroutine write_epsh(epsh1,epsh2,depsh,rxpsh)

!! writes the LOG DERIVATIVE ANALYSIS section to unit 6
!!
!! @param epsh1  lower limit of the log derivative analysis energy range (Ha)
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha)
!! @param depsh  energy step (Ha)
!! @param rxpsh  radius of the analysis, negative for the default

!Input variables
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh

!Local variables
 character(len=echo_len) :: text(echo_max)
 integer :: ii,ntext

 call epsh_lines(epsh1,epsh2,depsh,rxpsh,text,ntext)
 write(6,'(a)') (trim(text(ii)),ii=1,ntext)

 return
 end subroutine write_epsh

 subroutine write_input_echo(nkap,iformat,atsym,zz,nc,nv,iexc,psfile, &
&                            na,la,fa,lmax,rc0,ep,ncon,nbas,qcut, &
&                            lloc,lpopt,dvloc0,nproj,debl, &
&                            icmod,fcfact,rcfact,slack,fcgrid,rcgrid, &
&                            epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                            ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh)
!! writes the echo of the input data to unit 6, from the reference
!! configuration to the log mesh, as the psp8 and UPF files carry it
!!
!! @param nkap  number of kappa per l: 1, or 2 for the fully-relativistic code
!! @param iformat  1 for the psp8 file, 2 for the UPF file
!!
!! (the other arguments are the input data, as read by the input readers;
!! rc0 the input core radii, ep the energies for the first kappa and fcfact
!! the input value)

!Input variables
 integer, intent(in) :: nkap,iformat
 integer, intent(in) :: nc,nv,iexc,lmax,lloc,lpopt,icmod,ncnf
 integer, intent(in) :: na(30),la(30),ncon(6),nbas(6),nproj(6)
 integer, intent(in) :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp), intent(in) :: zz,dvloc0,fcfact,rcfact
 real(dp), intent(in) :: slack
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh,rlmax,drl
 real(dp), intent(in) :: fa(30),rc0(6),ep(6),qcut(6),debl(6),facnf(30,5)
 real(dp), intent(in) :: fcgrid(3),rcgrid(3),log_mesh(3)
 character(len=2), intent(in) :: atsym
 character(len=4), intent(in) :: psfile

!Local variables
 integer :: ii,jj,l1

 write(6,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
 write(6,'(a)') '# atsym  z   nc   nv     iexc    psfile'
 write(6,'(a,a,f6.2,2i5,i8,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&      '      ',trim(psfile)
!the energy column is named, though not written, except in the sr psp8 file
 if(iformat==1 .and. nkap==1) then
   write(6,'(a/a)') '#','#   n    l    f'
 else
   write(6,'(a/a)') '#','#   n    l    f        energy (Ha)'
 end if
 do ii=1,nc+nv
   write(6,'(2i5,f8.2)') na(ii),la(ii),fa(ii)
 end do

 write(6,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
 write(6,'(i5)')  lmax
 if(nkap==2) then
   write(6,'(a/a)') '#','#   l,   rc,     ep,   ncon, nbas, qcut'
 else
   write(6,'(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
 end if
 do l1=1,lmax+1
   write(6,'(i5,2f10.5,2i5,f10.5)') l1-1,rc0(l1),ep(l1),ncon(l1),&
&        nbas(l1),qcut(l1)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,f10.5,a,f10.5)') lloc,lpopt,rc0(5),'   ',dvloc0

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      '# l, nproj, debl'
 do l1=1,lmax+1
   write(6,'(2i5,f10.5)') l1-1,nproj(l1),debl(l1)
 end do

 call write_icmod(icmod,fcfact,rcfact,slack,fcgrid,rcgrid)

 call write_epsh(epsh1,epsh2,depsh,rxpsh)

 write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
 write(6,'(2f8.2)') rlmax,drl

 write(6,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(6,'(i5)') ncnf
 write(6,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,ncnf+1
   write(6,'(i5)') nvcnf(jj)
   do ii=nc+1,nc+nvcnf(jj)
     write(6,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
   end do
   write(6,'(a)') '#'
 end do

 if(log_mesh(1)>0.0d0) then
   write(6,'(a/a/a)') '#','# LOG MESH','# amesh, rr1, rmax'
   write(6,'(3es24.16)') log_mesh
 end if
 return
 end subroutine write_input_echo

 end module input_echo_m
