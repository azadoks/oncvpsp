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
! them echo as they always have. The *_lines routines return the lines of a
! section; write_icmod and write_epsh write them to unit 6.
! input_echo_lines returns the whole echo of the psp8, UPF and PSML files,
! and write_input_echo writes it to unit 6.
! fx and fe shape single values for the output files' fields.

 module input_echo_m

 use constants_m, only: dp, mxstate, mxcnf, mxlin, well_onset_default
 use relativity_m, only: c_light_default, c_light
 implicit none
 private
 public :: atom_lines, lmax_lines, icmod_lines, epsh_lines, mesh_lines
 public :: write_atom, write_lmax, write_icmod, write_epsh
 public :: write_mesh
 public :: write_input_echo, input_echo_lines
 public :: fx, fe, vkb_head, vkb_line

!> length of the lines returned by the *_lines routines
 integer, parameter, public :: echo_len=100
!> maximum number of lines returned by the *_lines routines
 integer, parameter, public :: echo_max=5

 contains

 function fx(xx,ww,dd) result(field)

!! xx as an fww.dd field with a leading blank: the field of fww.dd where the
!! value leaves one, so that the echo of an ordinary input is unchanged,
!! else widened until it does (an es field past 40 characters), so that
!! an extreme value never runs into the field before it
!!
!! @param xx  value
!! @param ww  width of the field
!! @param dd  digits after the decimal point

!Input variables
 real(dp), intent(in) :: xx
 integer, intent(in) :: ww,dd

!Output variables
 character(len=:), allocatable :: field

!Local variables
 character(len=48) :: buf
 character(len=16) :: fmt
 integer :: wt

 do wt=ww,40
   write(fmt,'(a,i0,a,i0,a)') '(f',wt,'.',dd,')'
   write(buf,fmt) xx
   if(buf(1:1)==' ' .and. index(buf(1:wt),'*')==0) then
     field=buf(1:wt)
     return
   end if
 end do
 write(buf,'(es24.16)') xx
 field=' '//trim(adjustl(buf))
 return
 end function fx

 elemental function fe(xx) result(yy)

!! xx for an E field that has no exponent width (e24.16, not e24.16e3): 0
!! where |xx| < 1e-99, whose 3-digit exponent would be written in place
!! of the E ("2.1869480064996404-119"), which many readers reject
!!
!! @param xx  value

!Input variables
 real(dp), intent(in) :: xx

!Output variables
 real(dp) :: yy

 yy=xx
 if(abs(xx)<1.0d-99) yy=0.0d0
 return
 end function fe

 subroutine atom_lines(atsym,zz,nc,nv,iexc,psfile,text,ntext)

!! ATOM AND REFERENCE CONFIGURATION heading and first line, with the speed
!! of light of relativity_m only when it is not the default
!!
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional
!! @param psfile  format of the pseudopotential file
!! @param text  the lines of the section
!! @param ntext  number of lines

!Input variables
 character(len=2), intent(in) :: atsym
 character(len=4), intent(in) :: psfile
 real(dp), intent(in) :: zz
 integer, intent(in) :: nc,nv,iexc

!Output variables
 character(len=echo_len), intent(out) :: text(echo_max)
 integer, intent(out) :: ntext

 text(1)='# ATOM AND REFERENCE CONFIGURATION'
 if(c_light/=c_light_default) then
   text(2)='# atsym  z   nc   nv     iexc    psfile  c_light'
   write(text(3),'(3a,2i5,i8,3a,es24.16)') '  ',trim(atsym),fx(zz,6,2), &
&        nc,nv,iexc,'      ',psfile,' ',c_light
 else
   text(2)='# atsym  z   nc   nv     iexc    psfile'
   write(text(3),'(3a,2i5,i8,2a)') '  ',trim(atsym),fx(zz,6,2),nc,nv, &
&        iexc,'      ',trim(psfile)
 end if
 ntext=3

 return
 end subroutine atom_lines

 subroutine lmax_lines(lmax,ep_as_given,well_onset,well_tail,wreg,well_norm, &
&                      text,ntext)

!! the head of the PSEUDOPOTENTIAL AND OPTIMIZATION section, to the lmax
!! line: ep_as_given (0 or 1), well_onset (0 for the default) and
!! well_tail (0 for "well", 1 for "free", 2 for "confined") are written
!! after lmax only when one of them is not its default, wave_regularization
!! after them only when it or well_norm is not 0, and well_norm (0 for that
!! of well_tail, 1 all space, 2 inside rc) after it only when not 0
!!
!! @param lmax  maximum angular momentum
!! @param ep_as_given  whether the well states take ep as given
!! @param well_onset  onset of the well as a multiple of rc
!! @param well_tail  the tail of the well states, 0, 1 or 2
!! @param wreg  the regularization of the pseudo-wave optimization
!! @param well_norm  the normalization of the well states, 0, 1 or 2
!! @param text  the lines of the section
!! @param ntext  number of lines

!Input variables
 integer, intent(in) :: lmax,well_tail,well_norm
 logical, intent(in) :: ep_as_given
 real(dp), intent(in) :: wreg
 real(dp), intent(in) :: well_onset

!Output variables
 character(len=echo_len), intent(out) :: text(echo_max)
 integer, intent(out) :: ntext

!Local variables
 real(dp) :: onset

 text(1)='#'
 text(2)='# PSEUDOPOTENTIAL AND OPTIMIZATION'
 onset=0.0d0
 if(well_onset/=well_onset_default) onset=well_onset
 if(well_norm/=0) then
   text(3)='# lmax, ep_as_given, well_onset, well_tail, '// &
&          'wave_regularization, well_norm'
   write(text(4),'(i5,i5,es24.16,i5,es24.16,i5)') lmax, &
&        merge(1,0,ep_as_given),onset,well_tail,wreg,well_norm
 else if(wreg/=0.0d0) then
   text(3)='# lmax, ep_as_given, well_onset, well_tail, wave_regularization'
   write(text(4),'(i5,i5,es24.16,i5,es24.16)') lmax, &
&        merge(1,0,ep_as_given),onset,well_tail,wreg
 else if(ep_as_given .or. well_tail/=0 .or. &
&        well_onset/=well_onset_default) then
   text(3)='# lmax, ep_as_given, well_onset, well_tail'
   write(text(4),'(i5,i5,es24.16,i5)') lmax,merge(1,0,ep_as_given), &
&        onset,well_tail
 else
   text(3)='# lmax'
   write(text(4),'(i5)') lmax
 end if
 ntext=4

 return
 end subroutine lmax_lines

 subroutine icmod_lines(icmod,fcfact,rcfact,slack,budget_cut,fcgrid, &
&                       rcgrid,text,ntext)

!! MODEL CORE CHARGE section
!! the icmod=4 search grid is written only when it differs from the default
!!
!! @param icmod  model core charge option
!! @param fcfact  prefactor for model amplitude
!! @param rcfact  prefactor for model scale
!! @param slack  share of the cutoff the model core must fit (icmod=6, 7)
!! @param budget_cut  the cutoff of the pseudo-waves (icmod=6, 7): 0 the
!!                    largest of the ground-state waves, 1 their weighted
!!                    mean, 2 the largest of all waves (budget_cutoff)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param text  the lines of the section
!! @param ntext  number of lines

!Input variables
 integer, intent(in) :: icmod
 real(dp), intent(in) :: fcfact,rcfact
 real(dp), intent(in) :: slack
 integer, intent(in) :: budget_cut
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Output variables
 character(len=echo_len), intent(out) :: text(echo_max)
 integer, intent(out) :: ntext

!Local variables
!defaults, as set by the input readers
 real(dp), parameter :: fcgrid_def(3)=[1.5d0,6.0d0,0.5d0]
 real(dp), parameter :: rcgrid_def(3)=[1.0d0,1.9d0,0.1d0]
 integer :: ii

 text(1)='#'
 text(2)='# MODEL CORE CHARGE'
 if(icmod==4 .and. (any(fcgrid/=fcgrid_def) &
&   .or. any(rcgrid/=rcgrid_def))) then
   text(3)='# icmod, fcfact, rcfact,'
   text(4)='#   fcfact_min, fcfact_max, fcfact_step, '// &
&   'rcfact_min, rcfact_max, rcfact_step'
   write(text(5),'(i5,8a)') icmod,fx(fcfact,10,5),fx(rcfact,10,5), &
&        (fx(fcgrid(ii),10,5),ii=1,3),(fx(rcgrid(ii),10,5),ii=1,3)
   ntext=5
 else if((icmod==6 .or. icmod==7) .and. budget_cut/=0) then
!icmod=6, 7 budget_cutoff is written if not "max" (0), with slack before it
   text(3)='# icmod, fcfact, rcfact, slack, budget_cutoff'
   write(text(4),'(i5,a,es12.4,a,i5)') icmod,fx(fcfact,10,5),rcfact, &
&        fx(slack,10,5),budget_cut
   ntext=4
 else if((icmod==6 .or. icmod==7) .and. slack/=1.0d0) then
!icmod=6, 7 rcfact is an error tolerance, slack is written if not 1
   text(3)='# icmod, fcfact, rcfact, slack'
   write(text(4),'(i5,a,es12.4,a)') icmod,fx(fcfact,10,5),rcfact, &
&        fx(slack,10,5)
   ntext=4
 else
   text(3)='# icmod, fcfact, rcfact'
   if(icmod==6 .or. icmod==7) then
     write(text(4),'(i5,a,es12.4)') icmod,fx(fcfact,10,5),rcfact
   else
     write(text(4),'(i5,2a)') icmod,fx(fcfact,10,5),fx(rcfact,10,5)
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
   write(text(4),'(4a)') fx(epsh1,8,2),fx(epsh2,8,2),fx(depsh,8,2), &
&        fx(rxpsh,10,5)
 else
   text(3)='# epsh1, epsh2, depsh'
   write(text(4),'(3a)') fx(epsh1,8,2),fx(epsh2,8,2),fx(depsh,8,2)
 end if
 ntext=4

 return
 end subroutine epsh_lines

 subroutine mesh_lines(log_mesh,text,ntext)

!! LOG MESH section, the optional last line of the input, written only
!! when the input sets it
!!
!! @param log_mesh  amesh, rr1, rmax of the log mesh, 0 unless the input
!!                  sets it
!! @param text  the lines of the sections
!! @param ntext  number of lines

!Input variables
 real(dp), intent(in) :: log_mesh(3)

!Output variables
 character(len=echo_len), intent(out) :: text(echo_max)
 integer, intent(out) :: ntext

 ntext=0
 if(log_mesh(1)>0.0d0) then
   text(1)='#'
   text(2)='# LOG MESH'
   text(3)='# amesh, rr1, rmax'
   write(text(4),'(3es24.16)') log_mesh
   ntext=4
 end if

 return
 end subroutine mesh_lines

 subroutine write_atom(atsym,zz,nc,nv,iexc,psfile)

!! writes the ATOM AND REFERENCE CONFIGURATION heading and first line to
!! unit 6 (see atom_lines)
!!
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional
!! @param psfile  format of the pseudopotential file

!Input variables
 character(len=2), intent(in) :: atsym
 character(len=4), intent(in) :: psfile
 real(dp), intent(in) :: zz
 integer, intent(in) :: nc,nv,iexc

!Local variables
 character(len=echo_len) :: text(echo_max)
 integer :: ii,ntext

 call atom_lines(atsym,zz,nc,nv,iexc,psfile,text,ntext)
 write(6,'(a)') (trim(text(ii)),ii=1,ntext)

 return
 end subroutine write_atom

 subroutine write_lmax(lmax,ep_as_given,well_onset,well_tail,wreg,well_norm)

!! writes the head of the PSEUDOPOTENTIAL AND OPTIMIZATION section, to the
!! lmax line, to unit 6
!!
!! @param lmax  maximum angular momentum
!! @param ep_as_given  whether the well states take ep as given
!! @param well_onset  onset of the well as a multiple of rc
!! @param well_tail  the tail of the well states, 0, 1 or 2
!! @param wreg  the regularization of the pseudo-wave optimization
!! @param well_norm  the normalization of the well states, 0, 1 or 2

!Input variables
 integer, intent(in) :: lmax,well_tail,well_norm
 logical, intent(in) :: ep_as_given
 real(dp), intent(in) :: wreg
 real(dp), intent(in) :: well_onset

!Local variables
 character(len=echo_len) :: text(echo_max)
 integer :: ii,ntext

 call lmax_lines(lmax,ep_as_given,well_onset,well_tail,wreg,well_norm,text, &
&                ntext)
 write(6,'(a)') (trim(text(ii)),ii=1,ntext)

 return
 end subroutine write_lmax

 character(len=echo_len) function vkb_head(lmax,well_n,ecut_wt)
!! the column line of the VANDERBILT-KLEINMAN-BYLANDER PROJECTORs section:
!! well_n is written only when the input sets it for some l, or sets
!! ecut_weight, which follows it
!!
!! @param lmax  maximum angular momentum
!! @param well_n  n of the first well state of each l (0 for the default)
!! @param ecut_wt  cutoff weight of each l (negative for the default)

!Input variables
 integer, intent(in) :: lmax,well_n(mxlin+1)
 real(dp), intent(in) :: ecut_wt(mxlin+1)

 if(any(ecut_wt(1:lmax+1)>=0.0d0)) then
   vkb_head='# l, nproj, debl, well_n, ecut_weight'
 else if(any(well_n(1:lmax+1)>0)) then
   vkb_head='# l, nproj, debl, well_n'
 else
   vkb_head='# l, nproj, debl'
 end if
 end function vkb_head

 character(len=echo_len) function vkb_line(lmax,l1,nproj,debl,well_n, &
&                                           ecut_wt)
!! the line of l = l1-1 of that section
!!
!! @param lmax  maximum angular momentum
!! @param l1  l+1
!! @param nproj  number of projectors of each l
!! @param debl  energy steps of each l
!! @param well_n  n of the first well state of each l (0 for the default)
!! @param ecut_wt  cutoff weight of each l (negative for the default)

!Input variables
 integer, intent(in) :: lmax,l1,nproj(mxlin+1),well_n(mxlin+1)
 real(dp), intent(in) :: debl(mxlin+1),ecut_wt(mxlin+1)

 if(any(ecut_wt(1:lmax+1)>=0.0d0)) then
   write(vkb_line,'(2i5,a,i5,a)') l1-1,nproj(l1),fx(debl(l1),10,5), &
&        well_n(l1),fx(ecut_wt(l1),10,5)
 else if(any(well_n(1:lmax+1)>0)) then
   write(vkb_line,'(2i5,a,i5)') l1-1,nproj(l1),fx(debl(l1),10,5),well_n(l1)
 else
   write(vkb_line,'(2i5,a)') l1-1,nproj(l1),fx(debl(l1),10,5)
 end if
 end function vkb_line

 subroutine write_icmod(icmod,fcfact,rcfact,slack,budget_cut,fcgrid,rcgrid)

!! writes the MODEL CORE CHARGE section to unit 6
!!
!! @param icmod  model core charge option
!! @param fcfact  prefactor for model amplitude
!! @param rcfact  prefactor for model scale
!! @param slack  share of the cutoff the model core must fit (icmod=6, 7)
!! @param budget_cut  the cutoff of the pseudo-waves (icmod=6, 7): 0 the
!!                    largest of the ground-state waves, 1 their weighted
!!                    mean, 2 the largest of all waves (budget_cutoff)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)

!Input variables
 integer, intent(in) :: icmod
 real(dp), intent(in) :: fcfact,rcfact
 real(dp), intent(in) :: slack
 integer, intent(in) :: budget_cut
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Local variables
 character(len=echo_len) :: text(echo_max)
 integer :: ii,ntext

 call icmod_lines(icmod,fcfact,rcfact,slack,budget_cut,fcgrid,rcgrid, &
&                 text,ntext)
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

 subroutine write_mesh(log_mesh)

!! writes the LOG MESH section to unit 6 (see mesh_lines)
!!
!! @param log_mesh  amesh, rr1, rmax of the log mesh, 0 unless the input
!!                  sets it

!Input variables
 real(dp), intent(in) :: log_mesh(3)

!Local variables
 character(len=echo_len) :: text(echo_max)
 integer :: ii,ntext

 call mesh_lines(log_mesh,text,ntext)
 if(ntext>0) write(6,'(a)') (trim(text(ii)),ii=1,ntext)

 return
 end subroutine write_mesh

 subroutine input_echo_lines(nkap,iformat,atsym,zz,nc,nv,iexc,psfile, &
&                            na,la,fa,lmax,rc0,ep,ncon,nbas,qcut,ep_as_given, &
&                            well_onset,well_tail,wreg,lloc,lpopt,dvloc0, &
&                            nproj,debl,well_n,well_norm,ecut_wt, &
&                            icmod,fcfact,rcfact,slack,budget_cut, &
&                            fcgrid,rcgrid, &
&                            epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                            ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh,upf_mesh, &
&                           text,ntext)
!! the echo of the input data, from the reference configuration to the log
!! mesh, as the psp8, UPF and PSML files carry it, line by line
!!
!! @param nkap  number of kappa per l: 1, or 2 for the fully-relativistic code
!! @param iformat  1 for the psp8 file (and PSML, with nkap 1), 2 for the
!!                 UPF file
!! @param text  the lines (output)
!! @param ntext  their number (output)
!!
!! (the other arguments are the input data, as read by the input readers;
!! rc0 the input core radii, ep the energies for the first kappa and fcfact
!! the input value)

!Input variables
 integer, intent(in) :: nkap,iformat
 integer, intent(in) :: nc,nv,iexc,lmax,lloc,lpopt,icmod,ncnf,upf_mesh
 integer, intent(in) :: na(mxstate),la(mxstate),ncon(mxlin+1),nbas(mxlin+1)
 integer, intent(in) :: nproj(mxlin+1)
 integer, intent(in) :: nvcnf(mxcnf),nacnf(mxstate,mxcnf),lacnf(mxstate,mxcnf)
 real(dp), intent(in) :: zz,dvloc0,fcfact,rcfact
 real(dp), intent(in) :: slack
 integer, intent(in) :: budget_cut
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh,rlmax,drl
 real(dp), intent(in) :: fa(mxstate),rc0(mxlin+1),ep(mxlin+1),qcut(mxlin+1)
 logical, intent(in) :: ep_as_given
 real(dp), intent(in) :: well_onset
 integer, intent(in) :: well_tail
 real(dp), intent(in) :: wreg
 real(dp), intent(in) :: debl(mxlin+1)
 integer, intent(in) :: well_n(mxlin+1)
 integer, intent(in) :: well_norm
 real(dp), intent(in) :: ecut_wt(mxlin+1)
 real(dp), intent(in) :: facnf(mxstate,mxcnf)
 real(dp), intent(in) :: fcgrid(3),rcgrid(3),log_mesh(3)
 character(len=2), intent(in) :: atsym
 character(len=4), intent(in) :: psfile

!Output variables
 character(len=echo_len), allocatable, intent(out) :: text(:)
 integer, intent(out) :: ntext

!Local variables
 character(len=echo_len) :: part(echo_max),line
 integer :: ii,jj,l1,npart

 allocate(text(64))
 ntext=0

 call atom_lines(atsym,zz,nc,nv,iexc,psfile,part,npart)
 call put_part()
!the energy column is named, though not written, except in the sr psp8 file
 call put('#')
 if(iformat==1 .and. nkap==1) then
   call put('#   n    l    f')
 else
   call put('#   n    l    f        energy (Ha)')
 end if
 do ii=1,nc+nv
   write(line,'(2i5,f8.2)') na(ii),la(ii),fa(ii)
   call put(line)
 end do

 call lmax_lines(lmax,ep_as_given,well_onset,well_tail,wreg,well_norm,part, &
&                npart)
 call put_part()
 call put('#')
 if(nkap==2) then
   call put('#   l,   rc,     ep,   ncon, nbas, qcut')
 else
   call put('#   l,   rc,      ep,       ncon, nbas, qcut')
 end if
 do l1=1,lmax+1
   write(line,'(i5,2a,2i5,a)') l1-1,fx(rc0(l1),10,5),fx(ep(l1),10,5), &
&        ncon(l1),nbas(l1),fx(qcut(l1),10,5)
   call put(line)
 end do

 call put('#')
 call put('# LOCAL POTENTIAL')
 call put('# lloc, lpopt,  rc(5),   dvloc0')
 write(line,'(2i5,3a)') lloc,lpopt,fx(rc0(5),10,5),'   ',fx(dvloc0,10,5)
 call put(line)

 call put('#')
 call put('# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs')
 call put(vkb_head(lmax,well_n,ecut_wt))
 do l1=1,lmax+1
   call put(vkb_line(lmax,l1,nproj,debl,well_n,ecut_wt))
 end do

 call icmod_lines(icmod,fcfact,rcfact,slack,budget_cut,fcgrid,rcgrid, &
&                 part,npart)
 call put_part()

 call epsh_lines(epsh1,epsh2,depsh,rxpsh,part,npart)
 call put_part()

 call put('#')
 call put('# OUTPUT GRID')
 if(upf_mesh==0) then
   call put('# rlmax, drl')
   write(line,'(2a)') fx(rlmax,8,2),fx(drl,8,2)
 else
   call put('# rlmax, drl, upf_mesh')
   write(line,'(2a,i8)') fx(rlmax,8,2),fx(drl,8,2),upf_mesh
 end if
 call put(line)

 call put('#')
 call put('# TEST CONFIGURATIONS')
 call put('# ncnf')
 write(line,'(i5)') ncnf
 call put(line)
 call put('# nvcnf')
 call put('#   n    l    f')
 do jj=2,ncnf+1
   write(line,'(i5)') nvcnf(jj)
   call put(line)
   do ii=nc+1,nc+nvcnf(jj)
     write(line,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
     call put(line)
   end do
   call put('#')
 end do

 call mesh_lines(log_mesh,part,npart)
 call put_part()
 return

 contains

 subroutine put(str)
!! appends a line
 character(len=*), intent(in) :: str
 character(len=echo_len), allocatable :: tmp(:)
 if(ntext==size(text)) then
   allocate(tmp(2*size(text)))
   tmp(1:ntext)=text(1:ntext)
   call move_alloc(tmp,text)
 end if
 ntext=ntext+1
 text(ntext)=str
 end subroutine put

 subroutine put_part()
!! appends the lines of part
 integer :: kk
 do kk=1,npart
   call put(part(kk))
 end do
 end subroutine put_part

 end subroutine input_echo_lines

 subroutine write_input_echo(nkap,iformat,atsym,zz,nc,nv,iexc,psfile, &
&                            na,la,fa,lmax,rc0,ep,ncon,nbas,qcut,ep_as_given, &
&                            well_onset,well_tail,wreg,lloc,lpopt,dvloc0, &
&                            nproj,debl,well_n,well_norm,ecut_wt, &
&                            icmod,fcfact,rcfact,slack,budget_cut, &
&                            fcgrid,rcgrid, &
&                            epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                            ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh,upf_mesh)
!! writes the echo of the input data to unit 6 (input_echo_lines)
!!
!! @param nkap  number of kappa per l: 1, or 2 for the fully-relativistic code
!! @param iformat  1 for the psp8 file, 2 for the UPF file
!!
!! (the other arguments are the input data, as read by the input readers)

!Input variables
 integer, intent(in) :: nkap,iformat
 integer, intent(in) :: nc,nv,iexc,lmax,lloc,lpopt,icmod,ncnf,upf_mesh
 integer, intent(in) :: na(mxstate),la(mxstate),ncon(mxlin+1),nbas(mxlin+1)
 integer, intent(in) :: nproj(mxlin+1)
 integer, intent(in) :: nvcnf(mxcnf),nacnf(mxstate,mxcnf),lacnf(mxstate,mxcnf)
 real(dp), intent(in) :: zz,dvloc0,fcfact,rcfact
 real(dp), intent(in) :: slack
 integer, intent(in) :: budget_cut
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh,rlmax,drl
 real(dp), intent(in) :: fa(mxstate),rc0(mxlin+1),ep(mxlin+1),qcut(mxlin+1)
 logical, intent(in) :: ep_as_given
 real(dp), intent(in) :: well_onset
 integer, intent(in) :: well_tail
 real(dp), intent(in) :: wreg
 real(dp), intent(in) :: debl(mxlin+1)
 integer, intent(in) :: well_n(mxlin+1)
 integer, intent(in) :: well_norm
 real(dp), intent(in) :: ecut_wt(mxlin+1)
 real(dp), intent(in) :: facnf(mxstate,mxcnf)
 real(dp), intent(in) :: fcgrid(3),rcgrid(3),log_mesh(3)
 character(len=2), intent(in) :: atsym
 character(len=4), intent(in) :: psfile

!Local variables
 character(len=echo_len), allocatable :: text(:)
 integer :: ii,ntext

 call input_echo_lines(nkap,iformat,atsym,zz,nc,nv,iexc,psfile, &
&                      na,la,fa,lmax,rc0,ep,ncon,nbas,qcut,ep_as_given, &
&                      well_onset,well_tail,wreg,lloc,lpopt,dvloc0, &
&                      nproj,debl,well_n,well_norm,ecut_wt, &
&                      icmod,fcfact,rcfact,slack,budget_cut,fcgrid,rcgrid, &
&                      epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                      ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh,upf_mesh, &
&                      text,ntext)
 do ii=1,ntext
   write(6,'(a)') trim(text(ii))
 end do
 return
 end subroutine write_input_echo

 end module input_echo_m
