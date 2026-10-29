!
! Copyright (c) 2025-2026 by Austin Zadoks
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
 module input_toml_m
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
! TOML Fortran library: documentation at toml-f.readthedocs.io
 implicit none
 private
 public :: read_input_toml
 contains

 subroutine read_input_toml(unit, inp, iout, ifail)

!! Reads the TOML input: the same data as read_input_text, with the same
!! defaults and the same checks against the array sizes (check_data
!! validates the values afterwards).  Every table and key is required unless
!! stated otherwise; a missing or mistyped one is a failure (ifail=1), with
!! a message naming it.
!!
!! scripts/input_text_to_toml.py converts a text input to this form.  The
!! tables and their keys are listed below.
!!
!! @param unit  the unit of the open TOML file
!! @param inp  the input values (input_m)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: mxstate, mxcnf, mxlin, well_onset_default
 use input_m, only: input_t
 use tomlf, only: toml_table, toml_array, toml_error, toml_stat, toml_load, &
&                  get_value, len

!   [oncvpsp]                  atsym, z, nc, nv, iexc, c_light (optional, 0
!                              or less for the default)
!   [reference_configuration]  n, l, f                (nc + nv each)
!   [pseudopotentials]         lmax, and l (0..lmax), rc, ep, ncon, nbas,
!                              qcut, ep_as_given (optional, false: the well
!                              states take ep as given, see wellstate),
!                              well_onset (optional, 1.006**5, at least 1:
!                              the onset of the well of the well states as
!                              a multiple of rc, see wellstate),
!                              well_tail (optional, "well", "free" or
!                              "confined": the tail of the well states
!                              in E_resid, see free_tail_m, and their
!                              well, see wellstate), wave_regularization
!                              (optional, 0, at least 0: its weight on
!                              E_resid at qcut/4 in the optimization, see
!                              optimize), well_norm (optional, "space" or
!                              "rc": the normalization of the well
!                              states; by default that of well_tail)
!   [local_potential]          lloc, lpopt, rcloc, dvloc0
!   [vkb_projectors]           l (0..lmax), nproj, debl, well_n
!                              (optional, with well_tail = "confined"
!                              only: the n of each l's first well state,
!                              0 for the default), ecut_weight (optional:
!                              the weight of the first wave of each l with
!                              no occupied state in the pseudo-waves'
!                              cutoff, negative for the default 2(2l+1))
!   [model_core_charge]        icmod, fcfact (icmod 1, 2, 3, 6, 7; else 0),
!                              rcfact (icmod 3; 1e-5 for icmod 6, 7; else
!                              0), slack (icmod 6, 7; 1), budget_cutoff
!                              (icmod 6, 7; "max", "weighted" or
!                              "max_all", "max"),
!                              fcfact_min, fcfact_max, fcfact_step, rcfact_min,
!                              rcfact_max, rcfact_step (optional, icmod 4 grid)
!   [log_derivative_analysis]  epsh1, epsh2, depsh, rxpsh (optional, -1 for
!                              auto)
!   [linear_mesh]              a (drl), rmax (rlmax)
!   [log_mesh]                 amesh, rr1, rmax (optional table and keys, each
!                              0 or less or missing for the built-in mesh),
!                              adams_order (optional, 4)
!   [pp_output]                psfile, upf_mesh (optional, 0)
!   [[test_configurations]]    n, l, f (valence states only; optional, up
!                              to 4)

!Input variables
!> The Fortran unit number to read from
 integer, intent(in) :: unit
!> unit of the text output
 integer, intent(in) :: iout

!Output variables
!> The input values
 type(input_t), intent(out) :: inp
!> 0, or 1 after a failure, whose message is written
 integer, intent(out) :: ifail

!Local variables
!> The TOML document
 type(toml_table), allocatable :: table
!> Parser error
 type(toml_error), allocatable :: error
!> A section of the document
 type(toml_table), pointer :: child
!> The array of test configurations
 type(toml_array), pointer :: arr
!> The n array of a test configuration
 type(toml_array), pointer :: narr
!> Status of a lookup
 integer :: stat
!> Configuration index
 integer :: jj
!> Angular momenta of a section's rows
 integer :: lt(mxlin+1)
!> Temporary string
 character(len=:), allocatable :: tmp_str
!> Message of an input error that names a limit
 character(len=80) :: msg

 ifail = 0
 associate(atsym=>inp%atsym, zz=>inp%zz, nc=>inp%nc, nv=>inp%nv, &
&          iexc=>inp%iexc, psfile=>inp%psfile, na=>inp%na, la=>inp%la, &
&          fa=>inp%fa, lmax=>inp%lmax, rc=>inp%rc, ep=>inp%ep, ncon=>inp%ncon, &
&          nbas=>inp%nbas, qcut=>inp%qcut, lloc=>inp%lloc, lpopt=>inp%lpopt, &
&          dvloc0=>inp%dvloc0, nproj=>inp%nproj, debl=>inp%debl, &
&          icmod=>inp%icmod, fcfact=>inp%fcfact, rcfact=>inp%rcfact, &
&          slack=>inp%slack, budget_cutoff=>inp%budget_cutoff, &
&          fcfact_min=>inp%fcfact_min, fcfact_max=>inp%fcfact_max, &
&          fcfact_step=>inp%fcfact_step, rcfact_min=>inp%rcfact_min, &
&          rcfact_max=>inp%rcfact_max, rcfact_step=>inp%rcfact_step, &
&          epsh1=>inp%epsh1, epsh2=>inp%epsh2, depsh=>inp%depsh, &
&          rxpsh=>inp%rxpsh, rlmax=>inp%rlmax, drl=>inp%drl, &
&          upf_mesh=>inp%upf_mesh, ncnf=>inp%ncnf, cnf=>inp%cnf, &
&          amesh=>inp%amesh, rr1=>inp%rr1, rmax=>inp%rmax, &
&          adams_order=>inp%adams_order, &
&          c_light=>inp%c_light, ep_as_given=>inp%ep_as_given, &
&          well_onset=>inp%well_onset, well_tail=>inp%well_tail, &
&          wave_regularization=>inp%wave_regularization, well_n=>inp%well_n, &
&          well_norm=>inp%well_norm, &
&          ecut_weight=>inp%ecut_weight)

! the values not read keep the defaults of input_t
   call toml_load(table, unit, error=error)
   if(allocated(error)) then
     call input_error('cannot parse the input: '//error%message, iout, ifail)
     if(ifail /= 0) return
   end if
   if(.not. allocated(table)) then
     call input_error('cannot parse the input', iout, ifail)
     if(ifail /= 0) return
   end if

! atom
   call get_section(table, 'oncvpsp', child, iout, ifail)
   if(ifail /= 0) return
   call get_string(child, 'oncvpsp', 'atsym', tmp_str, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   if(len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 2) then
     call input_error('atsym in [oncvpsp] must have 1 or 2 characters, '// &
&      'got "'//tmp_str//'"', iout, ifail)
     if(ifail /= 0) return
   end if
   atsym = trim(tmp_str)
   call get_real(child, 'oncvpsp', 'z', zz, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_int(child, 'oncvpsp', 'nc', nc, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_int(child, 'oncvpsp', 'nv', nv, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_real(child, 'oncvpsp', 'c_light', c_light, -1.0d0, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_int(child, 'oncvpsp', 'iexc', iexc, iout=iout, ifail=ifail)
   if(ifail /= 0) return
! counts beyond the array sizes would overrun before check_data runs
   if(nc < 0 .or. nv < 0 .or. nc + nv > mxstate) then
     write(msg, '(a,i0)') 'nc and nv in [oncvpsp] must be >= 0 with '// &
&      'nc + nv <= ', mxstate
     call input_error(trim(msg), iout, ifail)
     if(ifail /= 0) return
   end if

! reference configuration
   call get_section(table, 'reference_configuration', child, iout, ifail)
   if(ifail /= 0) return
   call get_int_array(child, 'reference_configuration', 'n', nc + nv, na, &
&        iout, ifail)
   if(ifail /= 0) return
   call get_int_array(child, 'reference_configuration', 'l', nc + nv, la, &
&        iout, ifail)
   if(ifail /= 0) return
   call get_real_array(child, 'reference_configuration', 'f', nc + nv, fa, &
&        iout, ifail)
   if(ifail /= 0) return

! pseudopotential and optimization
   call get_section(table, 'pseudopotentials', child, iout, ifail)
   if(ifail /= 0) return
   call get_int(child, 'pseudopotentials', 'lmax', lmax, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   if(lmax < 0 .or. lmax > mxlin) then
     write(msg, '(a,i0)') 'lmax in [pseudopotentials] must be in 0..', mxlin
     call input_error(trim(msg), iout, ifail)
     if(ifail /= 0) return
   end if
   call get_int_array(child, 'pseudopotentials', 'l', lmax + 1, lt, iout, ifail)
   if(ifail /= 0) return
   call check_l(lt, lmax, 'pseudopotentials', iout, ifail)
   if(ifail /= 0) return
   call get_real_array(child, 'pseudopotentials', 'rc', lmax + 1, rc, &
&        iout, ifail)
   if(ifail /= 0) return
   call get_real_array(child, 'pseudopotentials', 'ep', lmax + 1, ep(:, 1), &
&        iout, ifail)
   if(ifail /= 0) return
   ep(:, 2) = ep(:, 1)
   call get_int_array(child, 'pseudopotentials', 'ncon', lmax + 1, ncon, &
&        iout, ifail)
   if(ifail /= 0) return
   call get_int_array(child, 'pseudopotentials', 'nbas', lmax + 1, nbas, &
&        iout, ifail)
   if(ifail /= 0) return
   call get_real_array(child, 'pseudopotentials', 'qcut', lmax + 1, qcut, &
&        iout, ifail)
   if(ifail /= 0) return
   call get_logical(child, 'pseudopotentials', 'ep_as_given', ep_as_given, &
&                    .false., iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_real(child, 'pseudopotentials', 'well_onset', well_onset, &
&                well_onset_default, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   if(.not. (well_onset >= 1.0d0)) then
     call input_error('well_onset in [pseudopotentials] must be at least 1', &
&          iout,ifail)
     if(ifail /= 0) return
   end if
   call get_string(child, 'pseudopotentials', 'well_tail', tmp_str, 'well', &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   select case(tmp_str)
   case('well')
     well_tail = 0
   case('free')
     well_tail = 1
   case('confined')
     well_tail = 2
   case default
     call input_error('well_tail in [pseudopotentials] must be "well", '// &
&      '"free" or "confined", got "'//tmp_str//'"', iout, ifail)
     if(ifail /= 0) return
   end select

   call get_real(child, 'pseudopotentials', 'wave_regularization', &
&                wave_regularization, &
&                0.0d0, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   if(.not. (wave_regularization >= 0.0d0)) then
     call input_error('wave_regularization in [pseudopotentials] must be '// &
&      'at least 0', iout, ifail)
     if(ifail /= 0) return
   end if
   call get_string(child, 'pseudopotentials', 'well_norm', tmp_str, '-', &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   select case(tmp_str)
   case('-')
     well_norm = 0
   case('space')
     well_norm = 1
   case('rc')
     well_norm = 2
   case default
     call input_error('well_norm in [pseudopotentials] must be "space" or '// &
&      '"rc", got "'//tmp_str//'"', iout, ifail)
     if(ifail /= 0) return
   end select

! local potential
   call get_section(table, 'local_potential', child, iout, ifail)
   if(ifail /= 0) return
   call get_int(child, 'local_potential', 'lloc', lloc, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_int(child, 'local_potential', 'lpopt', lpopt, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_real(child, 'local_potential', 'rcloc', rc(5), &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_real(child, 'local_potential', 'dvloc0', dvloc0, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return

! Vanderbilt-Kleinman-Bylander projectors
   call get_section(table, 'vkb_projectors', child, iout, ifail)
   if(ifail /= 0) return
   call get_int_array(child, 'vkb_projectors', 'l', lmax + 1, lt, iout, ifail)
   if(ifail /= 0) return
   call check_l(lt, lmax, 'vkb_projectors', iout, ifail)
   if(ifail /= 0) return
   call get_int_array(child, 'vkb_projectors', 'nproj', lmax + 1, nproj, &
&        iout, ifail)
   if(ifail /= 0) return
   call get_real_array(child, 'vkb_projectors', 'debl', lmax + 1, debl, &
&        iout, ifail)
   if(ifail /= 0) return
   if(child%has_key('well_n')) then
     call get_int_array(child, 'vkb_projectors', 'well_n', lmax + 1, well_n, &
&          iout, ifail)
     if(ifail /= 0) return
     if(any(well_n(1:lmax + 1) < 0)) then
       call input_error('well_n in [vkb_projectors] must be at least 0', &
&            iout, ifail)
       if(ifail /= 0) return
     end if
   end if
   if(child%has_key('ecut_weight')) then
     call get_real_array(child, 'vkb_projectors', 'ecut_weight', lmax + 1, &
&                        ecut_weight, iout, ifail)
     if(ifail /= 0) return
   end if

! model core charge; check_data validates icmod and the values
   call get_section(table, 'model_core_charge', child, iout, ifail)
   if(ifail /= 0) return
   call get_int(child, 'model_core_charge', 'icmod', icmod, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   select case(icmod)
   case(1, 2, 6, 7)
     call get_real(child, 'model_core_charge', 'fcfact', fcfact, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
   case(3)
     call get_real(child, 'model_core_charge', 'fcfact', fcfact, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'model_core_charge', 'rcfact', rcfact, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
   case default
     call get_real(child, 'model_core_charge', 'fcfact', fcfact, 0.0d0, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
   end select
   select case(icmod)
   case(4, 5)
! fcfact and rcfact are ignored; the grid of the Teter parameters
     call get_real(child, 'model_core_charge', 'rcfact', rcfact, 0.0d0, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'model_core_charge', 'fcfact_min', fcfact_min, &
&      fcfact_min, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'model_core_charge', 'fcfact_max', fcfact_max, &
&      fcfact_max, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'model_core_charge', 'fcfact_step', fcfact_step, &
&      fcfact_step, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'model_core_charge', 'rcfact_min', rcfact_min, &
&      rcfact_min, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'model_core_charge', 'rcfact_max', rcfact_max, &
&      rcfact_max, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'model_core_charge', 'rcfact_step', rcfact_step, &
&      rcfact_step, iout=iout, ifail=ifail)
     if(ifail /= 0) return
   case(6, 7)
     call get_real(child, 'model_core_charge', 'rcfact', rcfact, 1.0d-5, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'model_core_charge', 'slack', slack, 1.0d0, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_string(child, 'model_core_charge', 'budget_cutoff', tmp_str, &
&      'max', iout=iout, ifail=ifail)
     if(ifail /= 0) return
     select case(tmp_str)
     case('max')
       budget_cutoff = 0
     case('weighted')
       budget_cutoff = 1
     case('max_all')
       budget_cutoff = 2
     case default
       call input_error('budget_cutoff in [model_core_charge] must be '// &
&        '"max", "weighted" or "max_all", got "'//tmp_str//'"', iout, ifail)
       if(ifail /= 0) return
     end select
   case default
! no further keys
   end select

! log derivative analysis
   call get_section(table, 'log_derivative_analysis', child, iout, ifail)
   if(ifail /= 0) return
   call get_real(child, 'log_derivative_analysis', 'epsh1', epsh1, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_real(child, 'log_derivative_analysis', 'epsh2', epsh2, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_real(child, 'log_derivative_analysis', 'depsh', depsh, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_real(child, 'log_derivative_analysis', 'rxpsh', rxpsh, -1.0d0, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return

! output grid
   call get_section(table, 'linear_mesh', child, iout, ifail)
   if(ifail /= 0) return
   call get_real(child, 'linear_mesh', 'a', drl, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   call get_real(child, 'linear_mesh', 'rmax', rlmax, iout=iout, ifail=ifail)
   if(ifail /= 0) return

! log mesh (optional)
   call get_value(table, 'log_mesh', child, requested=.false., stat=stat)
   if(associated(child)) then
     call get_real(child, 'log_mesh', 'amesh', amesh, -1.0d0, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'log_mesh', 'rr1', rr1, -1.0d0, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_real(child, 'log_mesh', 'rmax', rmax, -1.0d0, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
     call get_int(child, 'log_mesh', 'adams_order', adams_order, 4, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) return
   else if(stat /= toml_stat%success) then
     call input_error('[log_mesh] must be a table', iout, ifail)
     if(ifail /= 0) return
   end if

! pseudopotential output file
   call get_section(table, 'pp_output', child, iout, ifail)
   if(ifail /= 0) return
   call get_string(child, 'pp_output', 'psfile', tmp_str, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return
   if(len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 4) then
     call input_error('psfile in [pp_output] must have 1 to 4 characters, '// &
&      'got "'//tmp_str//'"', iout, ifail)
     if(ifail /= 0) return
   end if
   psfile = trim(tmp_str)
   call get_int(child, 'pp_output', 'upf_mesh', upf_mesh, 0, &
&        iout=iout, ifail=ifail)
   if(ifail /= 0) return

! test configurations (optional): the valence states, after the
! reference core states the driver copies in
   call get_value(table, 'test_configurations', arr, requested=.false., &
&    stat=stat)
   if(stat /= toml_stat%success) then
     call input_error('test_configurations must be an array of tables, '// &
&      '[[test_configurations]]', iout, ifail)
     if(ifail /= 0) return
   end if
   if(associated(arr)) then
     ncnf = len(arr)
     if(ncnf > mxcnf - 1) then
       write(msg, '(a,i0,a)') 'at most ', mxcnf - 1, &
&        ' [[test_configurations]] are allowed'
       call input_error(trim(msg), iout, ifail)
       if(ifail /= 0) return
     end if
     do jj = 2, ncnf + 1
       call get_value(arr, jj - 1, child, stat=stat)
       if(.not. associated(child)) then
         call input_error('test_configurations must be an array of '// &
&                         'tables, [[test_configurations]]', iout, ifail)
         if(ifail /= 0) return
       end if
! the number of valence states is that of n
       call get_value(child, 'n', narr, requested=.false., stat=stat)
       if(.not. associated(narr)) then
         call key_error('test_configurations', 'n', stat, 'an array', &
&              iout, ifail)
         if(ifail /= 0) return
       end if
       cnf(jj)%nv = len(narr)
       if(nc + cnf(jj)%nv > mxstate) then
         write(msg, '(a,i0)') 'nc plus the states of a '// &
&          '[[test_configurations]] must be <= ', mxstate
         call input_error(trim(msg), iout, ifail)
         if(ifail /= 0) return
       end if
       call get_int_array(child, 'test_configurations', 'n', cnf(jj)%nv, &
&        cnf(jj)%na(nc + 1:), iout, ifail)
       if(ifail /= 0) return
       call get_int_array(child, 'test_configurations', 'l', cnf(jj)%nv, &
&        cnf(jj)%la(nc + 1:), iout, ifail)
       if(ifail /= 0) return
       call get_real_array(child, 'test_configurations', 'f', cnf(jj)%nv, &
&        cnf(jj)%fa(nc + 1:), iout, ifail)
       if(ifail /= 0) return
     end do
   end if
 end associate
 return
 end subroutine read_input_toml

 subroutine input_error(message, iout, ifail)

!! report an input error in the output (unit 6) and on stderr, and fail
!!
!! @param message  description of the error
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input
 use tomlf, only: len

!Input variables
 character(len=*), intent(in) :: message
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

 write(iout,'(a)') 'TOML input ERROR: '//message
 write(iout,'(a)') 'Program will stop'
 write(stderr,'(a)') 'Error: '//message
 ifail=fail_input
 return
 end subroutine input_error

 subroutine get_section(table, name, child, iout, ifail)

!! the table [name] of the document
!!
!! @param table  the TOML document
!! @param name  name of the table
!! @param child  the table [name]; a failure if it is missing or not a table
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use tomlf, only: toml_table, toml_stat, get_value, len

!Input variables
 character(len=*), intent(in) :: name
 integer, intent(in) :: iout

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 type(toml_table), pointer, intent(out) :: child
 integer, intent(out) :: ifail

!Local variables
 integer :: stat

 ifail = 0
! a missing key leaves child null with success, a mistyped one with an error
 call get_value(table, name, child, requested=.false., stat=stat)
 if(.not. associated(child)) then
   if(stat == toml_stat%success) then
     call input_error('missing table ['//name//']', iout, ifail)
     if(ifail /= 0) return
   end if
   call input_error('['//name//'] must be a table', iout, ifail)
   if(ifail /= 0) return
 end if
 return
 end subroutine get_section

 subroutine key_error(section, key, stat, what, iout, ifail)

!! report a key that a lookup did not return: missing for missing_key, or
!! for success from a lookup with requested=.false., else mistyped
!!
!! @param section  name of the table holding the key
!! @param key  name of the key
!! @param stat  status of the lookup
!! @param what  the expected type, for the message (e.g. 'an integer')
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use tomlf, only: toml_stat, len

!Input variables
 character(len=*), intent(in) :: section, key, what
 integer, intent(in) :: stat
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

 ifail = 0
 if(stat == toml_stat%missing_key .or. stat == toml_stat%success) then
   call input_error('missing key '//key//' in ['//section//']', iout, ifail)
   if(ifail /= 0) return
 else
   call input_error(key//' in ['//section//'] must be '//what, iout, ifail)
   if(ifail /= 0) return
 end if
 return
 end subroutine key_error

 subroutine get_logical(table, section, key, val, default, iout, ifail)

!! a logical key, required unless default is given
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param val  the value of the key, or default if the key is missing
!! @param default  optional, the value for a missing key, which is then not
!!                 required
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use tomlf, only: toml_table, toml_keyval, toml_stat, get_value, len

!Input variables
 character(len=*), intent(in) :: section, key
 logical, intent(in), optional :: default
 integer, intent(in) :: iout

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 logical, intent(out) :: val
 integer, intent(out) :: ifail

!Local variables
 type(toml_keyval), pointer :: ptr
 integer :: stat

 ifail = 0
 if(present(default)) then
   call get_value(table, key, ptr, requested=.false., stat=stat)
   if(.not. associated(ptr) .and. stat == toml_stat%success) then
     val = default
     return
   end if
 end if
 call get_value(table, key, val, stat=stat)
 if(stat /= toml_stat%success) then
   call key_error(section, key, stat, 'a logical', iout, ifail)
   if(ifail /= 0) return
 end if
 return
 end subroutine get_logical

 subroutine get_int(table, section, key, val, default, iout, ifail)

!! an integer key, required unless default is given
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param val  the value of the key, or default if the key is missing
!! @param default  optional, the value for a missing key, which is then not
!!                 required
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use tomlf, only: toml_table, toml_keyval, toml_stat, get_value, len

!Input variables
 character(len=*), intent(in) :: section, key
 integer, intent(in), optional :: default
 integer, intent(in) :: iout

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 integer, intent(out) :: val
 integer, intent(out) :: ifail

!Local variables
 type(toml_keyval), pointer :: ptr
 integer :: stat

 ifail = 0
 if(present(default)) then
   call get_value(table, key, ptr, requested=.false., stat=stat)
   if(.not. associated(ptr) .and. stat == toml_stat%success) then
     val = default
     return
   end if
 end if
 call get_value(table, key, val, stat=stat)
 if(stat /= toml_stat%success) then
   call key_error(section, key, stat, 'an integer', iout, ifail)
   if(ifail /= 0) return
 end if
 return
 end subroutine get_int

 subroutine get_real(table, section, key, val, default, iout, ifail)

!! a real key (an integer is converted), required unless default is given;
!! nan and inf are errors
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param val  the value of the key, or default if the key is missing
!! @param default  optional, the value for a missing key, which is then not
!!                 required
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use input_text_m, only: finite
 use tomlf, only: toml_table, toml_keyval, toml_stat, get_value, len
 use constants_m, only: dp

!Input variables
 character(len=*), intent(in) :: section, key
 real(dp), intent(in), optional :: default
 integer, intent(in) :: iout

!Input/output variables
 type(toml_table), intent(inout) :: table
 real(dp), intent(inout) :: val

!Output variables
 integer, intent(out) :: ifail

!Local variables
 type(toml_keyval), pointer :: ptr
 integer :: stat

 ifail = 0
 call get_value(table, key, ptr, requested=.false., stat=stat)
 if(.not. associated(ptr)) then
   if(present(default) .and. stat == toml_stat%success) then
     val = default
     return
   end if
   call key_error(section, key, stat, 'a number', iout, ifail)
   if(ifail /= 0) return
 end if
 call get_value(ptr, val, stat=stat)
 if(stat /= toml_stat%success) then
   call key_error(section, key, stat, 'a number', iout, ifail)
   if(ifail /= 0) return
 end if
 if(.not. finite(val)) then
   call input_error(key//' in ['//section//'] must be finite', iout, ifail)
   if(ifail /= 0) return
 end if
 return
 end subroutine get_real

 subroutine get_string(table, section, key, val, default, iout, ifail)

!! a string key, required unless default is given
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param val  the value of the key, or default if the key is missing
!! @param default  optional, the value for a missing key, which is then not
!!                 required
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use tomlf, only: toml_table, toml_keyval, toml_stat, get_value, len

!Input variables
 character(len=*), intent(in) :: section, key
 character(len=*), intent(in), optional :: default
 integer, intent(in) :: iout

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 character(len=:), allocatable, intent(out) :: val
 integer, intent(out) :: ifail

!Local variables
 type(toml_keyval), pointer :: ptr
 integer :: stat

 ifail = 0
 if(present(default)) then
   call get_value(table, key, ptr, requested=.false., stat=stat)
   if(.not. associated(ptr) .and. stat == toml_stat%success) then
     val = default
     return
   end if
 end if
 call get_value(table, key, val, stat=stat)
 if(stat /= toml_stat%success .or. .not. allocated(val)) then
   call key_error(section, key, stat, 'a string', iout, ifail)
   if(ifail /= 0) return
 end if
 return
 end subroutine get_string

 subroutine get_array(table, section, key, nn, arr, iout, ifail)

!! the array key of a table, required, with nn elements
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param nn  required number of elements
!! @param arr  the array of the key
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use tomlf, only: toml_table, toml_array, get_value, len

!Input variables
 character(len=*), intent(in) :: section, key
 integer, intent(in) :: nn
 integer, intent(in) :: iout

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 type(toml_array), pointer, intent(out) :: arr
 integer, intent(out) :: ifail

!Local variables
 integer :: stat
 character(len=24) :: counts

 ifail = 0
 call get_value(table, key, arr, requested=.false., stat=stat)
 if(.not. associated(arr)) then
   call key_error(section, key, stat, 'an array', iout, ifail)
   if(ifail /= 0) return
 end if
 if(len(arr) /= nn) then
   write(counts, '(i0,a,i0)') len(arr), ' instead of ', nn
   call input_error(key//' in ['//section//'] has '//trim(counts)// &
&    ' elements', iout, ifail)
   if(ifail /= 0) return
 end if
 return
 end subroutine get_array

 subroutine get_int_array(table, section, key, nn, vals, iout, ifail)

!! an integer array key, required, with nn elements
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param nn  required number of elements
!! @param vals  the values, in vals(1:nn); the rest is unchanged
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use tomlf, only: toml_table, toml_array, toml_stat, get_value, len

!Input variables
 character(len=*), intent(in) :: section, key
 integer, intent(in) :: nn
 integer, intent(in) :: iout

!Input/output variables
 type(toml_table), intent(inout) :: table
 integer, intent(inout) :: vals(:)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 type(toml_array), pointer :: arr
 integer :: ii, stat

 ifail = 0
 call get_array(table, section, key, nn, arr, iout, ifail)
 if(ifail /= 0) return
 do ii = 1, nn
   call get_value(arr, ii, vals(ii), stat=stat)
   if(stat /= toml_stat%success) then
     call input_error(key//' in ['//section//'] must hold integers', &
&          iout, ifail)
     if(ifail /= 0) return
   end if
 end do
 return
 end subroutine get_int_array

 subroutine get_real_array(table, section, key, nn, vals, iout, ifail)

!! a real array key (integers are converted), required, with nn elements;
!! nan and inf are errors
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param nn  required number of elements
!! @param vals  the values, in vals(1:nn); the rest is unchanged
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use input_text_m, only: finite
 use tomlf, only: toml_table, toml_array, toml_stat, get_value, len
 use constants_m, only: dp

!Input variables
 character(len=*), intent(in) :: section, key
 integer, intent(in) :: nn
 integer, intent(in) :: iout

!Input/output variables
 type(toml_table), intent(inout) :: table
 real(dp), intent(inout) :: vals(:)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 type(toml_array), pointer :: arr
 integer :: ii, stat

 ifail = 0
 call get_array(table, section, key, nn, arr, iout, ifail)
 if(ifail /= 0) return
 do ii = 1, nn
   call get_value(arr, ii, vals(ii), stat=stat)
   if(stat /= toml_stat%success) then
     call input_error(key//' in ['//section//'] must hold numbers', iout, ifail)
     if(ifail /= 0) return
   end if
   if(.not. finite(vals(ii))) then
     call input_error(key//' in ['//section//'] must hold finite numbers', &
&      iout,ifail)
     if(ifail /= 0) return
   end if
 end do
 return
 end subroutine get_real_array

 subroutine check_l(lt, lmax, section, iout, ifail)

!! the angular momenta of the rows of a section must be 0..lmax in order
!!
!! @param lt  the angular momenta of the rows
!! @param lmax  maximum angular momentum
!! @param section  name of the table, for messages
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use tomlf, only: len

!Input variables
 integer, intent(in) :: lt(:), lmax
 character(len=*), intent(in) :: section
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: l1

 ifail = 0
 do l1 = 1, lmax + 1
   if(lt(l1) /= l1 - 1) then
     call input_error('l in ['//section//'] must be 0, 1, ..., lmax in order', &
&          iout, ifail)
     if(ifail /= 0) return
   end if
 end do
 return
 end subroutine check_l

 end module input_toml_m
