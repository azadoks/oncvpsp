 module input_toml_m
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
! TOML Fortran library: documentation at toml-f.readthedocs.io
 use tomlf, only: toml_table, toml_array, toml_keyval, toml_error, toml_stat, &
&                 toml_load, get_value, len
 use constants_m, only: mxstate, mxcnf, mxlin, mxkap, well_onset_default
 implicit none
 private
 public :: read_input_toml
 contains

 subroutine read_input_toml(unit, &
&                           atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
&                           lmax, rc, ep, ncon, nbas, qcut, &
&                           lloc, lpopt, dvloc0, nproj, debl, &
&                           icmod, fcfact, rcfact, slack, budget_weighted, &
&                           fcfact_min, fcfact_max, fcfact_step, &
&                           rcfact_min, rcfact_max, rcfact_step, &
&                           epsh1, epsh2, depsh, rxpsh, &
&                           rlmax, drl, upf_mesh, &
&                           ncnf, nvcnf, nacnf, lacnf, facnf, &
&                           amesh, rr1, rmax, c_light, ep_as_given, &
&                           well_onset, well_tail_free)

!! Reads the TOML input: the same data as read_input_text, with the same
!! defaults and the same checks against the array sizes (check_data
!! validates the values afterwards).  Every table and key is required unless
!! stated otherwise; a missing or mistyped one stops the program with a
!! message naming it.
!!
!! scripts/input_text_to_toml.py converts a text input to this form.  The
!! tables and their keys are listed below.

!   [oncvpsp]                  atsym, z, nc, nv, iexc, c_light (optional, 0
!                              or less for the default)
!   [reference_configuration]  n, l, f                (nc + nv each)
!   [pseudopotentials]         lmax, and l (0..lmax), rc, ep, ncon, nbas,
!                              qcut, ep_as_given (optional, false: the well
!                              states take ep as given, see wellstate),
!                              well_onset (optional, 1.006**5, at least 1:
!                              the onset of the well of the well states as
!                              a multiple of rc, see wellstate),
!                              well_tail (optional, "well" or "free":
!                              the tail of the well states in E_resid,
!                              see free_tail_m)
!   [local_potential]          lloc, lpopt, rcloc, dvloc0
!   [vkb_projectors]           l (0..lmax), nproj, debl
!   [model_core_charge]        icmod, fcfact (icmod 1, 2, 3, 6, 7; else 0),
!                              rcfact (icmod 3; 1e-5 for icmod 6, 7; else
!                              0), slack (icmod 6, 7; 1), budget_cutoff
!                              (icmod 6, 7; "max" or "weighted", "max"),
!                              fcfact_min, fcfact_max, fcfact_step, rcfact_min,
!                              rcfact_max, rcfact_step (optional, icmod 4 grid)
!   [log_derivative_analysis]  epsh1, epsh2, depsh, rxpsh (optional, -1 for
!                              auto)
!   [linear_mesh]              a (drl), rmax (rlmax)
!   [log_mesh]                 amesh, rr1, rmax (optional table and keys, each
!                              0 or less or missing for the built-in mesh)
!   [pp_output]                psfile, upf_mesh (optional, 0)
!   [[test_configurations]]    n, l, f (valence states only; optional, up
!                              to 4)

 implicit none
!Input variables
!> The Fortran unit number to read from
 integer, intent(in) :: unit

!Output variables
! [Atom and reference configuration]
!> Atomic symbol
 character(len=2), intent(out) :: atsym
!> Atomic number
 real(dp), intent(out) :: zz
!> Number of core states
 integer, intent(out) :: nc
!> Number of valence states
 integer, intent(out) :: nv
!> Exchange-correlation functional code
 integer, intent(out) :: iexc
!> Pseudopotential output file format flag
 character(len=4), intent(out) :: psfile
!> Principal quantum number array
 integer, intent(out) :: na(mxstate)
!> Orbital angular momentum quantum number array
 integer, intent(out) :: la(mxstate)
!> Occupation number array
 real(dp), intent(out) :: fa(mxstate)
! [Pseudopotential and optimization]
!> Maximum angular momentum
 integer, intent(out) :: lmax
!> Core radii for pseudopotentials (rc(5) is the local potential's)
 real(dp), intent(out) :: rc(mxlin+1)
!> Pseudopotential energies, the input value for both j=l+1/2 and l-1/2
!> (only the fully-relativistic code uses the second)
 real(dp), intent(out) :: ep(mxlin+1, mxkap)
!> Number of matching constraints for each pseudopotential
 integer, intent(out) :: ncon(mxlin+1)
!> Number of basis functions for each pseudopotential
 integer, intent(out) :: nbas(mxlin+1)
!> Residual-energy cutoff wave vector for each pseudopotential
 real(dp), intent(out) :: qcut(mxlin+1)
! [Local potential]
!> Angular momentum used for local potential
 integer, intent(out) :: lloc
!> Polynomial type of the local potential
 integer, intent(out) :: lpopt
!> Local potential offset at origin
 real(dp), intent(out) :: dvloc0
! [Vanderbilt-Kleinman-Bylander projectors]
!> Number of projectors for each angular momentum
 integer, intent(out) :: nproj(mxlin+1)
!> use ep as given for the well states (D-1), else wellstate's replacements
 logical, intent(out) :: ep_as_given
!> onset of the well of the well states as a multiple of rc (wellstate)
 real(dp), intent(out) :: well_onset
 logical, intent(out) :: well_tail_free
!> Energy added to ep for the second projector (Ha)
 real(dp), intent(out) :: debl(mxlin+1)
! [Model core charge]
!> Model core charge flag
 integer, intent(out) :: icmod
!> Scaling factor for core charge
 real(dp), intent(out) :: fcfact
!> Scaling factor range minimum for optimization
 real(dp), intent(out) :: fcfact_min
!> Scaling factor range maximum for optimization
 real(dp), intent(out) :: fcfact_max
!> Scaling factor step size for optimization
 real(dp), intent(out) :: fcfact_step
!> Core charge width factor
 real(dp), intent(out) :: rcfact
!> Share of the cutoff the model core must fit (icmod 6, 7)
 real(dp), intent(out) :: slack
!> Whether the cutoff of the pseudo-waves is that of their weighted mean
!> rather than the largest (icmod 6, 7; budget_cutoff "weighted")
 logical, intent(out) :: budget_weighted
!> Core charge width factor minimum for optimization
 real(dp), intent(out) :: rcfact_min
!> Core charge width factor maximum for optimization
 real(dp), intent(out) :: rcfact_max
!> Core charge width factor step size for optimization
 real(dp), intent(out) :: rcfact_step
! [Log derivative analysis]
!> Lower bound of energy range for log derivative analysis
 real(dp), intent(out) :: epsh1
!> Upper bound of energy range for log derivative analysis
 real(dp), intent(out) :: epsh2
!> Energy step size for log derivative analysis
 real(dp), intent(out) :: depsh
!> Radius for log derivative analysis
 real(dp), intent(out) :: rxpsh
! [Output grid]
!> Maximum radius for output grid
 real(dp), intent(out) :: rlmax
!> Grid spacing for output grid
 real(dp), intent(out) :: drl
!> Most points of the log mesh of the UPF file (0 for the linear output grid)
 integer, intent(out) :: upf_mesh
! [Test configurations]
!> Number of test configurations
 integer, intent(out) :: ncnf
!> Number of valence states in test configurations
 integer, intent(out) :: nvcnf(mxcnf)
!> Principal quantum number array for test configurations
 integer, intent(out) :: nacnf(mxstate, mxcnf)
!> Orbital angular momentum quantum number array for test configurations
 integer, intent(out) :: lacnf(mxstate, mxcnf)
!> Occupation number array for test configurations
 real(dp), intent(out) :: facnf(mxstate, mxcnf)
!> Ratio of successive points of the log radial mesh
!> (0 or less for the default)
 real(dp), intent(out) :: amesh
!> First point of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: rr1
!> Outer radius of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: rmax
!> Speed of light (0 or less for the default)
 real(dp), intent(out) :: c_light

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

! the defaults of read_input_text
 nproj(:) = 0
 fcfact = 0.0_dp
 rcfact = 0.0_dp
 slack = 1.0_dp
 budget_weighted = .false.
 rc(:) = 0.0_dp
 ep(:, :) = 0.0_dp
 rxpsh = -1.0_dp  ! negative value for auto setting of radius
 fcfact_min = 1.5_dp
 fcfact_max = 6.0_dp
 fcfact_step = 0.5_dp
 rcfact_min = 1.0_dp
 rcfact_max = 1.9_dp
 rcfact_step = 0.1_dp
 ncnf = 0
 amesh = -1.0_dp  ! negative values for the built-in mesh
 rr1 = -1.0_dp
 rmax = -1.0_dp
 atsym = ' '
 psfile = ' '

 call toml_load(table, unit, error=error)
 if(allocated(error)) then
   call input_error('cannot parse the input: '//error%message)
 end if
 if(.not. allocated(table)) call input_error('cannot parse the input')

! atom
 call get_section(table, 'oncvpsp', child)
 call get_string(child, 'oncvpsp', 'atsym', tmp_str)
 if(len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 2) then
   call input_error('atsym in [oncvpsp] must have 1 or 2 characters, '// &
&    'got "'//tmp_str//'"')
 end if
 atsym = trim(tmp_str)
 call get_real(child, 'oncvpsp', 'z', zz)
 call get_int(child, 'oncvpsp', 'nc', nc)
 call get_int(child, 'oncvpsp', 'nv', nv)
 call get_real(child, 'oncvpsp', 'c_light', c_light, -1.0_dp)
 call get_int(child, 'oncvpsp', 'iexc', iexc)
! counts beyond the array sizes would overrun before check_data runs
 if(nc < 0 .or. nv < 0 .or. nc + nv > mxstate) then
   write(msg, '(a,i0)') 'nc and nv in [oncvpsp] must be >= 0 with '// &
&    'nc + nv <= ', mxstate
   call input_error(trim(msg))
 end if

! reference configuration
 call get_section(table, 'reference_configuration', child)
 call get_int_array(child, 'reference_configuration', 'n', nc + nv, na)
 call get_int_array(child, 'reference_configuration', 'l', nc + nv, la)
 call get_real_array(child, 'reference_configuration', 'f', nc + nv, fa)

! pseudopotential and optimization
 call get_section(table, 'pseudopotentials', child)
 call get_int(child, 'pseudopotentials', 'lmax', lmax)
 if(lmax < 0 .or. lmax > mxlin) then
   write(msg, '(a,i0)') 'lmax in [pseudopotentials] must be in 0..', mxlin
   call input_error(trim(msg))
 end if
 call get_int_array(child, 'pseudopotentials', 'l', lmax + 1, lt)
 call check_l(lt, lmax, 'pseudopotentials')
 call get_real_array(child, 'pseudopotentials', 'rc', lmax + 1, rc)
 call get_real_array(child, 'pseudopotentials', 'ep', lmax + 1, ep(:, 1))
 ep(:, 2) = ep(:, 1)
 call get_int_array(child, 'pseudopotentials', 'ncon', lmax + 1, ncon)
 call get_int_array(child, 'pseudopotentials', 'nbas', lmax + 1, nbas)
 call get_real_array(child, 'pseudopotentials', 'qcut', lmax + 1, qcut)
 call get_logical(child, 'pseudopotentials', 'ep_as_given', ep_as_given, &
&                  .false.)
 call get_real(child, 'pseudopotentials', 'well_onset', well_onset, &
&              well_onset_default)
 if(.not. (well_onset >= 1.0_dp)) then
   call input_error('well_onset in [pseudopotentials] must be at least 1')
 end if
 call get_string(child, 'pseudopotentials', 'well_tail', tmp_str, 'well')
 select case(tmp_str)
 case('well')
   well_tail_free = .false.
 case('free')
   well_tail_free = .true.
 case default
   call input_error('well_tail in [pseudopotentials] must be "well" or '// &
&    '"free", got "'//tmp_str//'"')
 end select

! local potential
 call get_section(table, 'local_potential', child)
 call get_int(child, 'local_potential', 'lloc', lloc)
 call get_int(child, 'local_potential', 'lpopt', lpopt)
 call get_real(child, 'local_potential', 'rcloc', rc(5))
 call get_real(child, 'local_potential', 'dvloc0', dvloc0)

! Vanderbilt-Kleinman-Bylander projectors
 call get_section(table, 'vkb_projectors', child)
 call get_int_array(child, 'vkb_projectors', 'l', lmax + 1, lt)
 call check_l(lt, lmax, 'vkb_projectors')
 call get_int_array(child, 'vkb_projectors', 'nproj', lmax + 1, nproj)
 call get_real_array(child, 'vkb_projectors', 'debl', lmax + 1, debl)

! model core charge; check_data validates icmod and the values
 call get_section(table, 'model_core_charge', child)
 call get_int(child, 'model_core_charge', 'icmod', icmod)
 select case(icmod)
 case(1, 2, 6, 7)
   call get_real(child, 'model_core_charge', 'fcfact', fcfact)
 case(3)
   call get_real(child, 'model_core_charge', 'fcfact', fcfact)
   call get_real(child, 'model_core_charge', 'rcfact', rcfact)
 case default
   call get_real(child, 'model_core_charge', 'fcfact', fcfact, 0.0_dp)
 end select
 select case(icmod)
 case(4, 5)
! fcfact and rcfact are ignored; the grid of the Teter parameters
   call get_real(child, 'model_core_charge', 'rcfact', rcfact, 0.0_dp)
   call get_real(child, 'model_core_charge', 'fcfact_min', fcfact_min, &
&    fcfact_min)
   call get_real(child, 'model_core_charge', 'fcfact_max', fcfact_max, &
&    fcfact_max)
   call get_real(child, 'model_core_charge', 'fcfact_step', fcfact_step, &
&    fcfact_step)
   call get_real(child, 'model_core_charge', 'rcfact_min', rcfact_min, &
&    rcfact_min)
   call get_real(child, 'model_core_charge', 'rcfact_max', rcfact_max, &
&    rcfact_max)
   call get_real(child, 'model_core_charge', 'rcfact_step', rcfact_step, &
&    rcfact_step)
 case(6, 7)
   call get_real(child, 'model_core_charge', 'rcfact', rcfact, 1.0e-5_dp)
   call get_real(child, 'model_core_charge', 'slack', slack, 1.0_dp)
   call get_string(child, 'model_core_charge', 'budget_cutoff', tmp_str, &
&    'max')
   select case(tmp_str)
   case('max')
     budget_weighted = .false.
   case('weighted')
     budget_weighted = .true.
   case default
     call input_error('budget_cutoff in [model_core_charge] must be "max"'// &
&      ' or "weighted", got "'//tmp_str//'"')
   end select
 case default
! no further keys
 end select

! log derivative analysis
 call get_section(table, 'log_derivative_analysis', child)
 call get_real(child, 'log_derivative_analysis', 'epsh1', epsh1)
 call get_real(child, 'log_derivative_analysis', 'epsh2', epsh2)
 call get_real(child, 'log_derivative_analysis', 'depsh', depsh)
 call get_real(child, 'log_derivative_analysis', 'rxpsh', rxpsh, -1.0_dp)

! output grid
 call get_section(table, 'linear_mesh', child)
 call get_real(child, 'linear_mesh', 'a', drl)
 call get_real(child, 'linear_mesh', 'rmax', rlmax)

! log mesh (optional)
 call get_value(table, 'log_mesh', child, requested=.false., stat=stat)
 if(associated(child)) then
   call get_real(child, 'log_mesh', 'amesh', amesh, -1.0_dp)
   call get_real(child, 'log_mesh', 'rr1', rr1, -1.0_dp)
   call get_real(child, 'log_mesh', 'rmax', rmax, -1.0_dp)
 else if(stat /= toml_stat%success) then
   call input_error('[log_mesh] must be a table')
 end if

! pseudopotential output file
 call get_section(table, 'pp_output', child)
 call get_string(child, 'pp_output', 'psfile', tmp_str)
 if(len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 4) then
   call input_error('psfile in [pp_output] must have 1 to 4 characters, '// &
&    'got "'//tmp_str//'"')
 end if
 psfile = trim(tmp_str)
 call get_int(child, 'pp_output', 'upf_mesh', upf_mesh, 0)

! test configurations (optional): the valence states, after the
! reference core states the driver copies in
 call get_value(table, 'test_configurations', arr, requested=.false., &
&  stat=stat)
 if(stat /= toml_stat%success) then
   call input_error('test_configurations must be an array of tables, '// &
&    '[[test_configurations]]')
 end if
 if(associated(arr)) then
   ncnf = len(arr)
   if(ncnf > mxcnf - 1) then
     write(msg, '(a,i0,a)') 'at most ', mxcnf - 1, &
&      ' [[test_configurations]] are allowed'
     call input_error(trim(msg))
   end if
   do jj = 2, ncnf + 1
     call get_value(arr, jj - 1, child, stat=stat)
     if(.not. associated(child)) then
       call input_error('test_configurations must be an array of tables, ' // &
&                       '[[test_configurations]]')
     end if
! the number of valence states is that of n
     call get_value(child, 'n', narr, requested=.false., stat=stat)
     if(.not. associated(narr)) then
       call key_error('test_configurations', 'n', stat, 'an array')
     end if
     nvcnf(jj) = len(narr)
     if(nc + nvcnf(jj) > mxstate) then
       write(msg, '(a,i0)') 'nc plus the states of a '// &
&        '[[test_configurations]] must be <= ', mxstate
       call input_error(trim(msg))
     end if
     call get_int_array(child, 'test_configurations', 'n', nvcnf(jj), &
&      nacnf(nc + 1:, jj))
     call get_int_array(child, 'test_configurations', 'l', nvcnf(jj), &
&      lacnf(nc + 1:, jj))
     call get_real_array(child, 'test_configurations', 'f', nvcnf(jj), &
&      facnf(nc + 1:, jj))
   end do
 end if
 return
 end subroutine read_input_toml

 subroutine input_error(message)

!! report an input error in the output (unit 6) and on stderr, and stop
!!
!! @param message  description of the error

!Input variables
 character(len=*), intent(in) :: message

 write(6,'(a)') 'TOML input ERROR: '//message
 write(6,'(a)') 'Program will stop'
 write(stderr,'(a)') 'Error: '//message
 error stop 1
 end subroutine input_error

 subroutine get_section(table, name, child)

!! the table [name] of the document
!!
!! @param table  the TOML document
!! @param name  name of the table
!! @param child  the table [name]; the program stops if it is missing or not a
!!               table

!Input variables
 character(len=*), intent(in) :: name

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 type(toml_table), pointer, intent(out) :: child

!Local variables
 integer :: stat

! a missing key leaves child null with success, a mistyped one with an error
 call get_value(table, name, child, requested=.false., stat=stat)
 if(.not. associated(child)) then
   if(stat == toml_stat%success) call input_error('missing table ['//name//']')
   call input_error('['//name//'] must be a table')
 end if
 return
 end subroutine get_section

 subroutine key_error(section, key, stat, what)

!! report a key that a lookup did not return: missing for missing_key, or
!! for success from a lookup with requested=.false., else mistyped
!!
!! @param section  name of the table holding the key
!! @param key  name of the key
!! @param stat  status of the lookup
!! @param what  the expected type, for the message (e.g. 'an integer')

!Input variables
 character(len=*), intent(in) :: section, key, what
 integer, intent(in) :: stat

 if(stat == toml_stat%missing_key .or. stat == toml_stat%success) then
   call input_error('missing key '//key//' in ['//section//']')
 else
   call input_error(key//' in ['//section//'] must be '//what)
 end if
 return
 end subroutine key_error

 subroutine get_logical(table, section, key, val, default)

!! a logical key, required unless default is given
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param val  the value of the key, or default if the key is missing
!! @param default  optional, the value for a missing key, which is then not
!!                 required

!Input variables
 character(len=*), intent(in) :: section, key
 logical, intent(in), optional :: default

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 logical, intent(out) :: val

!Local variables
 type(toml_keyval), pointer :: ptr
 integer :: stat

 if(present(default)) then
   call get_value(table, key, ptr, requested=.false., stat=stat)
   if(.not. associated(ptr) .and. stat == toml_stat%success) then
     val = default
     return
   end if
 end if
 call get_value(table, key, val, stat=stat)
 if(stat /= toml_stat%success) then
   call key_error(section, key, stat, 'a logical')
 end if
 return
 end subroutine get_logical

 subroutine get_int(table, section, key, val, default)

!! an integer key, required unless default is given
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param val  the value of the key, or default if the key is missing
!! @param default  optional, the value for a missing key, which is then not
!!                 required

!Input variables
 character(len=*), intent(in) :: section, key
 integer, intent(in), optional :: default

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 integer, intent(out) :: val

!Local variables
 type(toml_keyval), pointer :: ptr
 integer :: stat

 if(present(default)) then
   call get_value(table, key, ptr, requested=.false., stat=stat)
   if(.not. associated(ptr) .and. stat == toml_stat%success) then
     val = default
     return
   end if
 end if
 call get_value(table, key, val, stat=stat)
 if(stat /= toml_stat%success) then
   call key_error(section, key, stat, 'an integer')
 end if
 return
 end subroutine get_int

 subroutine get_real(table, section, key, val, default)

!! a real key (an integer is converted), required unless default is given
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param val  the value of the key, or default if the key is missing
!! @param default  optional, the value for a missing key, which is then not
!!                 required

!Input variables
 character(len=*), intent(in) :: section, key
 real(dp), intent(in), optional :: default

!Input/output variables
 type(toml_table), intent(inout) :: table
 real(dp), intent(inout) :: val

!Local variables
 type(toml_keyval), pointer :: ptr
 integer :: stat

 call get_value(table, key, ptr, requested=.false., stat=stat)
 if(.not. associated(ptr)) then
   if(present(default) .and. stat == toml_stat%success) then
     val = default
     return
   end if
   call key_error(section, key, stat, 'a number')
 end if
 call get_value(ptr, val, stat=stat)
 if(stat /= toml_stat%success) call key_error(section, key, stat, 'a number')
 return
 end subroutine get_real

 subroutine get_string(table, section, key, val, default)

!! a string key, required unless default is given
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param val  the value of the key, or default if the key is missing
!! @param default  optional, the value for a missing key, which is then not
!!                 required

!Input variables
 character(len=*), intent(in) :: section, key
 character(len=*), intent(in), optional :: default

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 character(len=:), allocatable, intent(out) :: val

!Local variables
 type(toml_keyval), pointer :: ptr
 integer :: stat

 if(present(default)) then
   call get_value(table, key, ptr, requested=.false., stat=stat)
   if(.not. associated(ptr) .and. stat == toml_stat%success) then
     val = default
     return
   end if
 end if
 call get_value(table, key, val, stat=stat)
 if(stat /= toml_stat%success .or. .not. allocated(val)) then
   call key_error(section, key, stat, 'a string')
 end if
 return
 end subroutine get_string

 subroutine get_array(table, section, key, nn, arr)

!! the array key of a table, required, with nn elements
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param nn  required number of elements
!! @param arr  the array of the key

!Input variables
 character(len=*), intent(in) :: section, key
 integer, intent(in) :: nn

!Input/output variables
 type(toml_table), intent(inout) :: table

!Output variables
 type(toml_array), pointer, intent(out) :: arr

!Local variables
 integer :: stat
 character(len=24) :: counts

 call get_value(table, key, arr, requested=.false., stat=stat)
 if(.not. associated(arr)) call key_error(section, key, stat, 'an array')
 if(len(arr) /= nn) then
   write(counts, '(i0,a,i0)') len(arr), ' instead of ', nn
   call input_error(key//' in ['//section//'] has '//trim(counts)// &
&    ' elements')
 end if
 return
 end subroutine get_array

 subroutine get_int_array(table, section, key, nn, vals)

!! an integer array key, required, with nn elements
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param nn  required number of elements
!! @param vals  the values, in vals(1:nn); the rest is unchanged

!Input variables
 character(len=*), intent(in) :: section, key
 integer, intent(in) :: nn

!Input/output variables
 type(toml_table), intent(inout) :: table
 integer, intent(inout) :: vals(:)

!Local variables
 type(toml_array), pointer :: arr
 integer :: ii, stat

 call get_array(table, section, key, nn, arr)
 do ii = 1, nn
   call get_value(arr, ii, vals(ii), stat=stat)
   if(stat /= toml_stat%success) then
     call input_error(key//' in ['//section//'] must hold integers')
   end if
 end do
 return
 end subroutine get_int_array

 subroutine get_real_array(table, section, key, nn, vals)

!! a real array key (integers are converted), required, with nn elements
!!
!! @param table  the table holding the key
!! @param section  name of the table, for messages
!! @param key  name of the key
!! @param nn  required number of elements
!! @param vals  the values, in vals(1:nn); the rest is unchanged

!Input variables
 character(len=*), intent(in) :: section, key
 integer, intent(in) :: nn

!Input/output variables
 type(toml_table), intent(inout) :: table
 real(dp), intent(inout) :: vals(:)

!Local variables
 type(toml_array), pointer :: arr
 integer :: ii, stat

 call get_array(table, section, key, nn, arr)
 do ii = 1, nn
   call get_value(arr, ii, vals(ii), stat=stat)
   if(stat /= toml_stat%success) then
     call input_error(key//' in ['//section//'] must hold numbers')
   end if
 end do
 return
 end subroutine get_real_array

 subroutine check_l(lt, lmax, section)

!! the angular momenta of the rows of a section must be 0..lmax in order
!!
!! @param lt  the angular momenta of the rows
!! @param lmax  maximum angular momentum
!! @param section  name of the table, for messages

!Input variables
 integer, intent(in) :: lt(:), lmax
 character(len=*), intent(in) :: section

!Local variables
 integer :: l1

 do l1 = 1, lmax + 1
   if(lt(l1) /= l1 - 1) then
     call input_error('l in ['//section//'] must be 0, 1, ..., lmax in order')
   end if
 end do
 return
 end subroutine check_l

 end module input_toml_m
