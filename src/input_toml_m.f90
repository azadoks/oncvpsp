module input_toml_m
   use, intrinsic :: iso_fortran_env, only: dp => real64, stdout => output_unit, &
                                            stderr => error_unit
   ! TOML Fortran library: documentation at toml-f.readthedocs.io
   use tomlf, only: toml_table, toml_array, toml_keyval, toml_error, toml_stat, &
                    toml_load, get_value, len
   implicit none
   private
   public :: read_input_toml, read_input_toml_r
   contains

   ! Reads the TOML input of the scalar- and non-relativistic drivers: the
   ! same data as read_input_text, with the same defaults and the same checks
   ! against the array sizes (check_data validates the values afterwards).
   ! Every table and key is required unless stated otherwise; a missing or
   ! mistyped one stops the program with a message naming it.
   !
   !   [oncvpsp]                  atsym, z, nc, nv, iexc
   !   [reference_configuration]  n, l, f                (nc + nv each)
   !   [pseudopotentials]         lmax, and l (0..lmax), rc, ep, ncon, nbas, qcut
   !   [local_potential]          lloc, lpopt, rcloc, dvloc0
   !   [vkb_projectors]           l (0..lmax), nproj, debl
   !   [model_core_charge]        icmod, fcfact (icmod 1, 2, 3, 6; else 0),
   !                              rcfact (icmod 3; 1e-4 for icmod 6; else 0),
   !                              fcfact_min, fcfact_max, fcfact_step, rcfact_min,
   !                              rcfact_max, rcfact_step (optional, icmod 4 grid)
   !   [log_derivative_analysis]  epsh1, epsh2, depsh, rxpsh (optional, -1 for auto)
   !   [linear_mesh]              a (drl), rmax (rlmax)
   !   [log_mesh]                 amesh, rr1, rmax (optional table and keys, each
   !                              0 or less or missing for the built-in mesh)
   !   [pp_output]                psfile
   !   [[test_configurations]]    n, l, f (valence states only; optional, up to 4)
   !
   ! scripts/input_text_to_toml.py converts a text input to this form.
   subroutine read_input_toml(unit, &
                              atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                              lmax, rc, ep, ncon, nbas, qcut, &
                              lloc, lpopt, dvloc0, nproj, debl, &
                              icmod, fcfact, rcfact, &
                              fcfact_min, fcfact_max, fcfact_step, &
                              rcfact_min, rcfact_max, rcfact_step, &
                              epsh1, epsh2, depsh, rxpsh, &
                              rlmax, drl, &
                              ncnf, nvcnf, nacnf, lacnf, facnf, &
                              amesh, rr1, rmax)
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
      integer, intent(out) :: na(30)
      !> Orbital angular momentum quantum number array
      integer, intent(out) :: la(30)
      !> Occupation number array
      real(dp), intent(out) :: fa(30)
      ! [Pseudopotential and optimization]
      !> Maximum angular momentum
      integer, intent(out) :: lmax
      !> Core radii for pseudopotentials (rc(5) is the local potential's)
      real(dp), intent(out) :: rc(6)
      !> Pseudopotential energies
      real(dp), intent(out) :: ep(6)
      !> Number of matching constraints for each pseudopotential
      integer, intent(out) :: ncon(6)
      !> Number of basis functions for each pseudopotential
      integer, intent(out) :: nbas(6)
      !> Maximum wave number for each pseudopotential
      real(dp), intent(out) :: qcut(6)
      ! [Local potential]
      !> Angular momentum used for local potential
      integer, intent(out) :: lloc
      !> Polynomial type of the local potential
      integer, intent(out) :: lpopt
      !> Local potential offset at origin
      real(dp), intent(out) :: dvloc0
      ! [Vanderbilt-Kleinman-Bylander projectors]
      !> Number of projectors for each angular momentum
      integer, intent(out) :: nproj(6)
      !> Energy shift for unbound states
      real(dp), intent(out) :: debl(6)
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
      ! [Test configurations]
      !> Number of test configurations
      integer, intent(out) :: ncnf
      !> Number of valence states in test configurations
      integer, intent(out) :: nvcnf(5)
      !> Principal quantum number array for test configurations
      integer, intent(out) :: nacnf(30, 5)
      !> Orbital angular momentum quantum number array for test configurations
      integer, intent(out) :: lacnf(30, 5)
      !> Occupation number array for test configurations
      real(dp), intent(out) :: facnf(30, 5)
      !> Ratio of successive points of the log radial mesh (0 or less for the default)
      real(dp), intent(out) :: amesh
      !> First point of the log radial mesh (0 or less for the default)
      real(dp), intent(out) :: rr1
      !> Outer radius of the log radial mesh (0 or less for the default)
      real(dp), intent(out) :: rmax

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
      integer :: lt(6)
      !> Temporary string
      character(len=:), allocatable :: tmp_str

      ! the defaults of read_input_text
      nproj(:) = 0
      fcfact = 0.0_dp
      rcfact = 0.0_dp
      rc(:) = 0.0_dp
      ep(:) = 0.0_dp
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
      if (allocated(error)) call input_error('cannot parse the input: '//error%message)
      if (.not. allocated(table)) call input_error('cannot parse the input')

      ! atom
      call get_section(table, 'oncvpsp', child)
      call get_string(child, 'oncvpsp', 'atsym', tmp_str)
      if (len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 2) then
         call input_error('atsym in [oncvpsp] must have 1 or 2 characters, got "'//tmp_str//'"')
      end if
      atsym = trim(tmp_str)
      call get_real(child, 'oncvpsp', 'z', zz)
      call get_int(child, 'oncvpsp', 'nc', nc)
      call get_int(child, 'oncvpsp', 'nv', nv)
      call get_int(child, 'oncvpsp', 'iexc', iexc)
      ! counts beyond the array sizes would overrun before check_data runs
      if (nc < 0 .or. nv < 0 .or. nc + nv > 30) then
         call input_error('nc and nv in [oncvpsp] must be >= 0 with nc + nv <= 30')
      end if

      ! reference configuration
      call get_section(table, 'reference_configuration', child)
      call get_int_array(child, 'reference_configuration', 'n', nc + nv, na)
      call get_int_array(child, 'reference_configuration', 'l', nc + nv, la)
      call get_real_array(child, 'reference_configuration', 'f', nc + nv, fa)

      ! pseudopotential and optimization
      call get_section(table, 'pseudopotentials', child)
      call get_int(child, 'pseudopotentials', 'lmax', lmax)
      if (lmax < 0 .or. lmax > 5) call input_error('lmax in [pseudopotentials] must be in 0..5')
      call get_int_array(child, 'pseudopotentials', 'l', lmax + 1, lt)
      call check_l(lt, lmax, 'pseudopotentials')
      call get_real_array(child, 'pseudopotentials', 'rc', lmax + 1, rc)
      call get_real_array(child, 'pseudopotentials', 'ep', lmax + 1, ep)
      call get_int_array(child, 'pseudopotentials', 'ncon', lmax + 1, ncon)
      call get_int_array(child, 'pseudopotentials', 'nbas', lmax + 1, nbas)
      call get_real_array(child, 'pseudopotentials', 'qcut', lmax + 1, qcut)

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
      select case (icmod)
      case (1, 2, 6)
         call get_real(child, 'model_core_charge', 'fcfact', fcfact)
      case (3)
         call get_real(child, 'model_core_charge', 'fcfact', fcfact)
         call get_real(child, 'model_core_charge', 'rcfact', rcfact)
      case default
         call get_real(child, 'model_core_charge', 'fcfact', fcfact, 0.0_dp)
      end select
      select case (icmod)
      case (4, 5)
         ! fcfact and rcfact are ignored; the grid of the Teter parameters
         call get_real(child, 'model_core_charge', 'rcfact', rcfact, 0.0_dp)
         call get_real(child, 'model_core_charge', 'fcfact_min', fcfact_min, fcfact_min)
         call get_real(child, 'model_core_charge', 'fcfact_max', fcfact_max, fcfact_max)
         call get_real(child, 'model_core_charge', 'fcfact_step', fcfact_step, fcfact_step)
         call get_real(child, 'model_core_charge', 'rcfact_min', rcfact_min, rcfact_min)
         call get_real(child, 'model_core_charge', 'rcfact_max', rcfact_max, rcfact_max)
         call get_real(child, 'model_core_charge', 'rcfact_step', rcfact_step, rcfact_step)
      case (6)
         call get_real(child, 'model_core_charge', 'rcfact', rcfact, 1.0e-4_dp)
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
      if (associated(child)) then
         call get_real(child, 'log_mesh', 'amesh', amesh, -1.0_dp)
         call get_real(child, 'log_mesh', 'rr1', rr1, -1.0_dp)
         call get_real(child, 'log_mesh', 'rmax', rmax, -1.0_dp)
      else if (stat /= toml_stat%success) then
         call input_error('[log_mesh] must be a table')
      end if

      ! pseudopotential output file
      call get_section(table, 'pp_output', child)
      call get_string(child, 'pp_output', 'psfile', tmp_str)
      if (len_trim(tmp_str) < 1 .or. len_trim(tmp_str) > 4) then
         call input_error('psfile in [pp_output] must have 1 to 4 characters, got "'//tmp_str//'"')
      end if
      psfile = trim(tmp_str)

      ! test configurations (optional): the valence states, after the
      ! reference core states the driver copies in
      call get_value(table, 'test_configurations', arr, requested=.false., stat=stat)
      if (stat /= toml_stat%success) then
         call input_error('test_configurations must be an array of tables, [[test_configurations]]')
      end if
      if (associated(arr)) then
         ncnf = len(arr)
         if (ncnf > 4) call input_error('at most 4 [[test_configurations]] are allowed')
         do jj = 2, ncnf + 1
            call get_value(arr, jj - 1, child, stat=stat)
            if (.not. associated(child)) then
               call input_error('test_configurations must be an array of tables, ' // &
                                '[[test_configurations]]')
            end if
            ! the number of valence states is that of n
            call get_value(child, 'n', narr, requested=.false., stat=stat)
            if (.not. associated(narr)) call key_error('test_configurations', 'n', stat, 'an array')
            nvcnf(jj) = len(narr)
            if (nc + nvcnf(jj) > 30) then
               call input_error('nc plus the states of a [[test_configurations]] must be <= 30')
            end if
            call get_int_array(child, 'test_configurations', 'n', nvcnf(jj), nacnf(nc + 1:, jj))
            call get_int_array(child, 'test_configurations', 'l', nvcnf(jj), lacnf(nc + 1:, jj))
            call get_real_array(child, 'test_configurations', 'f', nvcnf(jj), facnf(nc + 1:, jj))
         end do
      end if
   end subroutine read_input_toml

   ! Reads the TOML input of the fully-relativistic driver: that of
   ! read_input_toml, with the same ep for both kappa, as read_input_text_r
   subroutine read_input_toml_r(unit, &
                                atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                                lmax, rc, ep, ncon, nbas, qcut, &
                                lloc, lpopt, dvloc0, nproj, debl, &
                                icmod, fcfact, rcfact, fcfact_min, fcfact_max, fcfact_step, &
                                rcfact_min, rcfact_max, rcfact_step, &
                                epsh1, epsh2, depsh, rxpsh, &
                                rlmax, drl, &
                                ncnf, nvcnf, nacnf, lacnf, facnf, &
                                amesh, rr1, rmax)
      implicit none
      integer, intent(in) :: unit
      character(len=2), intent(out) :: atsym
      real(dp), intent(out) :: zz
      integer, intent(out) :: nc, nv, iexc
      character(len=4), intent(out) :: psfile
      integer, intent(out) :: na(30), la(30)
      real(dp), intent(out) :: fa(30)
      integer, intent(out) :: lmax
      real(dp), intent(out) :: rc(6)
      !> Pseudopotential energies for both kappa
      real(dp), intent(out) :: ep(6, 2)
      integer, intent(out) :: ncon(6), nbas(6)
      real(dp), intent(out) :: qcut(6)
      integer, intent(out) :: lloc, lpopt
      real(dp), intent(out) :: dvloc0
      integer, intent(out) :: nproj(6)
      real(dp), intent(out) :: debl(6)
      integer, intent(out) :: icmod
      real(dp), intent(out) :: fcfact, rcfact, fcfact_min, fcfact_max, fcfact_step
      real(dp), intent(out) :: rcfact_min, rcfact_max, rcfact_step
      real(dp), intent(out) :: epsh1, epsh2, depsh, rxpsh
      real(dp), intent(out) :: rlmax, drl
      integer, intent(out) :: ncnf, nvcnf(5), nacnf(30, 5), lacnf(30, 5)
      real(dp), intent(out) :: facnf(30, 5)
      real(dp), intent(out) :: amesh, rr1, rmax

      real(dp) :: ep1(6)

      call read_input_toml(unit, &
                           atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                           lmax, rc, ep1, ncon, nbas, qcut, &
                           lloc, lpopt, dvloc0, nproj, debl, &
                           icmod, fcfact, rcfact, &
                           fcfact_min, fcfact_max, fcfact_step, &
                           rcfact_min, rcfact_max, rcfact_step, &
                           epsh1, epsh2, depsh, rxpsh, &
                           rlmax, drl, &
                           ncnf, nvcnf, nacnf, lacnf, facnf, &
                           amesh, rr1, rmax)
      ep(:, 1) = ep1(:)
      ep(:, 2) = ep1(:)
   end subroutine read_input_toml_r

   ! report an input error on stdout (the output file) and stderr, and stop
   subroutine input_error(message)
      character(len=*), intent(in) :: message
      write (stdout, '(a)') 'TOML input ERROR: '//message
      write (stdout, '(a)') 'Program will stop'
      write (stderr, '(a)') 'Error: '//message
      stop 1
   end subroutine input_error

   ! the table [name] of the document
   subroutine get_section(table, name, child)
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: name
      type(toml_table), pointer, intent(out) :: child
      integer :: stat
      ! a missing key leaves child null with success, a mistyped one with an error
      call get_value(table, name, child, requested=.false., stat=stat)
      if (.not. associated(child)) then
         if (stat == toml_stat%success) call input_error('missing table ['//name//']')
         call input_error('['//name//'] must be a table')
      end if
   end subroutine get_section

   ! report a key that a lookup did not return: missing for missing_key, or
   ! for success from a lookup with requested=.false., else mistyped
   subroutine key_error(section, key, stat, what)
      character(len=*), intent(in) :: section, key, what
      integer, intent(in) :: stat
      if (stat == toml_stat%missing_key .or. stat == toml_stat%success) then
         call input_error('missing key '//key//' in ['//section//']')
      else
         call input_error(key//' in ['//section//'] must be '//what)
      end if
   end subroutine key_error

   ! an integer key, required
   subroutine get_int(table, section, key, val)
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: section, key
      integer, intent(out) :: val
      integer :: stat
      call get_value(table, key, val, stat=stat)
      if (stat /= toml_stat%success) call key_error(section, key, stat, 'an integer')
   end subroutine get_int

   ! a real key (an integer is converted), required unless default is given
   subroutine get_real(table, section, key, val, default)
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: section, key
      real(dp), intent(inout) :: val
      real(dp), intent(in), optional :: default
      type(toml_keyval), pointer :: ptr
      integer :: stat
      call get_value(table, key, ptr, requested=.false., stat=stat)
      if (.not. associated(ptr)) then
         if (present(default) .and. stat == toml_stat%success) then
            val = default
            return
         end if
         call key_error(section, key, stat, 'a number')
      end if
      call get_value(ptr, val, stat=stat)
      if (stat /= toml_stat%success) call key_error(section, key, stat, 'a number')
   end subroutine get_real

   ! a string key, required
   subroutine get_string(table, section, key, val)
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: section, key
      character(len=:), allocatable, intent(out) :: val
      integer :: stat
      call get_value(table, key, val, stat=stat)
      if (stat /= toml_stat%success .or. .not. allocated(val)) then
         call key_error(section, key, stat, 'a string')
      end if
   end subroutine get_string

   ! the array key of a table, required, with n elements
   subroutine get_array(table, section, key, n, arr)
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: section, key
      integer, intent(in) :: n
      type(toml_array), pointer, intent(out) :: arr
      integer :: stat
      character(len=24) :: counts
      call get_value(table, key, arr, requested=.false., stat=stat)
      if (.not. associated(arr)) call key_error(section, key, stat, 'an array')
      if (len(arr) /= n) then
         write (counts, '(i0,a,i0)') len(arr), ' instead of ', n
         call input_error(key//' in ['//section//'] has '//trim(counts)//' elements')
      end if
   end subroutine get_array

   ! an integer array key, required, with n elements
   subroutine get_int_array(table, section, key, n, vals)
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: section, key
      integer, intent(in) :: n
      integer, intent(inout) :: vals(:)
      type(toml_array), pointer :: arr
      integer :: ii, stat
      call get_array(table, section, key, n, arr)
      do ii = 1, n
         call get_value(arr, ii, vals(ii), stat=stat)
         if (stat /= toml_stat%success) then
            call input_error(key//' in ['//section//'] must hold integers')
         end if
      end do
   end subroutine get_int_array

   ! a real array key (integers are converted), required, with n elements
   subroutine get_real_array(table, section, key, n, vals)
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: section, key
      integer, intent(in) :: n
      real(dp), intent(inout) :: vals(:)
      type(toml_array), pointer :: arr
      integer :: ii, stat
      call get_array(table, section, key, n, arr)
      do ii = 1, n
         call get_value(arr, ii, vals(ii), stat=stat)
         if (stat /= toml_stat%success) then
            call input_error(key//' in ['//section//'] must hold numbers')
         end if
      end do
   end subroutine get_real_array

   ! the angular momenta of the rows of a section must be 0..lmax in order
   subroutine check_l(lt, lmax, section)
      integer, intent(in) :: lt(:), lmax
      character(len=*), intent(in) :: section
      integer :: l1
      do l1 = 1, lmax + 1
         if (lt(l1) /= l1 - 1) then
            call input_error('l in ['//section//'] must be 0, 1, ..., lmax in order')
         end if
      end do
   end subroutine check_l

end module input_toml_m
