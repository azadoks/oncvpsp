program oncvatom
   use, intrinsic :: iso_fortran_env, only: stdin => input_unit, stdout => output_unit, stderr => error_unit
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use tomlf, only: toml_table, toml_array, toml_dump, set_value, add_array, destroy
   use input_text_m, only: read_input_text
   use input_toml_m, only: read_input_toml
   implicit none
   ! Constants
   integer, parameter :: INPUT_STDIN=1, INPUT_TEXT=2, INPUT_TOML=3
   real(dp), parameter :: eps = 1.0e-8_dp
   real(dp), parameter :: al = log(1.006_dp)
   integer, parameter :: maxprj = 5
   integer, parameter :: maxstates = 30
   integer, parameter :: maxl1 = 6
   integer, parameter :: maxcnf = 5
   real(dp), parameter :: rr1_max = 0.0005_dp / 10
   real(dp), parameter :: rrn_max = 45.0_dp
   logical, parameter :: srel = .true.

   ! Input variables
   character(len=1024) :: input_filename
   integer :: input_mode
   ! Control variables
   character(len=2) :: atsym
   real(dp) :: zz
   integer :: nc
   integer :: nv
   integer :: iexc
   character(len=4) :: psfile
   ! Reference configuration
   !> Principal quantum numbers
   integer :: na(maxstates)
   !> Angular momenta
   integer :: la(maxstates)
   !> Occupations
   real(dp) :: fa(maxstates)
   ! Pseudopotential and optimization
   integer :: lmax
   real(dp) :: rc(maxl1)
   real(dp) :: ep(maxl1)
   integer :: ncon(maxl1)
   integer :: nbas(maxl1)
   real(dp) :: qcut(maxl1)
   ! Local potential
   integer :: lloc
   integer :: lpopt
   real(dp) :: dvloc0
   ! VKB projectors
   integer :: nproj(maxl1)
   real(dp) :: debl(maxl1)
   ! Model core charge
   integer :: icmod
   real(dp) :: fcfact
   real(dp) :: rcfact
   real(dp) :: fcfact_min
   real(dp) :: fcfact_max
   integer :: fcfact_step
   real(dp) :: rcfact_min
   real(dp) :: rcfact_max
   integer :: rcfact_step
   ! Logarithmic derivative analysis
   real(dp) :: epsh1
   real(dp) :: epsh2
   real(dp) :: depsh
   real(dp) :: rxpsh
   ! Linear grid
   real(dp) :: rlmax
   real(dp) :: drl
   ! Test configurations
   integer :: ncnf
   integer :: nvcnf(maxcnf)
   integer :: nacnf(maxstates,maxcnf)
   integer :: lacnf(maxstates,maxcnf)
   real(dp) :: facnf(maxstates,maxcnf)
   ! Output files
   character(len=256) :: upffile
   character(len=16) :: upfgrid
   character(len=256) :: psp8file
   character(len=256) :: psmlfile

   ! Log radial grid
   real(dp) :: rr1
   integer :: mmax
   real(dp), allocatable :: rr(:)

   ! All-electron DFT outputs
   real(dp) :: ea(maxstates)
   real(dp) :: rpk(maxstates)
   integer :: imatch(maxstates)
   real(dp) :: rmatch(maxstates)
   integer :: it
   real(dp) :: etot
   real(dp) :: eeel
   real(dp) :: eexc
   integer :: ierr
   real(dp), allocatable :: rhoc(:)
   real(dp), allocatable :: rhov(:)
   real(dp), allocatable :: rhot(:)
   real(dp), allocatable :: vxc(:)
   real(dp), allocatable :: vha(:)
   real(dp), allocatable :: vfull(:)
   real(dp), allocatable :: uuaea(:,:)
   real(dp), allocatable :: upaea(:,:)

   ! TOML output variables
   type(toml_table) :: table
   type(toml_array), pointer :: array

   ! Local variables
   integer :: ii
   real(dp) :: e_in_out

   ! Read and setup input
   call parse_args(input_filename, input_mode)
   call read_input( &
      input_filename, input_mode, &
      atsym, zz, nc, nv, iexc, psfile, lmax, rc, ep, ncon, nbas, &
      qcut, lloc, lpopt, dvloc0, nproj, debl, icmod, fcfact, rcfact, epsh1, epsh2, depsh, &
      rlmax, drl, ncnf, nvcnf, nacnf, lacnf, facnf, upffile, upfgrid, psp8file, psmlfile
   )
   call set_first_test_config(nc, nv, na, la, fa, nvcnf, nacnf, lacnf, facnf)
   call set_model_core_params( &
      icmod, fcfact, rcfact, fcfact_min, fcfact_max, fcfact_step, rcfact_min, rcfact_max, &
      rcfact_step
   )
   call check_data( &
      atsym, zz, nc, nv, iexc, psfile, lmax, rc, ep, ncon, nbas, &
      qcut, lloc, lpopt, dvloc0, nproj, debl, icmod, fcfact, rcfact, epsh1, epsh2, depsh, &
      rlmax, drl, ncnf, nvcnf, nacnf, lacnf, facnf, upffile, upfgrid, psp8file, psmlfile
   )
   ! Echo input in .dat format
   call echo_input(
      atsym, zz, nc, nv, iexc, psfile, lmax, rc, ep, ncon, nbas, &
      qcut, lloc, lpopt, dvloc0, nproj, debl, icmod, fcfact, rcfact, epsh1, epsh2, depsh, &
      rlmax, drl, ncnf, nvcnf, nacnf, lacnf, facnf, upffile, upfgrid, psp8file, psmlfile
   )
   ! Set up logarithmic radial grid
   rr1 = min(rr1_max, 0.0005_dp / zz)
   mmax = int(log(rrn_max / rr1) / al)
   allocate(rr(mmax))
   do ii = 1, mmax
      rr(ii) = rr1 * exp(al * (ii - 1))
   end do
   ! Run all-electron atomic DFT calculation
   allocate(rhoc(mmax), rhov(mmax), rhot(mmax), vxc(mmax), vha(mmax), vfull(mmax))
   allocate(uuaea(mmax, nc + nv), upaea(mmax, nc + nv))
   ! Compute ea, rpk, it, rhoc, rhot, vfull, etot, and ierr
   call sratom(na, la, ea, fa, rpk, nc, nc + nv, it, rhoc, rhot, rr, vfull, zz, mmax, iexc, etot, ierr, srel)
   if (ierr /= 0) then
      write(stderr, '(a,i0)') 'Error in sratom, ierr = ', ierr
      stop 1
   end if
   ! Compute eexc, vxc, and eeel
   call vout(0, rhot, rhoc, vfull, vxc, sum(fa), eeel, eexc, rr, mmax, iexc)
   ! Compute valence charge density
   rhov(:) = rhot(:) - rhoc(:)
   ! Re-compute the wavefunctions in the full potential
   do ii = 1, nc + nv
      e_in_out = ea(ii)
      call lschfb(na(ii), la(ii), ierr, e_in_out, rr, vfull, uuaea(:, ii), upaea(:, ii), zz, mmax, imatch(ii), srel)
      if (ierr /= 0) then
         write(stderr, '(a,i0,a,i0)') 'Error in lschfb for state ', ii, ', ierr = ', ierr
         stop 1
      end if
      if (abs(e_in_out - ea(ii)) > eps) then
         write(stderr, '(a,i0)') 'Warning: energy mismatch in lschfb for state ', ii
         stop 1
      end if
      rmatch(ii) = rr(imatch(ii))
   end do

   ! Write output
   table = toml_table()
   call set_value(table, "it", it)

   call add_array(table, "rr", array)
   call set_value(array, "rr", rr)

   call add_array(table, "rhoc", array)
   call set_value(array, "rhoc", rhoc)

   call add_array(table, "rhov", array)
   call set_value(array, "rhov", rhov)

   call add_array(table, "rhot", array)
   call set_value(array, "rhot", rhot)

   call add_array(table, "vxc", array)
   call set_value(array, "vxc", vxc)

   call add_array(table, "vfull", array)
   call set_value(array, "vfull", vfull)

   call add_array(table, "etot", array)
   call set_value(array, "etot", etot)

   call add_array(table, "eeel", array)
   call set_value(array, "eeel", eeel)

   call add_array(table, "eexc", array)
   call set_value(array, "eexc", eexc)

   call add_array(table, "na", array)
   call set_value(array, "na", na)

   call add_array(table, "la", array)
   call set_value(array, "la", la)

   call add_array(table, "fa", array)
   call set_value(array, "fa", fa)

   call add_array(table, "ea", array)
   call set_value(array, "ea", ea)

   call add_array(table, "rpk", array)
   call set_value(array, "rpk", rpk)

   call add_array(table, "imatch", array)
   call set_value(array, "imatch", imatch)

   call add_array(table, "rmatch", array)
   call set_value(array, "rmatch", rmatch)

   call add_array(table, "uuaea", array)
   call set_value(array, "uuaea", uuaea)

   call add_array(table, "upaea", array)
   call set_value(array, "upaea", upaea)

   call toml_dump(table, stdout)

   deallocate(rr)
   deallocate(rhoc, rhov, rhot, vxc, vha, vfull)
   deallocate(uuaea, upaea)
   nullify(array)
   destroy(table)
end program oncvatom
