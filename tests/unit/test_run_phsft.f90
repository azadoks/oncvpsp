! Unit tests of the rxpsh check in run_phsft, non- and fully-relativistic,
! called directly as from the library: a radius inside the core radii must
! give a warning for each l, and one at or beyond them none. The programs
! never get the warning, because check_data keeps rxpsh beyond every rc. The
! atom is hydrogen-like (zz=1) with no projectors. Stops with a non-zero code
! at the first failure.
 program test_run_phsft
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use constants_m, only: mxl, mxlin, mxkap
 use pieces_m, only: piecewise_t, pw_alloc
 use run_phsft_m, only: run_phsft
 use relativity_m, only: c_light => c_light_default
 implicit none
 integer, parameter :: iout=6

!Local variables
 integer, parameter :: mmax = 950, mxprj = 5, lmax = 0, lloc = 0
 real(dp), parameter :: zz = 1.0_dp, rr1 = 0.0005_dp, amesh = 1.012_dp
 real(dp), parameter :: epsh2 = 0.5_dp, depsh = 0.1_dp
! the energies from epsh2 down to -0.5
 integer, parameter :: npsh = 10

 real(dp) :: rr(mmax), vfull(mmax), vp(mmax, mxl+2), vp_r(mmax, mxl+2, mxkap)
 real(dp) :: epa(mxprj, mxlin+1), epa_r(mxprj, mxlin+1, mxkap)
 real(dp) :: evkb(mxprj, mxl+1), evkb_r(mxprj, mxl+1, mxkap)
 real(dp) :: rphs(4), pshae(npsh, 4, mxkap), pshps(npsh, 4, mxkap)
 type(piecewise_t) :: vlocpw, vkbpw(mxprj, mxl+1), vkbpw_r(mxprj, mxl+1, mxkap)
 integer :: nproj(mxlin+1), ii, nfail = 0
 real(dp) :: rc(mxlin+1)

 do ii = 1, mmax
   rr(ii) = rr1*amesh**(ii - 1)
 end do
 vfull(:) = -zz/rr(:)
 vp(:, :) = spread(vfull, 2, 5)
 vp_r(:, :, 1) = vp
 vp_r(:, :, 2) = vp
 epa = 0.0_dp; epa_r = 0.0_dp
 evkb = 0.0_dp; evkb_r = 0.0_dp
 nproj(:) = 0
! the local potential is one piece over the whole mesh
 call pw_alloc(vlocpw, mmax, 1, [0.0_dp, rr(mmax)])
 vlocpw%vals(:, 1) = vfull
! every core radius at 2 bohr
 rc(:) = 2.0_dp

! l=0 is the local potential's l, and l=1-3 are beyond lmax: 1 + 3 warnings
! inside rc
 call check_case('run_phsft', 1.0_dp, 4)
 call check_case('run_phsft', rc(1), 0)
 call check_case('run_phsft', 3.0_dp, 0)
 call check_case('run_phsft_r', 1.0_dp, 4)
 call check_case('run_phsft_r', rc(1), 0)
 call check_case('run_phsft_r', 3.0_dp, 0)

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'run_phsft: all checks passed'

 contains

 subroutine check_case(name, rxpsh, nexpected)

!! runs name with rxpsh, its output going to a file, and counts its warnings
!!
!! @param name  run_phsft for the non-relativistic case, run_phsft_r for the
!!              fully-relativistic one, whose messages it names
!! @param rxpsh  radius for the log-derivative analysis (bohr)
!! @param nexpected  number of warnings it should give

!Input variables
 character(len=*), intent(in) :: name
 real(dp), intent(in) :: rxpsh
 integer, intent(in) :: nexpected

!Local variables
 integer :: ifail
 character(len=256) :: line
 integer :: unit, ios, nwarn

 open(unit=6, file=name//'.out', status='replace', action='write')
 if(name == 'run_phsft') then
   call run_phsft(lmax, lloc, nproj, epa, npsh, epsh2, depsh, rxpsh, &
&                 vlocpw, vkbpw, evkb, rr, vfull, vp, zz, mmax, mxprj, rc, &
&                 0, c_light, 1, rphs, pshae, pshps, iout, ifail)
   if(ifail /= 0) error stop 1
 else
   call run_phsft(lmax, lloc, nproj, epa_r, npsh, epsh2, depsh, rxpsh, &
&                 vlocpw, vkbpw_r, evkb_r, rr, vfull, vp_r, zz, mmax, mxprj, &
&                 rc, 2, c_light, 2, rphs, pshae, pshps, iout, ifail)
   if(ifail /= 0) error stop 1
 end if
 close(6)

 nwarn = 0
 open(newunit=unit, file=name//'.out', status='old', action='read')
 do
   read(unit, '(a)', iostat=ios) line
   if(ios /= 0) exit
   if(index(line, name//': WARNING rxpsh inside the') == 1) nwarn = nwarn + 1
 end do
 close(unit)
 if(nwarn /= nexpected) then
   write(stderr, '(3a,f8.4,2(a,i0))') 'FAILED: ', name, ' rxpsh=', rxpsh, &
&     ': warnings=', nwarn, ', expected ', nexpected
   nfail = nfail + 1
 end if
 return
 end subroutine check_case

 end program test_run_phsft
