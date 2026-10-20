! Unit tests of the rxpsh check in run_phsft, non- and fully-relativistic,
! called directly as from the library: a radius inside the core radii must
! give a warning for each l, and one at or beyond them none. The programs
! never get the warning, because check_data keeps rxpsh beyond every rc. The
! atom is hydrogen-like (zz=1) with no projectors. Stops with a non-zero code
! at the first failure.
 program test_run_phsft
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use constants_m, only: mxl, mxlin, mxkap
 implicit none
 external :: run_phsft

!Local variables
 integer, parameter :: mmax = 950, mxprj = 5, lmax = 0, lloc = 0
 real(dp), parameter :: zz = 1.0_dp, rr1 = 0.0005_dp, amesh = 1.012_dp
 real(dp), parameter :: epsh1 = -0.5_dp, epsh2 = 0.5_dp, depsh = 0.1_dp

 real(dp) :: rr(mmax), vfull(mmax), vp(mmax, mxl+2), vp_r(mmax, mxl+2, mxkap)
 real(dp) :: epa(mxprj, mxlin+1), epa_r(mxprj, mxlin+1, mxkap)
 real(dp) :: vkb(mmax, mxprj, mxl+1), vkb_r(mmax, mxprj, mxl+1, mxkap)
 real(dp) :: evkb(mxprj, mxl+1), evkb_r(mxprj, mxl+1, mxkap)
 integer :: nproj(mxlin+1), irc(mxlin+1), ii, nfail = 0

 do ii = 1, mmax
   rr(ii) = rr1*amesh**(ii - 1)
 end do
 vfull(:) = -zz/rr(:)
 vp(:, :) = spread(vfull, 2, 5)
 vp_r(:, :, 1) = vp
 vp_r(:, :, 2) = vp
 epa = 0.0_dp; epa_r = 0.0_dp
 vkb = 0.0_dp; vkb_r = 0.0_dp
 evkb = 0.0_dp; evkb_r = 0.0_dp
 nproj(:) = 0
! every core radius at the first mesh point beyond 2 bohr, as the programs
! place it
 irc(:) = findloc(rr > 2.0_dp, .true., dim=1)

! l=0 is the local potential's l, and l=1-3 are beyond lmax: 1 + 3 warnings
! inside rc
 call check_case('run_phsft', 1.0_dp, 4)
 call check_case('run_phsft', rr(irc(1)), 0)
 call check_case('run_phsft', 3.0_dp, 0)
 call check_case('run_phsft_r', 1.0_dp, 4)
 call check_case('run_phsft_r', rr(irc(1)), 0)
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
 character(len=256) :: line
 integer :: unit, ios, nwarn

 open(unit=6, file=name//'.out', status='replace', action='write')
 if(name == 'run_phsft') then
   call run_phsft(lmax, lloc, nproj, epa, epsh1, epsh2, depsh, rxpsh, &
&                 vkb, evkb, rr, vfull, vp, zz, mmax, mxprj, irc, 0, 1)
 else
   call run_phsft(lmax, lloc, nproj, epa_r, epsh1, epsh2, depsh, rxpsh, &
&                 vkb_r, evkb_r, rr, vfull, vp_r, zz, mmax, mxprj, irc, 2, 2)
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
