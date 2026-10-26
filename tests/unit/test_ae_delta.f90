! Unit tests of the relativistic term Delta = u'' - [l(l+1)/r**2 + 2(V-E)]u
! of an all-electron wave (ae_delta; D-4 and R-11 of doc/smoothness_plan.md),
! written before the code, on the program's default mesh (amesh=1.006,
! rr1=5e-5, rmax=45), hydrogen-like Z=50 (where the relativistic terms are
! large), the states solved by boundstate:
!   nr      Delta of the non-relativistic 2p is 0 exactly
!   sr      Delta of the scalar-relativistic 2p at the mesh points between
!           0.02 and 0.2 bohr against u'' from lm_derivs on uu (an 8th-order
!           stencil) less the Schroedinger term, to 1e-8 of max |u''|; and
!           Delta is not small: above 1e-3 of max |u''| somewhere there
!   fr      the same for the Dirac 2p3/2 large component
! Stops with a non-zero code if any check failed.
 program test_ae_delta
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use logmesh_m, only: lm_derivs
 use boundstate_m, only: boundstate
 implicit none
 external :: ae_delta

!Local variables
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: zz = 50.0_dp
 real(dp), allocatable :: rr(:), vv(:), uu(:, :), up(:, :), delta(:)
 real(dp) :: al, ee, dd(3), u2, err, big, scale
 integer :: nfail = 0, mmax, ii, irel, ierr, mch, i1, i2, kap
 character(len=2), parameter :: names(0:2) = ['nr', 'sr', 'fr']

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), vv(mmax), uu(mmax, 2), up(mmax, 2), delta(mmax))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do
 vv(:) = -zz/rr
 i1 = 1 + int(log(0.02_dp/rr1)/al)
 i2 = 1 + int(log(0.2_dp/rr1)/al)

 do irel = 0, 2
   kap = -2
   ee = -0.5_dp*zz**2/4.0_dp
   uu(:, :) = 0.0_dp
   up(:, :) = 0.0_dp
   mch = 0
   call boundstate(irel, 2, 1, kap, ierr, ee, rr, zz, vv, uu, up, mmax, mch)
   if(ierr /= 0) then
     write(stderr, '(a,a,a,i0)') 'FAILED: ', names(irel), ' boundstate ierr ', &
&      ierr
     nfail = nfail + 1
     cycle
   end if
   call ae_delta(irel, 1, kap, ee, rr, vv, uu(:, 1), up(:, 1), mmax, i1, i2, &
&                delta)
   err = 0.0_dp
   big = 0.0_dp
   scale = 0.0_dp
   do ii = i1, i2
     call lm_derivs(rr, uu(:, 1), mmax, rr(ii), 2, dd)
     u2 = dd(3) - (2.0_dp/rr(ii)**2 + 2.0_dp*(vv(ii) - ee))*uu(ii, 1)
     err = max(err, abs(delta(ii) - u2))
     big = max(big, abs(delta(ii)))
     scale = max(scale, abs(dd(3)))
   end do
   if(irel == 0) then
     call check(names(irel)//' Delta zero', big, 0.0_dp)
   else
     call check(names(irel)//' Delta vs stencil', err/scale, 1.0e-8_dp)
     if(.not. (big > 1.0e-3_dp*scale)) then
       write(stderr, '(a,a,a,es9.2)') 'FAILED: ', names(irel), &
&        ' Delta too small to test, relative ', big/scale
       nfail = nfail + 1
     end if
   end if
 end do

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'ae_delta: all checks passed'

 contains

 subroutine check(label, err, tol)

!! err against tol, printed and counted

 character(len=*), intent(in) :: label
 real(dp), intent(in) :: err, tol

 write(stderr, '(a,a,es9.2)') label, ': error ', err
 if(.not. (err <= tol)) then
   write(stderr, '(a,a,es10.3,a,es9.2)') 'FAILED: ', label, err, ' > ', tol
   nfail = nfail + 1
 end if
 end subroutine check

 end program test_ae_delta
