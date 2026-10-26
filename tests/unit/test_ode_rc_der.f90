! Unit tests of ode_rc_der, the derivatives of an all-electron wave function
! at any radius from the radial equation (phase 1 of
! doc/smoothness_plan.md), on hydrogen (Z=1) on the program's default mesh
! (amesh=1.006, rr1=5e-5, rmax=45), at radii off the mesh, from uu and
! d(uu)/dn as the solvers give them:
!   nr 1s, 2p      against the analytic u/r = 2 exp(-r), r exp(-r/2)
!   fr 1s          the Dirac ground state, whose large component is
!                  r**gamma exp(-r), gamma=sqrt(1-1/c**2), at E=c**2(gamma-1)
!   nr 1s solved   the state of boundstate against the analytic one: the
!                  derivatives carry the solver's own error in uu (1e-10),
!                  where the stencils of wf_rc_der reach 8e-7 in the 4th
!   sr 1s          the Koelling-Harmon state of boundstate, against the
!                  stencils of wf_rc_der at a mesh point (whose error sets
!                  the tolerance)
! The tolerance of the k-th derivative is tol(k) times the largest of the
! first five derivatives. Stops with a non-zero code if any check failed.
 program test_ode_rc_der
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use relativity_m, only: c_light
 use ode_rc_der_m, only: ode_rc_der
 use wf_rc_der_m, only: wf_rc_der
 implicit none
 external :: boundstate

!Local variables
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: radii(3) = [0.7_dp, 1.2345_dp, 2.718281828459045_dp]
 real(dp), parameter :: tol(0:4) = [1.0e-13_dp, 1.0e-12_dp, 1.0e-10_dp, &
&                                   1.0e-9_dp, 1.0e-8_dp]
 real(dp), allocatable :: rr(:), vv(:), uu(:, :), up(:, :)
 real(dp) :: al, uorder(5), exact(5), gam, ee, rc
 integer :: nfail = 0, mmax, ii, ierr, mch, irc

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), vv(mmax), uu(mmax, 2), up(mmax, 2))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do
 vv(:) = -1.0_dp/rr

! nr 1s: u/r = 2 exp(-r); 2p: u/r = r exp(-r/2)
 uu(:, 1) = 2.0_dp*rr*exp(-rr)
 up(:, 1) = al*rr*2.0_dp*(1.0_dp - rr)*exp(-rr)
 do ii = 1, 3
   rc = radii(ii)
   call ode_rc_der(0, 0, -1, -0.5_dp, rr, vv, uu(:, 1), up(:, 1), mmax, &
&                  rc, uorder)
   exact = 2.0_dp*exp(-rc)*[1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp]
   call compare('nr 1s', uorder, exact)
 end do
 uu(:, 1) = rr**2*exp(-0.5_dp*rr)
 up(:, 1) = al*rr*(2.0_dp*rr - 0.5_dp*rr**2)*exp(-0.5_dp*rr)
 do ii = 1, 3
   rc = radii(ii)
   call ode_rc_der(0, 1, -2, -0.125_dp, rr, vv, uu(:, 1), up(:, 1), mmax, &
&                  rc, uorder)
   exact = exp(-0.5_dp*rc)*[rc, 1.0_dp - 0.5_dp*rc, -1.0_dp + 0.25_dp*rc, &
&                            0.75_dp - 0.125_dp*rc, -0.5_dp + 0.0625_dp*rc]
   call compare('nr 2p', uorder, exact)
 end do

! fr 1s: g/r = r**(gamma-1) exp(-r)
 gam = sqrt(1.0_dp - 1.0_dp/c_light**2)
 ee = c_light**2*(gam - 1.0_dp)
 uu(:, 1) = rr**gam*exp(-rr)
 up(:, 1) = al*rr*(gam/rr - 1.0_dp)*uu(:, 1)
 do ii = 1, 3
   rc = radii(ii)
   call ode_rc_der(2, 0, -1, ee, rr, vv, uu(:, 1), up(:, 1), mmax, rc, uorder)
   call power_exp(gam - 1.0_dp, rc, exact)
   call compare('fr 1s', uorder, exact)
 end do

! nr 1s by the solver, against the analytic state
 call solve(0)
 if(ierr == 0) then
   do ii = 1, 3
     rc = radii(ii)
     call ode_rc_der(0, 0, -1, ee, rr, vv, uu(:, 1), up(:, 1), mmax, rc, &
&                    uorder)
     exact = 2.0_dp*exp(-rc)*[1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp]
     call compare_tol('nr 1s solved', uorder, exact, [(1.0e-9_dp, ii = 1, 5)])
   end do
 end if

! sr 1s by the solver, against the stencils at a mesh point
 call solve(1)
 if(ierr == 0) then
   irc = 0
   rc = 1.2345_dp
   call wf_rc_der(rr, uu(:, 1), al, rc, irc, mmax, exact)
   call ode_rc_der(1, 0, -1, ee, rr, vv, uu(:, 1), up(:, 1), mmax, rc, uorder)
   call compare_tol('sr 1s', uorder, exact, &
&                   [1.0e-13_dp, 1.0e-9_dp, 1.0e-8_dp, 1.0e-7_dp, 1.0e-5_dp])
 end if

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'ode_rc_der: all checks passed'

 contains

 subroutine solve(irel)

!! the hydrogen 1s state of boundstate, positive at large r, in uu and up;
!! counts a failure if there is none
!!
!! @param irel  0 non-, 1 scalar-relativistic

 integer, intent(in) :: irel

 ee = -0.5_dp
 uu(:, :) = 0.0_dp
 up(:, :) = 0.0_dp
 mch = 0
 call boundstate(irel, 1, 0, -1, ierr, ee, rr, 1.0_dp, vv, uu, up, mmax, mch)
 if(ierr /= 0) then
   write(stderr, '(a,i0,a,i0)') 'FAILED: irel ', irel, &
&    ' 1s: boundstate ierr ', ierr
   nfail = nfail + 1
 end if
 end subroutine solve

 subroutine compare(name, got, want)

!! counts a failure for each derivative off by more than tol
!!
!! @param name  what was checked
!! @param got  ode_rc_der's value and derivatives
!! @param want  the exact ones

 character(len=*), intent(in) :: name
 real(dp), intent(in) :: got(5), want(5)

 call compare_tol(name, got, want, tol)
 end subroutine compare

 subroutine compare_tol(name, got, want, tols)

!! counts a failure for each derivative k off by more than tols(k) times
!! the largest of want
!!
!! @param name  what was checked
!! @param got  ode_rc_der's value and derivatives
!! @param want  the reference
!! @param tols  relative tolerances of the orders 0..4

 character(len=*), intent(in) :: name
 real(dp), intent(in) :: got(5), want(5), tols(5)
 real(dp) :: scale, err
 integer :: kk

 scale = maxval(abs(want))
 do kk = 1, 5
   err = abs(got(kk) - want(kk))/scale
   if(.not. (err <= tols(kk))) then
     write(stderr, '(3a,i0,a,es10.3,a,es10.3,a,f7.4)') 'FAILED: ', name, &
&      ' order ', kk - 1, ': error ', err, ' > ', tols(kk), ' at r=', rc
     nfail = nfail + 1
   end if
 end do
 end subroutine compare_tol

 subroutine power_exp(pp, rq, dd)

!! value and four derivatives of r**p exp(-r) at rq, by Leibniz
!!
!! @param pp  the power p
!! @param rq  radius
!! @param dd  d^k/dr^k in dd(k+1)

 real(dp), intent(in) :: pp, rq
 real(dp), intent(out) :: dd(5)
 real(dp) :: fall(0:4), binom
 integer :: kk, jj

! fall(j) = d^j r**p / r**p = p(p-1)..(p-j+1) / r**j
 fall(0) = 1.0_dp
 do jj = 1, 4
   fall(jj) = fall(jj - 1)*(pp - jj + 1)/rq
 end do
 do kk = 0, 4
   dd(kk + 1) = 0.0_dp
   binom = 1.0_dp
   do jj = 0, kk
     dd(kk + 1) = dd(kk + 1) + binom*fall(jj)*(-1.0_dp)**(kk - jj)
     binom = binom*(kk - jj)/(jj + 1)
   end do
   dd(kk + 1) = dd(kk + 1)*rq**pp*exp(-rq)
 end do
 end subroutine power_exp

 end program test_ode_rc_der
