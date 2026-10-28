! Unit tests of the operations on a uniform grid (ugrid_m; Q-1 to Q-4 of
! doc/smoothness_plan.md), written before the code:
!   weights    the composite 8th-order weights integrate x**7 on [0,1] with
!              20 points exactly (1e-15), and exp(-x) on [0,4] with 33 and 65
!              points with an error ratio above 2**7 (8th order), the finer
!              below 1e-13
!   cumint     the cumulative integral from the top of exp(-x) on [0,4]
!              (65 points) equals exp(-x_i)-exp(-4) to 5e-13 at every point
!              (O(h**8); 1.4e-13 here, the shifted end stencils' error)
!   compensated  a sum of about 1 to which 2000 cells of 1e-16 are added
!              keeps them to 1e-2 (each would round away alone): the
!              compensation survives the compiler's flags (fast math)
!   interp     the interpolant of exp(-x) (81 points on [0,4]) at the
!              midpoints and at 0.3 of a cell is within 1e-11 (its error is
!              O(h**7), 6e-12 here), and at a node exact
!   root       the root of the interpolant of cos(x) (81 points on [0,4]) is
!              pi/2 to 1e-11
!   cross      the largest x where exp(-x) (81 points) falls to 0.1 is
!              -ln(0.1) to 1e-11; a level below the curve's end, 1e-3, is
!              not found
! Stops with a non-zero code if any check failed.
 program test_ugrid
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use ugrid_m, only: ug_weights, ug_cumint_down, ug_interp, ug_root, &
&                   ug_cross_down
 implicit none

!Local variables
 real(dp), parameter :: pi = 3.14159265358979323846_dp
 real(dp), allocatable :: xx(:), ff(:), ww(:), cum(:)
 real(dp) :: hh, ss, err1, err2, fq, tq, xq
 integer :: nfail = 0, nn, ii, ierr
 logical :: found

! weights
 call grid(20, 0.0_dp, 1.0_dp)
 ff(:) = xx**7
 call ug_weights(nn, hh, ww)
 ss = sum(ww*ff)
 call check('weights x**7', abs(ss - 0.125_dp), 1.0e-15_dp)
 call grid(33, 0.0_dp, 4.0_dp)
 ff(:) = exp(-xx)
 call ug_weights(nn, hh, ww)
 err1 = abs(sum(ww*ff) - (1.0_dp - exp(-4.0_dp)))
 call grid(65, 0.0_dp, 4.0_dp)
 ff(:) = exp(-xx)
 call ug_weights(nn, hh, ww)
 err2 = abs(sum(ww*ff) - (1.0_dp - exp(-4.0_dp)))
 write(stderr, '(a,2es10.2,a,f6.1)') 'weights exp(-x): errors', err1, err2, &
&  ' ratio', err1/max(err2, tiny(1.0_dp))
 call check('weights exp(-x) 65 points', err2, 1.0e-13_dp)
 if(.not. (err1/max(err2, tiny(1.0_dp)) >= 128.0_dp)) then
   write(stderr, '(a)') 'FAILED: weights: not 8th order'
   nfail = nfail + 1
 end if

! cumint
 call ug_cumint_down(ff, nn, hh, cum)
 err1 = 0.0_dp
 do ii = 1, nn
   err1 = max(err1, abs(cum(ii) - (exp(-xx(ii)) - exp(-4.0_dp))))
 end do
 call check('cumint exp(-x)', err1, 5.0e-13_dp)

! compensated: 1 in the last cells, then cells of 1e-16 each
 call grid(2101, 0.0_dp, 1.0_dp)
 ff(:) = 1.0e-16_dp/hh
 ff(nn) = 1.0_dp/hh
 call ug_cumint_down(ff, nn, hh, cum)
 call check('cumint compensated', abs((cum(1) - cum(nn - 100)) &
&           /(1.0e-16_dp*dble(nn - 101)) - 1.0_dp), 1.0e-2_dp)

! interp
 call grid(81, 0.0_dp, 4.0_dp)
 ff(:) = exp(-xx)
 err1 = 0.0_dp
 do ii = 1, nn - 1
   call ug_interp(ff, nn, dble(ii - 1) + 0.5_dp, fq)
   err1 = max(err1, abs(fq - exp(-(xx(ii) + 0.5_dp*hh))))
   call ug_interp(ff, nn, dble(ii - 1) + 0.3_dp, fq)
   err1 = max(err1, abs(fq - exp(-(xx(ii) + 0.3_dp*hh))))
 end do
 call check('interp exp(-x) between nodes', err1, 1.0e-11_dp)
 call ug_interp(ff, nn, 17.0_dp, fq)
 call check('interp at a node', abs(fq - ff(18)), 0.0_dp)

! root
 ff(:) = cos(xx)
 call ug_root(ff, nn, 20.0_dp, 40.0_dp, tq, ierr)
 call check('root of cos: ierr', dble(ierr), 0.0_dp)
 call check('root of cos', abs(tq*hh - 0.5_dp*pi), 1.0e-11_dp)

! cross
 ff(:) = exp(-xx)
 call ug_cross_down(ff, nn, hh, 0.0_dp, 0.1_dp, xq, found)
 if(.not. found) then
   write(stderr, '(a)') 'FAILED: cross: not found'
   nfail = nfail + 1
 else
   call check('cross of exp(-x) at 0.1', abs(xq + log(0.1_dp)), 1.0e-11_dp)
 end if
 call ug_cross_down(ff, nn, hh, 0.0_dp, 1.0e-3_dp, xq, found)
 if(found) then
   write(stderr, '(a)') 'FAILED: cross: a level below the curve was found'
   nfail = nfail + 1
 end if

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'ugrid: all checks passed'

 contains

 subroutine grid(npts, xa, xb)

!! the uniform grid of npts points on [xa,xb] in xx, and the work arrays

 integer, intent(in) :: npts
 real(dp), intent(in) :: xa, xb

 nn = npts
 hh = (xb - xa)/dble(nn - 1)
 if(allocated(xx)) deallocate(xx, ff, ww, cum)
 allocate(xx(nn), ff(nn), ww(nn), cum(nn))
 do ii = 1, nn
   xx(ii) = xa + dble(ii - 1)*hh
 end do
 end subroutine grid

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

 end program test_ugrid
