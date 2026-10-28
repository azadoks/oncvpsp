! Unit tests of the outward integration of the radial equation with
! piecewise sources and potential (outwf_m; R-10 of doc/smoothness_plan.md),
! written before the code, on the program's default mesh (amesh=1.006,
! rr1=5e-5, rmax=45), l=1, E=-0.8 Ha:
!   potential  -2/r on [0, bv], -2/r + 0.3 (r-bv)**2 beyond (C1 at bv)
!   sources    beta_1 = r**2 e**-r on [0, b1], plus 0.5 (r-b1)**2 r on
!              (b1, b2], 0 beyond (C1 at b1, a jump at b2); beta_2 =
!              r**2 (1+r) e**-2r on [0, b1], less (r-b1)**2 on (b1, b2]
!   with bv = 1.9, b1 = 1.3, b2 = 2.0, each piece held on the whole mesh:
!   fine     the log derivatives at r = 3 of the three solutions
!            (homogeneous, one per source) and the six overlaps int phi_i
!            beta_j dr over the solution's value at r = 3 (the amplitude
!            of each solution carries the error of its power-series start
!            at the first nodes, as in lschps) against the same integration
!            on the mesh with the step al/8 (the method is 4th order, its
!            truncation 2e-7 on the log derivatives at the default step,
!            falling 15.6, 15.4, 14.7 times on three halvings, so that is
!            the exact answer to 1e-10): to 3e-7
!   smooth   b1, then b2, then bv stepped over one mesh step in 41
!            values, the rest fixed: the same quantities fit a polynomial
!            of degree 6 to 1e-12 (the exact
!            answers are smooth in the bounds, and so must be the method's
!            error)
!   no kink  with every piece the extension of the first, the breakpoints
!            b1 and bv change nothing (b2 still ends the sources): moving
!            them moves the solutions by < 1e-10
!   count    the homogeneous solution of -2/r alone at E = -0.35 (between
!            the 2p and the 3p) has one node up to r = 20, counted across
!            the segments
! Stops with a non-zero code if any check failed.
 program test_outwf
 use constants_m, only: dp
 use unit_check_m, only: check, check_done
 use pieces_m, only: piecewise_t, pw_alloc, pw_free, pw_node_values
 use outwf_m, only: outwf_t, outwf_setup, outwf_integrate, outwf_overlap, &
&                   outwf_nodes, outwf_count, outwf_free
 use lapack_m, only: dgels
 implicit none
 integer, parameter :: iout=6

!Local variables
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: ee = -0.8_dp, r3 = 3.0_dp
 integer, parameter :: ll = 1, nfine = 8
 real(dp) :: bv = 1.9_dp, b1 = 1.3_dp, b2 = 2.0_dp
 real(dp) :: al, vals(3), ovl(3, 2), valf(3), ovlf(3, 2), err
 real(dp) :: ts(41), ys(41, 9), resid(41), base(9), scale
 integer :: mmax, i3, ii, jj, kk, iscan, node
 character(len=2), parameter :: names(3) = ['b1', 'b2', 'bv']

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 i3 = 1 + nint(log(r3/rr1)/al)

! fine
 call solve(amesh, mmax, i3, bv, b1, b2, .true., vals, ovl)
 call solve(amesh**(1.0_dp/nfine), nfine*(mmax - 1) + 1, nfine*(i3 - 1) + 1, &
&           bv, b1, b2, .true., valf, ovlf)
 err = 0.0_dp
 do ii = 1, 3
   err = max(err, abs(vals(ii) - valf(ii))/abs(valf(ii)))
   do jj = 1, 2
     err = max(err, abs(ovl(ii, jj) - ovlf(ii, jj))/abs(ovlf(ii, jj)))
   end do
 end do
 call check('fine mesh', err, 3.0e-7_dp)

! smooth
 do iscan = 1, 3
   do ii = 1, 41
     ts(ii) = 2.0_dp*(ii - 1)/40.0_dp - 1.0_dp
     scale = amesh**(0.5_dp*(ts(ii) + 1.0_dp))
     select case(iscan)
     case(1)
       call solve(amesh, mmax, i3, bv, b1*scale, b2, .true., vals, ovl)
     case(2)
       call solve(amesh, mmax, i3, bv, b1, b2*scale, .true., vals, ovl)
     case default
       call solve(amesh, mmax, i3, bv*scale, b1, b2, .true., vals, ovl)
     end select
     ys(ii, 1:3) = vals
     ys(ii, 4:9) = reshape(ovl, [6])
   end do
   err = 0.0_dp
   do kk = 1, 9
     call fit_residual(ts, ys(:, kk), 6, resid)
     err = max(err, maxval(abs(resid))/maxval(abs(ys(:, kk))))
   end do
   call check('smooth in '//names(iscan), err, 1.0e-12_dp)
 end do

! no kink
 call solve(amesh, mmax, i3, bv, b1, b2, .false., vals, ovl)
 base(1:3) = vals
 base(4:9) = reshape(ovl, [6])
 call solve(amesh, mmax, i3, bv*amesh**0.37_dp, b1*amesh**0.61_dp, b2, &
&           .false., vals, ovl)
 err = 0.0_dp
 do kk = 1, 3
   err = max(err, abs(vals(kk) - base(kk))/abs(base(kk)))
 end do
 err = max(err, maxval(abs(reshape(ovl, [6]) - base(4:9))/abs(base(4:9))))
 call check('no kink', err, 1.0e-10_dp)

! count
 call count_nodes(amesh, mmax, node)
 call check('node count', abs(dble(node - 1)), 0.0_dp)

 call check_done('outwf')

 contains

 subroutine solve(am, mm, i3q, bvq, b1q, b2q, kinks, uval, ovq)

!! the log derivatives of the three solutions at the node i3q (r3) and
!! their overlaps over their values there, on the mesh of ratio am and mm
!! points, with the bounds given; without kinks every piece is the first
!! one

 real(dp), intent(in) :: am, bvq, b1q, b2q
 integer, intent(in) :: mm, i3q
 logical, intent(in) :: kinks
 real(dp), intent(out) :: uval(3), ovq(3, 2)
 real(dp), allocatable :: rr(:), vv(:), uu(:), up(:)
 type(piecewise_t) :: vpw, spw(2)
 type(outwf_t) :: ow
 real(dp) :: alq, bnd(0:2)
 integer :: isol, jsrc, im
 integer :: ifail

 alq = log(am)
 allocate(rr(mm), vv(mm), uu(mm), up(mm))
 do im = 1, mm
   rr(im) = rr1*exp(alq*(im - 1))
 end do

 bnd = [0.0_dp, bvq, rr(mm)]
 call pw_alloc(vpw, mm, 2, bnd)
 vpw%vals(:, 1) = -2.0_dp/rr
 vpw%vals(:, 2) = -2.0_dp/rr
 if(kinks) vpw%vals(:, 2) = vpw%vals(:, 2) + 0.3_dp*(rr - bvq)**2
 call pw_node_values(rr, mm, vpw, vv)

 bnd = [0.0_dp, b1q, b2q]
 call pw_alloc(spw(1), mm, 2, bnd)
 call pw_alloc(spw(2), mm, 2, bnd)
 spw(1)%vals(:, 1) = rr**2*exp(-rr)
 spw(1)%vals(:, 2) = spw(1)%vals(:, 1)
 spw(2)%vals(:, 1) = rr**2*(1.0_dp + rr)*exp(-2.0_dp*rr)
 spw(2)%vals(:, 2) = spw(2)%vals(:, 1)
 if(kinks) then
   spw(1)%vals(:, 2) = spw(1)%vals(:, 2) + 0.5_dp*(rr - b1q)**2*rr
   spw(2)%vals(:, 2) = spw(2)%vals(:, 2) - (rr - b1q)**2
 end if

 call outwf_setup(ll, rr, mm, vv, vpw, 2, spw, ow, iout, ifail)
 if(ifail /= 0) error stop 1
 call outwf_integrate(ow, ee, i3q)
 do isol = 1, 3
   call outwf_nodes(ow, isol, i3q, uu, up)
   uval(isol) = up(i3q)/(alq*uu(i3q))
   do jsrc = 1, 2
     call outwf_overlap(ow, isol, jsrc, ovq(isol, jsrc))
     ovq(isol, jsrc) = ovq(isol, jsrc)/uu(i3q)
   end do
 end do
 call outwf_free(ow)
 call pw_free(vpw)
 call pw_free(spw(1))
 call pw_free(spw(2))
 deallocate(rr, vv, uu, up)
 end subroutine solve

 subroutine count_nodes(am, mm, nodeq)

!! nodes of the homogeneous solution of -2/r at E = -0.35 up to r = 20,
!! the potential split at bv

 real(dp), intent(in) :: am
 integer, intent(in) :: mm
 integer, intent(out) :: nodeq
 real(dp), allocatable :: rr(:), vv(:)
 type(piecewise_t) :: vpw, spw(1)
 type(outwf_t) :: ow
 real(dp) :: alq, bnd(0:2)
 integer :: i20, im
 integer :: ifail

 alq = log(am)
 allocate(rr(mm), vv(mm))
 do im = 1, mm
   rr(im) = rr1*exp(alq*(im - 1))
 end do
 i20 = 1 + nint(log(20.0_dp/rr1)/alq)
 bnd = [0.0_dp, bv, rr(mm)]
 call pw_alloc(vpw, mm, 2, bnd)
 vpw%vals(:, 1) = -2.0_dp/rr
 vpw%vals(:, 2) = -2.0_dp/rr
 call pw_node_values(rr, mm, vpw, vv)
 call outwf_setup(ll, rr, mm, vv, vpw, 0, spw, ow, iout, ifail)
 if(ifail /= 0) error stop 1
 call outwf_integrate(ow, -0.35_dp, i20)
 nodeq = outwf_count(ow, 1, 0.05_dp, i20)
 call outwf_free(ow)
 call pw_free(vpw)
 deallocate(rr, vv)
 end subroutine count_nodes

 subroutine fit_residual(tt, yy, ndeg, res)

!! residual of the least-squares polynomial of degree ndeg, by dgels

 real(dp), intent(in) :: tt(41), yy(41)
 integer, intent(in) :: ndeg
 real(dp), intent(out) :: res(41)
 real(dp) :: aa(41, ndeg + 1), bb(41), work(2048)
 integer :: j1, k1, info

 do j1 = 1, 41
   do k1 = 0, ndeg
     aa(j1, k1 + 1) = tt(j1)**k1
   end do
 end do
 bb(:) = yy(:)
 call dgels('N', 41, ndeg + 1, 1, aa, 41, bb, 41, work, size(work), info)
 do j1 = 1, 41
   res(j1) = yy(j1)
   do k1 = 0, ndeg
     res(j1) = res(j1) - bb(k1 + 1)*tt(j1)**k1
   end do
 end do
 end subroutine fit_residual

 end program test_outwf
