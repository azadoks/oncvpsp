! Unit tests of the scalar-relativistic and spin-orbit projectors from the
! j = l +/- 1/2 ones (sr_so_r; numerics audit C-14), on the program's
! default mesh, for l = 1 with two projectors per j, beta_1 = r**2 e**-r
! and beta_2 = r**3 e**-r/4 on [0, 2.5]:
!   equal   the j = l - 1/2 projectors and coefficients the same as the
!           j = l + 1/2 ones (the non-relativistic limit): the overlap of
!           the four is singular, its two null directions at round-off
!           (sqrt of a negative before the guard). The scalar kernel
!           sum_m esr_m |vsr_m><vsr_m| is sum_i e_i |beta_i><beta_i| and
!           the spin-orbit one 0, to 1e-13 of the largest kernel value
!   near    the j = l - 1/2 projectors beta_i (1 + 1e-5 r**2), and their
!           coefficients 1 + 1e-5 times as large: the overlap's smallest
!           eigenvalue about 1e-12 of its largest, above the floor, so
!           nothing is dropped and both kernels are the definitions'
!           sum_j fscl_j sum_i e_ij |beta_ij><beta_ij| and the same with
!           fso_j, to 1e-12
! Stops with a non-zero code if any check failed.
 program test_sr_so_r
 use constants_m, only: dp, mxl, mxlin, mxkap, mxprj
 use unit_check_m, only: check, check_done
 use pieces_m, only: piecewise_t, pw_alloc
 use sr_so_r_m, only: sr_so_r
 use failure_m, only: warnings_t
 use, intrinsic :: iso_fortran_env, only: output_unit
 implicit none

!Local variables
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: rc = 2.5_dp, delta = 1.0e-5_dp
 integer, parameter :: lmax = 1, l1 = 2, np = 2
 real(dp), allocatable :: rr(:), beta(:, :, :)
 real(dp), allocatable :: vsr(:, :, :), vso(:, :, :)
 real(dp) :: esr(2*mxprj, mxl + 1), eso(2*mxprj, mxl + 1)
 real(dp) :: evkb(mxprj, mxl + 1, mxkap), bound(0:1), ee(np)
 real(dp) :: fscl(mxkap), fso(mxkap)
 type(piecewise_t) :: vkbpw(mxprj, mxl + 1, mxkap)
 type(piecewise_t) :: vsrpw(2*mxprj, mxl + 1), vsopw(2*mxprj, mxl + 1)
 type(warnings_t) :: warns
 integer :: nproj(mxlin + 1), mmax, ii, ik, ifail
 real(dp) :: al

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), beta(mmax, np, mxkap))
 allocate(vsr(mmax, 2*mxprj, mxl + 1), vso(mmax, 2*mxprj, mxl + 1))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do
 nproj(:) = 0
 nproj(l1) = np
 ee = [-3.0_dp, 0.7_dp]
 fscl = [2.0_dp/3.0_dp, 1.0_dp/3.0_dp]
 fso = [2.0_dp/3.0_dp, -2.0_dp/3.0_dp]

! equal
 do ik = 1, mxkap
   beta(:, 1, ik) = rr**2*exp(-rr)
   beta(:, 2, ik) = 0.25_dp*rr**3*exp(-rr)
   evkb(1:np, l1, ik) = ee
 end do
 call run_case('equal', 1.0e-13_dp)

! near
 beta(:, 1, 2) = beta(:, 1, 1)*(1.0_dp + delta*rr**2)
 beta(:, 2, 2) = beta(:, 2, 1)*(1.0_dp + delta*rr**2)
 evkb(1:np, l1, 2) = ee*(1.0_dp + delta)
 call run_case('near', 1.0e-12_dp)

 call check_done('test_sr_so_r')

 contains

 subroutine run_case(label, tol)

!! sr_so_r on the projectors beta and coefficients evkb, and its kernels
!! against the definitions' on the nodes inside rc

 character(len=*), intent(in) :: label
 real(dp), intent(in) :: tol
 real(dp) :: ksr, kso, xsr, xso, scale, esr_err, eso_err
 integer :: ia, ib, ip, mm, nrc

 bound(0) = 0.0_dp
 bound(1) = rc
 do ik = 1, mxkap
   do ip = 1, np
     call pw_alloc(vkbpw(ip, l1, ik), mmax, 1, bound, ifail)
     vkbpw(ip, l1, ik)%vals(:, 1) = beta(:, ip, ik)
   end do
 end do
 vsr(:, :, :) = 0.0_dp
 vso(:, :, :) = 0.0_dp
 esr(:, :) = 0.0_dp
 eso(:, :) = 0.0_dp
 call sr_so_r(lmax, nproj, rr, mmax, evkb, vkbpw, vsr, esr, vso, eso, &
&             vsrpw, vsopw, warns, output_unit, ifail)
 call check(label//': ifail', real(ifail, dp), 0.0_dp)

 nrc = count(rr < rc)
 scale = 0.0_dp
 esr_err = 0.0_dp
 eso_err = 0.0_dp
 do ia = 1, nrc, 25
   do ib = 1, nrc, 25
     xsr = 0.0_dp
     xso = 0.0_dp
     do ik = 1, mxkap
       do ip = 1, np
         xsr = xsr + fscl(ik)*evkb(ip, l1, ik)*beta(ia, ip, ik)*beta(ib, ip, ik)
         xso = xso + fso(ik)*evkb(ip, l1, ik)*beta(ia, ip, ik)*beta(ib, ip, ik)
       end do
     end do
     ksr = 0.0_dp
     kso = 0.0_dp
     do mm = 1, 2*np
       ksr = ksr + esr(mm, l1)*vsr(ia, mm, l1)*vsr(ib, mm, l1)
       kso = kso + eso(mm, l1)*vso(ia, mm, l1)*vso(ib, mm, l1)
     end do
     scale = max(scale, abs(xsr))
     esr_err = max(esr_err, abs(ksr - xsr))
     eso_err = max(eso_err, abs(kso - xso))
   end do
 end do
 call check(label//': scalar kernel', esr_err/scale, tol)
 call check(label//': spin-orbit kernel', eso_err/scale, tol)
 end subroutine run_case

 end program test_sr_so_r
