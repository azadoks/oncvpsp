! Unit tests of the Adams-Moulton steps of order k = 4 to 8 (adams_m):
!   coefs    adams_init(4) is the default adams_t; the step of order k
!            integrates the polynomials of degree up to k exactly,
!            sum_i cm(i) d (1-i)**(d-1) = 1 for d = 1..k, to 1e-14 of the
!            sum of the terms' magnitudes
!   orders   the non-relativistic 2p and 3d of the hydrogen-like zz = 3
!            (lschfb, exact energies -zz**2/(2n**2)) on meshes of ratio
!            1.024 and 1.012 (rr1 = 5e-5, rmax = 45): the error of each
!            order falls by about 2**k from the first to the second, by at
!            least 2**(k-1) (unless it is at round-off, 1e-12), and at
!            1.012 the error of order 8 is below that of order 4 by 1e3
! Stops with a non-zero code if any check failed.
 program test_adams
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check, check_done, nfail
 use adams_m, only: adams_t, adams_init, adams_min, adams_max
 use lschfb_m, only: lschfb
 implicit none
 integer, parameter :: iout=6

!Local variables
 real(dp), parameter :: zz = 3.0_dp, rr1 = 5.0e-5_dp, rmax = 45.0_dp
 real(dp), parameter :: ameshes(2) = [1.024_dp, 1.012_dp]
 integer, parameter :: nstate = 2, nns(nstate) = [2, 3], lls(nstate) = [1, 2]
 type(adams_t) :: adams, adams4
 real(dp) :: err(adams_min:adams_max, 2, nstate), ss, sa, rate
 integer :: kk, id, ii, im, is

! coefs
 adams4 = adams_init(4)
 call check('adams_init(4) order', real(adams4%order - adams%order, dp), &
&           0.0_dp)
 call check('adams_init(4) coefficients', maxval(abs(adams4%cm - adams%cm)), &
&           0.0_dp)
 do kk = adams_min, adams_max
   adams = adams_init(kk)
   do id = 1, kk
     ss = 0.0_dp
     sa = 0.0_dp
     do ii = 0, kk - 1
       ss = ss + adams%cm(ii)*id*real(1 - ii, dp)**(id - 1)
       sa = sa + abs(adams%cm(ii)*id*real(1 - ii, dp)**(id - 1))
     end do
     call check('exact on polynomials', abs(ss - 1.0_dp)/sa, 1.0e-14_dp)
   end do
 end do

! orders
 do kk = adams_min, adams_max
   adams = adams_init(kk)
   do im = 1, 2
     do is = 1, nstate
       err(kk, im, is) = energy_error(ameshes(im), nns(is), lls(is))
     end do
   end do
 end do
 do is = 1, nstate
   write(iout, '(a,i0,a,i0)') 'n=', nns(is), ' l=', lls(is)
   do kk = adams_min, adams_max
     write(iout, '(a,i0,2es12.3)') '  order ', kk, err(kk, :, is)
     if(err(kk, 2, is) > 1.0e-12_dp) then
       rate = err(kk, 1, is)/err(kk, 2, is)
       call check('error falls by 2**(k-1)', &
&                 max(0.0_dp, 2.0_dp**(kk - 1) - rate), 0.0_dp)
     end if
   end do
   call check('order 8 against order 4', &
&             err(8, 2, is)/err(4, 2, is), 1.0e-3_dp)
 end do

 call check_done('test_adams')

 contains

 function energy_error(amesh, nn, ll) result(de)

!! |E - E_exact| of the state nn, ll on the mesh of ratio amesh, with the
!! steps of adams

 real(dp), intent(in) :: amesh
 integer, intent(in) :: nn, ll
 real(dp) :: de
 real(dp), allocatable :: rr(:), vv(:), uu(:), up(:)
 real(dp) :: al, ee
 integer :: mmax, jj, ierr, mch, ifail

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), vv(mmax), uu(mmax), up(mmax))
 do jj = 1, mmax
   rr(jj) = rr1*exp(al*(jj - 1))
 end do
 vv(:) = -zz/rr
 ee = -0.5_dp*zz**2/nn**2
 mch = 0
 call lschfb(0, 137.0_dp, nn, ll, ierr, ee, rr, zz, vv, uu, up, mmax, &
&            adams, mch, iout, ifail)
 if(ifail /= 0 .or. ierr /= 0) then
   write(stderr, '(a,i0,a,i0)') 'FAILED: lschfb ierr ', ierr, ' ifail ', &
&    ifail
   nfail = nfail + 1
 end if
 de = abs(ee + 0.5_dp*zz**2/nn**2)
 end function energy_error

 end program test_adams
