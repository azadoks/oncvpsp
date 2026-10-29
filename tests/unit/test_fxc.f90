! Unit tests of the analytic second derivatives of the exchange-correlation
! functionals (fxc_dispatch), against 4th-order differences of the potential
! (exc_dispatch) in the direction of a change of the density, rho + t drho
! with t = +-1e-3, +-2e-3, on a model atom (rho of a 1s, 2s, 2p and 3p
! shell of charges 2, 2, 6, 2, drho of the 3p shell):
!   lda      iexc 1, 2, 3, 6, 7 (built in), at the points where rho > 1e-10
!            (PZ, iexc 3, away from rs = 1, where its potential steps):
!            the largest difference below 1e-8 of the largest |dV|
!   libxc    an LDA (-001012) to the same 1e-8, and a GGA (-101130), the
!            derivatives on the mesh included, when Libxc has them. In the
!            core a GGA's potential is a difference of large terms, whose
!            round-off the differences raise to 2e-7 of the largest |dV|,
!            so for it the matrix element der2exc takes, the integral of
!            drho dV r**2, is checked, to 1e-8 relative
!   none     iexc 4 (the built-in PBE) has none: have is false
! Stops with a non-zero code if any check failed.
 program test_fxc
 use constants_m, only: dp, pi
 use unit_check_m, only: check, check_done
 use xc_dispatch_m, only: exc_dispatch, fxc_dispatch
 use exc_libxc_m, only: have_libxc
 implicit none
 integer, parameter :: iout=6

!Local variables
 integer, parameter :: mmax=1100
 real(dp), parameter :: c_light=137.035999084_dp
 real(dp), parameter :: rr1=1.0e-4_dp, tt=1.0e-3_dp
 real(dp) :: al,rr(mmax),rho(mmax),drho(mmax),dv(mmax),dfd(mmax)
 real(dp) :: vp1(mmax),vm1(mmax),vp2(mmax),vm2(mmax),exc(mmax),rs(mmax)
 logical :: have,mask(mmax)
 integer :: ii,kk
 integer, parameter :: nlda=5, iexcs(nlda)=[1, 2, 3, 6, 7]
 integer, parameter :: ilxc(2)=[-1012, -101130]

 al=log(1.012_dp)
 do ii=1,mmax
   rr(ii)=rr1*exp(al*(ii - 1))
 end do
! 4 pi n of hydrogen-like shells of exponents 14, 5 and 1.5, to 45 bohr
 rho(:)=2.0_dp*shell(14.0_dp, 0) + 8.0_dp*shell(5.0_dp, 1) &
&      + 2.0_dp*shell(1.5_dp, 2)
 drho(:)=shell(1.5_dp, 2)
 rs(:)=(3.0_dp/max(rho(:), 1.0e-30_dp))**(1.0_dp/3.0_dp)

! lda
 do kk=1,nlda
   mask(:)=rho(:)>1.0e-10_dp
   if(iexcs(kk)==3) mask(:)=mask(:) .and. abs(rs(:) - 1.0_dp)>1.0e-2_dp
   call compare(iexcs(kk))
 end do

! libxc
 if(have_libxc) then
   do kk=1,2
     mask(:)=rho(:)>1.0e-10_dp
     call compare(ilxc(kk))
   end do
 end if

! none
 call fxc_dispatch(4,c_light,al,rho,drho,dv,rr,mmax,have)
 call check('iexc 4 has no analytic fxc', .not. have)

 call check_done('test_fxc')

 contains

 function shell(zeta, nr) result(sh)

!! the normalized density 4 pi n of a hydrogen-like s shell of exponent
!! zeta with the power r**nr of its radial function

 real(dp), intent(in) :: zeta
 integer, intent(in) :: nr
 real(dp) :: sh(mmax)
 real(dp) :: norm

! integral of r**(2 nr + 2) exp(-2 zeta r) from 0 to infinity
 norm=gamma(real(2*nr + 3, dp))/(2.0_dp*zeta)**(2*nr + 3)
 sh(:)=rr(:)**(2*nr)*exp(-2.0_dp*zeta*rr(:))/norm
 end function shell

 subroutine compare(iexc)

!! fxc_dispatch against the differences of exc_dispatch for iexc

 integer, intent(in) :: iexc
 real(dp) :: err,scale,d2an,d2fd

 call fxc_dispatch(iexc,c_light,al,rho,drho,dv,rr,mmax,have)
 if(iexc<0 .and. .not. have) then
   write(iout,'(a,i0,a)') 'iexc ',iexc,': Libxc without fxc, skipped'
   return
 end if
 call check('iexc has an analytic fxc', have)
 call exc_dispatch(iexc,c_light,al,rho + tt*drho,vp1,exc,rr,mmax,'test')
 call exc_dispatch(iexc,c_light,al,rho - tt*drho,vm1,exc,rr,mmax,'test')
 call exc_dispatch(iexc,c_light,al,rho + 2.0_dp*tt*drho,vp2,exc,rr,mmax, &
&                  'test')
 call exc_dispatch(iexc,c_light,al,rho - 2.0_dp*tt*drho,vm2,exc,rr,mmax, &
&                  'test')
 dfd(:)=(8.0_dp*(vp1(:) - vm1(:)) - (vp2(:) - vm2(:)))/(12.0_dp*tt)
 scale=maxval(abs(dv), mask=mask)
 err=maxval(abs(dv - dfd), mask=mask)/scale
 write(iout,'(a,i0,a,es10.2,a,es10.2)') 'iexc ',iexc,': max |dV| ', &
&     scale,', relative difference ',err
! the matrix element, by the trapezoidal sum in ln r (the same for both)
 d2an=al*sum(drho*dv*rr**3, mask=mask)
 d2fd=al*sum(drho*dfd*rr**3, mask=mask)
 write(iout,'(a,es22.14,a,es10.2)') '  d2exc ',d2an,', relative difference ', &
&     abs(d2an - d2fd)/abs(d2an)
 call check('d2exc against differences', abs(d2an - d2fd)/abs(d2an), &
&           1.0e-8_dp)
 if(iexc /= -101130) call check('fxc against differences', err, 1.0e-8_dp)
 end subroutine compare

 end program test_fxc
