!
! Copyright (c) 2026 by Austin Zadoks
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!> Gauss-Legendre quadrature, for the integrals of the model core's budget
!> (modcore_budget) and of the confined atom's basis (ws_basis_m); and the
!> norm of a radial function on the log mesh of the bound and scattering
!> solvers; and the order-4 sum of the all-electron energies.
 module quad_m
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: quad_gauss_legendre, quad_log_norm, quad_log_norm_to
 public :: quad_gregory_sum

 contains

 subroutine quad_gauss_legendre(xa,xb,nn,xx,ww,ifail)

!! the nn-point Gauss-Legendre nodes and weights on [xa, xb], the roots of
!! P_nn by Newton's method
!!
!! @param xa  lower end
!! @param xb  upper end
!! @param nn  number of nodes
!! @param xx  nodes
!! @param ww  weights
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp, pi
 use failure_m, only: fail_convergence, failed, msg_gauss_legendre_node

!Input variables
 real(dp), intent(in) :: xa,xb
 integer, intent(in) :: nn

!Output variables
 real(dp), intent(out) :: xx(nn),ww(nn)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: zz,z1,p1,p2,p3,pp
 integer :: ii,jj,it
 ifail=0
 do ii=1,(nn+1)/2
   zz=cos(pi*(ii-0.25d0)/(nn+0.5d0))
   do it=1,100
     p1=1.0d0
     p2=0.0d0
     do jj=1,nn
       p3=p2
       p2=p1
       p1=((2*jj-1)*zz*p2 - (jj-1)*p3)/jj
     end do
     pp=nn*(zz*p1 - p2)/(zz**2 - 1.0d0)
     z1=zz
     zz=z1 - p1/pp
     if(abs(zz-z1)<1.0d-15) exit
     if(it==100) then
       write(error_unit,'(/a,i0,a,i0)') &
&            'quad_gauss_legendre: ERROR node not converged, nn=',nn, &
&            ', node=',ii
       ifail=failed(fail_convergence,msg_gauss_legendre_node)
       return
     end if
   end do
   xx(ii)=0.5d0*(xa+xb) - 0.5d0*(xb-xa)*zz
   xx(nn+1-ii)=0.5d0*(xa+xb) + 0.5d0*(xb-xa)*zz
   ww(ii)=(xb-xa)/((1.0d0 - zz**2)*pp**2)
   ww(nn+1-ii)=ww(ii)
 end do
 return
 end subroutine quad_gauss_legendre

 real(dp) function quad_log_norm(rr,ff,nn,ll,gam)

!! the integral of ff from 0 to rr(nn), ff going as r**(2 gam) at the
!! origin: ff(1) rr(1)/(2 gam+1) below rr(1), the Floater-Hormann
!! integral of lm_integral above it, from the first nn points alone (no
!! value beyond rr(nn) is read), as fpovlp's overlaps
!!
!! @param rr  log radial mesh (bohr)
!! @param ff  integrand: u**2, or g**2+f**2 of a Dirac solution
!! @param nn  the last point
!! @param ll  angular momentum (lm_integral's power at the origin, which
!!            the limit rr(1) leaves out)
!! @param gam  power of the solution at the origin, u ~ r**gam

 use constants_m, only: dp
 use logmesh_m, only: lm_integral

!Input variables
 integer, intent(in) :: nn,ll
 real(dp), intent(in) :: rr(nn),ff(nn),gam

!Local variables
 real(dp) :: si

 call lm_integral(rr,ff,nn,2*ll+2,rr(1),rr(nn),si)
 quad_log_norm=ff(1)*rr(1)/(2.0d0*gam+1.0d0)+si
 end function quad_log_norm

 real(dp) function quad_log_norm_to(rr,ff,nn,ll,gam,rq)

!! quad_log_norm to the radius rq instead of rr(nn), a smooth function of rq
!!
!! @param rr  log radial mesh (bohr)
!! @param ff  integrand: u**2, or g**2+f**2 of a Dirac solution
!! @param nn  the last point read, far enough past rq for the interpolant
!!            of lm_integral
!! @param ll  angular momentum (lm_integral's power at the origin, which
!!            the limit rr(1) leaves out)
!! @param gam  power of the solution at the origin, u ~ r**gam
!! @param rq  the upper limit, rr(1)<=rq<=rr(nn) (bohr)

 use constants_m, only: dp
 use logmesh_m, only: lm_integral

!Input variables
 integer, intent(in) :: nn,ll
 real(dp), intent(in) :: rr(nn),ff(nn),gam,rq

!Local variables
 real(dp) :: si

 call lm_integral(rr,ff,nn,2*ll+2,rr(1),rq,si)
 quad_log_norm_to=ff(1)*rr(1)/(2.0d0*gam+1.0d0)+si
 end function quad_log_norm_to

 real(dp) function quad_gregory_sum(ff,nn)

!! the sum of ff on the nodes of a uniform grid, the first three weighted
!! (9,28,23)/24 (Gregory's start): the integral from the first node over
!! the step, to O(h**4). The all-electron energies (vout, ae_energies)
!! still integrate by it (doc/numerics_audit.md, 6c item 2)
!!
!! @param ff  values on the nodes
!! @param nn  number of nodes, at least 3

 use constants_m, only: dp

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: ff(nn)

!Local variables
 integer :: ii

 quad_gregory_sum=(9.0d0*ff(1) + 28.0d0*ff(2) + 23.0d0*ff(3))/24.0d0
 do ii=4,nn
   quad_gregory_sum=quad_gregory_sum + ff(ii)
 end do
 return
 end function quad_gregory_sum

 end module quad_m
