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
!> solvers.
 module quad_m
 use constants_m, only: dp, pi
 use logmesh_m, only: lm_integral
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: gauss_legendre, log_norm, log_norm_to

 contains

 subroutine gauss_legendre(xa,xb,nn,xx,ww)

!! the nn-point Gauss-Legendre nodes and weights on [xa, xb], the roots of
!! P_nn by Newton's method
!!
!! @param xa  lower end
!! @param xb  upper end
!! @param nn  number of nodes
!! @param xx  nodes
!! @param ww  weights

!Input variables
 real(dp), intent(in) :: xa,xb
 integer, intent(in) :: nn

!Output variables
 real(dp), intent(out) :: xx(nn),ww(nn)

!Local variables
 real(dp) :: zz,z1,p1,p2,p3,pp
 integer :: ii,jj,it
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
       write(error_unit,'(/2a,i0,a,i0)') 'quad_m: ERROR gauss_legendre', &
&            ' node not converged, nn=',nn,', node=',ii
       error stop 1
     end if
   end do
   xx(ii)=0.5d0*(xa+xb) - 0.5d0*(xb-xa)*zz
   xx(nn+1-ii)=0.5d0*(xa+xb) + 0.5d0*(xb-xa)*zz
   ww(ii)=(xb-xa)/((1.0d0 - zz**2)*pp**2)
   ww(nn+1-ii)=ww(ii)
 end do
 return
 end subroutine gauss_legendre

 real(dp) function log_norm(rr,ff,nn,ll,gam)

!! the integral of ff from 0 to rr(nn), ff going as r**(2 gam) at the
!! origin: ff(1) rr(1)/(2 gam+1) below rr(1), the Floater-Hormann
!! integral of lm_integral above it, from the first nn points alone (no
!! value beyond rr(nn) is read), as fpovlp's overlaps
!!
!! @param rr  log radial mesh
!! @param ff  integrand: u**2, or g**2+f**2 of a Dirac solution
!! @param nn  the last point
!! @param ll  angular momentum (lm_integral's power at the origin, which
!!            the limit rr(1) leaves out)
!! @param gam  power of the solution at the origin, u ~ r**gam

!Input variables
 integer, intent(in) :: nn,ll
 real(dp), intent(in) :: rr(nn),ff(nn),gam

!Local variables
 real(dp) :: si

 call lm_integral(rr,ff,nn,2*ll+2,rr(1),rr(nn),si)
 log_norm=ff(1)*rr(1)/(2.0d0*gam+1.0d0)+si
 end function log_norm

 real(dp) function log_norm_to(rr,ff,nn,ll,gam,rq)

!! log_norm to the radius rq instead of rr(nn), a smooth function of rq
!!
!! @param rr  log radial mesh
!! @param ff  integrand: u**2, or g**2+f**2 of a Dirac solution
!! @param nn  the last point read, far enough past rq for the interpolant
!!            of lm_integral
!! @param ll  angular momentum (lm_integral's power at the origin, which
!!            the limit rr(1) leaves out)
!! @param gam  power of the solution at the origin, u ~ r**gam
!! @param rq  the upper limit, rr(1)<=rq<=rr(nn) (bohr)

!Input variables
 integer, intent(in) :: nn,ll
 real(dp), intent(in) :: rr(nn),ff(nn),gam,rq

!Local variables
 real(dp) :: si

 call lm_integral(rr,ff,nn,2*ll+2,rr(1),rq,si)
 log_norm_to=ff(1)*rr(1)/(2.0d0*gam+1.0d0)+si
 end function log_norm_to

 end module quad_m
