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
!> Values, derivatives, integrals and roots of a function given on the log
!> radial mesh, at any radius, not only at mesh points: the operations that
!> let a radius (a core radius, a matching or crossover radius) vary
!> continuously instead of being rounded to the mesh. Every result is a
!> smooth (C-infinity, to round-off) function of the radii it is evaluated
!> at, so a quantity built from them is as smooth in its inputs as the
!> mathematics allows. The method is in doc/smoothness_plan.md: the
!> Floater-Hormann rational interpolant of blending degree lm_blend in
!> x = ln(r), which is uniform on the mesh, of the values (lm_interp), of
!> derivative fields made with fixed stencils (lm_derivs) and of the
!> cumulative integral (lm_cumint, lm_integral).
!>
!> STUB: the interfaces only, for the tests written before the code
!> (tests/unit/test_logmesh.f90). Every routine returns huge() (lm_root
!> ierr=-1) until it is implemented.

 module logmesh_m
 use constants_m, only: dp
 implicit none
 private
 public :: lm_blend, lm_interp, lm_derivs, lm_cumint, lm_integral, lm_root

!> blending degree d of the Floater-Hormann interpolant: it reproduces
!> polynomials of degree d in ln(r), and its error is O(al**(d+1))
 integer, parameter :: lm_blend=6

 contains

 subroutine lm_interp(rr,ff,mmax,rq,fq)

!! value at the radius rq of the function given on the log mesh, from its
!! Floater-Hormann interpolant in ln(r)
!!
!! @param rr  log radial mesh
!! @param ff  function on the mesh
!! @param mmax  size of the log radial mesh
!! @param rq  radius, rr(1)<=rq<=rr(mmax) (bohr)
!! @param fq  value at rq

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),ff(mmax),rq

!Output variables
 real(dp), intent(out) :: fq

 fq=huge(1.0_dp)
 return
 end subroutine lm_interp

 subroutine lm_derivs(rr,ff,mmax,rq,nder,dd)

!! value and radial derivatives up to order nder at the radius rq of a smooth
!! function given on the log mesh: the derivative fields in ln(r) on the
!! mesh, made with fixed centred stencils, converted to d/dr and interpolated
!! to rq, so that the result is smooth in rq
!!
!! @param rr  log radial mesh
!! @param ff  function on the mesh
!! @param mmax  size of the log radial mesh
!! @param rq  radius, away from the ends of the mesh by the stencils' reach
!!            (bohr)
!! @param nder  highest derivative wanted, 0<=nder<=4
!! @param dd  d^k ff/dr^k at rq in dd(k+1), k=0..nder

!Input variables
 integer, intent(in) :: mmax,nder
 real(dp), intent(in) :: rr(mmax),ff(mmax),rq

!Output variables
 real(dp), intent(out) :: dd(nder+1)

 dd(:)=huge(1.0_dp)
 return
 end subroutine lm_derivs

 subroutine lm_cumint(rr,ff,mmax,mm,cum)

!! cumulative integral of a function given on the log mesh, from 0 to each
!! mesh point, the part below rr(1) from its power law r**mm
!!
!! @param rr  log radial mesh
!! @param ff  integrand on the mesh
!! @param mmax  size of the log radial mesh
!! @param mm  power of r with which ff goes to 0 at the origin
!! @param cum  int_0^rr(i) ff dr in cum(i)

!Input variables
 integer, intent(in) :: mmax,mm
 real(dp), intent(in) :: rr(mmax),ff(mmax)

!Output variables
 real(dp), intent(out) :: cum(mmax)

 cum(:)=huge(1.0_dp)
 return
 end subroutine lm_cumint

 subroutine lm_integral(rr,ff,mmax,mm,ra,rb,ss)

!! integral of a function given on the log mesh between any two radii: the
!! interpolated cumulative integral at rb less that at ra, so smooth in both
!!
!! @param rr  log radial mesh
!! @param ff  integrand on the mesh
!! @param mmax  size of the log radial mesh
!! @param mm  power of r with which ff goes to 0 at the origin
!! @param ra,rb  limits, 0<=ra, rb<=rr(mmax) (bohr)
!! @param ss  int_ra^rb ff dr

!Input variables
 integer, intent(in) :: mmax,mm
 real(dp), intent(in) :: rr(mmax),ff(mmax),ra,rb

!Output variables
 real(dp), intent(out) :: ss

 ss=huge(1.0_dp)
 return
 end subroutine lm_integral

 subroutine lm_root(rr,ff,mmax,ra,rb,rq,ierr)

!! root, between ra and rb, of the interpolant (lm_interp) of a function
!! given on the log mesh, converged to round-off
!!
!! @param rr  log radial mesh
!! @param ff  function on the mesh
!! @param mmax  size of the log radial mesh
!! @param ra,rb  bracket (bohr); the interpolant must change sign between
!!               them
!! @param rq  the root (bohr)
!! @param ierr  0, or 1 if the interpolant has the same sign at ra and rb

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),ff(mmax),ra,rb

!Output variables
 real(dp), intent(out) :: rq
 integer, intent(out) :: ierr

 rq=huge(1.0_dp)
 ierr=-1
 return
 end subroutine lm_root

 end module logmesh_m
