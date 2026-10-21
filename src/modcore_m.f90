!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
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
!> The pieces of the model core charges at continuous radii
!> (doc/smoothness_plan.md, phase 5): the Teter function and its
!> derivatives with no cutoff (M-5), the crossover radius where the core
!> charge equals fcfact times the pseudo valence charge (M-1), the
!> polynomial of icmod 1 with its constant in closed form (M-3), and the
!> scale of the Teter function of icmod 2 (M-4).
!> Stubs: written with tests/unit/test_modcore.f90, before the code.

 module modcore_m
 use constants_m, only: dp
 implicit none
 private
 public :: teter_derivs, core_crossover, poly_core_fit, teter_scale

 contains

 subroutine teter_derivs(xx,nder,gd)

!! the Teter function g(x) = (sin(2 pi x)/(2 pi x (1-4x**2)(1-x**2)))**2
!! and its derivatives at x >= 0, for every x (no cutoff)
!!
!! @param xx  argument, x >= 0
!! @param nder  highest derivative wanted, 0<=nder<=4
!! @param gd  d^k g/dx^k in gd(k), k=0..nder

!Input variables
 real(dp), intent(in) :: xx
 integer, intent(in) :: nder

!Output variables
 real(dp), intent(out) :: gd(0:nder)

 gd(:)=huge(1.0d0)
 if(xx<0.0d0) gd(:)=-huge(1.0d0)
 end subroutine teter_derivs

 subroutine core_crossover(rr,rhoc,rhotps,mmax,fcfact,rstar,ierr)

!! the crossover radius: the outermost radius where the core charge equals
!! fcfact times the pseudo valence charge, the root of
!! ln(rhoc) - ln(fcfact rhotps) on the log mesh, converged to round-off
!!
!! @param rr  log radial mesh
!! @param rhoc  core charge
!! @param rhotps  pseudo valence charge
!! @param mmax  size of the log radial mesh
!! @param fcfact  factor of the pseudo valence charge
!! @param rstar  the crossover radius (bohr)
!! @param ierr  0, or 1 if the core charge never exceeds fcfact rhotps

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rhoc(mmax),rhotps(mmax),fcfact

!Output variables
 real(dp), intent(out) :: rstar
 integer, intent(out) :: ierr

 rstar=huge(1.0d0)+0.0d0*(rr(1)+rhoc(1)+rhotps(1)+fcfact)
 ierr=0
 end subroutine core_crossover

 subroutine poly_core_fit(rstar,fmatch,a0,aco)

!! the polynomial model core of icmod 1, p(r) = a0 + sum_j aco(j) r**(j+2),
!! j=1..5, that matches fmatch at rstar and is barely monotonic on
!! [0, rstar]: a0 the least for which p' <= 0 there, in closed form
!!
!! @param rstar  matching radius (bohr)
!! @param fmatch  the core charge and its first 4 derivatives at rstar
!! @param a0  constant term
!! @param aco  coefficients of r**3..r**7

!Input variables
 real(dp), intent(in) :: rstar,fmatch(5)

!Output variables
 real(dp), intent(out) :: a0,aco(5)

 a0=huge(1.0d0)+0.0d0*(rstar+fmatch(1))
 aco(:)=huge(1.0d0)
 end subroutine poly_core_fit

 subroutine teter_scale(ymatch,xx,ierr)

!! the argument x of the Teter function, 0 < x < 1.5, at which its
!! logarithmic slope x g'(x)/g(x) is ymatch (< 0); monotone in x there
!!
!! @param ymatch  r rho'/rho of the core charge at the matching radius
!! @param xx  the argument
!! @param ierr  0, or 1 if ymatch is not in the range of x g'/g on (0, 1.5)

!Input variables
 real(dp), intent(in) :: ymatch

!Output variables
 real(dp), intent(out) :: xx
 integer, intent(out) :: ierr

 xx=huge(1.0d0)+0.0d0*ymatch
 ierr=0
 end subroutine teter_scale

 end module modcore_m
