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
!> Kind and physical and numerical constants shared by all routines

 module constants_m
 use, intrinsic :: iso_fortran_env, only: real64
 implicit none
 private
 public :: dp, pi, c_light, mxprj, dqout

!> kind of all real variables
 integer, parameter :: dp=real64

!> pi
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!> speed of light in atomic units (Hartree), 1/alpha.  The 4.0.1 value is
!> kept, rather than CODATA 137.035999..., so results are unchanged
 real(dp), parameter :: c_light=137.036_dp

!> maximum number of projectors for each l
 integer, parameter :: mxprj=5

!> spacing (bohr**-1) of the q mesh of the "Energy error per electron"
!> tables, as their interpolation assumes it (eresid moves the mesh points
!> onto its own q mesh, so the actual spacing differs slightly)
 real(dp), parameter :: dqout=0.5000001_dp

 end module constants_m
