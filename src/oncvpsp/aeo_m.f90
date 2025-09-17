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
! Adams extrapolation and interpolation formulas for
! outward and inward integration, Abramowitz and
! Stegun, p. 896
module aeo_m
    use constants_m, only: dp
    implicit none
    private
    public :: aeo, aio, aei, aii
contains
function aeo(xx, ii) result(yy)
    implicit none
    real(dp), dimension(:) :: xx
    real(dp) :: yy
    integer :: ii

    yy = (4.16666666667d-2) * (55.0d0 * xx(ii) - 59.0d0 * xx(ii - 1) + 37.0d0 * xx(ii - 2) - 9.0d0 * xx(ii - 3))
    return
end function aeo

function aio(xx, ii) result(yy)
    implicit none
    real(dp), dimension(:) :: xx
    real(dp) :: yy
    integer :: ii

    yy = (4.16666666667d-2) * (9.0d0 * xx(ii + 1) + 19.0d0 * xx(ii) - 5.0d0 * xx(ii - 1) + xx(ii - 2))
    return
end function aio

function aei(xx, ii) result(yy)
    implicit none
    real(dp), dimension(:) :: xx
    real(dp) :: yy
    integer :: ii

    yy = -(4.16666666667d-2) * (55.0d0 * xx(ii) - 59.0d0 * xx(ii + 1) + 37.0d0 * xx(ii + 2) - 9.0d0 * xx(ii + 3))
    return
end function aei

function aii(xx, ii) result(yy)
    implicit none
    real(dp), dimension(:) :: xx
    real(dp) :: yy
    integer :: ii

    yy = -(4.16666666667d-2) * (9.0d0 * xx(ii - 1) + 19.0d0 * xx(ii) - 5.0d0 * xx(ii + 1) + xx(ii + 2))
    return
end function aii
end module aeo_m
