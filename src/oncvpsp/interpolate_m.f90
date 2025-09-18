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
module interpolate_m
   use iso_fortran_env, only: stdout => output_unit
   use precision_m, only: dp
   implicit none
   private
   public :: interpolate

   !> Polynomial order for interpolation
   integer, parameter :: N_POLY = 7
contains

!> Local polynomial interpolation with special treatment of first point
!> Input grid x_in must be ordered in ascending order
subroutine interpolate(x_in, y_in, n_in, x_out, y_out, n_out)
   implicit none

   ! Input variables
   !> Number of points on input grid
   integer, intent(in) :: n_in
   !> Function values on input grid
   real(dp), intent(in) :: x_in(n_in)
   !> Input grid points
   real(dp), intent(in) :: y_in(n_in)
   !> Number of points on output grid
   integer, intent(in) :: n_out
   !> Output grid points
   real(dp), intent(in) :: x_out(n_out)

   ! Output variables
   !> Interpolated function values on output grid
   real(dp), intent(out) :: y_out(n_out)

   ! Local variables
   real(dp) :: sum
   real(dp) :: term
   real(dp) :: zz
   !> Loop index for input grid
   integer :: i
   !> Loop index for output grid
   integer :: j
   !> Loop index for interval-halving search
   integer :: k
   !> Index of lower bracketing point in input grid
   integer :: imin
   !> Index of upper bracketing point in input grid
   integer :: imax
   integer :: iprod
   integer :: iy
   integer :: istart
   !> Error message
   character(len=1024) :: error_msg

   if (n_in < N_POLY + 1) then
      write (error_msg, '(/a,i6,a,i4)') 'interpolate: ERROR - n_in=', n_in, '< N_POLY=', N_POLY
      error stop error_msg
   end if
   do i = 2, n_in
      if (x_in(i) <= x_in(i - 1)) then
         write (error_msg, '(/a)') 'interpolate: ERROR - x_in not in ascending order'
         error stop error_msg
      end if
   end do


   ! Output point 1 is skipped here because of special properties of pp data
   ! (i.e., x_in is a log grid).
   ! This point receives special treatment (see below).
   y_out(:) = 0.0_dp
   imin = 1
   do j = 2, n_out
      if (x_out(j) < x_in(1)) then
         write (error_msg, '(/a,i6)') 'interpolate: ERROR - x_out(', j, ') out of range'
         error stop error_msg
      end if
      if (x_out(j) > x_in(n_in)) then
         write (error_msg, '(/a,i6)') 'interpolate: ERROR - x_out(', j, ') out of range'
         error stop error_msg
      end if

      ! Interval-halving search for x_in(i) points bracketing x_out(j)
      imin = 1
      imax = n_in
      do k = 1, n_in
         i = (imin + imax) / 2
         if (x_out(j) > x_in(i)) then
            imin = i
         else
            imax = i
         end if
         if (imax - imin == 1) then
            exit
         end if
      end do

      zz = x_out(j)
      istart = min(imin - int(dble(N_POLY) / dble(2)), n_in - N_POLY)
      istart = max(istart, 1)
      sum = 0.0_dp
      do iy = istart, istart + N_POLY
         if (abs(y_in(iy)) < 1.0e-12_dp) cycle
         term = y_in(iy)
         do iprod = istart, istart + N_POLY
            if (iprod == iy) cycle
            term = term * (zz - x_in(iprod)) / (x_in(iy) - x_in(iprod))
         end do
         sum = sum + term
      end do
      y_out(j) = sum

      ! Special treatment for the origin
      ! Do order N_POLY extrapolation to the origin using the points inerpolated
      ! on the linear grid rather than the log grid, since this represents an
      ! extrapolation of 1 grid point.
      ! Extrapolation from the innermost points of the log grid would represent
      ! a huge extrapolation, and round-off error would not be acceptable.
      ! If the fitted function is a polynomial of order N_POLY or less, this is exact.
      if (abs(x_out(1)) < 1.0e-12_dp) then
         istart = 2
         sum = 0.0_dp
         do iy = istart, istart + N_POLY
            if (abs(y_out(iy)) < 1.0e-12_dp) cycle
            term = y_out(iy)
            do iprod = istart, istart + N_POLY
               if (iprod == iy) cycle
               term = term * (0.0_dp - x_out(iprod)) / (x_out(iy) - x_out(iprod))
            end do
            sum = sum + term
         end do
         y_out(1) = sum
      end if  ! iszero(x_out(1))
   end do  ! j
   return
end subroutine interpolate

end module interpolate_m
