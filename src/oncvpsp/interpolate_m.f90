module interpolate_m
   use precision_m, only: dp
   implicit none
   private
   public :: interpolate_poly

contains

subroutine interpolate_poly(x_src, y_src, x_dst, y_dst, order)
   ! Input variables
   real(dp), intent(in) :: x_src(:)
   real(dp), intent(in) :: y_src(:)
   real(dp), intent(in) :: x_dst(:)
   integer, intent(in), optional :: order

   ! Output variables
   real(dp), intent(out) :: y_dst(:)

   ! Local variables
   !> Number of points in the source domain
   integer :: np_src
   !> Number of points in the destination domain
   integer :: np_dst
   !> Loop indices
   integer :: i, j, k
   !> Low and high indices for bracketing
   integer :: lo, hi
   !> Starting index for interpolation
   integer :: start
   !> Actual interpolation order
   integer :: actual_order
   !> Temporary variable for destination value
   real(dp) :: y_dst_tmp

   if (.not. present(order)) then
      actual_order = 7
   else
      actual_order = order
   end if

   np_src = size(x_src)
   np_dst = size(x_dst)

   if (np_src < actual_order + 1) then
      error stop 'interpolate_poly: not enough source points for the specified order'
   end if

   y_dst = 0.0_dp
   do i = 2, size(y_dst)
      if ((x_dst(i) < x_src(1)) .or. (x_dst(i) > x_src(np_src))) then
         error stop 'interpolate_poly: interpolation ERROR - out of range'
      end if
      call bracket_index(x_src, x_dst(i), lo, hi)
      if ((mod(actual_order, 2) == 1) .or. (x_dst(i) - x_src(lo) < x_src(hi) - x_dst(i))) then
         start = lo - int(real(actual_order, dp) / 2.0_dp)
      else
         start = hi - int(real(actual_order, dp) / 2.0_dp)
      end if

      start = min(start, np_src - actual_order)
      start = max(start, 1)

      do j = start, start + actual_order
         if (abs(y_src(j)) < tiny(1.0_dp)) cycle
         y_dst_tmp = y_src(j)
         do k = start, start + actual_order
            if (k == j) cycle
            y_dst_tmp = y_dst_tmp * (x_dst(i) - x_src(k)) / (x_src(j) - x_src(k))
         end do  ! j = start, start + order
         y_dst(i) = y_dst(i) + y_dst_tmp
      end do  ! k = start, start + order

      ! special treatment for the origin
      ! Do order npoly EXTRAPOLATION to the origin using the points inerpolated
      ! on the linear grid rather than the log grid, since this represents an
      ! extrapolation of 1 grid point.  Exrapolation from the innermost points
      ! of the log grid would represent a huge extrapolation, and round-off
      ! error would not be acceptable If the fitted function is a polynomial
      ! of order npoly or less, this is exact.
      if (abs(x_dst(1)) < tiny(1.0_dp)) then
         y_dst(1) = 0.0_dp
         do j = 2, 2 + actual_order
            if (abs(y_dst(j)) < tiny(1.0_dp)) cycle
            y_dst_tmp = y_dst(j)
            do k = 2, 2 + actual_order
               if (k == j) cycle
               y_dst_tmp = y_dst_tmp * (x_dst(i) - x_dst(k)) / (x_dst(j) - x_dst(k))
            end do  ! k = 2, 2, + order
            y_dst(1) = y_dst(1) + y_dst_tmp
         end do  ! j = 2, 2 + order
      end if  ! x_dst(1) == 0
   end do  ! i = 2, size(dst)

   return
end subroutine interpolate_poly

subroutine bracket_index(x, y, lo, hi)
   ! Input variables
   real(dp), intent(in) :: x(:)
   real(dp), intent(in) :: y

   ! Output variables
   integer, intent(out) :: lo
   integer, intent(out) :: hi

   ! Local variables
   integer :: i, j
   integer :: np

   np = size(x)
   lo = 1
   hi = np
   do i = 1, np
      j = (lo + hi) / 2
      if (y > x(j)) then
         lo = j
      else
         hi = j
      end if
      if (hi - lo == 1) exit
   end do

   return
end subroutine bracket_index

end module interpolate_m
