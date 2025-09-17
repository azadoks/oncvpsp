module constants_m
   use precision_m, only: dp
   implicit none
   private
   !> π
   real(dp), parameter, public :: pi = 3.141592653589793238462643383279502884197_dp
   !> 2π
   real(dp), parameter, public :: twopi = 2.0_dp * pi
   !> 4π
   real(dp), parameter, public :: fourpi = 4.0_dp * pi
   !> 1/(4π)
   real(dp), parameter, public :: ifourpi = 1.0_dp / fourpi
   !> 1/3
   real(dp), parameter, public :: third = 1.0_dp / 3.0_dp
end module constants_m
