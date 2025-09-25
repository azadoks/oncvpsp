module orbital
   use precision_m, only: dp
   implicit none

   type orbital
      real(dp), allocatable :: rpsi
      real(dp), allocatable :: d_rpsi_dr
   end type orbital

contains

end module orbital
