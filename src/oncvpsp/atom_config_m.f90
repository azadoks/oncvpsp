module atom_config_m
   use precision_m, only: dp
   implicit none
   private
   type, public :: atom_config_t
      !> Atomic number
      real(dp) :: z
      !> State principal quantum numbers
      integer, allocatable :: ns(:)
      !> State angular momenta
      integer, allocatable :: ells(:)
      !> State occupancies
      real(dp), allocatable :: occs(:)
   end type atom_config_t
contains
end module atom_config_m
