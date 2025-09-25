module state_m
   use precision_m, only: dp
   implicit none
   private

   type, public :: state_t
      !> Quantum numbers
      type(qn_t) :: qn
      !> Angular momentum
      integer :: l
      !> Occupancy
      real(dp) :: occ
      !> Eigenvalue
      real(dp) :: eig
      !> Radial coordinate of the outermost peak of the wave function
      real(dp) :: peak
      !> Radial wavefunction on log mesh (r*ψ(r))
      real(dp), allocatable :: psi(:)
      !> Radial wavefunction derivative on log mesh (d(r*ψ(r))/dr)
      real(dp), allocatable :: dpsi(:)
      !> 1-electron charge density (ρ(r) / r**2)
      real(dp), allocatable :: rho(:)
   end type state_t

contains

end module state_m
