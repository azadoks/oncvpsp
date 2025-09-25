module atom_m
   use precision_m, only: dp
   use state_m, only: sr_state_t, init_sr_state, destroy_state
   implicit none
   private
   public :: atom_null
   public :: atom_end

   type, public :: atom_t
      !> Total energy
      real(dp) :: etot
      !> Sum of eigenvales
      real(dp) :: eeig
      !> Hartree energy
      real(dp) :: ehartree
      !> Exchange-correlation energy
      real(dp) :: exc
      !> Total charge density
      real(dp), allocatable :: rho_tot(:)
      !> All-electron total potential
      real(dp), allocatable :: vtot(:)
      !> All-electron hartree potential
      real(dp), allocatable :: vhartree(:)
      !> All-electron exchange-correlation potential
      real(dp), allocatable :: vxc(:)
      !> All-electron bound states
      type(sr_state_t), allocatable :: bound_states(:)
   end type atom_t

contains

subroutine atom_null(atom, n_states, mesh_len)
   ! Input variables
   !> Number of atomic states
   integer, intent(in) :: n_states
   !> Size of radial mesh
   integer, intent(in) :: mesh_len

   ! Output variables
   type(atom_t), intent(out) :: atom

   ! Local variables
   integer :: i

   atom%etot = 0.0_dp
   atom%ehartree = 0.0_dp
   atom%exc = 0.0_dp
   atom%eeig = 0.0_dp
   allocate(atom%rho_tot(mesh_len)); atom%rho_tot = 0.0_dp
   allocate(atom%vtot(mesh_len)); atom%vtot = 0.0_dp
   allocate(atom%vhartree(mesh_len)); atom%vhartree = 0.0_dp
   allocate(atom%vxc(mesh_len)); atom%vxc = 0.0_dp
   allocate(atom%bound_states(n_states))
   do i = 1, n_states
      call init_sr_state(atom%bound_states(i), 0, 0, 0.0_dp, 0.0_dp, 0.0_dp, mesh_len)
   end do

   return
end subroutine atom_null

subroutine atom_end(atom)
   ! Input/Output variables
   type(atom_t), intent(in out) :: atom

   integer :: i

   atom%etot = 0.0_dp
   atom%ehartree = 0.0_dp
   atom%exc = 0.0_dp
   atom%eeig = 0.0_dp
   if (allocated(atom%rho_tot)) deallocate(atom%rho_tot)
   if (allocated(atom%vtot)) deallocate(atom%vtot)
   if (allocated(atom%vhartree)) deallocate(atom%vhartree)
   if (allocated(atom%vxc)) deallocate(atom%vxc)
   if (allocated(atom%bound_states)) then
      do i = 1, size(atom%bound_states)
         call destroy_state(atom%bound_states(i))
      end do
      deallocate(atom%bound_states)
   end if

   return
end subroutine atom_end

end module atom_m
