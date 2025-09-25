module mesh_m
   use precision_m, only: dp
   implicit none
   private

   interface assignment (=)
      module procedure mesh_copy
   end interface assignment (=)

   interface operator (==)
      module procedure mesh_equal
   end interface operator (==)

   type mesh_t
      integer, private :: type_
      real(dp), public :: a
      real(dp), public :: b
      integer, public :: np
      real(dp), allocatable, public :: r(:)
      real(dp), allocatable, public :: dr(:)
   end type mesh_t

   ! Mesh types
   integer, parameter :: MESH_LINEAR = 1
   integer, parameter :: MESH_LOG    = 2

contains

subroutine mesh_null(mesh)
   type(mesh_t), intent(out) :: mesh

   mesh%type_ = 0
   mesh%a = 0.0_dp
   mesh%b = 0.0_dp
   mesh%np = 0
   if (allocated(mesh%r)) deallocate(mesh%r)
   if (allocated(mesh%dr)) deallocate(mesh%dr)

   return
end subroutine mesh_null

subroutine mesh_init(mesh, type_, rmin, rmax, np, a)
   type(mesh_t), intent(in out) :: mesh
   integer, intent(in) :: type_
   real(dp), intent(in) :: rmin
   real(dp), intent(in), optional :: rmax
   integer, intent(in), optional :: np
   real(dp), intent(in), optional :: a

   if ((type_ == MESH_LOG) .and. (present(rmax) .and. present(a))) then
      call generate_log_mesh(mesh, rmin, rmax, a)
   else if ((type_ == MESH_LINEAR) .and. (present(rmax) .and. present(np))) then
      call generate_lin_mesh(mesh, rmin, a, np)
   else
      error stop 'mesh_init: must specify (LOG_MESH, rmax, a) or (LIN_MESH, a, np)'
   end if

   return
end subroutine mesh_init

subroutine generate_log_mesh(mesh, rmin, rmax, a)
   type(mesh_t), intent(out) :: mesh
   real(dp), intent(in) :: rmin
   real(dp), intent(in) :: rmax
   real(dp), intent(in) :: a

   integer :: i

   mesh%type_ = MESH_LOG
   mesh%a = a
   mesh%b = rmax
   mesh%np = int(log(rmax / rmin) / a) + 1
   allocate(mesh%r(mesh%np))
   allocate(mesh%dr(mesh%np))
   do i = 1, mesh%np
      mesh%r(i) = rmin * exp(a * real(i - 1, dp))
      mesh%dr(i) = a * mesh%r(i)
   end do

   return
end subroutine generate_log_mesh

subroutine generate_lin_mesh(mesh, rmin, a, np)
   type(mesh_t), intent(out) :: mesh
   real(dp), intent(in) :: rmin
   real(dp), intent(in) :: a
   integer, intent(in) :: np

   integer :: i

   mesh%type_ = MESH_LINEAR
   mesh%a = a
   mesh%b = rmin + a * real(np - 1, dp)
   mesh%np = np
   allocate(mesh%r(mesh%np))
   allocate(mesh%dr(mesh%np))
   do i = 1, mesh%np
      mesh%r(i) = rmin + a * real(i - 1, dp)
      mesh%dr(i) = a
   end do

   return
end subroutine generate_lin_mesh

function mesh_equal(mesh_lhs, mesh_rhs) result(is_equal)
   type(mesh_t), intent(in) :: mesh_lhs
   type(mesh_t), intent(in) :: mesh_rhs
   logical :: is_equal

   is_equal = .false.
   if (mesh_lhs%type_ /= mesh_rhs%type_) return
   if (abs(mesh_lhs%a - mesh_rhs%a) > tiny(1.0_dp)) return
   if (abs(mesh_lhs%b - mesh_rhs%b) > tiny(1.0_dp)) return
   if (mesh_lhs%np /= mesh_rhs%np) return
   if (.not. allocated(mesh_lhs%r) .or. .not. allocated(mesh_rhs%r)) return
   is_equal = .true.

   return
end function mesh_equal

subroutine mesh_copy(mesh_lhs, mesh_rhs)
   type(mesh_t), intent(out) :: mesh_lhs
   type(mesh_t), intent(in) :: mesh_rhs

   if (.not. allocated(mesh_rhs%r) .or. .not. allocated(mesh_rhs%dr)) then
      error stop 'mesh_copy: rhs mesh not allocated'
   end if

   mesh_lhs%type_ = mesh_rhs%type_
   mesh_lhs%a = mesh_rhs%a
   mesh_lhs%b = mesh_rhs%b
   mesh_lhs%np = mesh_rhs%np

   if (allocated(mesh_lhs%r)) deallocate(mesh_lhs%r)
   allocate(mesh_lhs%r(mesh_rhs%np))
   mesh_lhs%r(:) = mesh_rhs%r(:)

   if (allocated(mesh_lhs%dr)) deallocate(mesh_lhs%dr)
   allocate(mesh_lhs%dr(mesh_rhs%np))
   mesh_lhs%dr(:) = mesh_rhs%dr(:)

   return
end subroutine mesh_copy

subroutine mesh_transfer(mesh_src, f_src, mesh_dst, f_dst)
   type(mesh_t), intent(in) :: mesh_src
   real(dp), intent(in) :: f_src(mesh_src%np)
   type(mesh_t), intent(in) :: mesh_dst
   real(dp), intent(out) :: f_dst(mesh_dst%np)

   if (mesh_src%type_ /= MESH_LOG) then
      error stop 'mesh_transfer: only log mesh supported for source'
   end if
   if (mesh_dst%type_ /= MESH_LINEAR) then
      error stop 'mesh_transfer: only linear mesh supported for destination'
   end if

   call interpolate_poly(mesh_src%r, f_src, mesh_src%np, mesh_dst%r, f_dst, mesh_dst%np, order=7)

   return

end subroutine mesh_transfer

subroutine mesh_end(mesh)
   type(mesh_t), intent(in out) :: mesh

   mesh%type_ = 0
   mesh%a = 0.0_dp
   mesh%b = 0.0_dp
   mesh%np = 0
   if (allocated(mesh%r)) deallocate(mesh%r)
   if (allocated(mesh%dr)) deallocate(mesh%dr)

   return
end subroutine mesh_end

end module mesh_m
