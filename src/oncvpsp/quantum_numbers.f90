module quantum_numbers_m
   use precision_m, only: dp
   implicit none
   private
   public :: qn_t, &
      qn_null, &
      qn_init, &
      qn_equal, &
      qn_label, &
      qn_max_occ, &
      qn_wf_dim, &
      operator(==)

   interface operator (==)
      module procedure qn_equal
   end interface operator (==)

   type :: qn_t
      !> Principal quantum number
      integer :: n
      !> Angular momentum quantum number
      integer :: l
      !> Spin quantum number
      real(dp) :: s
      !> Total angular momentum quantum number
      real(dp) :: j
      !> Magnetic quantum number
      real(dp) :: m
      !> Spin projection quantum number
      real(dp) :: sigma
      !> Kappa quantum number
      integer :: kappa
   end type qn_t

contains

!> Set quantum numbers to null (zero) values.
subroutine qn_null(qn)
   type(qn_t), intent(out) :: qn

   qn%n = 0
   qn%l = 0
   qn%j = 0.0_dp
   qn%kappa = 0
   return
end subroutine qn_null

!> Initialize a set of quantum numbers.
subroutine qn_init(qn, n, l, j)
   ! Input variables
   integer, intent(in) :: n
   integer, intent(in) :: l
   real(dp), intent(in), optional :: j

   ! Output variables
   type(qn_t), intent(out) :: qn

   qn%n = n
   qn%l = l
   if (present(j)) then
      qn%j = j
   else
      qn%j = 0.0_dp
   end if
   ! Determine kappa
   if (abs(qn%j) > tiny(1.0_dp)) then
      qn%kappa = -int(2.0_dp * (qn%j - dble(qn%l)) * (qn%j + 0.5_dp))
   else
      qn%kappa = 0
   end if

   ! Check consistency
   if (n <= 0) error stop 'qn_init: n must be positive'
   if (l < 0) error stop 'qn_init: l < 0'
   if (l > 3) error stop 'qn_init: l > 3 not supported'
   if (l > n - 1) error stop 'qn_init: l > n-1'
   if (abs(qn%j) > tiny(1.0_dp)) then
      if((abs(qn%j - dble(l) - 0.5_dp) > tiny(1.0_dp)) .and. &
        (abs(qn%j - dble(l) + 0.5_dp) > tiny(1.0_dp))) then
         error stop 'qn_init: j inconsistent with l'
      end if
   end if

   return
end subroutine qn_init

!> Whether two sets of quantum numbers are equal.
elemental function qn_equal(qn_lhs, qn_rhs) result(is_equal)
   type(qn_t), intent(in) :: qn_lhs
   type(qn_t), intent(in) :: qn_rhs
   logical :: is_equal

   is_equal = (qn_lhs%n == qn_rhs%n) .and. &
      (qn_lhs%l == qn_rhs%l) .and. &
      (abs(qn_lhs%j - qn_rhs%j) < tiny(1.0_dp)) .and. &
      (qn_lhs%kappa == qn_rhs%kappa)

   return
end function qn_equal

!> String representation of the quantum numbers.
function qn_label(qn) result(label)
   type(qn_t), intent(in) :: qn

   character(len=10) :: label

   select case (qn%l)
    case (0)
      write(label,'(i1,a)') qn%n, 's'
    case (1)
      write(label,'(i1,a)') qn%n, 'p'
    case (2)
      write(label,'(i1,a)') qn%n, 'd'
    case (3)
      write(label,'(i1,a)') qn%n, 'f'
    case default
      error stop 'qn_label: l > 3 not supported'
   end select

   if (abs(qn%j) > tiny(1.0_dp)) then
      write(label(3:5), '(i1,".5")') int(qn%j - 0.5_dp)
   end if

   return
end function qn_label

!> Maximum occupancy of the state with given quantum numbers.
elemental function qn_max_occ(qn) result(max_occ)
   type(qn_t), intent(in) :: qn
   real(dp) :: max_occ

   max_occ =  2.0_dp * qn%j + 1.0_dp

   return
end function qn_max_occ

!> Number of dimensions of the wavefunction.
!> 2 if relativistic, 1 if scalar- or non-relativistic.
elemental function qn_wf_dim(qn) result(dim)
   type(qn_t), intent(in) :: qn
   integer :: dim

   if (abs(qn%j) > tiny(1.0_dp)) then
      dim = 2
   else
      dim = 1
   end if

   return
end function qn_wf_dim


end module quantum_numbers_m
