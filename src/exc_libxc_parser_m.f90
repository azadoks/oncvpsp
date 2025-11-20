module exc_libxc_parser_m
   implicit none
   interface parse_variable
      module procedure real_parser ! interface for a module
      module procedure logical_parser ! procedure is implicit
   end interface parse_variable
contains
subroutine real_parser(varval, var)
   real(8), intent(in) :: varval
   real(8), intent(out) :: var
   var = varval
end subroutine real_parser
subroutine logical_parser(varval, var)
   logical, intent(in) :: varval
   logical, intent(out) :: var
   var = varval
end subroutine logical_parser
end module exc_libxc_parser_m
