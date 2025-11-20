module exc_libxc_xc_messages_m
   implicit none
   ! Messages
   character(1200) :: message(2)
contains
subroutine messages_info(lmess, iunit)
   integer, intent(in) :: lmess, iunit
   integer :: i
   do i=1,lmess
      write(iunit,'("LXC ",A)') trim(message(i))
   end do
end subroutine messages_info

subroutine messages_fatal(lmess)
   integer, intent(in) :: lmess
   integer :: i
   do i=1,lmess
      write(*,'("Fatal ",A)') trim(message(i))
   end do
   stop
end subroutine messages_fatal

subroutine messages_warning(lmess)
   integer, intent(in) :: lmess
   integer :: i
   do i=1,lmess
      write(*,'("Warning ",A)') trim(message(i))
   end do
end subroutine messages_warning

subroutine messages_input_error(str1,str2)
   character(*), intent(in) :: str1,str2
   write(*,'("Input error ",A,A)') trim(str1),trim(str2)
end subroutine messages_input_error

subroutine messages_experimental(str)
   character(*), intent(in) :: str
   write(*,'("Experimental xc ",A)') trim(str)
end subroutine messages_experimental

subroutine messages_not_implemented(str)
   character(*), intent(in) :: str
   write(*,'("Not implemented ",A)') trim(str)
end subroutine messages_not_implemented
end module exc_libxc_xc_messages_m
