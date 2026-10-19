! Unit tests of the text input reader (input_text_m): next_line, which skips
! the blank and comment lines and returns lines of any length with their line
! numbers, and read_input_text on an input laid out in all the ways the reader
! must accept. Stops with a non-zero code at the first failure.
 program test_input_text
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use input_text_m, only: next_line, read_input_text
 implicit none

!Local variables
 character(len=*), parameter :: tab = achar(9)
 integer :: nfail = 0

 call test_next_line()
 call test_read_input_text(final_newline=.true., with_mesh=.false.)
 call test_read_input_text(final_newline=.false., with_mesh=.false.)
 call test_read_input_text(final_newline=.true., with_mesh=.true.)
 call test_read_input_text(final_newline=.false., with_mesh=.true.)

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(*, '(a)') 'input_text_m: all checks passed'

 contains

 subroutine check(ok, what)

!counts and reports a failed check

!ok  result of the check
!what  what was checked

!Input variables
 logical, intent(in) :: ok
 character(len=*), intent(in) :: what

 if(.not. ok) then
   write(stderr, '(a)') 'FAILED: '//what
   nfail = nfail + 1
 end if
 return
 end subroutine check

 integer function scratch_input(lines, final_newline) result(unit)

!opens a scratch file holding lines, with or without a newline after the
!last, and returns its unit

!lines  the lines of the file
!final_newline  .true. to end the last line with a newline

!Input variables
 character(len=*), intent(in) :: lines(:)
 logical, intent(in) :: final_newline

!Local variables
 integer :: ii

 open(newunit=unit, status='scratch', action='readwrite', access='stream', &
&     form='formatted')
 do ii = 1, size(lines)
   if(ii < size(lines) .or. final_newline) then
     write(unit, '(a)') trim(lines(ii))
   else
     write(unit, '(a)', advance='no') trim(lines(ii))
   end if
 end do
 rewind(unit)
 return
 end function scratch_input

 subroutine test_next_line()

!checks next_line on comments, blank, long and unterminated lines

!Local variables
 character(len=700) :: lines(8)
 character(len=:), allocatable :: line
 integer :: unit, inline
 logical :: eof

 lines(1) = '# comment'
 lines(2) = ''
 lines(3) = '   # indented comment'
 lines(4) = tab//'# comment after a tab'
 lines(5) = '  first data'
 lines(6) = '#'
 lines(7) = repeat('x', 600)//' long'
 lines(8) = 'last'
 unit = scratch_input(lines, final_newline=.false.)
 inline = 0
 call next_line(unit, line, inline)
 call check(line == '  first data' .and. inline == 5, &
&           'next_line skips blank and comment lines')
 call next_line(unit, line, inline)
 call check(line == repeat('x', 600)//' long' .and. inline == 7, &
&           'next_line reads a long line whole')
 call next_line(unit, line, inline)
 call check(line == 'last' .and. inline == 8, &
&           'next_line reads a last line without a newline')
 call next_line(unit, line, inline, eof)
 call check(eof .and. inline == 8, &
&           'next_line reports the end of the input when asked')
 close(unit)
 return
 end subroutine test_next_line

 subroutine test_read_input_text(final_newline, with_mesh)

!checks read_input_text on 14_Si.dat, with blank lines and comments inside
!the blocks, a long line and comments after the data, and with or without
!the optional log mesh line

!final_newline  .true. to end the last line with a newline
!with_mesh  .true. to add the log mesh line

!Input variables
 logical, intent(in) :: final_newline, with_mesh

!Local variables
 character(len=700) :: lines(33)
 character(len=700), allocatable :: all_lines(:)
 integer :: unit, inline
 character(len=2) :: atsym
 character(len=4) :: psfile
 real(dp) :: zz, fa(30), rc(6), ep(6), qcut(6), dvloc0, debl(6)
 real(dp) :: fcfact, rcfact, fcfact_min, fcfact_max, fcfact_step
 real(dp) :: rcfact_min, rcfact_max, rcfact_step
 real(dp) :: epsh1, epsh2, depsh, rxpsh, rlmax, drl, facnf(30, 5)
 real(dp) :: amesh, rr1, rmax
 integer :: nc, nv, iexc, na(30), la(30), lmax, ncon(6), nbas(6)
 integer :: lloc, lpopt, nproj(6)
 integer :: icmod, ncnf, nvcnf(5), nacnf(30, 5), lacnf(30, 5)

 lines = [character(len=700) :: &
&   '# atsym, z, nc, nv, iexc   psfile', &
&   '    Si  14     3   2   3   psp8', &
&   '', &
&   '# n, l, f  (nc+nv lines)', &
&   '    1    0    2.0', &
&   '    2    0    2.0', &
&   '  # a comment inside the block', &
&   '    2    1    6.0', &
&   '', &
&   '    3    0    2.0', &
&   tab//'3'//tab//'1'//tab//'2.0', &
&   '# lmax', &
&   '    2', &
&   '    0    2.70    0.00    5    8    5.00'//repeat(' ', 500)// &
&   'text after the data', &
&   '    1    2.70    0.00    5    8    5.00', &
&   '    2    2.70    0.10    5    8    5.00', &
&   '#', &
&   '    4    5    2.7  0.0', &
&   '    0    2    1.0', &
&   '    1    2    0.75', &
&   '    2    2    0.75', &
&   '# icmod, fcfact', &
&   '   1    0.25', &
&   '   -2.0  2.0  0.02', &
&   '    6.0  0.01', &
&   '# ncnf', &
&   '    1', &
&   '#', &
&   '    2', &
&   '    3    0    2.0', &
&   '', &
&   '    3    1    1.0', &
&   '# a comment after the data']
 if(with_mesh) then
   all_lines = [character(len=700) :: lines, '# amesh, rr1, rmax', &
&               '  1.003  -1  60.0', '# a comment after the mesh']
 else
   all_lines = lines
 end if
 unit = scratch_input(all_lines, final_newline)
 call read_input_text(unit, inline, &
&                     atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
&                     lmax, rc, ep, ncon, nbas, qcut, &
&                     lloc, lpopt, dvloc0, nproj, debl, &
&                     icmod, fcfact, rcfact, &
&                     fcfact_min, fcfact_max, fcfact_step, &
&                     rcfact_min, rcfact_max, rcfact_step, &
&                     epsh1, epsh2, depsh, rxpsh, &
&                     rlmax, drl, &
&                     ncnf, nvcnf, nacnf, lacnf, facnf, &
&                     amesh, rr1, rmax)
 close(unit)

 call check(atsym == 'Si' .and. zz == 14 .and. nc == 3 .and. nv == 2 &
&           .and. iexc == 3 .and. psfile == 'psp8', 'atom')
 call check(all(na(1:5) == [1, 2, 2, 3, 3]) &
&           .and. all(la(1:5) == [0, 0, 1, 0, 1]) &
&           .and. all(fa(1:5) == [2, 2, 6, 2, 2]), 'reference configuration')
 call check(lmax == 2 .and. all(rc(1:3) == 2.7_dp) &
&           .and. all(ep(1:3) == [0.0_dp, 0.0_dp, 0.1_dp]) &
&           .and. all(ncon(1:3) == 5) .and. all(nbas(1:3) == 8) &
&           .and. all(qcut(1:3) == 5), &
&           'pseudopotential and optimization')
 call check(lloc == 4 .and. lpopt == 5 .and. rc(5) == 2.7_dp &
&           .and. dvloc0 == 0, &
&           'local potential')
 call check(all(nproj(1:3) == 2) &
&           .and. all(debl(1:3) == [1.0_dp, 0.75_dp, 0.75_dp]), &
&           'projectors')
 call check(icmod == 1 .and. fcfact == 0.25_dp, 'model core charge')
 call check(epsh1 == -2 .and. epsh2 == 2 .and. depsh == 0.02_dp &
&           .and. rxpsh < 0, &
&           'log derivative analysis')
 call check(rlmax == 6 .and. drl == 0.01_dp, 'output grid')
 call check(ncnf == 1 .and. nvcnf(2) == 2 .and. all(nacnf(4:5, 2) == 3) &
&           .and. all(lacnf(4:5, 2) == [0, 1]) &
&           .and. all(facnf(4:5, 2) == [2, 1]), &
&           'test configurations')
 if(with_mesh) then
   call check(amesh == 1.003_dp .and. rr1 < 0 .and. rmax == 60, 'log mesh')
   call check(inline == 35, 'line number of the last data line')
 else
   call check(amesh < 0 .and. rr1 < 0 .and. rmax < 0, 'log mesh defaults')
   call check(inline == 32, 'line number of the last data line')
 end if
 return
 end subroutine test_read_input_text

 end program test_input_text
