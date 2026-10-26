! Unit tests of the text input reader (input_text_m): next_line, which skips
! the blank and comment lines and returns lines of any length with their line
! numbers, and read_input_text on an input laid out in all the ways the reader
! must accept. Stops with a non-zero code at the first failure.
 program test_input_text
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use input_text_m, only: next_line, read_input_text
 use input_m, only: input_t
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

!! counts and reports a failed check
!!
!! @param ok  result of the check
!! @param what  what was checked

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

!! opens a scratch file holding lines, with or without a newline after the
!! last, and returns its unit
!!
!! @param lines  the lines of the file
!! @param final_newline  .true. to end the last line with a newline

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

!! checks next_line on comments, blank, long and unterminated lines

!Local variables
 integer :: ifail
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
 call next_line(unit, line, inline, ifail=ifail)
 if(ifail /= 0) error stop 1
 call check(line == '  first data' .and. inline == 5, &
&           'next_line skips blank and comment lines')
 call next_line(unit, line, inline, ifail=ifail)
 if(ifail /= 0) error stop 1
 call check(line == repeat('x', 600)//' long' .and. inline == 7, &
&           'next_line reads a long line whole')
 call next_line(unit, line, inline, ifail=ifail)
 if(ifail /= 0) error stop 1
 call check(line == 'last' .and. inline == 8, &
&           'next_line reads a last line without a newline')
 call next_line(unit, line, inline, eof, ifail=ifail)
 if(ifail /= 0) error stop 1
 call check(eof .and. inline == 8, &
&           'next_line reports the end of the input when asked')
 close(unit)
 return
 end subroutine test_next_line

 subroutine test_read_input_text(final_newline, with_mesh)

!! checks read_input_text on 14_Si.dat, with blank lines and comments inside
!! the blocks, a long line and comments after the data, and with or without
!! the optional log mesh line and UPF mesh size
!!
!! @param final_newline  .true. to end the last line with a newline
!! @param with_mesh  .true. to add the log mesh line, upf_mesh to the
!!                   output grid line, ep_as_given, well_onset,
!!                   well_tail and wave_regularization to the lmax line,
!!                   and icmod 6 with slack and budget_cutoff

!Input variables
 logical, intent(in) :: final_newline, with_mesh

!Local variables
 integer :: ifail
 character(len=700) :: lines(33)
 character(len=700), allocatable :: all_lines(:)
 integer :: unit, inline
 type(input_t) :: inp

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
   lines(25) = '    6.0  0.01  2000'
   lines(13) = '    2    1    1.5    2    2.5e-9    2'
   lines(20) = '    1    2    0.75    5'
   lines(21) = '    2    2    0.75    0    0.3'
   lines(23) = '   6    0.0    1.0e-6    0.8    1'
   lines(2) = '    Si  14     3   2   3   psp8   137.0359895'
   all_lines = [character(len=700) :: lines, '# amesh, rr1, rmax', &
&               '  1.003  -1  60.0', &
&               '# a comment after the mesh']
 else
   all_lines = lines
 end if
 unit = scratch_input(all_lines, final_newline)
 call read_input_text(unit, inline, inp, ifail)
 if(ifail /= 0) error stop 1
 close(unit)

 associate(atsym=>inp%atsym, zz=>inp%zz, nc=>inp%nc, nv=>inp%nv, &
&          iexc=>inp%iexc, psfile=>inp%psfile, na=>inp%na, la=>inp%la, &
&          fa=>inp%fa, lmax=>inp%lmax, rc=>inp%rc, ep=>inp%ep, ncon=>inp%ncon, &
&          nbas=>inp%nbas, qcut=>inp%qcut, lloc=>inp%lloc, lpopt=>inp%lpopt, &
&          dvloc0=>inp%dvloc0, nproj=>inp%nproj, debl=>inp%debl, &
&          icmod=>inp%icmod, fcfact=>inp%fcfact, rcfact=>inp%rcfact, &
&          slack=>inp%slack, budget_cut=>inp%budget_cut, &
&          fcfact_min=>inp%fcfact_min, fcfact_max=>inp%fcfact_max, &
&          fcfact_step=>inp%fcfact_step, rcfact_min=>inp%rcfact_min, &
&          rcfact_max=>inp%rcfact_max, rcfact_step=>inp%rcfact_step, &
&          epsh1=>inp%epsh1, epsh2=>inp%epsh2, depsh=>inp%depsh, &
&          rxpsh=>inp%rxpsh, rlmax=>inp%rlmax, drl=>inp%drl, &
&          upf_mesh=>inp%upf_mesh, ncnf=>inp%ncnf, nvcnf=>inp%nvcnf, &
&          nacnf=>inp%nacnf, lacnf=>inp%lacnf, facnf=>inp%facnf, &
&          amesh=>inp%amesh, rr1=>inp%rr1, rmax=>inp%rmax, &
&          c_light=>inp%c_light, ep_as_given=>inp%ep_as_given, &
&          well_onset=>inp%well_onset, well_tail=>inp%well_tail, &
&          wreg=>inp%wreg, well_n=>inp%well_n, well_norm=>inp%well_norm, &
&          ecut_wt=>inp%ecut_wt)

   call check(atsym == 'Si' .and. zz == 14 .and. nc == 3 .and. nv == 2 &
&             .and. iexc == 3 .and. psfile == 'psp8', 'atom')
   call check(all(na(1:5) == [1, 2, 2, 3, 3]) &
&             .and. all(la(1:5) == [0, 0, 1, 0, 1]) &
&             .and. all(fa(1:5) == [2, 2, 6, 2, 2]), 'reference configuration')
   call check(lmax == 2 .and. all(rc(1:3) == 2.7_dp) &
&             .and. all(ep(1:3, 1) == [0.0_dp, 0.0_dp, 0.1_dp]) &
&             .and. all(ep(1:3, 2) == ep(1:3, 1)) &
&             .and. all(ncon(1:3) == 5) .and. all(nbas(1:3) == 8) &
&             .and. all(qcut(1:3) == 5), &
&             'pseudopotential and optimization')
   call check(lloc == 4 .and. lpopt == 5 .and. rc(5) == 2.7_dp &
&             .and. dvloc0 == 0, &
&             'local potential')
   call check(all(nproj(1:3) == 2) &
&             .and. all(debl(1:3) == [1.0_dp, 0.75_dp, 0.75_dp]), &
&             'projectors')
   if(with_mesh) then
     call check(icmod == 6 .and. fcfact == 0 .and. rcfact == 1.0e-6_dp &
&               .and. slack == 0.8_dp .and. budget_cut == 1, &
&               'model core charge, slack and budget_cutoff')
   else
     call check(icmod == 1 .and. fcfact == 0.25_dp .and. slack == 1 &
&               .and. budget_cut == 0, 'model core charge')
   end if
   call check(epsh1 == -2 .and. epsh2 == 2 .and. depsh == 0.02_dp &
&             .and. rxpsh < 0, &
&             'log derivative analysis')
   call check(rlmax == 6 .and. drl == 0.01_dp, 'output grid')
   call check(ncnf == 1 .and. nvcnf(2) == 2 .and. all(nacnf(4:5, 2) == 3) &
&             .and. all(lacnf(4:5, 2) == [0, 1]) &
&             .and. all(facnf(4:5, 2) == [2, 1]), &
&             'test configurations')
   if(with_mesh) then
     call check(upf_mesh == 2000, 'UPF mesh size')
     call check(amesh == 1.003_dp .and. rr1 < 0 .and. rmax == 60, 'log mesh')
     call check(c_light == 137.0359895_dp, 'speed of light')
     call check(ep_as_given .and. well_onset == 1.5_dp .and. well_tail == 2 &
&               .and. wreg == 2.5e-9_dp, &
&               'ep_as_given, well_onset, well_tail, wave_regularization')
     call check(all(well_n(1:3) == [0, 5, 0]), 'well_n')
     call check(well_norm == 2, 'well_norm')
     call check(all(ecut_wt(1:3) == [-1.0_dp, -1.0_dp, 0.3_dp]), 'ecut_weight')
     call check(inline == 35, 'line number of the last data line')
   else
     call check(upf_mesh == 0, 'UPF mesh size default')
     call check(amesh < 0 .and. rr1 < 0 .and. rmax < 0, 'log mesh defaults')
     call check(c_light < 0, 'speed of light default')
     call check(.not. ep_as_given .and. well_onset == 1.006_dp**5 &
&               .and. well_tail == 0 .and. wreg == 0, &
&               'ep_as_given, well_onset, well_tail, wave_regularization '// &
&               'defaults')
     call check(all(well_n(1:3) == 0), 'well_n default')
     call check(well_norm == 0, 'well_norm default')
     call check(all(ecut_wt(1:3) < 0), 'ecut_weight default')
     call check(inline == 32, 'line number of the last data line')
   end if
 end associate
 return
 end subroutine test_read_input_text

 end program test_input_text
