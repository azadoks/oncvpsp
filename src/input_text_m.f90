 module input_text_m
 use, intrinsic :: iso_fortran_env, only: dp => real64, stdout => output_unit, &
&                                         stderr => error_unit
 implicit none
 private
 public :: read_input_text, read_input_text_r, next_line
 contains

 subroutine read_input_text(unit, inline, &
&                           atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
&                           lmax, rc, ep, ncon, nbas, qcut, &
&                           lloc, lpopt, dvloc0, nproj, debl, &
&                           icmod, fcfact, rcfact, &
&                           fcfact_min, fcfact_max, fcfact_step, &
&                           rcfact_min, rcfact_max, rcfact_step, &
&                           epsh1, epsh2, depsh, rxpsh, &
&                           rlmax, drl, &
&                           ncnf, nvcnf, nacnf, lacnf, facnf, &
&                           amesh, rr1, rmax)
!> File unit number for input
 integer, intent(in) :: unit
!> Current line number in input file
 integer, intent(out) :: inline
! [Atom and reference configuration]
!> Atomic symbol
 character(len=2), intent(out) :: atsym
!> Atomic number
 real(dp), intent(out) :: zz
!> Number of core electrons
 integer, intent(out) :: nc
!> Number of valence electrons
 integer, intent(out) :: nv
!> Exchange-correlation functional code
 integer, intent(out) :: iexc
!> Pseudopotential output file format flag
 character(len=4), intent(out) :: psfile
!> Principal quantum number array
 integer, intent(out) :: na(30)
!> Orbital angular momentum quantum number array
 integer, intent(out) :: la(30)
!> Occupation number array
 real(dp), intent(out) :: fa(30)
! [Pseudopotential and optimization]
!> Maximum angular momentum
 integer, intent(out) :: lmax
!> Core radii for pseudopotentials
 real(dp), intent(out) :: rc(6)
!> Pseudopotential energies
 real(dp), intent(out) :: ep(6)
!> Number of matching constraints for each pseudopotential
 integer, intent(out) :: ncon(6)
!> Number of basis functions for each pseudopotential
 integer, intent(out) :: nbas(6)
!> Maximum wave number for each pseudopotential
 real(dp), intent(out) :: qcut(6)
! [Local potential]
!> Angular momentum used for local potential
 integer, intent(out) :: lloc
!>
 integer, intent(out) :: lpopt
!> Local potential offset at origin
 real(dp), intent(out) :: dvloc0
! [Vanderbilt-Kleinman-Bylander projectors]
!> Number of projectors for each angular momentum
 integer, intent(out) :: nproj(6)
!> Energy shift for unbound states
 real(dp), intent(out) :: debl(6)
! [Model core charge]
!> Model core charge flag
 integer, intent(out) :: icmod
!> Scaling factor for core charge
 real(dp), intent(out) :: fcfact
!> Scaling factor range minimum for optimization
 real(dp), intent(out) :: fcfact_min
!> Scaling factor range maximum for optimization
 real(dp), intent(out) :: fcfact_max
!> Scaling factor step size for optimization
 real(dp), intent(out) :: fcfact_step
!> Core charge width factor
 real(dp), intent(out) :: rcfact
!> Core charge width factor minimum for optimization
 real(dp), intent(out) :: rcfact_min
!> Core charge width factor maximum for optimization
 real(dp), intent(out) :: rcfact_max
!> Core charge width factor step size for optimization
 real(dp), intent(out) :: rcfact_step
! [Log derivative analysis]
!> Lower bound of energy range for log derivative analysis
 real(dp), intent(out) :: epsh1
!> Upper bound of energy range for log derivative analysis
 real(dp), intent(out) :: epsh2
!> Energy step size for log derivative analysis
 real(dp), intent(out) :: depsh
!> Radius for log derivative analysis
 real(dp), intent(out) :: rxpsh
! [Output grid]
!> Maximum radius for output grid
 real(dp), intent(out) :: rlmax
!> Grid spacing for output grid
 real(dp), intent(out) :: drl
! [Test configurations]
!> Number of test configurations
 integer, intent(out) :: ncnf
!> Number of valence states in test configurations
 integer, intent(out) :: nvcnf(5)
!> Principal quantum number array for test configurations
 integer, intent(out) :: nacnf(30, 5)
!> Orbital angular momentum quantum number array for test configurations
 integer, intent(out) :: lacnf(30, 5)
!> Occupation number array for test configurations
 real(dp), intent(out) :: facnf(30, 5)
! [Log mesh]
!> Ratio of successive points of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: amesh
!> First point of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: rr1
!> Outer radius of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: rmax

!Local variables
!> I/O status variable
 integer :: ios
!> Loop indices
 integer :: ii, jj
!> Angular momentum index (l+1)
 integer :: l1
!> Temporary variable for reading angular momentum
 integer :: lt
!> The input line being read
 character(len=:), allocatable :: line
!> Whether the input has ended
 logical :: eof

 nproj(:) = 0
 fcfact = 0.d0
 rcfact = 0.d0
 rc(:) = 0.d0
 ep(:) = 0.d0
 rxpsh = -1.d0  ! negative value for auto setting of radius
! default values for icmod=4 optimization ranges
! overridden by input values if given with icmod=4
 fcfact_min = 1.5d0
 fcfact_max = 6.0d0
 fcfact_step = 0.5d0
 rcfact_min = 1.0d0
 rcfact_max = 1.9d0
 rcfact_step = 0.1d0

 amesh = -1.d0  ! negative values for the built-in mesh
 rr1 = -1.d0
 rmax = -1.d0

! read input data
 inline = 0

! atom and reference configuration
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) atsym, zz, nc, nv, iexc, psfile
! counts beyond the array sizes would overrun before check_data runs
 if (ios == 0 .and. (nc < 0 .or. nv < 0 .or. nc + nv > 30)) ios = 999
 call read_error(ios, inline)

 do ii = 1, nc + nv
   call next_line(unit, line, inline)
   read (line, *, iostat=ios) na(ii), la(ii), fa(ii)
   call read_error(ios, inline)
 end do

! pseudopotential and optimization
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) lmax
 if (ios == 0 .and. lmax > 5) ios = 999
 call read_error(ios, inline)

 do l1 = 1, lmax + 1
   call next_line(unit, line, inline)
   read (line, *, iostat=ios) lt, rc(l1), ep(l1), ncon(l1), nbas(l1), qcut(l1)
   if (ios == 0 .and. lt /= l1 - 1) ios = 999
   call read_error(ios, inline)
 end do

! local potential
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) lloc, lpopt, rc(5), dvloc0
 call read_error(ios, inline)

! Vanderbilt-Kleinman-Bylander projectors
 do l1 = 1, lmax + 1
   call next_line(unit, line, inline)
   read (line, *, iostat=ios) lt, nproj(l1), debl(l1)
   if (ios == 0 .and. lt /= l1 - 1) ios = 999
   call read_error(ios, inline)
 end do

! model core charge
 call next_line(unit, line, inline)
! icmod=0,1,2,4[default grid] : icmod, fcfact
 read (line, *, iostat=ios) icmod, fcfact
! icmod==3 : icmod, fcfact, rcfact
 if (ios == 0 .and. icmod == 3) then
   read (line, *, iostat=ios) icmod, fcfact, rcfact
 end if
! icmod==6 : icmod, fcfact (cutoff, Ha; <= 0 for -fcfact, or 1 if 0, times the
!   cutoff of the pseudo-waves), rcfact (xc energy error, Ha; default 1e-4)
 if (ios == 0 .and. icmod == 6) then
   read (line, *, iostat=ios) icmod, fcfact, rcfact
   if (ios /= 0) then
     read (line, *, iostat=ios) icmod, fcfact
     rcfact = 1.d-4
   end if
 end if
! icmod==4[explicit grid]: icmod, fcfact(ignored), rcfact(ignored),
!   fcfact_min, fcfact_max, fcfact_step,
!   rcfact_min, rcfact_max, rcfact_step
 if (ios == 0 .and. (icmod == 4 .or. icmod == 5)) then
   read (line, *, iostat=ios) icmod, fcfact, rcfact, &
&     fcfact_min, fcfact_max, fcfact_step, &
&     rcfact_min, rcfact_max, rcfact_step
! icmod==4[implicit grid]: *_min, *_max, *_step not provided, fall back to defaults
   if (ios /= 0) then
     read (line, *, iostat=ios) icmod, fcfact
     rcfact = 0.d0
     fcfact_min = 1.5d0
     fcfact_max = 6.0d0
     fcfact_step = 0.5d0
     rcfact_min = 1.0d0
     rcfact_max = 1.9d0
     rcfact_step = 0.1d0
   end if
 end if
 call read_error(ios, inline)

! log derivative analysis
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) epsh1, epsh2, depsh, rxpsh
 if (ios /= 0) then
   read (line, *, iostat=ios) epsh1, epsh2, depsh
   rxpsh = -1.d0  ! negative value for auto setting of radius
 end if
 call read_error(ios, inline)

! output grid
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) rlmax, drl
 call read_error(ios, inline)

! test configurations
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) ncnf
 if (ios == 0 .and. ncnf > 4) ios = 999
 call read_error(ios, inline)

 do jj = 2, ncnf + 1
   call next_line(unit, line, inline)
   read (line, *, iostat=ios) nvcnf(jj)
   if (ios == 0 .and. (nvcnf(jj) < 0 .or. nc + nvcnf(jj) > 30)) ios = 999
   call read_error(ios, inline)
   do ii = nc + 1, nc + nvcnf(jj)
     call next_line(unit, line, inline)
     read (line, *, iostat=ios) nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
     call read_error(ios, inline)
   end do
 end do

! log mesh (optional, the last line): amesh, rr1, rmax
 call next_line(unit, line, inline, eof)
 if (.not. eof) then
   read (line, *, iostat=ios) amesh, rr1, rmax
   call read_error(ios, inline)
 end if

! end of reading input data
 end subroutine read_input_text

 subroutine read_input_text_r(unit, inline, &
&                             atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
&                             lmax, rc, ep, ncon, nbas, qcut, &
&                             lloc, lpopt, dvloc0, nproj, debl, &
&                             icmod, fcfact, rcfact, fcfact_min, fcfact_max, fcfact_step, &
&                             rcfact_min, rcfact_max, rcfact_step, &
&                             epsh1, epsh2, depsh, rxpsh, &
&                             rlmax, drl, &
&                             ncnf, nvcnf, nacnf, lacnf, facnf, &
&                             amesh, rr1, rmax)
!> File unit number for input
 integer, intent(in) :: unit
!> Current line number in input file
 integer, intent(out) :: inline
!> Atomic symbol
 character(len=2), intent(out) :: atsym
!> Atomic number
 real(dp), intent(out) :: zz
!> Number of core electrons
 integer, intent(out) :: nc
!> Number of valence electrons
 integer, intent(out) :: nv
!> Exchange-correlation functional code
 integer, intent(out) :: iexc
!> Pseudopotential output file format flag
 character(len=4), intent(out) :: psfile
!> Principal quantum number array
 integer, intent(out) :: na(30)
!> Orbital angular momentum quantum number array
 integer, intent(out) :: la(30)
!> Occupation number array
 real(dp), intent(out) :: fa(30)
!> Maximum angular momentum
 integer, intent(out) :: lmax
!> Pseudopotential cutoff radii
 real(dp), intent(out) :: rc(6)
!> Pseudopotential energies
 real(dp), intent(out) :: ep(6, 2)
!> Number of matching constraints for each pseudopotential
 integer, intent(out) :: ncon(6)
!> Number of basis functions for each pseudopotential
 integer, intent(out) :: nbas(6)
!> Maximum wave number for each pseudopotential
 real(dp), intent(out) :: qcut(6)
!> Angular momentum used for local potential
 integer, intent(out) :: lloc
!>
 integer, intent(out) :: lpopt
!> Local potential offset at origin
 real(dp), intent(out) :: dvloc0
!> Number of projectors for each angular momentum
 integer, intent(out) :: nproj(6)
!> Energy shift for unbound states
 real(dp), intent(out) :: debl(6)
!> Model core charge flag
 integer, intent(out) :: icmod
!> Scaling factor for core charge
 real(dp), intent(out) :: fcfact
!> Scaling factor range minimum for optimization
 real(dp), intent(out) :: fcfact_min
!> Scaling factor range maximum for optimization
 real(dp), intent(out) :: fcfact_max
!> Scaling factor step size for optimization
 real(dp), intent(out) :: fcfact_step
!> Core charge width factor
 real(dp), intent(out) :: rcfact
!> Core charge width factor minimum for optimization
 real(dp), intent(out) :: rcfact_min
!> Core charge width factor maximum for optimization
 real(dp), intent(out) :: rcfact_max
!> Core charge width factor step size for optimization
 real(dp), intent(out) :: rcfact_step
!> Lower bound of energy range for log derivative analysis
 real(dp), intent(out) :: epsh1
!> Upper bound of energy range for log derivative analysis
 real(dp), intent(out) :: epsh2
!> Energy step size for log derivative analysis
 real(dp), intent(out) :: depsh
!> Radius for log derivative analysis
 real(dp), intent(out) :: rxpsh
!> Maximum radius for output grid
 real(dp), intent(out) :: rlmax
!> Grid spacing for output grid
 real(dp), intent(out) :: drl
!> Number of test configurations
 integer, intent(out) :: ncnf
!> Number of valence states in test configurations
 integer, intent(out) :: nvcnf(5)
!> Principal quantum number array for test configurations
 integer, intent(out) :: nacnf(30, 5)
!> Orbital angular momentum quantum number array for test configurations
 integer, intent(out) :: lacnf(30, 5)
!> Occupation number array for test configurations
 real(dp), intent(out) :: facnf(30, 5)
! [Log mesh]
!> Ratio of successive points of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: amesh
!> First point of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: rr1
!> Outer radius of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: rmax

!Local variables
!> I/O status variable
 integer :: ios
!> Loop indices
 integer :: ii, jj
!> Angular momentum index (l+1)
 integer :: l1
!> Temporary variable for reading angular momentum
 integer :: lt
!> The input line being read
 character(len=:), allocatable :: line
!> Whether the input has ended
 logical :: eof

 nproj(:) = 0
 fcfact = 0.d0
 rcfact = 0.d0
 rc(:) = 0.d0
 ep(:, :) = 0.d0
 rxpsh = -1.d0  ! negative value for auto setting of radius
! default values for icmod=4 optimization ranges
 fcfact_min = 1.5d0
 fcfact_max = 6.0d0
 fcfact_step = 0.5d0
 rcfact_min = 1.0d0
 rcfact_max = 1.9d0
 rcfact_step = 0.1d0

 amesh = -1.d0  ! negative values for the built-in mesh
 rr1 = -1.d0
 rmax = -1.d0

! read input data
 inline = 0

! atom and reference configuration
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) atsym, zz, nc, nv, iexc, psfile
! counts beyond the array sizes would overrun before check_data runs
 if (ios == 0 .and. (nc < 0 .or. nv < 0 .or. nc + nv > 30)) ios = 999
 call read_error(ios, inline)

 do ii = 1, nc + nv
   call next_line(unit, line, inline)
   read (line, *, iostat=ios) na(ii), la(ii), fa(ii)
   call read_error(ios, inline)
 end do

! pseudopotential and optimization
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) lmax
 if (ios == 0 .and. lmax > 5) ios = 999
 call read_error(ios, inline)

 do l1 = 1, lmax + 1
   call next_line(unit, line, inline)
   read (line, *, iostat=ios) lt, rc(l1), ep(l1, 1), ncon(l1), nbas(l1), qcut(l1)
   if (ios == 0 .and. lt /= l1 - 1) ios = 999
   call read_error(ios, inline)
   ep(l1, 2) = ep(l1, 1)
 end do

! local potential
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) lloc, lpopt, rc(5), dvloc0
 call read_error(ios, inline)

! Vanderbilt-Kleinman-Bylander projectors
 do l1 = 1, lmax + 1
   call next_line(unit, line, inline)
   read (line, *, iostat=ios) lt, nproj(l1), debl(l1)
   if (ios == 0 .and. lt /= l1 - 1) ios = 999
   call read_error(ios, inline)
 end do

! model core charge
 call next_line(unit, line, inline)
! icmod=0,1,2,4[default grid] : icmod, fcfact
 read (line, *, iostat=ios) icmod, fcfact
! icmod==3 : icmod, fcfact, rcfact
 if (ios == 0 .and. icmod == 3) then
   read (line, *, iostat=ios) icmod, fcfact, rcfact
 end if
! icmod==6 : icmod, fcfact (cutoff, Ha; <= 0 for -fcfact, or 1 if 0, times the
!   cutoff of the pseudo-waves), rcfact (xc energy error, Ha; default 1e-4)
 if (ios == 0 .and. icmod == 6) then
   read (line, *, iostat=ios) icmod, fcfact, rcfact
   if (ios /= 0) then
     read (line, *, iostat=ios) icmod, fcfact
     rcfact = 1.d-4
   end if
 end if
! icmod==4[explicit grid]: icmod, fcfact(ignored), rcfact(ignored),
!   fcfact_min, fcfact_max, fcfact_step,
!   rcfact_min, rcfact_max, rcfact_step
 if (ios == 0 .and. (icmod == 4 .or. icmod == 5)) then
   read (line, *, iostat=ios) icmod, fcfact, rcfact, &
&     fcfact_min, fcfact_max, fcfact_step, &
&     rcfact_min, rcfact_max, rcfact_step
! icmod==4[implicit grid]: *_min, *_max, *_step not provided, fall back to defaults
   if (ios /= 0) then
     read (line, *, iostat=ios) icmod, fcfact
     rcfact = 0.d0
     fcfact_min = 1.5d0
     fcfact_max = 6.0d0
     fcfact_step = 0.5d0
     rcfact_min = 1.0d0
     rcfact_max = 1.9d0
     rcfact_step = 0.1d0
   end if
 end if
 call read_error(ios, inline)

! log derivative analysis
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) epsh1, epsh2, depsh, rxpsh
 if (ios /= 0) then
   read (line, *, iostat=ios) epsh1, epsh2, depsh
   rxpsh = -1.d0  ! negative value for auto setting of radius
 end if
 call read_error(ios, inline)

! output grid
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) rlmax, drl
 call read_error(ios, inline)

! test configurations
 call next_line(unit, line, inline)
 read (line, *, iostat=ios) ncnf
 if (ios == 0 .and. ncnf > 4) ios = 999
 call read_error(ios, inline)

 do jj = 2, ncnf + 1
   call next_line(unit, line, inline)
   read (line, *, iostat=ios) nvcnf(jj)
   if (ios == 0 .and. (nvcnf(jj) < 0 .or. nc + nvcnf(jj) > 30)) ios = 999
   call read_error(ios, inline)
   do ii = nc + 1, nc + nvcnf(jj)
     call next_line(unit, line, inline)
     read (line, *, iostat=ios) nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
     call read_error(ios, inline)
   end do
 end do

! log mesh (optional, the last line): amesh, rr1, rmax
 call next_line(unit, line, inline, eof)
 if (.not. eof) then
   read (line, *, iostat=ios) amesh, rr1, rmax
   call read_error(ios, inline)
 end if

! end of reading input data
 end subroutine read_input_text_r

 subroutine next_line(unit, line, inline, eof)
! reads into line the next line of the input that is neither blank nor a comment (one whose
! first nonblank character is #), and sets inline to its line number. The data are then read
! from line rather than from unit, so the input is never backspaced, and can come from a pipe.
! The end of the input is an error, unless eof is present: it is then set instead, and inline
! left at the line it had on entry

!Input variables
 integer, intent(in) :: unit

!Output variables
 character(len=:), allocatable, intent(out) :: line

!Input/output variables
 integer, intent(inout) :: inline

!Optional output variable
 logical, intent(out), optional :: eof

!Local variables
 character(len=256) :: chunk
 integer :: ios, nread, first, inline0

 if (present(eof)) eof = .false.
 inline0 = inline
 do
! a line of any length, read in chunks
   line = ''
   do
     read (unit, '(a)', advance='no', size=nread, iostat=ios) chunk
     line = line//chunk(1:nread)
     if (ios /= 0) exit
   end do
   inline = inline + 1
   if (present(eof) .and. is_iostat_end(ios) .and. len(line) == 0) then
     eof = .true.
     inline = inline0
     return
   end if
   if (.not. is_iostat_eor(ios)) call read_error(ios, inline)
   first = verify(line, ' '//achar(9))
   if (first > 0) then
     if (line(first:first) /= '#') exit
   end if
 end do

 return
 end subroutine next_line

 subroutine read_error(ios, inline)
! report data read error and stop

!Input variables
 integer, intent(in) :: ios
 integer, intent(inout) :: inline

 if (ios /= 0) then
   write (6, '(a,i4)') 'Read ERROR, input data file line', inline
   write (6, '(a)') 'Program will stop'
   error stop 1
 end if

 return
 end subroutine read_error

 end module input_text_m
