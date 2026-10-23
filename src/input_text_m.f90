 module input_text_m
 use, intrinsic :: iso_fortran_env, only: dp => real64, int64
 use constants_m, only: mxstate, mxcnf, mxlin, mxkap, well_onset_default
 implicit none
 private
 public :: read_input_text, next_line, finite
 contains

 subroutine read_input_text(unit, inline, &
&                           atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
&                           lmax, rc, ep, ncon, nbas, qcut, &
&                           lloc, lpopt, dvloc0, nproj, debl, &
&                           icmod, fcfact, rcfact, slack, budget_weighted, &
&                           fcfact_min, fcfact_max, fcfact_step, &
&                           rcfact_min, rcfact_max, rcfact_step, &
&                           epsh1, epsh2, depsh, rxpsh, &
&                           rlmax, drl, upf_mesh, &
&                           ncnf, nvcnf, nacnf, lacnf, facnf, &
&                           amesh, rr1, rmax, c_light, ep_as_given, &
&                           well_onset, well_tail, wreg, well_n)

!! reads the text input file (see doc/32_Ge_annotated.dat), stopping on a
!! read error

!Input variables
!> File unit number for input
 integer, intent(in) :: unit

!Output variables
!> Line number in the input file of the last line read
 integer, intent(out) :: inline
! [Atom and reference configuration]
!> Atomic symbol
 character(len=2), intent(out) :: atsym
!> Atomic number
 real(dp), intent(out) :: zz
!> Number of core states
 integer, intent(out) :: nc
!> Number of valence states
 integer, intent(out) :: nv
!> Exchange-correlation functional code
 integer, intent(out) :: iexc
!> Pseudopotential output file format flag
 character(len=4), intent(out) :: psfile
!> Principal quantum number array
 integer, intent(out) :: na(mxstate)
!> Orbital angular momentum quantum number array
 integer, intent(out) :: la(mxstate)
!> Occupation number array
 real(dp), intent(out) :: fa(mxstate)
! [Pseudopotential and optimization]
!> Maximum angular momentum
 integer, intent(out) :: lmax
!> Whether the well states take ep as given (wellstate; optional on the
!> lmax line, 0 or 1)
 logical, intent(out) :: ep_as_given
!> Onset of the well of the well states as a multiple of rc (optional on
!> the lmax line, 0 or less for the default)
 real(dp), intent(out) :: well_onset
!> The tail of the well states (optional on the lmax line): 0 "well",
!> Hamann's well; 1 "free", the free tail (free_tail_m); 2 "confined",
!> a wall that rises without bound, with the state's n as given
 integer, intent(out) :: well_tail
 real(dp), intent(out) :: wreg
!> Core radii for pseudopotentials
 real(dp), intent(out) :: rc(mxlin+1)
!> Pseudopotential energies, the input value for both j=l+1/2 and l-1/2
!> (only the fully-relativistic code uses the second)
 real(dp), intent(out) :: ep(mxlin+1, mxkap)
!> Number of matching constraints for each pseudopotential
 integer, intent(out) :: ncon(mxlin+1)
!> Number of basis functions for each pseudopotential
 integer, intent(out) :: nbas(mxlin+1)
!> Residual-energy cutoff wave vector for each pseudopotential
 real(dp), intent(out) :: qcut(mxlin+1)
! [Local potential]
!> Angular momentum used for local potential
 integer, intent(out) :: lloc
!> Type of polynomial continuation of the lloc=4 local potential
 integer, intent(out) :: lpopt
!> Local potential offset at origin
 real(dp), intent(out) :: dvloc0
! [Vanderbilt-Kleinman-Bylander projectors]
!> Number of projectors for each angular momentum
 integer, intent(out) :: nproj(mxlin+1)
!> Energy added to ep for the second projector (Ha)
 real(dp), intent(out) :: debl(mxlin+1)
!> n of the first well state of each l with well_tail 2 (optional on the
!> projector line, 0 for the default)
 integer, intent(out) :: well_n(mxlin+1)
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
!> Share of the cutoff the model core must fit (icmod 6, 7)
 real(dp), intent(out) :: slack
!> Whether the cutoff of the pseudo-waves is that of their weighted mean
!> rather than the largest (icmod 6, 7; optional after slack, 0 "max", 1
!> "weighted")
 logical, intent(out) :: budget_weighted
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
!> Most points of the log mesh of the UPF file (0 for the linear output grid)
 integer, intent(out) :: upf_mesh
! [Test configurations]
!> Number of test configurations
 integer, intent(out) :: ncnf
!> Number of valence states in test configurations
 integer, intent(out) :: nvcnf(mxcnf)
!> Principal quantum number array for test configurations
 integer, intent(out) :: nacnf(mxstate, mxcnf)
!> Orbital angular momentum quantum number array for test configurations
 integer, intent(out) :: lacnf(mxstate, mxcnf)
!> Occupation number array for test configurations
 real(dp), intent(out) :: facnf(mxstate, mxcnf)
! [Log mesh]
!> Ratio of successive points of the log radial mesh
!> (0 or less for the default)
 real(dp), intent(out) :: amesh
!> First point of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: rr1
!> Outer radius of the log radial mesh (0 or less for the default)
 real(dp), intent(out) :: rmax
! [Atom, optional]
!> Speed of light (0 or less for the default)
 real(dp), intent(out) :: c_light

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
!> The optional values of the lmax line
 integer :: igiven, itail
!> The optional budget_cutoff of the icmod line
 integer :: ibudget

 nproj(:) = 0
 fcfact = 0.d0
 rcfact = 0.d0
 slack = 1.d0
 budget_weighted = .false.
 rc(:) = 0.d0
 ep(:, :) = 0.d0
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
 c_light = -1.d0  ! negative value for the default

! read input data
 inline = 0

! atom and reference configuration, and optionally the speed of light
 call next_line(unit, line, inline)
 read(line, *, iostat=ios) atsym, zz, nc, nv, iexc, psfile, c_light
! a line of six values ends before c_light; a bad value is an error
 if(ios < 0) then
   c_light = -1.d0
   read(line, *, iostat=ios) atsym, zz, nc, nv, iexc, psfile
 end if
! counts beyond the array sizes would overrun before check_data runs
 if(ios == 0 .and. (nc < 0 .or. nv < 0 .or. nc + nv > mxstate)) ios = 999
 call check_finite(ios, [zz, c_light])
 call read_error(ios, inline)

 do ii = 1, nc + nv
   call next_line(unit, line, inline)
   read(line, *, iostat=ios) na(ii), la(ii), fa(ii)
   call check_finite(ios, [fa(ii)])
   call read_error(ios, inline)
 end do

! pseudopotential and optimization: lmax, and optionally ep_as_given (0 or
! 1), well_onset (0 or less for the default), well_tail (0 for "well",
! 1 for "free", 2 for "confined") and wave_regularization (at least 0); a
! line that ends early leaves the rest at their defaults, a bad value is
! an error
 call next_line(unit, line, inline)
 itail = 0
 wreg = 0.d0
 read(line, *, iostat=ios) lmax, igiven, well_onset, itail, wreg
 if(ios < 0) then
   wreg = 0.d0
   read(line, *, iostat=ios) lmax, igiven, well_onset, itail
 end if
 if(ios < 0) then
   itail = 0
   read(line, *, iostat=ios) lmax, igiven, well_onset
 end if
 if(ios < 0) then
   well_onset = 0.d0
   read(line, *, iostat=ios) lmax, igiven
 end if
 if(ios < 0) then
   igiven = 0
   read(line, *, iostat=ios) lmax
 end if
 if(ios == 0) then
   if(lmax > mxlin .or. igiven < 0 .or. igiven > 1 .or. itail < 0 &
&     .or. itail > 2 .or. (well_onset > 0.d0 .and. well_onset < 1.d0) &
&     .or. .not. (wreg >= 0.d0)) then
     ios = 999
   end if
 end if
 call check_finite(ios, [well_onset, wreg])
 call read_error(ios, inline)
 ep_as_given = igiven == 1
 well_tail = itail
 if(well_onset <= 0.d0) well_onset = well_onset_default

 do l1 = 1, lmax + 1
   call next_line(unit, line, inline)
   read(line, *, iostat=ios) lt, rc(l1), ep(l1, 1), ncon(l1), nbas(l1), &
&    qcut(l1)
   if(ios == 0 .and. lt /= l1 - 1) ios = 999
   call check_finite(ios, [rc(l1), ep(l1, 1), qcut(l1)])
   call read_error(ios, inline)
   ep(l1, 2) = ep(l1, 1)
 end do

! local potential
 call next_line(unit, line, inline)
 read(line, *, iostat=ios) lloc, lpopt, rc(5), dvloc0
 call check_finite(ios, [rc(5), dvloc0])
 call read_error(ios, inline)

! Vanderbilt-Kleinman-Bylander projectors: l, nproj, debl, and optionally
! well_n (0 for the default, at least 0)
 well_n(:) = 0
 do l1 = 1, lmax + 1
   call next_line(unit, line, inline)
   read(line, *, iostat=ios) lt, nproj(l1), debl(l1), well_n(l1)
   if(ios < 0) then
     well_n(l1) = 0
     read(line, *, iostat=ios) lt, nproj(l1), debl(l1)
   end if
   if(ios == 0 .and. (lt /= l1 - 1 .or. well_n(l1) < 0)) ios = 999
   call check_finite(ios, [debl(l1)])
   call read_error(ios, inline)
 end do

! model core charge
 call next_line(unit, line, inline)
! icmod=0,1,2,4[default grid] : icmod, fcfact
 read(line, *, iostat=ios) icmod, fcfact
! icmod==3 : icmod, fcfact, rcfact
 if(ios == 0 .and. icmod == 3) then
   read(line, *, iostat=ios) icmod, fcfact, rcfact
 end if
! icmod==6, 7 : icmod, fcfact (cutoff, Ha; <= 0 for -fcfact, or 1 if 0,
!   times the cutoff of the pseudo-waves), rcfact (xc energy error, Ha;
!   default 1e-5), slack (share of the cutoff the model must fit; default 1),
!   budget_cutoff (0 for the largest cutoff of the pseudo-waves, 1 for that
!   of their weighted mean; default 0, a bad value is an error)
 if(ios == 0 .and. (icmod == 6 .or. icmod == 7)) then
   ibudget = 0
   read(line, *, iostat=ios) icmod, fcfact, rcfact, slack, ibudget
   if(ios < 0) then
     ibudget = 0
     read(line, *, iostat=ios) icmod, fcfact, rcfact, slack
     if(ios /= 0) then
       slack = 1.d0
       read(line, *, iostat=ios) icmod, fcfact, rcfact
     end if
     if(ios /= 0) then
       read(line, *, iostat=ios) icmod, fcfact
       rcfact = 1.d-5
     end if
   end if
   if(ios == 0 .and. (ibudget < 0 .or. ibudget > 1)) ios = 999
   budget_weighted = ibudget == 1
 end if
! icmod==4[explicit grid]: icmod, fcfact(ignored), rcfact(ignored),
!   fcfact_min, fcfact_max, fcfact_step,
!   rcfact_min, rcfact_max, rcfact_step
 if(ios == 0 .and. (icmod == 4 .or. icmod == 5)) then
   read(line, *, iostat=ios) icmod, fcfact, rcfact, &
&     fcfact_min, fcfact_max, fcfact_step, &
&     rcfact_min, rcfact_max, rcfact_step
! icmod==4[implicit grid]: *_min, *_max, *_step not provided, fall back to
!   defaults
   if(ios /= 0) then
     read(line, *, iostat=ios) icmod, fcfact
     rcfact = 0.d0
     fcfact_min = 1.5d0
     fcfact_max = 6.0d0
     fcfact_step = 0.5d0
     rcfact_min = 1.0d0
     rcfact_max = 1.9d0
     rcfact_step = 0.1d0
   end if
 end if
 call check_finite(ios, [fcfact, rcfact, slack, fcfact_min, fcfact_max, &
&  fcfact_step, rcfact_min, rcfact_max, rcfact_step])
 call read_error(ios, inline)

! log derivative analysis
 call next_line(unit, line, inline)
 read(line, *, iostat=ios) epsh1, epsh2, depsh, rxpsh
 if(ios /= 0) then
   read(line, *, iostat=ios) epsh1, epsh2, depsh
   rxpsh = -1.d0  ! negative value for auto setting of radius
 end if
 call check_finite(ios, [epsh1, epsh2, depsh, rxpsh])
 call read_error(ios, inline)

! output grid: rlmax, drl, and optionally upf_mesh
 call next_line(unit, line, inline)
 read(line, *, iostat=ios) rlmax, drl, upf_mesh
! a line of two values ends before upf_mesh; a bad value is an error
 if(ios < 0) then
   read(line, *, iostat=ios) rlmax, drl
   upf_mesh = 0  ! the UPF file on the linear output grid
 end if
 call check_finite(ios, [rlmax, drl])
 call read_error(ios, inline)

! test configurations
 call next_line(unit, line, inline)
 read(line, *, iostat=ios) ncnf
 if(ios == 0 .and. ncnf > mxcnf - 1) ios = 999
 call read_error(ios, inline)

 do jj = 2, ncnf + 1
   call next_line(unit, line, inline)
   read(line, *, iostat=ios) nvcnf(jj)
   if(ios == 0 .and. (nvcnf(jj) < 0 .or. nc + nvcnf(jj) > mxstate)) ios = 999
   call read_error(ios, inline)
   do ii = nc + 1, nc + nvcnf(jj)
     call next_line(unit, line, inline)
     read(line, *, iostat=ios) nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
     call check_finite(ios, [facnf(ii, jj)])
     call read_error(ios, inline)
   end do
 end do

! log mesh (optional, the last line): amesh, rr1, rmax
 call next_line(unit, line, inline, eof)
 if(.not. eof) then
   read(line, *, iostat=ios) amesh, rr1, rmax
   call check_finite(ios, [amesh, rr1, rmax])
   call read_error(ios, inline)
 end if
 return

! end of reading input data
 end subroutine read_input_text

 subroutine next_line(unit, line, inline, eof)
!! reads into line the next line of the input that is neither blank nor a
!! comment (one whose first nonblank character is #), and sets inline to its
!! line number. The data are then read from line rather than from unit, so
!! the input is never backspaced, and can come from a pipe. The end of the
!! input is an error, unless eof is present: it is then set instead, and
!! inline left at the line it had on entry
!!
!! @param unit  input file unit
!! @param line  the line read, of any length
!! @param inline  line number of the last line read, updated on return
!! @param eof  optional, .true. if the input ended before another data line

!Input variables
 integer, intent(in) :: unit

!Input/output variables
 integer, intent(inout) :: inline

!Output variables
 character(len=:), allocatable, intent(out) :: line
 logical, intent(out), optional :: eof

!Local variables
 character(len=256) :: chunk
 integer :: ios, nread, first, inline0

 if(present(eof)) eof = .false.
 inline0 = inline
 do
! a line of any length, read in chunks
   line = ''
   do
     read(unit, '(a)', advance='no', size=nread, iostat=ios) chunk
     line = line//chunk(1:nread)
     if(ios /= 0) exit
   end do
   inline = inline + 1
   if(present(eof) .and. is_iostat_end(ios) .and. len(line) == 0) then
     eof = .true.
     inline = inline0
     return
   end if
   if(.not. is_iostat_eor(ios)) call read_error(ios, inline)
   first = verify(line, ' '//achar(9))
   if(first > 0) then
     if(line(first:first) /= '#') exit
   end if
 end do

 return
 end subroutine next_line

 subroutine check_finite(ios, vals)
!! makes a read that succeeded (ios == 0) an error if it gave a value that
!! is not finite: list-directed input takes NaN and Inf as reals
!!
!! @param ios  iostat of the read, set to 999 for a value not finite
!! @param vals  the reals the read gave

!Input variables
 real(dp), intent(in) :: vals(:)

!Input/output variables
 integer, intent(inout) :: ios

 if(ios == 0 .and. .not. all(finite(vals))) ios = 999

 return
 end subroutine check_finite

 elemental function finite(xx) result(ok)
!! whether xx is finite (neither NaN nor Inf), from its exponent bits:
!! under -ffast-math (-ffinite-math-only) a compiler may take ieee_is_nan
!! and x /= x to be false
!!
!! @param xx  the value

!Input variables
 real(dp), intent(in) :: xx

!Output variables
 logical :: ok

 ok = iand(shiftr(transfer(xx, 0_int64), 52), 2047_int64) /= 2047_int64

 return
 end function finite

 subroutine read_error(ios, inline)
!! report data read error and stop
!!
!! @param ios  iostat of the read, non-zero for an error
!! @param inline  line number of the input file line being read

!Input variables
 integer, intent(in) :: ios

!Input/output variables
 integer, intent(inout) :: inline

 if(ios /= 0) then
   write(6, '(a,i4)') 'Read ERROR, input data file line', inline
   write(6, '(a)') 'Program will stop'
   error stop 1
 end if

 return
 end subroutine read_error

 end module input_text_m
