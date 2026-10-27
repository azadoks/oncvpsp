 module input_text_m
 use, intrinsic :: iso_fortran_env, only: dp => real64, int64
 use constants_m, only: mxstate, mxcnf, mxlin, mxkap, well_onset_default
 use input_m, only: input_t
 implicit none
 private
 public :: read_input_text, next_line, finite
 contains

 subroutine read_input_text(unit, inline, inp, iout, ifail)

!! reads the text input file (see doc/32_Ge_annotated.dat), failing on a
!! read error

!Input variables
!> File unit number for input
 integer, intent(in) :: unit

!Output variables
!> Line number in the input file of the last line read
 integer, intent(out) :: inline
!> The input values
 type(input_t), intent(out) :: inp
!> unit of the text output
 integer, intent(in) :: iout
!> 0, or 1 after a failure, whose message is written
 integer, intent(out) :: ifail

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

 ifail = 0
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

   nproj(:) = 0
   fcfact = 0.d0
   rcfact = 0.d0
   slack = 1.d0
   budget_cut = 0
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
   call next_line(unit, line, inline, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   read(line, *, iostat=ios) atsym, zz, nc, nv, iexc, psfile, c_light
! a line of six values ends before c_light; a bad value is an error
   if(ios < 0) then
     c_light = -1.d0
     read(line, *, iostat=ios) atsym, zz, nc, nv, iexc, psfile
   end if
! counts beyond the array sizes would overrun before check_data runs
   if(ios == 0 .and. (nc < 0 .or. nv < 0 .or. nc + nv > mxstate)) ios = 999
   call check_finite(ios, [zz, c_light])
   call read_error(ios, inline, iout, ifail)
   if(ifail /= 0) return

   do ii = 1, nc + nv
     call next_line(unit, line, inline, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     read(line, *, iostat=ios) na(ii), la(ii), fa(ii)
     call check_finite(ios, [fa(ii)])
     call read_error(ios, inline, iout, ifail)
     if(ifail /= 0) return
   end do

! pseudopotential and optimization: lmax, and optionally ep_as_given (0 or
! 1), well_onset (0 or less for the default), well_tail (0 for "well",
! 1 for "free", 2 for "confined"), wave_regularization (at least 0) and
! well_norm (0, 1 or 2); a line that ends early leaves the rest at their
! defaults, a bad value is an error
   call next_line(unit, line, inline, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   itail = 0
   wreg = 0.d0
   well_norm = 0
   read(line, *, iostat=ios) lmax, igiven, well_onset, itail, wreg, well_norm
   if(ios < 0) then
     well_norm = 0
     read(line, *, iostat=ios) lmax, igiven, well_onset, itail, wreg
   end if
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
&       .or. itail > 2 .or. (well_onset > 0.d0 .and. well_onset < 1.d0) &
&       .or. .not. (wreg >= 0.d0) .or. well_norm < 0 .or. well_norm > 2) then
       ios = 999
     end if
   end if
   call check_finite(ios, [well_onset, wreg])
   call read_error(ios, inline, iout, ifail)
   if(ifail /= 0) return
   ep_as_given = igiven == 1
   well_tail = itail
   if(well_onset <= 0.d0) well_onset = well_onset_default

   do l1 = 1, lmax + 1
     call next_line(unit, line, inline, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     read(line, *, iostat=ios) lt, rc(l1), ep(l1, 1), ncon(l1), nbas(l1), &
&      qcut(l1)
     if(ios == 0 .and. lt /= l1 - 1) ios = 999
     call check_finite(ios, [rc(l1), ep(l1, 1), qcut(l1)])
     call read_error(ios, inline, iout, ifail)
     if(ifail /= 0) return
     ep(l1, 2) = ep(l1, 1)
   end do

! local potential
   call next_line(unit, line, inline, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   read(line, *, iostat=ios) lloc, lpopt, rc(5), dvloc0
   call check_finite(ios, [rc(5), dvloc0])
   call read_error(ios, inline, iout, ifail)
   if(ifail /= 0) return

! Vanderbilt-Kleinman-Bylander projectors: l, nproj, debl, and optionally
! well_n (0 for the default, at least 0) and ecut_weight (negative for the
! default)
   well_n(:) = 0
   ecut_wt(:) = -1.d0
   do l1 = 1, lmax + 1
     call next_line(unit, line, inline, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     read(line, *, iostat=ios) lt, nproj(l1), debl(l1), well_n(l1), ecut_wt(l1)
     if(ios < 0) then
       ecut_wt(l1) = -1.d0
       read(line, *, iostat=ios) lt, nproj(l1), debl(l1), well_n(l1)
     end if
     if(ios < 0) then
       well_n(l1) = 0
       read(line, *, iostat=ios) lt, nproj(l1), debl(l1)
     end if
     if(ios == 0 .and. (lt /= l1 - 1 .or. well_n(l1) < 0)) ios = 999
     call check_finite(ios, [debl(l1), ecut_wt(l1)])
     call read_error(ios, inline, iout, ifail)
     if(ifail /= 0) return
   end do

! model core charge
   call next_line(unit, line, inline, iout=iout, ifail=ifail)
   if(ifail /= 0) return
! icmod=0,1,2,4[default grid] : icmod, fcfact
   read(line, *, iostat=ios) icmod, fcfact
! icmod==3 : icmod, fcfact, rcfact
   if(ios == 0 .and. icmod == 3) then
     read(line, *, iostat=ios) icmod, fcfact, rcfact
   end if
! icmod==6, 7 : icmod, fcfact (cutoff, Ha; <= 0 for -fcfact, or 1 if 0,
!   times the cutoff of the pseudo-waves), rcfact (xc energy error, Ha;
!   default 1e-5), slack (share of the cutoff the model must fit; default 1),
!   budget_cutoff (0 for the largest cutoff of the ground-state
!   pseudo-waves, 1 for that of their weighted mean, 2 for the largest of
!   all; default 0, a bad value is an error)
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
     if(ios == 0 .and. (ibudget < 0 .or. ibudget > 2)) ios = 999
     budget_cut = ibudget
   end if
! icmod==4[explicit grid]: icmod, fcfact(ignored), rcfact(ignored),
!   fcfact_min, fcfact_max, fcfact_step,
!   rcfact_min, rcfact_max, rcfact_step
   if(ios == 0 .and. (icmod == 4 .or. icmod == 5)) then
     read(line, *, iostat=ios) icmod, fcfact, rcfact, &
&       fcfact_min, fcfact_max, fcfact_step, &
&       rcfact_min, rcfact_max, rcfact_step
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
&    fcfact_step, rcfact_min, rcfact_max, rcfact_step])
   call read_error(ios, inline, iout, ifail)
   if(ifail /= 0) return

! log derivative analysis
   call next_line(unit, line, inline, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   read(line, *, iostat=ios) epsh1, epsh2, depsh, rxpsh
   if(ios /= 0) then
     read(line, *, iostat=ios) epsh1, epsh2, depsh
     rxpsh = -1.d0  ! negative value for auto setting of radius
   end if
   call check_finite(ios, [epsh1, epsh2, depsh, rxpsh])
   call read_error(ios, inline, iout, ifail)
   if(ifail /= 0) return

! output grid: rlmax, drl, and optionally upf_mesh
   call next_line(unit, line, inline, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   read(line, *, iostat=ios) rlmax, drl, upf_mesh
! a line of two values ends before upf_mesh; a bad value is an error
   if(ios < 0) then
     read(line, *, iostat=ios) rlmax, drl
     upf_mesh = 0  ! the UPF file on the linear output grid
   end if
   call check_finite(ios, [rlmax, drl])
   call read_error(ios, inline, iout, ifail)
   if(ifail /= 0) return

! test configurations
   call next_line(unit, line, inline, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   read(line, *, iostat=ios) ncnf
   if(ios == 0 .and. ncnf > mxcnf - 1) ios = 999
   call read_error(ios, inline, iout, ifail)
   if(ifail /= 0) return

   do jj = 2, ncnf + 1
     call next_line(unit, line, inline, iout=iout, ifail=ifail)
     if(ifail /= 0) return
     read(line, *, iostat=ios) nvcnf(jj)
     if(ios == 0 .and. (nvcnf(jj) < 0 .or. nc + nvcnf(jj) > mxstate)) ios = 999
     call read_error(ios, inline, iout, ifail)
     if(ifail /= 0) return
     do ii = nc + 1, nc + nvcnf(jj)
       call next_line(unit, line, inline, iout=iout, ifail=ifail)
       if(ifail /= 0) return
       read(line, *, iostat=ios) nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
       call check_finite(ios, [facnf(ii, jj)])
       call read_error(ios, inline, iout, ifail)
       if(ifail /= 0) return
     end do
   end do

! log mesh (optional, the last line): amesh, rr1, rmax
   call next_line(unit, line, inline, eof, iout=iout, ifail=ifail)
   if(ifail /= 0) return
   if(.not. eof) then
     read(line, *, iostat=ios) amesh, rr1, rmax
     call check_finite(ios, [amesh, rr1, rmax])
     call read_error(ios, inline, iout, ifail)
     if(ifail /= 0) return
   end if
 end associate
 return

! end of reading input data
 end subroutine read_input_text

 subroutine next_line(unit, line, inline, eof, iout, ifail)
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
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: unit

!Input/output variables
 integer, intent(inout) :: inline

!Output variables
 character(len=:), allocatable, intent(out) :: line
 logical, intent(out), optional :: eof
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 character(len=256) :: chunk
 integer :: ios, nread, first, inline0

 ifail = 0
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
   if(.not. is_iostat_eor(ios)) then
     call read_error(ios, inline, iout, ifail)
     if(ifail /= 0) return
   end if
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

 subroutine read_error(ios, inline, iout, ifail)
!! report a data read error and fail
!!
!! @param ios  iostat of the read, non-zero for an error
!! @param inline  line number of the input file line being read
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: ios

!Input/output variables
 integer, intent(inout) :: inline
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

 ifail = 0
 if(ios /= 0) then
   write(iout, '(a,i4)') 'Read ERROR, input data file line', inline
   write(iout, '(a)') 'Program will stop'
   ifail = 1
   return
 end if

 return
 end subroutine read_error

 end module input_text_m
