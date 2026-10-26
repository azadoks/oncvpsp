!
! Copyright (c) 2026 by Austin Zadoks
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!> The last line of oncvpsp_confine.x, confine_summary, from the sphere
!> (confine_ws_m) and the walls (confine_wall_m).
 module confine_summary_m
 use constants_m, only: dp
 use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
&                                         ieee_positive_inf
 implicit none
 private
 public :: write_confine_summary

!> the counts of the ws_missing, ws_extra, ws_paired_off and wall_failed
!> lines, the largest |e_ps - e_ae| and cutoffs (1e-4, 1e-5 Ha; then those
!> of the top edges alone) of the edges below 0, whether one of them does
!> not reach it, and the largest |de_ps - de_ae| of the walls; -1 for a
!> part left out
 type, public :: confine_summary_t
   integer :: nmiss=-1,nextra=-1,noff=-1,nwfail=-1
   real(dp) :: smax_de=-1.0d0,scut(4)=-1.0d0,sdde=-1.0d0
   logical :: sinf(4)=.false.
 end type confine_summary_t

 contains

 subroutine write_confine_summary(sm)
!! the header and the line of the summary
 type(confine_summary_t), intent(in) :: sm
 write(6,'(3a)') '# confine_summary max_de(Ha) ecut4(Ha) ecut5(Ha)', &
&     ' n_missing n_extra n_paired_off n_wall_failed max_dde(Ha)', &
&     ' ecut4_top(Ha) ecut5_top(Ha)'
 write(6,'(a,1p,3e24.16,4i6,3e24.16)') 'confine_summary', &
&     summ(sm%smax_de,.false.),summ(sm%scut(1),sm%sinf(1)), &
&     summ(sm%scut(2),sm%sinf(2)),sm%nmiss,sm%nextra,sm%noff,sm%nwfail, &
&     summ(sm%sdde,.false.),summ(sm%scut(3),sm%sinf(3)), &
&     summ(sm%scut(4),sm%sinf(4))
 end subroutine write_confine_summary

 real(dp) function summ(xx,inf)
!! a value of the summary: Infinity if inf, NaN if xx < 0 (none)
 real(dp), intent(in) :: xx
 logical, intent(in) :: inf
 if(inf) then
   summ=ieee_value(1.0d0,ieee_positive_inf)
 else if(xx<0.0d0) then
   summ=ieee_value(1.0d0,ieee_quiet_nan)
 else
   summ=xx
 end if
 end function summ

 end module confine_summary_m
