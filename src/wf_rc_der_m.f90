!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
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

 module wf_rc_der_m

 implicit none
 private
 public :: wf_rc_der

 contains

 subroutine wf_rc_der(rr,uu,al,rc,irc,mmax,uorder)

!! calculate derivatives of radial wave functin at core radius rc
!!
!! @param rr  log radial mesh
!! @param uu  rr*radial wavefunction
!! @param al  log of radial mesh factor
!! @param mmax  last point used on radial mesh
!! @param irc  index of radial mesh point used as rc
!! @param rc  input core radius.  If irc = 0 on input, irc is calculated and rc
!!            is exact mesh point value on output (next larger point)
!! @param uorder  value and four derivatives of uu/rr at final rc

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),uu(mmax)
 real(dp), intent(in) :: al

!Input/output variables
 real(dp), intent(inout) :: rc
 integer, intent(inout) :: irc

!Output variables
 real(dp), intent(out) :: uorder(5)

!Local variables
 integer :: nder
 integer :: ii,jj
 real(dp),allocatable :: work(:,:)

!set rc to exact mesh value if not specified by irc
 if(irc == 0) then
   do ii=1,mmax
     if(irc == 0 .and. rr(ii)>=rc) then
       irc=ii
       rc=rr(ii)
       exit
     end if
   end do
 else
   rc=rr(irc)
 end if

 allocate(work(mmax,5))

 do ii=1,mmax
   work(ii,1)=uu(ii)/rr(ii)
 end do

 do jj=2,5

! 7-point numerical first derivatives applied successively
   do ii=irc-25+3*jj,irc+25-3*jj
     work(ii,jj)=(-work(ii-3,jj-1)+ 9.d0*work(ii-2,jj-1)&
&     -45.d0*work(ii-1,jj-1)+45.d0*work(ii+1,jj-1)&
&     -9.d0*work(ii+2,jj-1)+work(ii+3,jj-1))&
&     /(60.d0*al*rr(ii))
   end do
 end do

 do jj=1,5
   uorder(jj)=work(irc,jj)
 end do

 deallocate(work)

 return
 end subroutine wf_rc_der

 end module wf_rc_der_m
