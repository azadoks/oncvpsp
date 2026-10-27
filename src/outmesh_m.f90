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
!> Functions of the log mesh on the output meshes of the psp8 and UPF files
!> (doc/smoothness_plan.md, O-1): at each output radius the interpolant
!> (lm_interp) of a smooth function, or of the piece of the radius's
!> segment of a piecewise one, 0 beyond its last bound, so no interpolant
!> spans a kink and a projector is exactly 0 beyond its rc. An output point
!> at the origin, below the log mesh, takes the limit dpnint took: the
!> polynomial of degree 7 through the next 8 output points.

 module outmesh_m
 use constants_m, only: dp
 use logmesh_m, only: lm_interp
 use pieces_m, only: piecewise_t, pw_node_values, pw_segment
 implicit none
 private
 public :: om_values, om_piece_values, om_local_values

!> degree of the extrapolation to the origin
 integer, parameter :: npoly=7

 contains

 subroutine om_values(rr,ff,mmax,rl,fl,nrl,iout,ifail)

!! a smooth function on the output mesh, 0 beyond its last nonzero node
!! (the end of a model core's support, phase 5)
!!
!! @param rr  log radial mesh
!! @param ff  the function on the log mesh
!! @param mmax  size of the log radial mesh
!! @param rl  output mesh, ascending, rl(1) 0 or on the log mesh
!! @param fl  the function on the output mesh
!! @param nrl  size of the output mesh
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),ff(mmax),rl(nrl)

!Output variables
 real(dp), intent(out) :: fl(nrl)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,nlast

 ifail=0
 call check_range(rr,mmax,rl,nrl,iout,ifail)
 if(ifail/=0) return
 nlast=mmax
 do while(nlast>1 .and. ff(nlast)==0.0d0)
   nlast=nlast-1
 end do
 do ii=1,nrl
   if(rl(ii)<rr(1)) cycle
   if(rl(ii)>rr(nlast)) then
     fl(ii)=0.0d0
   else
     call lm_interp(rr,ff,mmax,rl(ii),fl(ii))
   end if
 end do
 call origin(rr,mmax,rl,fl,nrl,iout,ifail)
 if(ifail/=0) return
 return
 end subroutine om_values

 subroutine om_piece_values(rr,mmax,pw,rl,fl,nrl,iout,ifail)

!! a piecewise function on the output mesh: at each radius the interpolant
!! of the piece of its segment, 0 beyond the last bound
!!
!! @param rr  log radial mesh
!! @param mmax  size of the log radial mesh
!! @param pw  the function
!! @param rl  output mesh, ascending, rl(1) 0 or on the log mesh
!! @param fl  the function on the output mesh
!! @param nrl  size of the output mesh
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),rl(nrl)
 type(piecewise_t), intent(in) :: pw

!Output variables
 real(dp), intent(out) :: fl(nrl)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,kk

 ifail=0
 call check_range(rr,mmax,rl,nrl,iout,ifail)
 if(ifail/=0) return
 do ii=1,nrl
   if(rl(ii)<rr(1)) cycle
   kk=pw_segment(pw,rl(ii))
   if(kk>0) then
     call lm_interp(rr,pw%vals(:,kk),mmax,rl(ii),fl(ii))
   else
     fl(ii)=0.0d0
   end if
 end do
 call origin(rr,mmax,rl,fl,nrl,iout,ifail)
 if(ifail/=0) return
 return
 end subroutine om_piece_values

 subroutine om_local_values(rr,vv,mmax,vlocpw,rl,fl,nrl,iout,ifail)

!! a local potential on the output mesh: vv less the node values of the
!! pieces vlocpw is smooth (the screening and the ionic tail), and is added
!! to the piece of each radius's segment (as in outwf_setup)
!!
!! @param rr  log radial mesh
!! @param vv  the potential on the log mesh
!! @param mmax  size of the log radial mesh
!! @param vlocpw  the pieces, the last bound the mesh end
!! @param rl  output mesh, ascending, rl(1) 0 or on the log mesh
!! @param fl  the potential on the output mesh
!! @param nrl  size of the output mesh
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),vv(mmax),rl(nrl)
 type(piecewise_t), intent(in) :: vlocpw

!Output variables
 real(dp), intent(out) :: fl(nrl)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: vs(:)
 integer :: ii,kk

 ifail=0
 call check_range(rr,mmax,rl,nrl,iout,ifail)
 if(ifail/=0) return
 allocate(vs(mmax))
 call pw_node_values(rr,mmax,vlocpw,vs)
 vs(:)=vv(:)-vs(:)
 do ii=1,nrl
   if(rl(ii)<rr(1)) cycle
   kk=pw_segment(vlocpw,rl(ii))
   if(kk==0) then
     write(iout,'(a)') 'om_local_values: ERROR the pieces end inside the mesh'
     ifail=1
     return
   end if
   call lm_interp(rr,vlocpw%vals(:,kk)+vs(:),mmax,rl(ii),fl(ii))
 end do
 call origin(rr,mmax,rl,fl,nrl,iout,ifail)
 if(ifail/=0) return
 deallocate(vs)
 return
 end subroutine om_local_values

 subroutine check_range(rr,mmax,rl,nrl,iout,ifail)

!! the output mesh within the log mesh, but for a first point at the
!! origin

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),rl(nrl)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

 ifail=0
 if(rl(nrl)>rr(mmax)) then
   write(iout,'(a,f10.4,a,f10.4)') 'outmesh: ERROR output mesh to', &
&       rl(nrl),' bohr, beyond the log mesh, ',rr(mmax)
   ifail=1
   return
 end if
 if(nrl>1) then
   if(rl(2)<rr(1)) then
     write(iout,'(a)') &
&         'outmesh: ERROR more than one output point below the log mesh'
     ifail=1
     return
   end if
 end if
 return
 end subroutine check_range

 subroutine origin(rr,mmax,rl,fl,nrl,iout,ifail)

!! a first output point below the log mesh (the origin): the polynomial
!! of degree npoly through fl(2..npoly+2), at rl(1)

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),rl(nrl)

!Input/output variables
 real(dp), intent(inout) :: fl(nrl)

!Output variables
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: term
 integer :: iy,ip

 ifail=0
 if(rl(1)>=rr(1)) return
 if(nrl<npoly+2) then
   write(iout,'(a,i0)') &
&       'outmesh: ERROR too few points to extrapolate to the origin, ',nrl
   ifail=1
   return
 end if
 fl(1)=0.0d0
 do iy=2,npoly+2
   term=fl(iy)
   do ip=2,npoly+2
     if(ip==iy) cycle
     term=term*(rl(1)-rl(ip))/(rl(iy)-rl(ip))
   end do
   fl(1)=fl(1)+term
 end do
 return
 end subroutine origin

 end module outmesh_m
