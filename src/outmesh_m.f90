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
!> spans a kink and a projector is exactly 0 beyond its rc; a piece with a
!> closed form (the model core's) is evaluated instead. An output point
!> at the origin, below the log mesh, takes the value there of the
!> interpolant in r of the log mesh's values (lm_origin), of the piece of
!> the first segment for a piecewise function.

 module outmesh_m
 implicit none
 private
 public :: om_values, om_piece_values, om_local_values

 contains

 subroutine om_values(rr,ff,mmax,rl,fl,nrl,iout,ifail)

!! a smooth function on the output mesh, 0 beyond its last nonzero node
!! (the end of a model core's support, phase 5)
!!
!! @param rr  log radial mesh (bohr)
!! @param ff  the function on the log mesh
!! @param mmax  size of the log radial mesh
!! @param rl  output mesh (bohr), ascending, rl(1) 0 or on the log mesh
!! @param fl  the function on the output mesh
!! @param nrl  size of the output mesh
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use logmesh_m, only: lm_interp

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),ff(mmax),rl(nrl)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: fl(nrl)
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
 call origin(rr,ff,mmax,rl,fl,nrl,iout,ifail)
 if(ifail/=0) return
 return
 end subroutine om_values

 subroutine om_piece_values(rr,mmax,pw,rl,fl,nrl,iout,ifail)

!! a piecewise function on the output mesh: at each radius the interpolant
!! of the piece of its segment, 0 beyond the last bound
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param pw  the function
!! @param rl  output mesh (bohr), ascending, rl(1) 0 or on the log mesh
!! @param fl  the function on the output mesh
!! @param nrl  size of the output mesh
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use logmesh_m, only: lm_interp
 use pieces_m, only: piecewise_t, pw_segment, pw_form_value

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),rl(nrl)
 type(piecewise_t), intent(in) :: pw
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: fl(nrl)
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,kk

 ifail=0
 call check_range(rr,mmax,rl,nrl,iout,ifail)
 if(ifail/=0) return
 do ii=1,nrl
   if(rl(ii)<rr(1)) cycle
   kk=pw_segment(pw,rl(ii))
   if(pw_form_value(pw,kk,rl(ii),fl(ii))) cycle
   if(kk>0) then
     call lm_interp(rr,pw%vals(:,kk),mmax,rl(ii),fl(ii))
   else
     fl(ii)=0.0d0
   end if
 end do
 kk=pw_segment(pw,rr(1))
 if(rl(1)<rr(1)) then
   if(pw_form_value(pw,kk,0.0d0,fl(1))) return
 end if
 if(kk>0) then
   call origin(rr,pw%vals(:,kk),mmax,rl,fl,nrl,iout,ifail)
   if(ifail/=0) return
 else if(rl(1)<rr(1)) then
   fl(1)=0.0d0
 end if
 return
 end subroutine om_piece_values

 subroutine om_local_values(rr,vv,mmax,vlocpw,rl,fl,nrl,iout,ifail)

!! a local potential on the output mesh: vv less the node values of the
!! pieces vlocpw is smooth (the screening and the ionic tail), and is added
!! to the piece of each radius's segment (as in outwf_setup)
!!
!! @param rr  log radial mesh (bohr)
!! @param vv  the potential on the log mesh
!! @param mmax  size of the log radial mesh
!! @param vlocpw  the pieces, the last bound the mesh end
!! @param rl  output mesh (bohr), ascending, rl(1) 0 or on the log mesh
!! @param fl  the potential on the output mesh
!! @param nrl  size of the output mesh
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use logmesh_m, only: lm_interp
 use pieces_m, only: piecewise_t, pw_node_values, pw_segment

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),vv(mmax),rl(nrl)
 type(piecewise_t), intent(in) :: vlocpw
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: fl(nrl)
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
     write(iout,'(/a)') 'om_local_values: ERROR the pieces end inside the mesh'
     ifail=1
     return
   end if
   call lm_interp(rr,vlocpw%vals(:,kk)+vs(:),mmax,rl(ii),fl(ii))
 end do
 kk=pw_segment(vlocpw,rr(1))
 call origin(rr,vlocpw%vals(:,kk)+vs(:),mmax,rl,fl,nrl,iout,ifail)
 if(ifail/=0) return
 deallocate(vs)
 return
 end subroutine om_local_values

 subroutine check_range(rr,mmax,rl,nrl,iout,ifail)

!! the output mesh within the log mesh, but for a first point at the
!! origin
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param rl  output mesh (bohr)
!! @param nrl  size of the output mesh
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),rl(nrl)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

 ifail=0
 if(rl(nrl)>rr(mmax)) then
   write(iout,'(/a,f12.6,a,f12.6)') 'check_range: ERROR output mesh to', &
&       rl(nrl),' bohr, beyond the log mesh, ',rr(mmax)
   ifail=1
   return
 end if
 if(nrl>1) then
   if(rl(2)<rr(1)) then
     write(iout,'(/a)') &
&         'check_range: ERROR more than one output point below the log mesh'
     ifail=1
     return
   end if
 end if
 return
 end subroutine check_range

 subroutine origin(rr,ff,mmax,rl,fl,nrl,iout,ifail)

!! a first output point below the log mesh (the origin): the value at
!! r = 0 of the interpolant in r of the function on the log mesh
!! (lm_origin)
!!
!! @param rr  log radial mesh (bohr)
!! @param ff  the function on the log mesh
!! @param mmax  size of the log radial mesh
!! @param rl  output mesh (bohr)
!! @param fl  the function on the output mesh, fl(1) set
!! @param nrl  size of the output mesh
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use logmesh_m, only: lm_origin

!Input variables
 integer, intent(in) :: mmax,nrl
 real(dp), intent(in) :: rr(mmax),ff(mmax),rl(nrl)
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: fl(nrl)

!Output variables
 integer, intent(out) :: ifail

 ifail=0
 if(rl(1)>=rr(1)) return
 call lm_origin(rr,ff,mmax,fl(1),ifail)
 if(ifail/=0) then
   write(iout,'(/a)') &
&       'origin: ERROR log mesh too short for the value at the origin'
 end if
 return
 end subroutine origin

 end module outmesh_m
