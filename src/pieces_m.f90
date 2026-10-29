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
!> Functions built from pieces on the log mesh (doc/smoothness_plan.md,
!> R-6): a projector is a smooth function on [0, 0.95 rc], another on
!> [0.95 rc, rc] (the taper), zero beyond, with a kink at each boundary;
!> the local potential a polynomial inside rc(5) and the all-electron
!> potential beyond. Each piece is kept as a smooth function on the whole
!> mesh (its extension past its segment, which the piece's formula gives,
!> at least pw_reach nodes past each boundary, 0 far away), and every
!> operation is done per piece on that extension: an integral is the sum
!> over the segments of lm_integral between the segment's bounds, so no
!> stencil or interpolant spans a kink, and the result is smooth in the
!> bounds. A piece may also have a closed form (the model core's fill and
!> Teter function), which then gives its values between the nodes, and at
!> the origin, in place of the interpolant (pw_value, outmesh_m).

 module pieces_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: piecewise_t, pw_maxseg, pw_reach, pw_alloc, pw_free
 public :: pw_integral, pw_product_integral, pw_node_values, pw_segment
 public :: pw_value, pw_resample, pw_form, pw_set_form, pw_form_value

!> the most segments a function has
 integer, parameter :: pw_maxseg=4
!> nodes on each side of a segment on which its piece must hold the piece's
!> smooth extension: the reach of the interpolant and of the cell rule of
!> lm_integral at the segment's bounds, and of the grids of outwf_m, which
!> run og_ext points past a segment's end and interpolate the piece there
!> (og_ext + lm_fh_reach = 59)
 integer, parameter :: pw_reach=64

 abstract interface
!> a closed form of a piece: the function and its first 4 derivatives at
!> rad, from the parameters par
   subroutine pw_form(par,rad,dd)
   import :: dp
   implicit none
   real(dp), intent(in) :: par(:),rad
   real(dp), intent(out) :: dd(5)
   end subroutine pw_form
 end interface

 type piecewise_t
   integer :: nseg=0  !< number of segments
!> bound(k-1) < r <= bound(k) is segment k; bound(0)=0
   real(dp) :: bound(0:pw_maxseg)=0.0d0
!> vals(:,k) the k-th piece on the mesh
   real(dp), allocatable :: vals(:,:)
!> the segment whose piece has a closed form (0 none) (pw_set_form)
   integer :: iform=0
!> the element of the form's result the piece is (1 the function, 2..5
!> its derivatives)
   integer :: icol=1
!> the parameters of the form
   real(dp), allocatable :: par(:)
!> the closed form of the piece of segment iform
   procedure(pw_form), pointer, nopass :: form=>null()
 end type piecewise_t

 contains

 subroutine pw_alloc(pw,mmax,nseg,bound,ifail)

!! a piecewise function of nseg segments with the bounds given, its pieces
!! zero
!!
!! @param pw  the function
!! @param mmax  size of the log radial mesh
!! @param nseg  number of segments, 1..pw_maxseg
!! @param bound  bound(0..nseg), nondecreasing, bound(0)=0
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_internal, failed, msg_pw_alloc_order, &
&                    msg_pw_alloc_segments

!Input variables
 integer, intent(in) :: mmax,nseg
 real(dp), intent(in) :: bound(0:nseg)

!Output variables
 type(piecewise_t), intent(out) :: pw
 integer, intent(out) :: ifail

!Local variables
 integer :: kk

 ifail=0
 if(nseg<1 .or. nseg>pw_maxseg) then
   write(error_unit,'(/a,i0)') 'pw_alloc: ERROR number of segments ',nseg
   ifail=failed(fail_internal,msg_pw_alloc_segments)
   return
 end if
 do kk=1,nseg
   if(bound(kk)<bound(kk-1)) then
     write(error_unit,'(/a,i0)') &
&         'pw_alloc: ERROR bounds out of order at segment ',kk
     ifail=failed(fail_internal,msg_pw_alloc_order)
     return
   end if
 end do
 pw%nseg=nseg
 pw%bound(:)=0.0d0
 pw%bound(0:nseg)=bound(0:nseg)
 allocate(pw%vals(mmax,nseg))
 pw%vals(:,:)=0.0d0
 pw%iform=0
 return
 end subroutine pw_alloc

 subroutine pw_set_form(pw,iseg,form,par,icol)

!! a closed form for the piece of segment iseg: element icol of what form
!! returns (the function, icol=1, or its derivative icol-1)
!!
!! @param pw  the piecewise function
!! @param iseg  the segment
!! @param form  the closed form
!! @param par  its parameters
!! @param icol  the element of its result, 1..5

!Input variables
 integer, intent(in) :: iseg,icol
 real(dp), intent(in) :: par(:)
 procedure(pw_form) :: form

!Input/output variables
 type(piecewise_t), intent(inout) :: pw

 pw%iform=iseg
 pw%icol=icol
 pw%par=par
 pw%form=>form
 return
 end subroutine pw_set_form

 subroutine pw_form_value(pw,kk,rq,fq,found)

!! whether the piece of segment kk has a closed form, and if so its value
!! at rq in fq
!!
!! @param pw  the piecewise function
!! @param kk  the segment
!! @param rq  radius, 0 or more (bohr)
!! @param fq  the value, when the piece has a closed form (0 otherwise)
!! @param found  whether the piece has a closed form

!Input variables
 type(piecewise_t), intent(in) :: pw
 integer, intent(in) :: kk
 real(dp), intent(in) :: rq

!Output variables
 real(dp), intent(out) :: fq
 logical, intent(out) :: found

!Local variables
 real(dp) :: dd(5)

 fq=0.0d0
 found=(kk>0 .and. kk==pw%iform)
 if(.not. found) return
 call pw%form(pw%par,rq,dd)
 fq=dd(pw%icol)
 return
 end subroutine pw_form_value

 subroutine pw_free(pw)

!! releases the pieces
!!
!! @param pw  the piecewise function; its arrays deallocated

!Input/output variables
 type(piecewise_t), intent(inout) :: pw

 if(allocated(pw%vals)) deallocate(pw%vals)
 if(allocated(pw%par)) deallocate(pw%par)
 pw%nseg=0
 pw%iform=0
 return
 end subroutine pw_free

 subroutine pw_integral(rr,mmax,mm,ff,pw,ss)

!! int_0^bound(nseg) ff(r) pw(r) dr for a smooth ff given on the mesh:
!! the sum over the segments of lm_integral of ff times the piece
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param mm  power of r with which ff*pw goes to 0 at the origin
!! @param ff  smooth function on the mesh
!! @param pw  the piecewise function
!! @param ss  the integral

 use logmesh_m, only: lm_integral

!Input variables
 integer, intent(in) :: mmax,mm
 real(dp), intent(in) :: rr(mmax),ff(mmax)
 type(piecewise_t), intent(in) :: pw

!Output variables
 real(dp), intent(out) :: ss

!Local variables
 real(dp) :: sk
 integer :: kk

 ss=0.0d0
 do kk=1,pw%nseg
   if(pw%bound(kk)<=pw%bound(kk-1)) cycle
   call lm_integral(rr,ff*pw%vals(:,kk),mmax,mm,pw%bound(kk-1), &
&                   pw%bound(kk),sk)
   ss=ss+sk
 end do
 return
 end subroutine pw_integral

 subroutine pw_product_integral(rr,mmax,mm,pwa,pwb,ss,ifail)

!! int_0^bound(nseg) pwa(r) pwb(r) dr for two piecewise functions with the
!! same segments
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param mm  power of r with which the product goes to 0 at the origin
!! @param pwa,pwb  the functions
!! @param ss  the integral
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_internal, failed, msg_pw_product_segments
 use logmesh_m, only: lm_integral

!Input variables
 integer, intent(in) :: mmax,mm
 real(dp), intent(in) :: rr(mmax)
 type(piecewise_t), intent(in) :: pwa,pwb

!Output variables
 real(dp), intent(out) :: ss
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: sk
 integer :: kk

 ifail=0
 ss=0.0d0
 if(pwa%nseg/=pwb%nseg .or. any(pwa%bound(0:pwa%nseg)/=pwb%bound(0:pwa%nseg))) &
&  then
   write(error_unit,'(/a)') 'pw_product_integral: ERROR the segments differ'
   ifail=failed(fail_internal,msg_pw_product_segments)
   return
 end if
 ss=0.0d0
 do kk=1,pwa%nseg
   if(pwa%bound(kk)<=pwa%bound(kk-1)) cycle
   call lm_integral(rr,pwa%vals(:,kk)*pwb%vals(:,kk),mmax,mm, &
&                   pwa%bound(kk-1),pwa%bound(kk),sk)
   ss=ss+sk
 end do
 return
 end subroutine pw_product_integral

 subroutine pw_node_values(rr,mmax,pw,ff)

!! the function on the mesh: at each node the piece of the segment the
!! node is in, 0 beyond the last bound
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param pw  the piecewise function
!! @param ff  its values at the nodes

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax)
 type(piecewise_t), intent(in) :: pw

!Output variables
 real(dp), intent(out) :: ff(mmax)

!Local variables
 integer :: ii,kk

 do ii=1,mmax
   kk=pw_segment(pw,rr(ii))
   if(kk>0) then
     ff(ii)=pw%vals(ii,kk)
   else
     ff(ii)=0.0d0
   end if
 end do
 return
 end subroutine pw_node_values

 integer function pw_segment(pw,rq)

!! the segment that holds the radius rq, bound(k-1) < rq <= bound(k)
!! (segment 1 from the origin), or 0 beyond the last bound
!!
!! @param pw  the piecewise function
!! @param rq  radius (bohr)

!Input variables
 type(piecewise_t), intent(in) :: pw
 real(dp), intent(in) :: rq

!Local variables
 integer :: kk

 pw_segment=0
 do kk=1,pw%nseg
   if(rq<=pw%bound(kk)) then
     pw_segment=kk
     return
   end if
 end do
 return
 end function pw_segment

 real(dp) function pw_value(rr,mmax,pw,rq)

!! the function at a radius: the closed form or the interpolant
!! (lm_interp) of the piece of the segment that holds it, 0 beyond the
!! last bound
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param pw  the piecewise function
!! @param rq  radius, rr(1) <= rq <= rr(mmax) (bohr)

 use logmesh_m, only: lm_interp

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rq
 type(piecewise_t), intent(in) :: pw

!Local variables
 integer :: kk
 logical :: found

 pw_value=0.0d0
 kk=pw_segment(pw,rq)
 call pw_form_value(pw,kk,rq,pw_value,found)
 if(found) return
 if(kk>0) call lm_interp(rr,pw%vals(:,kk),mmax,rq,pw_value)
 return
 end function pw_value

 subroutine pw_resample(rr,mmax,pw,rb,nb,pwb,ifail)

!! the function on another log mesh of the same step, rb(i) = s rr(i),
!! s >= 1: the same bounds, each piece the interpolant (lm_interp) of the
!! piece on rr at the nodes of rb up to rr(mmax), 0 beyond
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial meshes
!! @param pw  the piecewise function on rr
!! @param rb  the other mesh
!! @param nb  the nodes of rb at which the pieces are needed, rb(nb)<=rr(mmax)
!! @param pwb  the function on rb
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use logmesh_m, only: lm_interp

!Input variables
 integer, intent(in) :: mmax,nb
 real(dp), intent(in) :: rr(mmax),rb(mmax)
 type(piecewise_t), intent(in) :: pw

!Output variables
 type(piecewise_t), intent(out) :: pwb
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,kk

 call pw_alloc(pwb,mmax,pw%nseg,pw%bound(0:pw%nseg),ifail)
 if(ifail/=0) return
 do kk=1,pw%nseg
   do ii=1,nb
     call lm_interp(rr,pw%vals(:,kk),mmax,rb(ii),pwb%vals(ii,kk))
   end do
 end do
 return
 end subroutine pw_resample

 end module pieces_m
