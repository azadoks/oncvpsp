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
!> The log mesh of a well state, continued beyond rr(mmax) as far as the
!> state needs. A well state's tail beyond the well's onset decays only as
!> fast as the well rises, and a wide well (a small ep) held it out to and
!> past rmax: the mesh's end was then the box, and the state, its energy
!> and its tail in E_resid depended on rmax. The mesh is continued,
!> rr(i) = rr(1) amesh**(i-1), until the WKB decay exp(-int kappa dr)
!> beyond the well's outer turning point at the target energy has reached
!> machine precision, the start inward_start gives the solvers; the
!> all-electron potential beyond rr(mmax) is vfull(mmax) rr(mmax)/r, as in
!> free_tail_m. A state that decays within the mesh gets no extension and
!> is solved as before.

 module well_mesh_m
 use constants_m, only: dp
 implicit none
 private
 public :: well_mesh_t, wm_build, wm_short

!> a well state on its own mesh: the mesh, the well's potential on it, and
!> the state (the large and small components for irel=2) and its
!> derivative in the mesh index; mx = 0 for a state on the log mesh
 type :: well_mesh_t
   integer :: mx=0
   real(dp), allocatable :: rr(:),vv(:),uu(:,:),up(:,:)
 end type well_mesh_t

!> the most the mesh is continued to (bohr)
 real(dp), parameter :: rlimit=1.0d4

 contains

 subroutine wm_build(rr,vfull,mmax,mx,rx,vx)
!! the log mesh continued to mx >= mmax points, and the all-electron
!! potential on it
!!
!! @param rr  log radial mesh
!! @param vfull  all-electron potential
!! @param mmax  size of the log mesh
!! @param mx  size of the continued mesh
!! @param rx  the continued mesh
!! @param vx  vfull on it, vfull(mmax) rr(mmax)/r beyond rr(mmax)

 use logmesh_m, only: lm_step

!Input variables
 integer, intent(in) :: mmax,mx
 real(dp), intent(in) :: rr(mmax),vfull(mmax)
!Output variables
 real(dp), allocatable, intent(out) :: rx(:),vx(:)
!Local variables
 real(dp) :: al
 integer :: ii

 al=lm_step(rr,mmax)
 allocate(rx(mx),vx(mx))
 rx(1:mmax)=rr(:)
 vx(1:mmax)=vfull(:)
 do ii=mmax+1,mx
   rx(ii)=rr(1)*exp(al*dble(ii-1))
   vx(ii)=vfull(mmax)*rr(mmax)/rx(ii)
 end do
 return
 end subroutine wm_build

 subroutine wm_short(ll,ee,rx,vx,mx,mneed,iout,ifail)
!! whether a state at ee in the potential vx needs a longer mesh: mneed is
!! 0 if the WKB decay beyond the potential's outer turning point reaches
!! machine precision before mx-4 (inward_start's budget), else a size
!! that adds a quarter of the radius, or ERROR beyond rlimit
!!
!! @param ll  angular momentum
!! @param ee  energy
!! @param rx  log mesh
!! @param vx  potential on it
!! @param mx  its size
!! @param mneed  0, or the size to continue the mesh to
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: ll,mx
 real(dp), intent(in) :: ee
 real(dp), intent(in) :: rx(mx),vx(mx)
 integer, intent(in) :: iout
!Output variables
 integer, intent(out) :: mneed
 integer, intent(out) :: ifail
!Local variables
 real(dp) :: budget,sk,xk,xk1,sls,al
 integer :: ii,itp

 ifail=0
 budget=-log(epsilon(1.0d0))
 sls=ll*(ll+1)
!the outer turning point: the last point classically allowed
 itp=1
 do ii=mx,2,-1
   if(sls/rx(ii)**2 + 2.0d0*(vx(ii)-ee)<0.0d0) then
     itp=ii
     exit
   end if
 end do
 mneed=0
 sk=0.0d0
 xk1=0.0d0
 do ii=itp+1,mx-4
   xk=sqrt(max(sls/rx(ii)**2 + 2.0d0*(vx(ii)-ee),0.0d0))
   sk=sk + 0.5d0*(xk + xk1)*(rx(ii) - rx(ii-1))
   xk1=xk
   if(sk>budget) return
 end do
 al=log(rx(mx)/rx(1))/dble(mx-1)
 mneed=mx+max(10,ceiling(log(1.25d0)/al))
 if(rx(1)*exp(al*dble(mneed-1))>rlimit) then
   write(iout,'(/a,f8.4,a,i2,a,f10.1,a)') &
&       'wm_short: ERROR the well state at',ee,' Ha, l=',ll, &
&       ' does not decay within',rlimit,' bohr'
   ifail=1
   return
 end if
 return
 end subroutine wm_short

 end module well_mesh_m
