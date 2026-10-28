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

 module inward_start_m

 implicit none
 private
 public :: inward_start

 contains

 subroutine inward_start(ll,ee,rr,vv,mch,mmax,nin)

!! start of the inward integration of a bound state: the first point
!! beyond the classical turning point where the WKB decay
!! exp(-int kappa dr) has reached machine precision, or mmax-4 if it
!! does not on the mesh (the practical infinity of RadialAtomicSolvers.jl)
!!
!! The wave function is set to zero beyond nin.  A state barely bound in
!! a high-l channel, an open 4f shell, decays slowly beyond its turning
!! point: at 10 times its radius, the start used before, it had only
!! fallen by exp(-5), and the step in the density made the all-electron
!! energy converge to first order in the mesh spacing.
!!
!! @param ll  angular momentum
!! @param ee  energy
!! @param rr  log radial mesh
!! @param vv  local atomic potential
!! @param mch  classical turning point (matching point)
!! @param mmax  size of log grid
!! @param nin  start of the inward integration (output)

 use constants_m, only: dp

!Input variables
 integer, intent(in) :: ll,mch,mmax
 real(dp), intent(in) :: ee
 real(dp), intent(in) :: rr(mmax),vv(mmax)

!Output variables
 integer, intent(out) :: nin

!Local variables
 real(dp) :: budget,sk,xk,xk1,sls
 integer :: ii

 budget=-log(epsilon(1.0d0))
 sls=ll*(ll+1)

 nin=mmax-4
 sk=0.0d0
 xk1=0.0d0
 do ii=mch+1,mmax-4
   xk=sqrt(max(sls/rr(ii)**2 + 2.0d0*(vv(ii)-ee),0.0d0))
   sk=sk + 0.5d0*(xk + xk1)*(rr(ii) - rr(ii-1))
   xk1=xk
   if(sk>budget) then
     nin=ii
     exit
   end if
 end do

 return
 end subroutine inward_start

 end module inward_start_m
