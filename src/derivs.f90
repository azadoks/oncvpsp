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
 subroutine derivs(mmax, rho, al, rr, dpr, dppr, dlap)
!! first and second radial derivatives and Laplacian of rho on the log
!! mesh, by finite differences, for excggc, exc_libxc and functionals_m
!! NB: presumes incoming rho is properly scaled without 4pi or r**2 factors.
!!
!! @param mmax  dimension of log radial grid
!! @param rho  charge density
!! @param al  log of the mesh ratio amesh of the log radial mesh
!! @param rr  log radial mesh
!! @param dpr  first radial derivative of rho
!! @param dppr  second radial derivative of rho
!! @param dlap  Laplacian of rho

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: rho(mmax)

!Output variables
 real(dp), intent(out) :: dpr(mmax), dppr(mmax), dlap(mmax)

!Local variables
 integer :: ii
 real(dp) :: dpn(mmax), dppn(mmax)
 real(dp) :: c11,c12,c13,c14,c15
 real(dp) :: c21,c22,c23,c24,c25

 c11 =   2.0d0 / 24.0d0
 c12 = -16.0d0 / 24.0d0
 c13 =   0.0d0 / 24.0d0
 c14 =  16.0d0 / 24.0d0
 c15 =  -2.0d0 / 24.0d0
!
 c21 =   -1.0d0 / 12.0d0
 c22 =   16.0d0 / 12.0d0
 c23 =  -30.0d0 / 12.0d0
 c24 =   16.0d0 / 12.0d0
 c25 =   -1.0d0 / 12.0d0
!
! n derivatives of d
!
 ii=1
 dpn(ii) = -25.d0/12.d0*rho(ii) +4.d0*rho(ii+1) -3.d0*rho(ii+2) &
&         +4.d0/3.d0*rho(ii+3) -1.d0/4.d0*rho(ii+4)
 dppn(ii) = 15.d0/4.d0*rho(ii) -77.d0/6.d0*rho(ii+1) +107.d0/6.d0*rho(ii+2) &
&         -13.d0*rho(ii+3) +61.d0/12.d0*rho(ii+4) -5.d0/6.d0*rho(ii+5)
 ii=2
 dpn(ii) = -25.d0/12.d0*rho(ii) +4.d0*rho(ii+1) -3.d0*rho(ii+2)  &
&         +4.d0/3.d0*rho(ii+3) -1.d0/4.d0*rho(ii+4)
 dppn(ii) = 15.d0/4.d0*rho(ii) -77.d0/6.d0*rho(ii+1) +107.d0/6.d0*rho(ii+2) &
&         -13.d0*rho(ii+3) +61.d0/12.d0*rho(ii+4) -5.d0/6.d0*rho(ii+5)

 do ii = 3, mmax - 2
   dpn(ii) =  c11*rho(ii-2) + c12*rho(ii-1) + c14*rho(ii+1) + c15*rho(ii+2)
   dppn(ii) = c21*rho(ii-2) + c22*rho(ii-1) + c23*rho(ii)   + c24*rho(ii+1) &
&           +c25*rho(ii+2)
 end do

 ii=mmax-1
 dpn(ii) = +25.d0/12.d0*rho(ii) -4.d0*rho(ii-1) +3.d0*rho(ii-2) &
&         -4.d0/3.d0*rho(ii-3) +1.d0/4.d0*rho(ii-4)
 dppn(ii) = -15.d0/4.d0*rho(ii) +77.d0/6.d0*rho(ii-1) -107.d0/6.d0*rho(ii-2) &
&          +13.d0*rho(ii-3) -61.d0/12.d0*rho(ii-4) +5.d0/6.d0*rho(ii-5)
 ii=mmax
 dpn(ii) = +25.d0/12.d0*rho(ii) -4.d0*rho(ii-1) +3.d0*rho(ii-2) &
&         -4.d0/3.d0*rho(ii-3) +1.d0/4.d0*rho(ii-4)
 dppn(ii) = -15.d0/4.d0*rho(ii) +77.d0/6.d0*rho(ii-1) -107.d0/6.d0*rho(ii-2) &
&          +13.d0*rho(ii-3) -61.d0/12.d0*rho(ii-4) +5.d0/6.d0*rho(ii-5)

!
! r derivatives of d
!
 do ii = 1, mmax
   dpr(ii) = dpn(ii) / (al * rr(ii))
   dppr(ii) = (dppn(ii) - al * dpn(ii)) / (al * rr(ii))**2
   dlap(ii) = (dppn(ii) + al * dpn(ii)) / (al * rr(ii))**2
 end do
 return

 end subroutine derivs
