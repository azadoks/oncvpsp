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
module exc_m
   use precision_m, only: dp
   use constants_m, only: pi, fourpi, ifourpi, third
   implicit none
   private
   public :: exchdl
   public :: excpzca
   public :: excwig
contains

!> Calculates Hedin-Lundquist exchange-correlation potential and energy density
subroutine exchdl(rho, vxc, exc, mmax)

   !Inpupt variables
   integer :: mmax
   real(dp) :: rho(mmax)

   !Output variables
   real(dp) :: vxc(mmax), exc(mmax)

   !Local variables
   integer :: ii
   real(dp) :: conrs, convx, conex
   real(dp) :: rs, ecp, aln, xx, rh

   conrs = (3.d0 / fourpi)**third
   convx = (1.5d0 / pi)**(2.0d0 / 3.0d0)
   conex = 0.75d0 * convx

   ! Hedin-Lundqvist correlation
   do ii = 1, mmax
      if (rho(ii) > 1.0d-20) then
         rh = rho(ii) * ifourpi
         rs = conrs / rh**third
         xx = rs / 21.0d0
         aln = log(1.0d0 + 1.0d0 / xx)
         ecp = aln + (xx**3 * aln - xx * xx) + 0.5d0 * xx - third
         exc(ii) = -conex / rs - 0.0225d0 * ecp
         vxc(ii) = -convx / rs - 0.0225d0 * aln
      else
         vxc(ii) = 0.0d0; exc(ii) = 0.0d0
      end if
   end do

   return
end subroutine exchdl

!> Wigner ingerpolation formula for exchange-correlation potential and energy density
subroutine excwig(rho, vxc, exc, mmax)
   !rho  charge density
   !vxc  exchange-correlation potential
   !exc  exchange-correlation energy density
   !mmax  dimension of log radial grid

   !Input variables
   real(dp) :: rho(mmax)
   integer :: mmax

   !Output variables
   real(dp) :: vxc(mmax), exc(mmax)

   !Local vafriables
   real(dp) :: rh, rs
   real(dp) :: conrs, convx, conex
   integer :: ii

   conrs = (3.d0 / fourpi)**third
   convx = (1.5d0 / pi)**(2.0d0 / 3.0d0)
   conex = 0.75d0 * convx

   ! Wigner interpolation formula, E. Wigner, Phys. Rev. 46, 1002 (1934).

   do ii = 1, mmax
      if (rho(ii) >= 1.0d-15) then
         rh = rho(ii) * ifourpi
         rs = conrs * rh**(-third)
         vxc(ii) = -convx / rs - 0.44d0 * (4.0d0 * third * rs + 7.8d0) / (rs + 7.8d0)**2
         exc(ii) = -conex / rs - 0.44d0 / (rs + 7.8d0)
      end if
   end do
   return
end subroutine excwig

!> Calculates Perdew-Zunger-Ceperly_Alder exchange-correlation potential
!> and energy density J. P. Perdew and A. Zunger, Phys. Rev. B23, 5048 (1981)
subroutine excpzca(rho, vxc, exc, mmax)
   !rho  charge density
   !vxc  exchange-correlation potential
   !exc  exchange-correlation energy density
   !mmax  dimension of log radial grid

   !Input variables
   integer :: mmax
   real(dp) :: rho(mmax)

   !Output variables
   real(dp) :: vxc(mmax), exc(mmax)

   !Local variables
   integer :: ii
   real(dp) :: rs, rh, sqrs, rsl, den

   real(dp), parameter :: conrs = (3.d0 / fourpi)**third
   real(dp), parameter :: convx = (1.5d0 / pi)**(2.0d0 / 3.0d0)
   real(dp), parameter :: conex = 0.75d0 * convx

   !rs>1, Ceperley-Alder fit to QMC
   real(dp), parameter :: gam = -0.1423d0
   real(dp), parameter :: bt1 = 1.0529d0
   real(dp), parameter :: bt2 = 0.3334d0

   !rs<1, Perdew-Zunger leading terms of RPA high-density expansion
   !AA and BB known from Gell-Mann & Bruckner
   !PZ adjust CC and DD for contiunity of exc and vxc at rs=1
   real(dp), parameter :: AA = 0.0311d0
   real(dp), parameter :: BB = -0.0480d0
   real(dp), parameter :: CC = 0.0020d0
   real(dp), parameter :: DD = -0.0116d0

   !PZ modified coeffients from libxc
   !real(dp), parameter :: CC =   0.0020191519406228d0
   !real(dp), parameter :: DD =  -0.0116320663789130d0

   do ii = 1, mmax
      if (rho(ii) > 1.0d-20) then
         rh = rho(ii) / fourpi
         rs = conrs * rh**(-third)
         if (rs > 1.0d0) then
            sqrs = dsqrt(rs)
            den = 1.0d0 + bt1 * sqrs + bt2 * rs
            exc(ii) = -conex / rs + gam / den
            vxc(ii) = -convx / rs + gam * (1.0d0 + (3.5d0 * bt1 * sqrs &
               &            + 4.0d0 * bt2 * rs) * third) / den**2
         else
            rsl = dlog(rs)
            exc(ii) = -conex / rs + AA * rsl + BB + CC * rs * rsl + DD * rs
            vxc(ii) = -convx / rs + AA * rsl + (BB - third * AA) + 2.0d0 * third * CC * rs * rsl &
               &            + third * (2.0d0 * DD - CC) * rs
         end if
      else
         vxc(ii) = 0.0d0; exc(ii) = 0.0d0
      end if
   end do

   return
end subroutine excpzca

end module exc_m
