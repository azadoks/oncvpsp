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
subroutine excwig(rho, vxc, exc, mmax)

    ! Wigner ingerpolation formula for exchange-correlation potential and
    ! energy density

    !rho  charge density
    !vxc  exchange-correlation potential
    !exc  exchange-correlation energy density
    !mmax  dimension of log radial grid
    use m_constants, only: dp, pi, fourpi, ifourpi, third
    implicit none

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
