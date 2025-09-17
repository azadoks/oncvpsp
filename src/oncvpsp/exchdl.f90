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
subroutine exchdl(rho, vxc, exc, mmax)

    !calculates Hedin-Lundquist exchange-correlation potential and energy
    !density
    use m_constants, only: dp, pi, fourpi, ifourpi, third
    implicit none

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
            aln = dlog(1.0d0 + 1.0d0 / xx)
            ecp = aln + (xx**3 * aln - xx * xx) + 0.5d0 * xx - third
            !    ecp = aln+(xx**3*aln-xx*xx)+xx/2-1.0d0/3.0d0
            !    exc(ii)=-0.458175d0/rs - 0.0225d0*ecp
            !    vxc(ii)=-0.6109d0/rs - 0.0225d0*aln
            exc(ii) = -conex / rs - 0.0225d0 * ecp
            vxc(ii) = -convx / rs - 0.0225d0 * aln
        else
            vxc(ii) = 0.0d0; exc(ii) = 0.0d0
        end if
    end do

    return
end subroutine exchdl
