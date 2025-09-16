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
subroutine run_diag_r(lmax, npa, epa, lloc, irc, &
    &                    vkb, evkb, nproj, rr, vfull, vp, zz, mmax, mxprj)

!diagnostics for semi-local and Vanderbilt-Kleinman-bylander pseudopotentials
!checks bound-state energies, norms, and slopes, and pseudo-bound state
!quantities for positive-energy projectors by matching all-electron
!log derivatives at rc.
!Should test as essentially exact for non-relativitic calculations
!Error for relativistic is B matrix Hermiticity error

!lmax  maximum angular momentum
!npa  principal quantum number for corresponding all-electron state
!epa  bound-state or scattering state reference energies for vkb potentials
!lloc  l for local potential
!irc  core radii indices
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vfull  all-electron potential
!vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
!zz  atomic number
!mmax  size of radial grid
!mxprj  dimension of number of projectors

    implicit none
    integer, parameter :: dp = kind(1.0d0)

!Input variables
    integer :: lmax, lloc, mmax, mxprj
    integer :: npa(mxprj, 6), irc(6), nproj(6)
    real(dp) :: zz
    real(dp) :: rr(mmax), vp(mmax, 5, 2), vfull(mmax), vkb(mmax, mxprj, 4, 2)
    real(dp) :: epa(mxprj, 6, 2), evkb(mxprj, 4, 2)

!Output variables - printing only

!Local variables
    integer :: ll, l1, ikap, kap, mkap, ierr, mch, mchf
    integer :: iprj, nnae, nnp, npr
    real(dp) :: al, emax, emin, etest, umch, upmch, uldf, gam, gpr, cnorm
    real(dp), allocatable :: uu(:), up(:), ur(:, :), urp(:, :)

    allocate (uu(mmax), up(mmax), ur(mmax, 2), urp(mmax, 2))

    al = 0.01d0 * dlog(rr(101) / rr(1))

! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors
!
    write (6, '(/a/a)') &
        & 'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
        & '  relativistic with J = L +/- 1/2 projectors'
    write (6, '(/2a)') '   l  kap   rcore       rmatch      e in        ', &
        &   'delta e    norm test   slope test'
!
    do l1 = 1, lmax + 1
        write (6, '(a)') ''
        ll = l1 - 1
        nnp = l1
        do iprj = 1, max(1, nproj(l1))
            if (l1 == 1) then
                mkap = 1
            else
                mkap = 2
            end if
!  loop on J = ll +/- 1/2
            do ikap = 1, mkap
                if (ikap == 1) then
                    kap = -(ll + 1)
                else
                    kap = ll
                end if

                npr = nproj(l1)

!   find cutoff radius for projectors
                mchf = max(irc(l1), irc(lloc + 1)) + 5

                etest = epa(iprj, l1, ikap)
                if (epa(iprj, l1, 1) < 0.0d0) then
                    call ldiracfb(npa(iprj, l1), ll, kap, ierr, etest, &
                        &                     rr, zz, vfull, ur, urp, mmax, mch)
                    if (ierr /= 0) then
                        write (6, '(/a,4i4)') 'run_diag_r: ldiracfb convergence ERROR &
                            &         n,l,,kap,ierr=', &
                            &        npa(iprj, l1), ll, kap, ierr
                        stop
                    end if
                else
                    call ldiracfs(nnae, ll, kap, ierr, etest, &
                        &                  rr, zz, vfull, ur, urp, mmax, mchf)
                end if
                call renorm_r(ur, rr, ll, kap, zz, mmax, cnorm)
                umch = ur(mchf, 1)
                upmch = cnorm * urp(mchf, 1)
                uldf = upmch / umch

                if (l1 == lloc + 1) npr = 0

                etest = epa(iprj, l1, ikap)
                emax = etest + 0.1d0 * dabs(etest) + 0.05d0
                emin = etest - 0.1d0 * dabs(etest) - 0.05d0

                if (epa(iprj, l1, 1) < 0.0d0) then

                    emax = dmin1(emax, 0.0d0)

                    call lschvkbb(ll + iprj, ll, npr, ierr, etest, emin, emax, &
                        &                   rr, vp(1, lloc + 1, ikap), vkb(1, 1, l1, ikap), evkb(1, l1, ikap), &
                        &                   uu, up, mmax, mch)
                    if (ierr /= 0) then
                        write (6, '(a,4i4,1p,2e16.8)') 'run_diag_r: lschvkbb ERROR', &
                            &             iprj, ll, kap, ierr, epa(iprj, l1, ikap), etest
                    end if

                    nnp = nnp + 1
                else
!     calculate effective pseudo wave function principal quantum number
!     from all-electron node count

                    nnp = nnae - npa(1, l1) + ll + 1
                    call lschvkbbe(nnp, ll, npr, ierr, etest, uldf, emin, emax, &
                        &                    rr, vp(1, lloc + 1, ikap), vkb(1, 1, l1, ikap), evkb(1, l1, ikap), &
                        &                    uu, up, mmax, mchf)

                    if (ierr /= 0) then
                        write (6, '(a,4i4,1p,2e16.8)') 'run_diag_r: lschvkbbe ERROR', &
                            &             iprj, ll, kap, ierr, epa(iprj, l1, ikap), etest
                    end if
                end if !epa<0 (bound or not)

                gam = dabs(umch / uu(mchf))
                gpr = dabs(upmch / up(mchf))
!
                write (6, '(2i4,6f12.7)') ll, kap, rr(irc(l1)), rr(mchf), epa(iprj, l1, ikap), &
                    &         etest - epa(iprj, l1, ikap), gam, gpr

            end do  !ikap
        end do  !iprj
    end do !l1

    deallocate (uu, up)
    return
end subroutine run_diag_r
