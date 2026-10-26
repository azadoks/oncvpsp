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
! Creates monotonic polynomial model core charge matching all-electron
! core charge and 4 derivatives at the "crossover" radius r*, where the core
! charge equals fcfact times the pseudo valence charge.
! Polynomial is 7th-order with no linear or quadratic term.
! r* is the root of ln(rhoc) - ln(fcfact rhotps), the derivatives at r* are
! interpolated (lm_derivs), and the constant term is in closed form
! (modcore_m), so the model is smooth in the inputs
! (doc/smoothness_plan.md, M-1 to M-3).

! Performs analysis and based on "hardness" criterion described in
! Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as

 subroutine modcore(icmod,rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
&                   corepw,fcfact,rcfact,rcmax,mmax,rr,nc,nv,la,zion, &
&                   iexc)

!! @param icmod  (not used)
!! @param rhops  state-by-state pseudocharge density
!! @param rhotps  total pseudocharge density
!! @param rhoc  core-charge density
!! @param rhoae  state-by-state all-electron valence charge density
!! @param rhotae  total all-electron valence charge density
!! @param rhomod  model core density and 4 derivatives
!! @param corepw  its pieces, split at r* (modcore_m, core_pieces)
!! @param fcfact  the model is matched where rhoc equals fcfact times the total
!!                pseudocharge
!! @param rcfact  (not used)
!! @param rcmax  maximum rc
!! @param mmax  dimension of log grid
!! @param rr  log radial grid
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular-momenta
!! @param zion  ion charge
!! @param iexc  exchange-correlation function to be used

 use constants_m, only: dp, mxstate, nrhomod
 use logmesh_m, only: lm_derivs
 use modcore_m, only: core_crossover, poly_core_fit, core_pieces, &
&                     d2exc_baseline, d2exc_report
 use pieces_m, only: piecewise_t
 implicit none

!Input variables
 integer, intent(in) :: icmod,nv,nc,iexc,mmax
 integer, intent(in) :: la(mxstate)
 real(dp), intent(in) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: zion,fcfact,rcfact,rcmax

!Output variables
 real(dp), intent(out) :: rhomod(mmax,nrhomod)
 type(piecewise_t), intent(out) :: corepw(nrhomod)

!Local variables
 real(dp) :: a0,d2mdiff,rstar,rcut
 real(dp) :: aco(5),fmatch(5)
 real(dp), allocatable :: d2excae(:,:),d2excps(:,:)
 integer :: ierr

 allocate(d2excae(nv,nv),d2excps(nv,nv))

 d2excae(:,:)=0.0d0
 d2excps(:,:)=0.0d0

 call core_crossover(rr,rhoc,rhotps,mmax,fcfact,rstar,ierr)
 if(ierr/=0) then
   write(6,'(/a)') 'rhomod: ERROR ircc (core-valence charge crossover) &
&        not found'
   error stop 1
 end if

!set limit for d2Exc calculation to radius beyond which rhomod=rhoc
 rcut=max(rcmax,rstar)

 call d2exc_baseline(rhotae,rhoc,rhoae,rhotps,rhops,rr,zion,iexc,nc,nv, &
&                    la,rcut,mmax,rhomod,d2excae,d2excps,d2mdiff)

! core charge and its derivatives at r*, and the barely monotonic
! polynomial that matches them

 call lm_derivs(rr,rhoc,mmax,rstar,4,fmatch)
 call poly_core_fit(rstar,fmatch,a0,aco)

! fill in model and derivatives: the polynomial inside r*, the core charge
! outside
 call core_pieces(rr,rhoc,mmax,rstar,poly_fill,corepw,rhomod)


!test model
 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
&                   zion,iexc,nc,nv,la,rcut,mmax)

 write(6,'(/a/)') 'Polynomial model core charge'
 write(6,'(a,f12.6)') '  model radius r*',rstar
 call d2exc_report('d2excps - pseudofunction derivatives with core '// &
&                  'correction',d2excps,nv,d2mdiff)

 deallocate(d2excae,d2excps)

 return

 contains

 subroutine poly_fill(rad,dd)

!! the polynomial a0 + sum_j aco(j) r**(j+2) and its first 4 derivatives
!!
!! @param rad  radius (bohr)
!! @param dd  the polynomial and its derivatives at rad

!Input variables
 real(dp), intent(in) :: rad

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 real(dp) :: fac
 integer :: ii,jj,ll

 do ii=1,5
   if(ii==1) then
     dd(ii)=a0
   else
     dd(ii)=0.0d0
   end if
   do jj=1,5
     fac=1.0d0
     do ll=0,ii-2
       fac=fac*(jj+2-ll)
     end do
     dd(ii)=dd(ii)+fac*aco(jj)*rad**(jj+3-ii)
   end do
 end do
 end subroutine poly_fill

 end subroutine modcore
