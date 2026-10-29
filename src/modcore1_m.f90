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

 module modcore1_m

 implicit none
 private
 public :: modcore

 contains

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
&                   corepw,fcfact,rcfact,rcmax,mmax,adams,rr,nc,nv,la,zion, &
&                   iexc,c_light,iout,ifail)

!! @param icmod  (not used)
!! @param rhops  state-by-state pseudocharge density
!! @param rhotps  pseudo valence charge density
!! @param rhoc  core charge density
!! @param rhoae  state-by-state all-electron valence charge density
!! @param rhotae  all-electron valence charge density
!! @param rhomod  model core charge and its first 4 derivatives
!! @param corepw  its pieces, split at r* (modcore_m, core_pieces)
!! @param fcfact  the model is matched where rhoc equals fcfact times the total
!!                pseudocharge
!! @param rcfact  (not used)
!! @param rcmax  largest core radius (bohr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param rr  log radial mesh (bohr)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular momenta
!! @param zion  ion charge
!! @param iexc  exchange-correlation functional
!! @param c_light  speed of light (Hartree atomic units)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_not_found
 use constants_m, only: dp, mxstate, nrhomod
 use logmesh_m, only: lm_derivs
 use modcore_m, only: core_crossover, poly_core_fit, core_pieces, poly_form, &
&                     d2exc_baseline, d2exc_report
 use pieces_m, only: piecewise_t
 use der2exc_m, only: der2exc
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: icmod,nv,nc,iexc,mmax
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 integer, intent(in) :: la(mxstate)
 real(dp), intent(in) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: zion,fcfact,rcfact,rcmax
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: rhomod(mmax,nrhomod)
 type(piecewise_t), intent(out) :: corepw(nrhomod)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: a0,d2mdiff,rstar,rcut
 real(dp) :: aco(5),fmatch(5)
 real(dp), allocatable :: d2excae(:,:),d2excps(:,:)
 integer :: ierr

 ifail=0
 allocate(d2excae(nv,nv),d2excps(nv,nv))

 d2excae(:,:)=0.0d0
 d2excps(:,:)=0.0d0

 call core_crossover(rr,rhoc,rhotps,mmax,fcfact,rstar,ierr,ifail)
 if(ifail/=0) return
 if(ierr/=0) then
   write(iout,'(/a)') 'rhomod: ERROR ircc (core-valence charge crossover) &
&        not found'
   ifail=fail_not_found
   return
 end if

!set limit for d2Exc calculation to radius beyond which rhomod=rhoc
 rcut=max(rcmax,rstar)

 call d2exc_baseline(rhotae,rhoc,rhoae,rhotps,rhops,rr,zion,iexc,c_light,nc, &
&                    nv,la,rcut,mmax,adams,rhomod,d2excae,d2excps,d2mdiff, &
&                    iout,ifail)
 if(ifail/=0) return

! core charge and its derivatives at r*, and the barely monotonic
! polynomial that matches them

 call lm_derivs(rr,rhoc,mmax,rstar,4,fmatch)
 call poly_core_fit(rstar,fmatch,a0,aco,iout,ifail)
 if(ifail/=0) return

! fill in model and derivatives: the polynomial inside r*, the core charge
! outside
 call core_pieces(rr,rhoc,mmax,rstar,poly_form,[a0,aco],corepw,rhomod,ifail)
 if(ifail/=0) return


!test model
 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
&                   zion,iexc,c_light,nc,nv,la,rcut,mmax,adams,ifail)
 if(ifail/=0) return

 write(iout,'(/a/)') 'Polynomial model core charge'
 write(iout,'(a,f12.6)') '  model radius r*',rstar
 call d2exc_report('d2excps - pseudofunction derivatives with core '// &
&                  'correction',d2excps,nv,d2mdiff,iout)

 deallocate(d2excae,d2excps)

 return

 end subroutine modcore

 end module modcore1_m
