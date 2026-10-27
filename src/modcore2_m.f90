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

 module modcore2_m

 implicit none
 private
 public :: modcore2

 contains

! Creates a Teter-function model core charge whose amplitude and scale fit
! the value and slope of the all-electron core charge at the radius r* where
! it equals fcfact times the valence pseudocharge (icmod=2).
! r* is the root of ln(rhoc) - ln(fcfact rhotps), the value and slope at r*
! are interpolated (lm_derivs), the scale is found by Newton's method, and
! the Teter function and its derivatives are analytic everywhere, tapered
! to 0 by a C-infinity switch between x = 3 and 4 instead of cut at 3
! (modcore_m, doc/smoothness_plan.md M-1, M-2, M-4, M-5).

! Performs analysis and based on "hardness" criterion described in
! Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as

 subroutine modcore2(icmod,rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
&                   fcfact,rcfact,mmax,rr,nc,nv,la,zion,iexc,c_light,iout,ifail)

!! @param icmod  (not used)
!! @param rhops  state-by-state pseudocharge density
!! @param rhotps  total pseudocharge density
!! @param rhoc  core-charge density
!! @param rhoae  state-by-state all-electron valence charge density
!! @param rhotae  total all-electron valence charge density
!! @param rhomod  model core density and 4 derivatives
!! @param fcfact  the model is matched where rhoc equals fcfact times the total
!!                pseudocharge
!! @param rcfact  (not used)
!! @param mmax  dimension of log grid
!! @param rr  log radial grid
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular-momenta
!! @param zion  ion charge
!! @param iexc  exchange-correlation function to be used
!! @param c_light  speed of light (Hartree atomic units)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxstate, nrhomod
 use logmesh_m, only: lm_derivs
 use modcore_m, only: core_crossover, teter_derivs, teter_core_derivs, &
&                     teter_scale, d2exc_baseline, d2exc_report
 use der2exc_m, only: der2exc
 implicit none

!Input variables
 integer, intent(in) :: icmod,nv,nc,iexc,mmax
 real(dp), intent(in) :: c_light
 integer, intent(in) :: la(mxstate)
 real(dp), intent(in) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: zion,fcfact,rcfact

!Output variables
 real(dp), intent(out) :: rhomod(mmax,nrhomod)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: d2mdiff,rmatch,rhocmatch,rstar
 real(dp) :: xx,a0,b0,ymatch
 real(dp) :: fmatch(2),gd(0:4)
 real(dp), allocatable :: d2excae(:,:),d2excps(:,:)
 integer :: ii,ierr,kk

 ifail=0
 allocate(d2excae(nv,nv),d2excps(nv,nv))

 d2excae(:,:)=0.0d0
 d2excps(:,:)=0.0d0

! find valence pseudocharge - core charge crossover
! this is needed for compatability with icmod=3 option
 call core_crossover(rr,rhoc,rhotps,mmax,1.0d0,rmatch,ierr)
 if(ierr/=0) then
   write(iout,'(/2a)') 'modcore2: ERROR core charge never exceeds', &
&        ' valence pseudocharge'
   ifail=1
   return
 end if
 call lm_derivs(rr,rhoc,mmax,rmatch,0,fmatch)
 rhocmatch=fmatch(1)

! find scaled valence pseudocharge - core charge crossover
 call core_crossover(rr,rhoc,rhotps,mmax,fcfact,rstar,ierr)
 if(ierr/=0) then
   write(iout,'(/a)') 'modcore2: ERROR ircc (core-valence charge crossover) &
&        not found'
   ifail=1
   return
 end if

 call d2exc_baseline(rhotae,rhoc,rhoae,rhotps,rhops,rr,zion,iexc,c_light,nc, &
&                    nv,la,rr(mmax),mmax,rhomod,d2excae,d2excps,d2mdiff,iout)

! Fit Teter function to value and slope at r*

 call lm_derivs(rr,rhoc,mmax,rstar,1,fmatch)
 ymatch=rstar*fmatch(2)/fmatch(1)

 call teter_scale(ymatch,xx,ierr,iout,ifail)
 if(ifail/=0) return
 if(ierr/=0) then
   write(iout,'(/a,1p,e12.4)') 'modcore2: ERROR no Teter scale for the slope', &
&        ymatch
   ifail=1
   return
 end if
 call teter_derivs(xx,0,gd)

 b0=xx/rstar
 a0=fmatch(1)/gd(0)

 write(iout,'(/a)') 'amplitude prefactor, scale prefactor'
 write(iout,'(2f10.4)') a0/rhocmatch,1.0d0/(b0*rmatch)


 do ii=1,mmax
   call teter_core_derivs(b0*rr(ii),4,gd)
   do kk=0,4
     rhomod(ii,kk+1)=a0*gd(kk)*b0**kk
   end do
 end do

! test model
 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
&                   zion,iexc,c_light,nc,nv,la,rr(mmax),mmax)

 call d2exc_report('d2excps - pseudofunction derivatives with core '// &
&                  'correction',d2excps,nv,d2mdiff,iout)

 deallocate(d2excae,d2excps)

 return
 end subroutine modcore2

 end module modcore2_m
