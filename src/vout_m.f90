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

 module vout_m

 implicit none
 private
 public :: vout

 contains

! calculates Hartree and exchange-correlation potentials and total-energy
! term to add to eigenvalue sum

 subroutine vout(add_core,rho,rhoc,vo,vxc,zion,eeel,eexc, &
&                rr,mmax,adams,iexc,c_light)

!! @param add_core  whether the exchange-correlation is of rho+rhoc, or of
!!                  rho alone
!! @param rho  total charge density or valence/pseudovalence charge density
!! @param rhoc  core charge density, full or model
!! @param vo  output total potential
!! @param vxc  exchange-correlation potential (Ha), output
!! @param zion  ion charge
!! @param eeel  electron-electron interaction energy
!! @param eexc  exchange-correlation correction to eigenvalue sum for total
!!              energy
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton quadrature of
!!               the Hartree potential
!! @param iexc  exchange-correlation functional (1-4, 6, 7 built in, <0 libxc)
!! @param c_light  speed of light (Hartree atomic units)

 use constants_m, only: dp
 use logmesh_m, only: lm_head, lm_step
 use xc_dispatch_m, only: exc_dispatch
 use quad_m, only: quad_gregory_sum
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: iexc,mmax
 logical, intent(in) :: add_core
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: rho(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: zion
 type(adams_t), intent(in) :: adams

!Input/output variables
 real(dp), intent(inout) :: vxc(mmax)

!Output variables
 real(dp), intent(out) :: vo(mmax)
 real(dp), intent(out) :: eeel,eexc

!Local variables
 integer :: ii,i0
 real(dp) :: al,tv

 real(dp), allocatable :: rvp(:),rv(:),excca(:),difxc(:)
 real(dp), allocatable :: exca(:),vxcd(:),rhot(:)
 real(dp), allocatable :: foutv(:),foutx(:),foutfc(:)

 allocate(rvp(mmax),rv(mmax),excca(mmax),difxc(mmax))
 allocate(exca(mmax),vxcd(mmax),rhot(mmax))
 allocate(foutv(mmax),foutx(mmax),foutfc(mmax))


 al=lm_step(rr,mmax)

 foutfc(:)=0.0d0

! integration for electrostatic potential
 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**3
 end do

! the cumulative Adams-Moulton rule of order k inward from the last k-1
! points
 i0=mmax-adams%order+2
 rv(i0:mmax)=zion

 do ii=i0,2,-1
   rv(ii-1)=rv(ii)+am_in(adams,rvp,ii)-adams%cm(0)*rvp(ii-1)
 end do

 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**2
 end do

 tv=0.0d0
 do ii=i0,2,-1
   tv=tv+am_in(adams,rvp,ii)-adams%cm(0)*rvp(ii-1)
   rv(ii-1)=rv(ii-1)-rr(ii-1)*tv
 end do

 do ii=1,mmax
   vo(ii)=rv(ii)/rr(ii)
 end do

! electron-electron interaction for total energy
 do ii=1,mmax
   foutv(ii)=rho(ii)*vo(ii)*rr(ii)**3
 end do

 eeel=al*quad_gregory_sum(foutv,mmax) + lm_head(foutv(1),rho(1),rho(2),al,3.0d0)

 if(add_core) then
   do ii = 1, mmax
     rhot(ii) = rho(ii) + rhoc(ii)
   end do
 else
   do ii = 1, mmax
     rhot(ii) = rho(ii)
   end do
 end if

! exchange-correlation potential added

 call exc_dispatch(iexc,c_light,al,rhot,vxc,exca,rr,mmax,'vout')

 do ii = 1, mmax
   difxc(ii)=exca(ii) - vxc(ii)
   vo(ii)=vo(ii) + vxc(ii)
 end do

! exchange-correlation correction for total energy

 do ii=1,mmax
   foutx(ii)=rho(ii)*difxc(ii)*rr(ii)**3
 end do
 eexc=quad_gregory_sum(foutx,mmax)

 if(add_core .and. rhoc(1) /= 0.0d0) then
   call exc_dispatch(iexc,c_light,al,rhoc,vxcd,excca,rr,mmax,'vout')

   do ii=1,mmax
     foutfc(ii)=rhoc(ii)*(exca(ii) - excca(ii))*rr(ii)**3
   end do
   eexc=eexc + quad_gregory_sum(foutfc,mmax)
 end if

 eexc=al*eexc + lm_head(foutx(1),rho(1),rho(2),al,3.0d0)


 if(add_core) then
   eexc=eexc + lm_head(foutfc(1),rhoc(1),rhoc(2),al,3.0d0)
 end if

 deallocate(rvp,rv,excca,difxc)
 deallocate(exca,vxcd,rhot)
 deallocate(foutv,foutx,foutfc)
 return

 contains

 include 'adams.inc'

 end subroutine vout

 end module vout_m
