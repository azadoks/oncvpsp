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
! calculates Hartree and exchange-correlation potentials and total-energy
! term to add to eigenvalue sum

 subroutine vout(mode,rho,rhoc,vo,vxc,zion,eeel,eexc, &
&                rr,mmax,iexc)

!! @param mode  1=> add rhoc to rho, 0 => don't
!! @param rho  total charge density or valence/pseudovalence charge density
!! @param rhoc  core charge density, full or model
!! @param vo  output total potential
!! @param vxc  output exchange-correlation potential
!! @param zion  charge of screened ion
!! @param eeel  electron-electron interaction energy
!! @param eexc  exchange-correlation correction to eigenvalue sum for total
!!              energy
!! @param rr  log radial mesh
!! @param mmax  dimension of log radial grid
!! @param iexc  exchange-correlation functional (1-4, 6, 7 built in, <0 libxc)

 use constants_m, only: dp, pi
 use logmesh_m, only: lm_head
 use xc_dispatch_m, only: exc_dispatch
 implicit none

!Input variables
 integer, intent(in) :: iexc,mmax,mode
 real(dp), intent(in) :: rho(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: zion

!Input/output variables
 real(dp), intent(inout) :: vxc(mmax)

!Output variables
 real(dp), intent(out) :: vo(mmax)
 real(dp), intent(out) :: eeel,eexc

!Local variables
 real(dp), parameter :: pi4=4.0_dp*pi
 integer :: ii
 real(dp) :: al,tv

 real(dp), allocatable :: rvp(:),rv(:),excca(:),difxc(:)
 real(dp), allocatable :: exca(:),vxcd(:),rhot(:)
 real(dp), allocatable :: foutv(:),foutx(:),foutfc(:)

 allocate(rvp(mmax),rv(mmax),excca(mmax),difxc(mmax))
 allocate(exca(mmax),vxcd(mmax),rhot(mmax))
 allocate(foutv(mmax),foutx(mmax),foutfc(mmax))


 al = 0.01d0 * log(rr(101) / rr(1))

 foutfc(:)=0.0d0

! integration for electrostatic potential
 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**3
 end do

 rv(mmax)=zion
 rv(mmax-1)=zion
 rv(mmax-2)=zion

 do ii=mmax-2,2,-1
   rv(ii-1)=rv(ii)+aii(rvp,ii)
 end do

 do ii=1,mmax
   rvp(ii)=rho(ii)*al*rr(ii)**2
 end do

 tv=0.0d0
 do ii=mmax-2,2,-1
   tv=tv+aii(rvp,ii)
   rv(ii-1)=rv(ii-1)-rr(ii-1)*tv
 end do

 do ii=1,mmax
   vo(ii)=rv(ii)/rr(ii)
 end do

! electron-electron interaction for total energy
 do ii=1,mmax
   foutv(ii)=rho(ii)*vo(ii)*rr(ii)**3
 end do

 eeel=(9.0d0*foutv(1) + 28.0d0*foutv(2) &
&   + 23.0d0*foutv(3))/24.0d0
 do ii=4,mmax
   eeel=eeel + foutv(ii)
 end do
 eeel=al*eeel + lm_head(foutv(1),rho(1),rho(2),al,3.0d0)

 if(mode == 0) then
   do ii = 1, mmax
     rhot(ii) = rho(ii)
   end do
 else if(mode == 1) then
   do ii = 1, mmax
     rhot(ii) = rho(ii) + rhoc(ii)
   end do
 else
   write(6,'(/a,i4)') 'vout: ERROR bad input mode =',mode
   error stop 1
 end if

! exchange-correlation potential added

 call exc_dispatch(iexc,al,rhot,vxc,exca,rr,mmax,'vout')

 do ii = 1, mmax
   difxc(ii)=exca(ii) - vxc(ii)
   vo(ii)=vo(ii) + vxc(ii)
 end do

! exchange-correlation correction for total energy

 do ii=1,mmax
   foutx(ii)=rho(ii)*difxc(ii)*rr(ii)**3
 end do
 eexc=(9.0d0*foutx(1) + 28.0d0*foutx(2) &
&   + 23.0d0*foutx(3))/24.0d0
 do ii=4,mmax
   eexc=eexc + foutx(ii)
 end do

 if(mode == 1 .and. rhoc(1) /= 0.0d0) then
   call exc_dispatch(iexc,al,rhoc,vxcd,excca,rr,mmax,'vout')

   do ii=1,mmax
     foutfc(ii)=rhoc(ii)*(exca(ii) - excca(ii))*rr(ii)**3
   end do
   eexc=eexc + (9.0d0*foutfc(1) + 28.0d0*foutfc(2) &
&     + 23.0d0*foutfc(3))/24.0d0
   do ii=4,mmax
     eexc=eexc + foutfc(ii)
   end do
 end if

 eexc=al*eexc + lm_head(foutx(1),rho(1),rho(2),al,3.0d0)


 if(mode == 1) then
   eexc=eexc + lm_head(foutfc(1),rhoc(1),rhoc(2),al,3.0d0)
 end if

 deallocate(rvp,rv,excca,difxc)
 deallocate(exca,vxcd,rhot)
 deallocate(foutv,foutx,foutfc)
 return

 contains

 include 'adams.inc'

 end subroutine vout
