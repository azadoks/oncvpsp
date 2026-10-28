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

 module der2exc_m

 implicit none
 private
 public :: der2exc

 contains

! computes the 2nd derivative of the contribution to the exchange-
! correlation energy from the region where the valence pseudo wave functions
! differ from the all-electron wave functions
!
 subroutine der2exc(rhotot,rhoc,rho,rr,d2exc,d2ref,d2mdiff, &
&                   zion,iexc,c_light,nc,nv,la,rcut,mmax)

!! @param rhotot  total valence charge, all-electron of pseudo
!! @param rhoc  core charge density, all-electron or model
!! @param rho  valence state-by-state charge (one-electron)
!! @param rr  log radial mesh (bohr)
!! @param d2exc  Exc 2nd-derivative matrix
!! @param d2ref  reference matrix
!! @param d2mdiff  root-mean-squared differrenc between d2exc and d2ref
!! @param zion  ion charge
!! @param iexc  exchange-correlation functional
!! @param c_light  speed of light (Hartree atomic units)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular momenta of the all-electron atom
!! @param rcut  radius to which the matrix elements are integrated, at least
!!              that beyond which all-electron and pseudo charges agree (bohr)
!! @param mmax  size of the log radial mesh

 use constants_m, only: dp
 use logmesh_m, only: lm_integral
 use vout_m, only: vout

!Input variables
 integer, intent(in) :: iexc,nc,nv,mmax
 real(dp), intent(in) :: rhotot(mmax),rhoc(mmax),rho(mmax,nv),rr(mmax)
 real(dp), intent(in) :: d2ref(nv,nv)
 real(dp), intent(in) :: zion,rcut
 integer, intent(in) :: la(nv+nc)
 real(dp), intent(in) :: c_light

!Output variables
 real(dp), intent(out) :: d2exc(nv,nv)
 real(dp), intent(out) :: d2mdiff

!Local variables
 real(dp) :: hh,eeel,eexc,ss
 real(dp), allocatable :: vo(:),vxct(:),rhot(:),dvxc(:,:)
 integer :: ii,jj,kk,l1

 allocate(vo(mmax),vxct(mmax),rhot(mmax),dvxc(mmax,nv))

 hh=0.15d0

 do kk=1,nv
   dvxc(:,kk)=0.0d0

   rhot(:)=rhotot(:)-2.0d0*hh*rho(:,kk)
   call vout(.true.,rhot,rhoc,vo,vxct,zion,eeel,eexc,rr,mmax,iexc,c_light)
   dvxc(:,kk)=dvxc(:,kk)+( 2.0d0/(24.0d0*hh))*vxct(:)

   rhot(:)=rhotot(:)     -hh*rho(:,kk)
   call vout(.true.,rhot,rhoc,vo,vxct,zion,eeel,eexc,rr,mmax,iexc,c_light)
   dvxc(:,kk)=dvxc(:,kk)+(-16.0d0/(24.0d0*hh))*vxct(:)

   rhot(:)=rhotot(:)     +hh*rho(:,kk)
   call vout(.true.,rhot,rhoc,vo,vxct,zion,eeel,eexc,rr,mmax,iexc,c_light)
   dvxc(:,kk)=dvxc(:,kk)+( 16.0d0/(24.0d0*hh))*vxct(:)

   rhot(:)=rhotot(:)+2.0d0*hh*rho(:,kk)
   call vout(.true.,rhot,rhoc,vo,vxct,zion,eeel,eexc,rr,mmax,iexc,c_light)
   dvxc(:,kk)=dvxc(:,kk)+(-2.0d0/(24.0d0*hh))*vxct(:)

 end do !kk

! compute Exc 2nd-derivative wrt occupation numbers matrix, integrated to
! the radius rcut (smooth in it)

 do kk=1,nv
   l1=la(nc+kk)+1
   do jj=1,nv
     rhot(:)=rho(:,kk)*rr(:)**2*dvxc(:,jj)
     call lm_integral(rr,rhot,mmax,2*l1,0.0d0,rcut,d2exc(kk,jj))
   end do !jj
 end do !kk

 ss=0.0d0
 do kk=1,nv
   do jj=1,nv
     ss=ss+(d2exc(jj,kk)-d2ref(jj,kk))**2
   end do
 end do
 d2mdiff=sqrt(ss/nv**2)

 deallocate(vo,vxct,rhot,dvxc)
 return
 end subroutine der2exc

 end module der2exc_m
