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
 subroutine lschfs(nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,irel)

!! integrates radial Pauli-type scalar-relativistic equation
!! on a logarithmic mesh
!! modified routine to be used in finding norm-conserving
!! pseudopotential
!!
!! The Koelling-Harmon equation is integrated outward in its first-order
!! flux form (sr_flux_m), as for the bound states (lschfb)
!!
!! @param nn  effective principal quantum number based on nodes inside mch
!!            (output)
!! @param ll  angular-momentum quantum number
!! @param ierr  non-zero return if error
!! @param ee  bound-state energy, input guess and output calculated value
!! @param rr  log radial mesh
!! @param vv  local atomic potential
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param zz  atomic number
!! @param mmax  size of log grid
!! @param mch  matching mesh point for inward-outward integrations
!! @param irel  0 non-relativistic, 1 scalar-relativistic

 use constants_m, only: dp
 use relativity_m, only: c_light
 use sr_flux_m, only: sr_gamma, sr_flux_coefs, sr_outward
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: ll,mch
 integer, intent(in) :: irel
 real(dp), intent(in) :: ee

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ierr,nn


!Local variables
 real(dp) :: amesh,al
 real(dp) :: cn
 real(dp) :: fss, gamma, ro
 real(dp) :: sn
 integer :: ii,node

 real(dp), allocatable :: qq(:),qp(:),pm(:),qm(:),xm(:)

 allocate(qq(mmax),qp(mmax),pm(mmax),qm(mmax),xm(mmax))


 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

 ierr = 0

! relativistic - non-relativistic switch
 if(irel/=0) then
   fss=(1.0d0/c_light)**2
 else
   fss=0.0d0
 end if

 gamma=sr_gamma(ll,zz,fss)

! null arrays to remove leftover garbage

 uu(:)=0.0d0
 up(:)=0.0d0
 qq(:)=0.0d0
 qp(:)=0.0d0

 call sr_flux_coefs(ll,ee,fss,al,rr,vv,mmax,xm,pm,qm)
 call sr_outward(ll,ee,zz,fss,al,rr,vv,pm,qm,mmax,mch, &
&                uu,up,qq,qp,node)

!perform normalization sum

 ro=rr(1)/sqrt(amesh)
 sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)

 do ii=1,mch-3
   sn=sn+al*rr(ii)*uu(ii)**2
 end do

 sn =sn + al*(23.0d0*rr(mch-2)*uu(mch-2)**2 &
&           + 28.0d0*rr(mch-1)*uu(mch-1)**2 &
&          +  9.0d0*rr(mch  )*uu(mch  )**2)/24.0d0

!normalize u

 cn=1.0d0/sqrt(sn)

 do ii=1,mch
   up(ii)=cn*up(ii)
   uu(ii)=cn*uu(ii)
 end do
 do ii=mch+1,mmax
   uu(ii)=0.0d0
 end do

!calculate effective principal quantum number as if this were a bound
!state with a barrier at mch
 nn=node+ll+1

 deallocate(qq,qp,pm,qm,xm)

 return
 end  subroutine lschfs
