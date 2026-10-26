!
! Copyright (c) 1989-2012 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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

 module ldiracfs_m

 implicit none
 private
 public :: ldiracfs

 contains

 subroutine ldiracfs(nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,mch)

!! Finds relativistic scattering states of an al-electron potential
!!
!! @param nn  effective principal quantum number based on nodes inside mch
!!            (output)
!! @param ll  angular-momentum quantum number
!! @param kap  =l, -(l+1) for j=l -/+ 1/2
!! @param ierr  non-zero return if error
!! @param ee  bound-state energy, input guess and output calculated value
!! @param rr  log radial mesh
!! @param zz  atomic number
!! @param vv  local psp
!! @param uu(mmax,jj)  output radial wave functions (*rr) jj=1,2 for large,
!!                     small
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of log grid
!! @param mch  matching mesh point for inward-outward integrations

 use constants_m, only: dp
 use relativity_m, only: c_light
 use dirac_m, only: dirac_outward
 use quad_m, only: log_norm
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: ll,kap,mmax
 real(dp), intent(in) :: ee
 integer, intent(in) :: mch

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 integer, intent(inout) :: nn

!Output variables
 integer, intent(out) :: ierr

!Local variables

 real(dp), allocatable :: gu(:),fu(:),gup(:),fup(:),cf(:)

 real(dp) :: cc,cci,gam,cof
 real(dp) :: sn,cn,uout,upout
 real(dp) :: amesh,al
 integer :: ii,node

 cc=c_light
 cci=1.0d0/cc

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

 ierr = 60

! check arguments
 if(kap/=ll .and. kap/=-(ll+1)) then
   write(6,'(/a,i4,a,i4)') 'ldiracfs: ERROR kap =',kap,' ll =',ll
   ierr=2
   return
 end if
 if(zz<1.0d0) then
   write(6,'(/a,f12.8)') 'ldiracfs: ERROR zz =',zz
   ierr=3
   return
 end if

 allocate(gu(mmax),fu(mmax),gup(mmax),fup(mmax),cf(mmax))

! null arrays
 gu(:)=0.0d0; fu(:)=0.0d0; gup(:)=0.0d0; fup(:)=0.0d0; cf(:)=0.0d0

 node=0
 ierr=0

 gam=sqrt(kap**2-(zz*cci)**2)

 cof=(kap+gam)*(cc/zz)

! outward from the series at the origin (dirac_m)
 call dirac_outward(ll,kap,ee,zz,al,rr,vv,mmax,mch,gu,fu,gup,fup,node)

 uout=gu(mch)
 upout=gup(mch)


! perform normalization sum

 sn=log_norm(rr,gu(1:mch)**2+fu(1:mch)**2,mch,ll,gam)

! normalize u

 cn=1.0d0/sqrt(sn)
 uout=cn*uout
 upout=cn*upout

 do ii=1,mch
   gup(ii)=cn*gup(ii)
   gu(ii)=cn*gu(ii)
   fup(ii)=cn*fup(ii)
   fu(ii)=cn*fu(ii)
 end do
 do ii=mch+1,mmax
   gu(ii)=0.0d0
   fu(ii)=0.0d0
 end do

! copy local arrays for output
 uu(:,1)=gu(:)
 up(:,1)=gup(:)
 uu(:,2)=fu(:)
 up(:,2)=fup(:)

!calculate effective principal quantum number as if this were a bound
!state with a barrier at mch
 nn=node+ll+1

 deallocate(gu,fu,gup,fup)
 deallocate(cf)
 return

 end subroutine ldiracfs

 end module ldiracfs_m
