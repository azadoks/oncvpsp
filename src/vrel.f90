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
 subroutine vrel(ll,ee,rr,vv,vr,uu,up,zz,mmax,irc,irel)

! Finds an effective potential representing the Pauli-type scalar-relativistic
! operators in the radial Schroedinger equation

!ll  angular-momentum quantum number
!ee  bound-state energy
!rr  log radial mesh
!vv  local atomic potential
!vr  effective scalar-relativistic potential, set for
!    0.5*rr(irc)<=rr<=1.5*rr(irc) and zero elsewhere
!uu  radial wave function (*rr)
!up  d(uu)/dr
!zz  atomic number
!mmax  size of log grid
!irc  mesh index of the radius around which vr is calculated
!irel  0 non-relativistic (vr is then not set), 1 scalar-, 2 fully-
!      relativistic

 use constants_m, only: dp, c_light
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: uu(mmax),up(mmax)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: ee
 integer, intent(in) :: ll,irc
 integer, intent(in) :: mmax
 integer, intent(in) :: irel

!Output variables
 real(dp), intent(out) :: vr(mmax)

!Local variables
 integer :: nn,ierr

 real(dp) :: als
 real(dp) :: de,emax,emin
 real(dp) :: eps,fss,tfss,gamma,ro,sc
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap
 real(dp) :: amesh,al,xx
 integer :: ii,it,nint,node,nin

 real(dp), allocatable :: cf(:),dv(:),fr(:),frp(:)


 allocate(cf(mmax),dv(mmax),fr(mmax),frp(mmax))

 vr(:)=0.0d0

! relativistic - non-relativistic switch
 if(irel==0) return

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
 eps=1.0d-10

 vr(:)=0.0d0

 fss=(1.0d0/c_light)**2

 if(ll==0) gamma=sqrt(1.0d0-fss*zz**2)
 if(ll>0) gamma=(ll*sqrt(ll**2-fss*zz**2) + &
& (ll+1)*sqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

 sls=ll*(ll+1)

 als=al**2

! coefficient array for u in differential eq.
 do ii=1,mmax
   cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
 end do

! calculate dv/dr for darwin correction
 dv(:)=0.0d0

 dv(1)=(-50.d0*vv(1)+96.d0*vv(2)-72.d0*vv(3)+32.d0*vv(4) &
&       -6.d0*vv(5))/(24.d0*al*rr(1))
 dv(2)=(-6.d0*vv(1)-20.d0*vv(2)+36.d0*vv(3)-12.d0*vv(4) &
&       +2.d0*vv(5))/(24.d0*al*rr(2))

 do ii=3,mmax-2
   dv(ii)=(2.d0*vv(ii-2)-16.d0*vv(ii-1)+16.d0*vv(ii+1) &
&        -2.d0*vv(ii+2))/(24.d0*al*rr(ii))
 end do

!  relativistic coefficient arrays for u (fr) and up (frp).
 do ii=1,mmax
   tfss=fss
   fr(ii)=als*(rr(ii)**2)*(-tfss*(vv(ii)-ee)**2 + 0.5d0*tfss*dv(ii)/ &
&   (rr(ii)*(1.0d0+0.5d0*tfss*(ee-vv(ii)))))
   frp(ii)=-al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0+0.5d0*tfss*(ee-vv(ii)))
 end do

 do ii=1,mmax
   if(rr(ii)<0.5d0*rr(irc) .or. rr(ii)>1.5d0*rr(irc)) cycle
   if(abs(uu(ii))>eps) then
     vr(ii)=0.5d0*(fr(ii) + frp(ii)*up(ii)/uu(ii))/(als*rr(ii)**2)
   end if
 end do

 do ii=1,mmax
   if(rr(ii)<0.5d0*rr(irc) .or. rr(ii)>1.5d0*rr(irc)) cycle
   if(vr(ii)==0.0d0) vr(ii)=0.5d0*(vr(ii-1)+vr(ii+1))
 end do

 deallocate(cf,dv,fr,frp)
 return

 end subroutine vrel
