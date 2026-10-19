!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University, and M. Verstratte, University of Liege
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

 subroutine exc_libxc(iexc,al,rho,vxc,exc,rr,mmax)
!
!  interface to libxc library for potential and energy density calculation
!
!iexc  exchange-correlation function to be used: -(1000*id1+id2) for
!      the libxc functionals id1 and id2
!al  log of the mesh ratio amesh of the log radial mesh
!rho  charge density
!vxc  exchange-correlation potential
!exc  exchange-correlation energy density
!rr  log radial mesh
!mmax  dimension of log radial grid

 use functionals_m, only: xc_functl_t, xc_functl_init_functl, &
&                         xc_functl_get_vxc, xc_functl_write_info, &
&                         xc_functl_end, XC_DERIV_NUMERICAL
 use constants_m, only: dp, pi
 implicit none

!Input variables
 integer, intent(in) :: mmax, iexc
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rho(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: vxc(mmax),exc(mmax)

!Local variables
 integer :: ii, unit
 integer ::  id(2), ifunc
 integer, save :: functprinted = 0
 type(xc_functl_t) :: functls(2)
 integer :: nspin, irel, deriv_method
 real(dp),allocatable :: ip(:)

 real(dp), allocatable :: vv(:), ee(:)
 real(dp), allocatable :: dpr(:), dppr(:), dlap(:)
 real(dp), allocatable :: tau_dummy(:), vtau_dummy(:)
 real(dp), allocatable :: rho_loc(:)

 allocate(vv(mmax), ee(mmax))
 allocate(dpr(mmax), dppr(mmax), dlap(mmax))
 allocate(tau_dummy(mmax), vtau_dummy(mmax))
 allocate(rho_loc(mmax))

 rho_loc = rho(1:mmax) /4.0d0 / pi

 call derivs(mmax, rho_loc, al, rr, dpr, dppr, dlap)

!
! determine xc to be used from input code (negative by convention)

 id(1) = abs(iexc)/1000
 id(2) = abs(iexc) - id(1)*1000
 nspin = 1
 irel = 0 !XC_NON_RELATIVISTIC
!numerical derivatives (XC_DERIV_ANALYTICAL is the alternative)
 deriv_method = XC_DERIV_NUMERICAL
 call xc_functl_init_functl(functls(1),id(1),3,0.d0,nspin,deriv_method)
 call xc_functl_init_functl(functls(2),id(2),3,0.d0,nspin,deriv_method)

 if(functprinted==0) then
   call xc_functl_write_info(functls(1), 6)
   call xc_functl_write_info(functls(2), 6)
   write(*,*)
   functprinted = 1
 end if

! no use of vtau kin en density optional argument
! ip is an ionization potential parameter for certain xc functionals, not
! supported here. Set to 0.5 Ha
 allocate(ip(nspin))
 ip = 0.5_dp
 tau_dummy = 0.0d0
 vtau_dummy = 0.0d0
 vxc(1:mmax) = 0.0d0
 exc(1:mmax) = 0.0d0
 do ifunc = 1, 2
   vv = 0.0d0
   ee = 0.0d0
   call xc_functl_get_vxc(functls(ifunc), mmax, al, rr, rho_loc, dpr, dlap, &
&         tau_dummy, ip, vv, ee, vtau_dummy)
   vxc(1:mmax) = vxc(1:mmax) + vv(1:mmax)
   exc(1:mmax) = exc(1:mmax) + ee(1:mmax)
 end do

! convert to the convention of D. R. Hamann's built-in routines (excwig
! etc.).  Why this is needed is not understood: in principle libxc
! outputs the energy per particle, yet dividing by rho once more here
! gives an LDA exc identical to that of the built-in routines.
 do ii=1,mmax
   exc(ii) = exc(ii) / max(rho_loc(ii),1.d-20)
 end do

! finish
 call xc_functl_end(functls(1))
 call xc_functl_end(functls(2))
 deallocate(ip)
 deallocate(vv, ee)
 deallocate(dpr, dppr, dlap)
 deallocate(tau_dummy, vtau_dummy)
 deallocate(rho_loc)
 return

 end subroutine exc_libxc


 subroutine derivs(mmax, rho, al, rr, dpr, dppr, dlap)
! first and second radial derivatives and Laplacian of rho on the log
! mesh, by finite differences, for exc_libxc
! NB: presumes incoming rho is properly scaled without 4pi or r**2 factors.

!mmax  dimension of log radial grid
!rho  charge density
!al  log of the mesh ratio amesh of the log radial mesh
!rr  log radial mesh
!dpr  first radial derivative of rho
!dppr  second radial derivative of rho
!dlap  Laplacian of rho

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: rho(mmax)

!Output variables
 real(dp), intent(out) :: dpr(mmax), dppr(mmax), dlap(mmax)

!Local variables
 integer :: ii
 real(dp) :: dpn(mmax), dppn(mmax)
 real(dp) :: c11,c12,c13,c14,c15
 real(dp) :: c21,c22,c23,c24,c25

 c11 =   2.0d0 / 24.0d0
 c12 = -16.0d0 / 24.0d0
 c13 =   0.0d0 / 24.0d0
 c14 =  16.0d0 / 24.0d0
 c15 =  -2.0d0 / 24.0d0
!
 c21 =   -1.0d0 / 12.0d0
 c22 =   16.0d0 / 12.0d0
 c23 =  -30.0d0 / 12.0d0
 c24 =   16.0d0 / 12.0d0
 c25 =   -1.0d0 / 12.0d0
!
! n derivatives of d
!
 ii=1
 dpn(ii) = -25.d0/12.d0*rho(ii) +4.d0*rho(ii+1) -3.d0*rho(ii+2) &
&         +4.d0/3.d0*rho(ii+3) -1.d0/4.d0*rho(ii+4)
 dppn(ii) = 15.d0/4.d0*rho(ii) -77.d0/6.d0*rho(ii+1) +107.d0/6.d0*rho(ii+2) &
&         -13.d0*rho(ii+3) +61.d0/12.d0*rho(ii+4) -5.d0/6.d0*rho(ii+5)
 ii=2
 dpn(ii) = -25.d0/12.d0*rho(ii) +4.d0*rho(ii+1) -3.d0*rho(ii+2)  &
&         +4.d0/3.d0*rho(ii+3) -1.d0/4.d0*rho(ii+4)
 dppn(ii) = 15.d0/4.d0*rho(ii) -77.d0/6.d0*rho(ii+1) +107.d0/6.d0*rho(ii+2) &
&         -13.d0*rho(ii+3) +61.d0/12.d0*rho(ii+4) -5.d0/6.d0*rho(ii+5)

 do ii = 3, mmax - 2
   dpn(ii) =  c11*rho(ii-2) + c12*rho(ii-1) + c14*rho(ii+1) + c15*rho(ii+2)
   dppn(ii) = c21*rho(ii-2) + c22*rho(ii-1) + c23*rho(ii)   + c24*rho(ii+1) &
&           +c25*rho(ii+2)
 end do

 ii=mmax-1
 dpn(ii) = +25.d0/12.d0*rho(ii) -4.d0*rho(ii-1) +3.d0*rho(ii-2) &
&         -4.d0/3.d0*rho(ii-3) +1.d0/4.d0*rho(ii-4)
 dppn(ii) = -15.d0/4.d0*rho(ii) +77.d0/6.d0*rho(ii-1) -107.d0/6.d0*rho(ii-2) &
&          +13.d0*rho(ii-3) -61.d0/12.d0*rho(ii-4) +5.d0/6.d0*rho(ii-5)
 ii=mmax
 dpn(ii) = +25.d0/12.d0*rho(ii) -4.d0*rho(ii-1) +3.d0*rho(ii-2) &
&         -4.d0/3.d0*rho(ii-3) +1.d0/4.d0*rho(ii-4)
 dppn(ii) = -15.d0/4.d0*rho(ii) +77.d0/6.d0*rho(ii-1) -107.d0/6.d0*rho(ii-2) &
&          +13.d0*rho(ii-3) -61.d0/12.d0*rho(ii-4) +5.d0/6.d0*rho(ii-5)

!
! r derivatives of d
!
 do ii = 1, mmax
   dpr(ii) = dpn(ii) / (al * rr(ii))
   dppr(ii) = (dppn(ii) - al * dpn(ii)) / (al * rr(ii))**2
   dlap(ii) = (dppn(ii) + al * dpn(ii)) / (al * rr(ii))**2
 end do
 return

 end subroutine derivs
