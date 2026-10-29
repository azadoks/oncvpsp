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

 module exc_libxc_m

 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: exc_libxc, exc_libxc_fxc, exc_libxc_check, exc_libxc_info
 public :: have_libxc

!> whether this build has libxc, for check_data
 logical, parameter :: have_libxc=.true.

 contains

 subroutine exc_libxc_check(iexc,ierr)
!! whether Libxc can evaluate the functionals of iexc: ierr 0, or that
!! of lxc_init, ierr_lxc_unknown, ierr_lxc_no_vxc or ierr_lxc_family
!! (failure_m), for the first it cannot (check_data)
!!
!! @param iexc  exchange-correlation functional: -(1000*id1+id2) for the libxc
!!              functionals id1 and id2
!! @param ierr  0 on success; the non-zero values differ by routine:
!!              the reasons above

 use libxc_api_m, only: lxc_t, lxc_init, lxc_end

!Input variables
 integer, intent(in) :: iexc

!Output variables
 integer, intent(out) :: ierr

!Local variables
 type(lxc_t) :: fn
 integer :: id(2),ii

 call exc_libxc_ids(iexc,id)
 ierr=0
 do ii=1,2
   if(id(ii)==0) cycle
   call lxc_init(id(ii),fn,ierr)
   if(ierr/=0) return
   call lxc_end(fn)
 end do
 return
 end subroutine exc_libxc_check

 subroutine exc_libxc_ids(iexc,id)
!! the two Libxc numbers of iexc, 0 for none
!!
!! @param iexc  exchange-correlation functional: -(1000*id1+id2) for the libxc
!!              functionals id1 and id2
!! @param id  id1 and id2 (output)

!Input variables
 integer, intent(in) :: iexc

!Output variables
 integer, intent(out) :: id(2)

 id(1)=abs(iexc)/1000
 id(2)=abs(iexc) - id(1)*1000
 return
 end subroutine exc_libxc_ids

 subroutine exc_libxc_init(iexc,fns,ifail)
!! the functionals of iexc; check_data has rejected those Libxc cannot
!! evaluate (exc_libxc_check)
!!
!! @param iexc  exchange-correlation functional: -(1000*id1+id2) for the libxc
!!              functionals id1 and id2
!! @param fns  the functionals, family 0 for none (output)
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_libxc_init
 use libxc_api_m, only: lxc_t, lxc_init, lxc_end

!Input variables
 integer, intent(in) :: iexc

!Output variables
 type(lxc_t), intent(out) :: fns(2)
 integer, intent(out) :: ifail

!Local variables
 integer :: id(2),ii,ierr

 ifail=0
 call exc_libxc_ids(iexc,id)
 do ii=1,2
   if(id(ii)==0) cycle
   call lxc_init(id(ii),fns(ii),ierr)
   if(ierr/=0) then
     write(error_unit,'(/a,i0)') &
&         'exc_libxc_init: ERROR Libxc cannot evaluate functional ',id(ii)
! the first, if the second fails, is freed
     if(ii==2 .and. id(1)/=0) call lxc_end(fns(1))
     ifail=failed(fail_input,msg_libxc_init)
     return
   end if
 end do
 return
 end subroutine exc_libxc_init

 subroutine exc_libxc(iexc,al,rho,vxc,exc,rr,mmax,ifail)
!! the exchange-correlation potential and energy per electron of the
!! Libxc functionals of iexc. A GGA's potential is d e/d n less the
!! divergence of d e/d grad(n), the radial derivatives by derivs
!!
!! @param iexc  exchange-correlation functional: -(1000*id1+id2) for the libxc
!!              functionals id1 and id2
!! @param al  log of the mesh ratio amesh of the log radial mesh
!! @param rho  charge density
!! @param vxc  exchange-correlation potential (Ha)
!! @param exc  exchange-correlation energy per electron (Ha)
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp, pi
 use derivs_m, only: derivs
 use libxc_api_m, only: lxc_t, lxc_end, lxc_eval

!Input variables
 integer, intent(in) :: mmax, iexc
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rho(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: vxc(mmax),exc(mmax)
 integer, intent(out) :: ifail

!Local variables
 type(lxc_t) :: fns(2)
 integer :: ii
 real(dp), allocatable :: nn(:),dpr(:),dppr(:),dlap(:),sigma(:)
 real(dp), allocatable :: zk(:),vrho(:),vsigma(:),gg(:),dgg(:)

 allocate(nn(mmax),dpr(mmax),dppr(mmax),dlap(mmax),sigma(mmax))
 allocate(zk(mmax),vrho(mmax),vsigma(mmax),gg(mmax),dgg(mmax))

 call exc_libxc_init(iexc,fns,ifail)
 if(ifail/=0) return

! the electron density, without the 4 pi of rho, and its squared gradient
 nn(:)=rho(:)/4.0d0/pi
 call derivs(mmax,nn,al,rr,dpr,dppr,dlap)
 sigma(:)=dpr(:)**2

 vxc(:)=0.0d0
 exc(:)=0.0d0
 do ii=1,2
   if(fns(ii)%family==0) cycle
   call lxc_eval(fns(ii),mmax,nn,sigma,zk,vrho,vsigma)
   if(fns(ii)%family==2) then
     gg(:)=2.0d0*vsigma(:)*dpr(:)
     call derivs(mmax,gg,al,rr,dgg,dppr,dlap)
     vrho(:)=vrho(:) - (2.0d0/rr(:)*gg(:) + dgg(:))
   end if
   vxc(:)=vxc(:) + vrho(:)
   exc(:)=exc(:) + zk(:)
   call lxc_end(fns(ii))
 end do

 deallocate(nn,dpr,dppr,dlap,sigma)
 deallocate(zk,vrho,vsigma,gg,dgg)
 return
 end subroutine exc_libxc

 subroutine exc_libxc_fxc(iexc,al,rho,drho,dvxc,rr,mmax,have,ifail)
!! the change of exc_libxc's potential at rho for the change drho of the
!! density, from the functionals' second derivatives: the derivative of
!! the potential on the mesh, derivs included, to which its differences
!! converge. have is false, and dvxc zero, when a functional of iexc
!! lacks them (or Libxc was built without them)
!!
!! @param iexc  exchange-correlation functional: -(1000*id1+id2) for the libxc
!!              functionals id1 and id2
!! @param al  log of the mesh ratio amesh of the log radial mesh
!! @param rho  charge density
!! @param drho  change of the charge density
!! @param dvxc  change of the exchange-correlation potential (Ha)
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param have  whether dvxc was evaluated
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp, pi
 use derivs_m, only: derivs
 use libxc_api_m, only: lxc_t, lxc_end, lxc_eval, lxc_fxc

!Input variables
 integer, intent(in) :: mmax, iexc
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rho(mmax),drho(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: dvxc(mmax)
 logical, intent(out) :: have
 integer, intent(out) :: ifail

!Local variables
 type(lxc_t) :: fns(2)
 integer :: ii
 real(dp), allocatable :: nn(:),dn(:),dpr(:),ddpr(:),dppr(:),dlap(:)
 real(dp), allocatable :: sigma(:),dsigma(:),zk(:),vrho(:),vsigma(:)
 real(dp), allocatable :: v2rho2(:),v2rhosigma(:),v2sigma2(:)
 real(dp), allocatable :: dg(:),ddg(:),dv(:)

 dvxc(:)=0.0d0
 have=.false.
 call exc_libxc_init(iexc,fns,ifail)
 if(ifail/=0) return
 have=all(fns(:)%family==0 .or. fns(:)%has_fxc)
 if(.not. have) then
   call lxc_end(fns(1))
   call lxc_end(fns(2))
   return
 end if

 allocate(nn(mmax),dn(mmax),dpr(mmax),ddpr(mmax),dppr(mmax),dlap(mmax))
 allocate(sigma(mmax),dsigma(mmax),zk(mmax),vrho(mmax),vsigma(mmax))
 allocate(v2rho2(mmax),v2rhosigma(mmax),v2sigma2(mmax))
 allocate(dg(mmax),ddg(mmax),dv(mmax))

 nn(:)=rho(:)/4.0d0/pi
 dn(:)=drho(:)/4.0d0/pi
 call derivs(mmax,nn,al,rr,dpr,dppr,dlap)
 sigma(:)=dpr(:)**2
 call derivs(mmax,dn,al,rr,ddpr,dppr,dlap)
 dsigma(:)=2.0d0*dpr(:)*ddpr(:)

 do ii=1,2
   if(fns(ii)%family==0) cycle
   call lxc_fxc(fns(ii),mmax,nn,sigma,v2rho2,v2rhosigma,v2sigma2)
   if(fns(ii)%family==1) then
     dv(:)=v2rho2(:)*dn(:)
   else
! the change of d e/d n, less the divergence of the change of
! d e/d grad(n) = 2 vsigma n'
     call lxc_eval(fns(ii),mmax,nn,sigma,zk,vrho,vsigma)
     dv(:)=v2rho2(:)*dn(:) + v2rhosigma(:)*dsigma(:)
     dg(:)=2.0d0*((v2rhosigma(:)*dn(:) + v2sigma2(:)*dsigma(:))*dpr(:) &
&          + vsigma(:)*ddpr(:))
     call derivs(mmax,dg,al,rr,ddg,dppr,dlap)
     dv(:)=dv(:) - (2.0d0/rr(:)*dg(:) + ddg(:))
   end if
   dvxc(:)=dvxc(:) + dv(:)
   call lxc_end(fns(ii))
 end do

 deallocate(nn,dn,dpr,ddpr,dppr,dlap)
 deallocate(sigma,dsigma,zk,vrho,vsigma)
 deallocate(v2rho2,v2rhosigma,v2sigma2)
 deallocate(dg,ddg,dv)
 return
 end subroutine exc_libxc_fxc

 subroutine exc_libxc_info(iexc,iout,ifail)
!! the description of the libxc functionals of iexc, written once per run
!! (oncvpsp_run)
!!
!! @param iexc  exchange-correlation functional: -(1000*id1+id2) for the libxc
!!              functionals id1 and id2
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use libxc_api_m, only: lxc_t, lxc_end, lxc_info

!Input variables
 integer, intent(in) :: iexc,iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 type(lxc_t) :: fns(2)
 integer :: ii,jj,nref
 character(len=32) :: kind,family
 character(len=1024) :: name
 character(len=1024) :: refs(9)

 call exc_libxc_init(iexc,fns,ifail)
 if(ifail/=0) return
 do ii=1,2
   if(fns(ii)%family==0) cycle
   call lxc_info(fns(ii),kind,name,family,nref,refs)
   write(iout,'(a,2x,a)') 'LXC ',trim(kind)
   write(iout,'(a,4x,4a)') 'LXC ',trim(name),' (',trim(family),')'
   do jj=1,nref
     write(iout,'(a,4x,a,i1,2a)') 'LXC ','[',jj,'] ',trim(refs(jj))
   end do
   if(.not. fns(ii)%has_exc) then
     write(iout,'(a,4x,a)') 'LXC ', &
&         'no energy per electron: its part of E_xc is left out'
   end if
   call lxc_end(fns(ii))
 end do
 write(iout,'(a)') ''
 return
 end subroutine exc_libxc_info

 end module exc_libxc_m
