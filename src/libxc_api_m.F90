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
!> the calls of Libxc that ONCVPSP makes, spin-unpolarized, and the only
!> place that knows which Fortran interface Libxc has: xc_f90 before
!> version 5, xc_f03 from it. ONCVPSP_LIBXC_FXC is defined (by CMake)
!> when the library has the second derivatives (fxc) compiled in.
 module libxc_api_m

#include "xc_version.h"

#if XC_MAJOR_VERSION<=4
 use xc_f90_types_m, only: xc_f90_pointer_t
#else
 use xc_f03_lib_m, only: xc_f03_func_t, xc_f03_func_info_t
#endif
 implicit none
 private
 public :: lxc_t, lxc_init, lxc_end, lxc_eval, lxc_fxc, lxc_info

!> one Libxc functional
 type lxc_t
!> its Libxc number
   integer :: id=0
!> family: 1 LDA, 2 GGA, any other value not usable
   integer :: family=0
!> whether it gives the energy per electron
   logical :: has_exc=.false.
!> whether it gives the second derivatives (and they are compiled in)
   logical :: has_fxc=.false.
#if XC_MAJOR_VERSION<=4
!> Libxc's handle of the functional
   type(xc_f90_pointer_t) :: conf
!> Libxc's handle of its description
   type(xc_f90_pointer_t) :: info
#else
!> Libxc's handle of the functional
   type(xc_f03_func_t) :: conf
!> Libxc's handle of its description
   type(xc_f03_func_info_t) :: info
#endif
 end type lxc_t

 contains

 subroutine lxc_init(id,fn,ierr)
!! the functional number id, unpolarized; ierr ierr_lxc_unknown for a
!! number Libxc does not know, ierr_lxc_no_vxc for one that gives no
!! potential, ierr_lxc_family for a family other than LDA and GGA
!! (meta-GGA, hybrid, ...) (failure_m)
!!
!! @param id  the Libxc number of the functional
!! @param fn  the functional (output)
!! @param ierr  0 on success; the non-zero values differ by routine:
!!              the reasons above

 use failure_m, only: ierr_lxc_family, ierr_lxc_no_vxc, ierr_lxc_unknown

#if XC_MAJOR_VERSION<=4
 use xc_f90_lib_m, only: xc_f90_family_from_id, xc_f90_func_init, &
&                        xc_f90_info_flags, XC_UNPOLARIZED, &
&                        XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FLAGS_HAVE_EXC, &
&                        XC_FLAGS_HAVE_VXC
#ifdef ONCVPSP_LIBXC_FXC
 use xc_f90_lib_m, only: XC_FLAGS_HAVE_FXC
#endif
#else
 use xc_f03_lib_m, only: xc_f03_family_from_id, xc_f03_func_init, &
&                        xc_f03_func_get_info, xc_f03_func_info_get_flags, &
&                        XC_UNPOLARIZED, XC_FAMILY_LDA, XC_FAMILY_GGA, &
&                        XC_FLAGS_HAVE_EXC, XC_FLAGS_HAVE_VXC
#ifdef ONCVPSP_LIBXC_FXC
 use xc_f03_lib_m, only: XC_FLAGS_HAVE_FXC
#endif
#endif

!Input variables
 integer, intent(in) :: id

!Output variables
 type(lxc_t), intent(out) :: fn
 integer, intent(out) :: ierr

!Local variables
 integer :: family,flags

 ierr=0
 fn%id=id
#if XC_MAJOR_VERSION<=4
 family=xc_f90_family_from_id(id)
#else
 family=xc_f03_family_from_id(id)
#endif
 if(family==XC_FAMILY_LDA) then
   fn%family=1
 else if(family==XC_FAMILY_GGA) then
   fn%family=2
 else if(family<=0) then
   ierr=ierr_lxc_unknown
   return
 else
   ierr=ierr_lxc_family
   return
 end if

#if XC_MAJOR_VERSION<=4
 call xc_f90_func_init(fn%conf,fn%info,id,XC_UNPOLARIZED)
 flags=xc_f90_info_flags(fn%info)
#else
 call xc_f03_func_init(fn%conf,id,XC_UNPOLARIZED)
 fn%info=xc_f03_func_get_info(fn%conf)
 flags=xc_f03_func_info_get_flags(fn%info)
#endif
 fn%has_exc=iand(flags,XC_FLAGS_HAVE_EXC)/=0
#ifdef ONCVPSP_LIBXC_FXC
 fn%has_fxc=iand(flags,XC_FLAGS_HAVE_FXC)/=0
#else
 fn%has_fxc=.false.
#endif
 if(iand(flags,XC_FLAGS_HAVE_VXC)==0) then
   call lxc_end(fn)
   ierr=ierr_lxc_no_vxc
 end if
 return
 end subroutine lxc_init

 subroutine lxc_end(fn)
!! free the functional of lxc_init
!!
!! @param fn  the functional

#if XC_MAJOR_VERSION<=4
 use xc_f90_lib_m, only: xc_f90_func_end
#else
 use xc_f03_lib_m, only: xc_f03_func_end
#endif

!Input/output variables
 type(lxc_t), intent(inout) :: fn

 if(fn%family==0) return
#if XC_MAJOR_VERSION<=4
 call xc_f90_func_end(fn%conf)
#else
 call xc_f03_func_end(fn%conf)
#endif
 fn%family=0
 return
 end subroutine lxc_end

 subroutine lxc_eval(fn,mmax,nn,sigma,zk,vrho,vsigma)
!! the energy per electron and the first derivatives of the energy
!! density at mmax points; zk is zero for a functional without the
!! energy, and sigma and vsigma are for a GGA only
!!
!! @param fn  the functional
!! @param mmax  size of the log radial mesh
!! @param nn  electron density (bohr^-3)
!! @param sigma  squared gradient of nn
!! @param zk  energy per electron (Ha)
!! @param vrho  derivative of the energy density by nn
!! @param vsigma  derivative of the energy density by sigma

#if XC_MAJOR_VERSION<=4
 use xc_f90_lib_m, only: xc_f90_lda_exc_vxc, xc_f90_lda_vxc, &
&                        xc_f90_gga_exc_vxc, xc_f90_gga_vxc
#else
 use, intrinsic :: iso_c_binding, only: c_size_t
 use xc_f03_lib_m, only: xc_f03_lda_exc_vxc, xc_f03_lda_vxc, &
&                        xc_f03_gga_exc_vxc, xc_f03_gga_vxc
#endif
 use constants_m, only: dp

!Input variables
 type(lxc_t), intent(in) :: fn
 integer, intent(in) :: mmax
 real(dp), intent(in) :: nn(mmax),sigma(mmax)

!Output variables
 real(dp), intent(out) :: zk(mmax),vrho(mmax),vsigma(mmax)

!Local variables
#if XC_MAJOR_VERSION<=4
 integer :: ii
#else
 integer(c_size_t) :: np
#endif

 zk(:)=0.0d0; vrho(:)=0.0d0; vsigma(:)=0.0d0
#if XC_MAJOR_VERSION<=4
! one point a call: the xc_f90 interfaces take scalars
 do ii=1,mmax
   if(fn%family==1 .and. fn%has_exc) then
     call xc_f90_lda_exc_vxc(fn%conf,1,nn(ii),zk(ii),vrho(ii))
   else if(fn%family==1) then
     call xc_f90_lda_vxc(fn%conf,1,nn(ii),vrho(ii))
   else if(fn%has_exc) then
     call xc_f90_gga_exc_vxc(fn%conf,1,nn(ii),sigma(ii),zk(ii),vrho(ii), &
&                            vsigma(ii))
   else
     call xc_f90_gga_vxc(fn%conf,1,nn(ii),sigma(ii),vrho(ii),vsigma(ii))
   end if
 end do
#else
 np=int(mmax,c_size_t)
 if(fn%family==1 .and. fn%has_exc) then
   call xc_f03_lda_exc_vxc(fn%conf,np,nn,zk,vrho)
 else if(fn%family==1) then
   call xc_f03_lda_vxc(fn%conf,np,nn,vrho)
 else if(fn%has_exc) then
   call xc_f03_gga_exc_vxc(fn%conf,np,nn,sigma,zk,vrho,vsigma)
 else
   call xc_f03_gga_vxc(fn%conf,np,nn,sigma,vrho,vsigma)
 end if
#endif
 return
 end subroutine lxc_eval

 subroutine lxc_fxc(fn,mmax,nn,sigma,v2rho2,v2rhosigma,v2sigma2)
!! the second derivatives of the energy density at mmax points, for a
!! functional with has_fxc; v2rhosigma and v2sigma2 are for a GGA only
!!
!! @param fn  the functional
!! @param mmax  size of the log radial mesh
!! @param nn  electron density (bohr^-3)
!! @param sigma  squared gradient of nn
!! @param v2rho2  second derivative of the energy density by nn
!! @param v2rhosigma  its derivative by nn and sigma
!! @param v2sigma2  its second derivative by sigma

#ifdef ONCVPSP_LIBXC_FXC
#if XC_MAJOR_VERSION<=4
 use xc_f90_lib_m, only: xc_f90_lda_fxc, xc_f90_gga_fxc
#else
 use, intrinsic :: iso_c_binding, only: c_size_t
 use xc_f03_lib_m, only: xc_f03_lda_fxc, xc_f03_gga_fxc
#endif
#endif
 use constants_m, only: dp

!Input variables
 type(lxc_t), intent(in) :: fn
 integer, intent(in) :: mmax
 real(dp), intent(in) :: nn(mmax),sigma(mmax)

!Output variables
 real(dp), intent(out) :: v2rho2(mmax),v2rhosigma(mmax),v2sigma2(mmax)

!Local variables
#ifdef ONCVPSP_LIBXC_FXC
#if XC_MAJOR_VERSION<=4
 integer :: ii
#else
 integer(c_size_t) :: np
#endif
#endif

 v2rho2(:)=0.0d0; v2rhosigma(:)=0.0d0; v2sigma2(:)=0.0d0
 if(.not. fn%has_fxc) return
#ifdef ONCVPSP_LIBXC_FXC
#if XC_MAJOR_VERSION<=4
 do ii=1,mmax
   if(fn%family==1) then
     call xc_f90_lda_fxc(fn%conf,1,nn(ii),v2rho2(ii))
   else
     call xc_f90_gga_fxc(fn%conf,1,nn(ii),sigma(ii),v2rho2(ii), &
&                        v2rhosigma(ii),v2sigma2(ii))
   end if
 end do
#else
 np=int(mmax,c_size_t)
 if(fn%family==1) then
   call xc_f03_lda_fxc(fn%conf,np,nn,v2rho2)
 else
   call xc_f03_gga_fxc(fn%conf,np,nn,sigma,v2rho2,v2rhosigma,v2sigma2)
 end if
#endif
#endif
 return
 end subroutine lxc_fxc

 subroutine lxc_info(fn,kind,name,family,nref,refs)
!! the description of the functional: its kind (Exchange, Correlation or
!! Exchange-correlation, blank for any other), name, family (LDA, GGA)
!! and its first nref references (at most size(refs))
!!
!! @param fn  the functional
!! @param kind  Exchange, Correlation or Exchange-correlation (output)
!! @param name  its name (output)
!! @param family  LDA or GGA (output)
!! @param nref  number of references (output)
!! @param refs  the references (output)

#if XC_MAJOR_VERSION<=4
 use xc_f90_lib_m, only: xc_f90_info_kind, xc_f90_info_name, &
&                        xc_f90_info_refs, XC_EXCHANGE, XC_CORRELATION, &
&                        XC_EXCHANGE_CORRELATION
#else
 use xc_f03_lib_m, only: xc_f03_func_info_get_kind, &
&                        xc_f03_func_info_get_name, &
&                        xc_f03_func_info_get_references, &
&                        xc_f03_func_reference_get_ref, &
&                        xc_f03_func_reference_t, XC_EXCHANGE, &
&                        XC_CORRELATION, XC_EXCHANGE_CORRELATION
#endif

!Input variables
 type(lxc_t), intent(in) :: fn

!Output variables
 character(len=32), intent(out) :: kind,family
 character(len=1024), intent(out) :: name
 character(len=1024), intent(out) :: refs(:)
 integer, intent(out) :: nref

!Local variables
 integer :: ikind,ii
 character(len=1024) :: text
#if XC_MAJOR_VERSION>4
 type(xc_f03_func_reference_t) :: xc_ref
#endif

#if XC_MAJOR_VERSION<=4
 ikind=xc_f90_info_kind(fn%info)
 call xc_f90_info_name(fn%info,text)
#else
 ikind=xc_f03_func_info_get_kind(fn%info)
 text=xc_f03_func_info_get_name(fn%info)
#endif
 name=text
 if(ikind==XC_EXCHANGE) then
   kind='Exchange'
 else if(ikind==XC_CORRELATION) then
   kind='Correlation'
 else if(ikind==XC_EXCHANGE_CORRELATION) then
   kind='Exchange-correlation'
 else
   kind=''
 end if
 if(fn%family==1) then
   family='LDA'
 else
   family='GGA'
 end if

! each call gives one reference and the number of the next, -1 after
! the last
 nref=0
 ii=0
#if XC_MAJOR_VERSION<=4
 call xc_f90_info_refs(fn%info,ii,text)
#else
 xc_ref=xc_f03_func_info_get_references(fn%info,ii)
 text=xc_f03_func_reference_get_ref(xc_ref)
#endif
 do while(ii>=0 .and. nref<size(refs))
   nref=nref+1
   refs(nref)=text
#if XC_MAJOR_VERSION<=4
   call xc_f90_info_refs(fn%info,ii,text)
#else
   xc_ref=xc_f03_func_info_get_references(fn%info,ii)
   text=xc_f03_func_reference_get_ref(xc_ref)
#endif
 end do
 return
 end subroutine lxc_info

 end module libxc_api_m
