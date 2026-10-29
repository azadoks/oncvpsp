!
! Copyright (c) 2026 by Austin Zadoks
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
!> The last line of oncvpsp_confine.x, confine_summary, from the sphere
!> The confined atom of oncvpsp_confine.x and oncv_confine, from the
!> all-electron atom and the pseudopotential data file's contents: the
!> Wigner-Seitz sphere (confine_ws_m), the walls (confine_wall_m) and
!> the summary (confine_summary_m), its lines to the text output.
 module confine_atom_m
 use constants_m, only: dp
 implicit none
 private
 public :: confine_atom, rsdef, descan, kdef

!> the default radii, in units of the largest core radius
 real(dp), parameter :: rsdef(3)=[1.25d0,1.5d0,2.0d0]
!> the default step of the scan (Ha): it only brackets the roots of
!> monotonic phases, which no step can miss
 real(dp), parameter :: descan=0.1d0
!> the default strengths of the wall (Ha/bohr**2)
 real(dp), parameter :: kdef(3)=[0.01d0,0.03d0,0.1d0]

 contains

 subroutine confine_atom(pf,irel,c_light,zz,mmax,rr,vfull,nc,nv,iexc,na, &
&                        la,fa,ea,chan,sphere,nq,rq,ee,ecutmax,table,nk, &
&                        kwall,r0,sm,iout,ifail)

!! the confined atom: the libxc functional's description, the sphere at
!! the radii rq (by default rsdef times the largest core radius), each
!! beyond every core radius, in the window ee (by default the input's
!! window of the log derivatives, with the step descan), then the atom in
!! the walls kwall(1:nk) beyond r0 (by default the largest core radius),
!! and the summary line
!!
!! @param pf  the pseudopotential data file's contents (ps_file_m)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional
!! @param na  principal quantum numbers
!! @param la  angular momenta
!! @param fa  occupancies
!! @param ea  eigenvalue of each state and kappa (Ha)
!! @param chan  whether to treat each l, 0..3
!! @param sphere  whether to make the sphere
!! @param nq  number of radii, 0 for the default
!! @param rq  the radii (bohr), nq of them
!! @param ee  the window and step of the scan (Ha), ee(3) 0 or less for the
!!            default
!! @param ecutmax  the largest cutoff of the basis (Ha)
!! @param table  whether to write the error at every basis size
!! @param nk  number of wall strengths, 0 for none
!! @param kwall  the wall strengths (Ha/bohr**2), nk of them
!! @param r0  the onset of the wall (bohr), 0 or less for the default
!! @param sm  the summary (confine_summary_m)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: mxstate, mxkap
 use ps_file_m, only: ps_file_t
 use confine_ws_m, only: confine_sphere
 use confine_wall_m, only: confine_walls
 use confine_summary_m, only: confine_summary_t, write_confine_summary
 use exc_libxc_m, only: exc_libxc_info
 use adams_m, only: adams_t
 use failure_m, only: fail_input, failed, msg_confine_radius

!Input variables
 type(ps_file_t), intent(in) :: pf
 integer, intent(in) :: irel,mmax,nc,nv,iexc,nq,nk
 integer, intent(in) :: na(mxstate),la(mxstate)
 real(dp), intent(in) :: c_light,zz,rr(mmax),vfull(mmax)
 real(dp), intent(in) :: fa(mxstate),ea(mxstate,mxkap)
 real(dp), intent(in) :: rq(:),ee(3),ecutmax,kwall(:),r0
 logical, intent(in) :: chan(0:3),sphere,table
 integer, intent(in) :: iout

!Output variables
 type(confine_summary_t), intent(out) :: sm
 integer, intent(out) :: ifail

!Local variables
!> the Adams-Moulton steps, of the default order
 type(adams_t) :: adams
 real(dp), allocatable :: rs(:),kw(:)
 real(dp) :: rcmax,emin,emax,de,rw
 integer :: iq,lc,ns

 ifail=0
 if(iexc<0) then
   call exc_libxc_info(iexc,iout,ifail)
   if(ifail/=0) return
 end if

!the radii: beyond every core radius, where the two atoms' potentials are
!the same and the phases comparable
 lc=pf%pot%lloc+1
 rcmax=max(maxval(pf%rc(1:pf%pot%lmax+1)),pf%rc(lc))
 if(nq>0) then
   rs=rq(1:nq)
 else
   rs=rsdef(:)*rcmax
 end if
 ns=size(rs)
 do iq=1,ns
   if(rs(iq)<=rcmax) then
     write(iout,'(/a,f12.6,a,f12.6,a)') 'confine_atom: ERROR the radius', &
&         rs(iq),' is not beyond the largest core radius',rcmax,' bohr'
     ifail=failed(fail_input,msg_confine_radius)
     return
   end if
 end do

!the window of the edges, that of the log derivatives
 if(ee(3)>0.0d0) then
   emin=ee(1)
   emax=ee(2)
   de=ee(3)
 else
   emin=pf%epsh1
   emax=pf%epsh2
   de=descan
 end if

 if(sphere) then
   call confine_sphere(pf,irel,c_light,zz,mmax,adams,rr,vfull,nc,na,la,chan, &
&                      ns,rs,emin,emax,de,ecutmax,table,sm,iout,ifail)
   if(ifail/=0) return
 end if

!the atom in walls: k=0 the free atom
 if(nk>0) then
   rw=r0
   if(rw<=0.0d0) rw=rcmax
   allocate(kw(0:nk))
   kw(1:nk)=kwall(1:nk)
   call confine_walls(pf,irel,c_light,zz,mmax,adams,rr,vfull,nc,nv,iexc,na,la, &
&                     fa,ea,nk,kw,rw,sm,iout,ifail)
   if(ifail/=0) return
 end if

 call write_confine_summary(sm,iout)
 return
 end subroutine confine_atom

 end module confine_atom_m
