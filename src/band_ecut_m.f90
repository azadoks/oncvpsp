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
!> The cutoff of a channel at the energies of the valence bands. The
!> convergence tables are those of the pseudo-waves, at the projectors'
!> energies: a channel with no occupied state (Si d) has only well states,
!> at energies the user chose and with a normalization that is a convention.
!> In a solid, the charge of such a channel comes from hybridization, at the
!> energies of the occupied bands. Here the channel's pseudo scattering
!> state at an energy e (the nonlocal equation, outward, lschvkbs) is
!> normalized to 1 inside rc, continued beyond rc as the windowed solution
!> of free_tail_m, and its residual kinetic energy E(q) above each q is
!> that of eresid on one linear mesh from r = 0: the error per electron
!> inside rc of a cutoff q**2/2, at that energy.
!> Beyond rc the state is in a solid's interstitial, where band states lie
!> above the potential and oscillate. The isolated atom's potential tends
!> to 0, so at e < 0 its solution would grow, and under free_tail_m's
!> window that growth swamps E(q). The tail's potential is therefore the
!> local potential going over to a flat interstitial level vint below the
!> window, by vloc + s (vint - vloc), s rising from 0 at rc as an erf step
!> of width w = 8/qcut centred 3 w beyond rc (its spectrum, exp(-(q w)**2
!> /4), below exp(-16) above qcut); with the energy and potential shifted
!> by vint the free tail is that of a positive energy e - vint, with no
!> growth, and its content lies at k = sqrt(2 (e - vint)), far below qcut.
!> The energies are the valence window of the reference configuration:
!> from its highest occupied valence eigenvalue down band_width, or to its
!> lowest valence eigenvalue if that is higher (the semicore states, flat
!> bands, are left out) (main).

 module band_ecut_m
 use constants_m, only: dp, pi
 use pieces_m, only: piecewise_t, pw_reach
 use logmesh_m, only: lm_integral, lm_fh_reach
 use free_tail_m, only: free_tail_whole, free_tail_k
 use ugrid_m, only: ug_cross_down
 use, intrinsic :: iso_fortran_env, only: qp => real128
 implicit none
 private
 public :: band_ecut

!> the width of the valence window (Ha) below the highest occupied valence
!> eigenvalue, about a valence band's (Si: 12 eV)
 real(dp), parameter, public :: band_width=0.5d0
!> the interstitial level below the lowest energy of the window (Ha)
 real(dp), parameter, public :: band_vint=0.2d0
!> the error levels of the cutoffs (Ha per electron), as the tables'
 integer, parameter, public :: nband_lev=4
 real(dp), parameter, public :: band_lev(nband_lev)= &
&  [1.0d-2,1.0d-3,1.0d-4,1.0d-5]

 contains

 subroutine band_ecut(ll,nvkb,ee,vint,rr,vloc,vlocpw,vkbpw,evkb,mmax,rc, &
&                     irc,qcut,qmax,eqcut,ecut)
!! the residual kinetic energy at qcut, and the cutoffs at the error levels
!! band_lev, of the pseudo scattering state of l at the energy ee,
!! normalized to 1 inside rc (above)
!!
!! @param ll  angular momentum
!! @param nvkb  number of projectors of l (0 for the local channel)
!! @param ee  energy (Ha)
!! @param vint  the interstitial level of the tail's potential (Ha)
!! @param rr  log radial mesh
!! @param vloc  the screened local pseudopotential
!! @param vlocpw  its pieces
!! @param vkbpw  the projectors of l as pieces
!! @param evkb  their coefficients
!! @param mmax  size of the log mesh
!! @param rc  core radius of l
!! @param irc  last mesh point at or inside rc
!! @param qcut  cutoff wave vector of l
!! @param qmax  q taken as infinity in the optimization of l
!! @param eqcut  E_resid above qcut (Ha per electron inside rc); -1 where
!!               the free tail's wave vector is not below qcut, which
!!               leaves it undefined (free_tail_width)
!! @param ecut  the cutoff (Ha) at each level, 0 for a level not reached

 external :: eresid, lschvkbs

!Input variables
 integer, intent(in) :: ll,nvkb,mmax,irc
 real(dp), intent(in) :: ee,vint,rc,qcut,qmax
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)

!Output variables
 real(dp), intent(out) :: eqcut,ecut(nband_lev)

!Local variables
 real(dp), allocatable :: uu(:),up(:),ulin(:),e0(:),ed(:,:),em(:,:,:)
 real(dp), allocatable :: vt(:)
 real(dp) :: hh,dq,q0,sn,qerr,qout(1),eq(1),edum(0,1),emdum(0,0,1)
 real(dp) :: qroot(0),ww,ss,e3
 real(qp) :: ps0(0),psn(0,0)
 integer :: mout,nq,nr,ii
 logical :: found

!undefined where the tail oscillates as fast as the cutoff
 if(free_tail_k(ee-vint)>=qcut) then
   eqcut=-1.0d0
   ecut(:)=0.0d0
   return
 end if

!the state, past the reach of the interpolants at rc
 mout=min(mmax,irc+pw_reach+16+lm_fh_reach()+2)
 allocate(uu(mmax),up(mmax))
 uu(:)=0.0d0
 up(:)=0.0d0
 call lschvkbs(ll,nvkb,ee,rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,mout)

!normalized to 1 inside rc
 call lm_integral(rr,uu**2,mmax,mout,0.0d0,rc,sn)
 uu(:)=uu(:)/sqrt(sn)
 up(:)=up(:)/sqrt(sn)

!the tail's potential, going over to the interstitial level, less that
!level (above)
 allocate(vt(mmax))
 ww=8.0d0/qcut
 e3=erfc(3.0d0)
 do ii=1,mmax
   ss=0.0d0
   if(rr(ii)>rc) then
     ss=(erfc(-(rr(ii)-rc-3.0d0*ww)/ww)-e3)/(2.0d0-e3)
   end if
   vt(ii)=vloc(ii)+ss*(vint-vloc(ii))-vint
 end do

!on the linear mesh from 0 of eresid's spacing, with its windowed tail
!beyond rc; the q grid of run_optimize
 hh=2.0d0*pi/(75.0d0*qmax)
 call free_tail_whole(0,ll,0,ee-vint,rr,vt,uu,up,mmax,rc,qcut,hh,nr,ulin)
 dq=2.0d0*pi/(150.0d0*rc)
 nq=1+int(qmax/dq)
 q0=qmax-dble(nq-1)*dq
 qout(1)=qcut
 allocate(e0(nq),ed(nq,0),em(nq,0,0))
 call eresid(ll,0.0d0,0,0,nq,dq,1,qout,qmax,qroot,hh,nr,ulin,ps0,psn, &
&            e0,ed,em,eq,edum,emdum)
 eqcut=eq(1)
 do ii=1,nband_lev
   call ug_cross_down(e0,nq,dq,q0,band_lev(ii),qerr,found)
   ecut(ii)=0.0d0
   if(found) ecut(ii)=0.5d0*qerr**2
 end do
 deallocate(uu,up,ulin,e0,ed,em,vt)
 return
 end subroutine band_ecut

 end module band_ecut_m
