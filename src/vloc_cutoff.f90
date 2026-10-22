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
! The Fourier convergence of the local potential, for the DIAGNOSTICS
! block: the cutoffs at which a plane-wave code's local energy of the
! pseudo valence density is converged to the levels of the projectors'
! tables

 subroutine vloc_cutoff(rr,mmax,rho,vloc,zion,cvl,echk)

!! The local energy E_loc = int n V_loc d^3r of the reference
!! configuration's pseudo valence density n is, in q space,
!! (1/2 pi^2) int q^2 n(q) V(q) dq. A plane-wave code sums it over the
!! density's sphere, |G| < Q = 2 sqrt(2 Ecut) for the wave cutoff Ecut,
!! so its truncation per electron (N = int n) is bounded by
!!
!!   T(Q) = (1/(2 pi^2 N)) [int_Q q^2 n(q)^2 dq int_Q q^2 V(q)^2 dq]^(1/2)
!!
!! (Cauchy-Schwarz), which falls monotonically; its crossings at the
!! levels of the tables (1e-2 Ha down by half decades to 1e-5) are
!! returned as cutoffs Ecut = Q^2/8 (a level not reached below qend has
!! none). The signed truncation oscillates in Q and has no single
!! crossing; the tail of |n(q) V(q)|, whose cutoffs are up to 35 % below
!! the bound's (14_Si, 07_N_opt, 60_Nd_GHOST, 73_Ta), has a kink at each
!! zero of n V, which move with the inputs (its cutoffs were smooth to
!! 5e-7 only, and jumped at 1e-2 in 07_N_opt's rcloc). A spike of V_loc
!! near the origin, where n is finite for s states, shows as a slow tail
!! of V(q).
!!
!! The transforms are taken on a uniform r grid (16 points per period of
!! qend) to rr(mmax), the functions interpolated from the log mesh, with
!! V(q) = 4 pi int (V_loc + zion/r) j0(qr) r^2 dr - 4 pi zion/q^2 (V_loc is
!! -zion/r far out; V(q) falls faster than 1/q^2, as V_loc is finite at
!! the origin). qend = 60 bohr^-1 (Ecut 450 Ha): the levels below 1e-5
!! would be the quadrature's (the reference outputs' 1e-7 levels moved by
!! a factor of 3 with its points per period).
!!
!! @param rr  log radial mesh
!! @param mmax  size of the log radial mesh
!! @param rho  the pseudo valence density, 4 pi n (int rho r^2 dr = N)
!! @param vloc  the local potential (Ha)
!! @param zion  the ionic charge
!! @param cvl  (Ecut, err) per level, Ecut 0 for a level not reached
!! @param echk  E_loc per electron from the q grid less that from r space
!!              (Ha): the transforms' error and the part above qend

 use constants_m, only: dp, pi, nerrlv
 use logmesh_m, only: lm_interp, lm_integral
 use ugrid_m, only: ug_weights, ug_cumint_down, ug_cross_down
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rho(mmax),vloc(mmax),zion

!Output variables
 real(dp), intent(out) :: cvl(2,nerrlv),echk

!Local variables
 real(dp), parameter :: qend=60.0d0,dq=0.02d0,npp=16.0d0
 real(dp), allocatable :: xx(:),wx(:),fn(:),fv(:),gg(:),tt(:),wq(:)
 real(dp), allocatable :: work(:),an(:),av(:)
 real(dp) :: hr,qq,sn,sv,sq,sp,sm,cs,nel,eloc_r,eloc_q,err,qerr,vq
 integer :: ii,kk,nr,nq
 logical :: found

 allocate(work(mmax))
 work(:)=rho(:)*rr(:)**2
 call lm_integral(rr,work,mmax,2,0.0d0,rr(mmax),nel)
 work(:)=rho(:)*vloc(:)*rr(:)**2
 call lm_integral(rr,work,mmax,2,0.0d0,rr(mmax),eloc_r)
 deallocate(work)

!the r grid from 0, the integrands r n(r) and r (V_loc + zion/r): below
!the mesh's first point, n and V_loc at it
 hr=2.0d0*pi/(npp*qend)
 nr=int(rr(mmax)/hr)+1
 allocate(xx(nr),wx(nr),fn(nr),fv(nr))
 call ug_weights(nr,hr,wx)
 do ii=1,nr
   xx(ii)=hr*dble(ii-1)
   if(xx(ii)<rr(1)) then
     sn=rho(1)
     sv=vloc(1)
   else
     call lm_interp(rr,rho,mmax,xx(ii),sn)
     call lm_interp(rr,vloc,mmax,xx(ii),sv)
   end if
   fn(ii)=wx(ii)*xx(ii)*sn
   fv(ii)=wx(ii)*(xx(ii)*sv+zion)
 end do

!the integrand q^2 n(q) V(q) on the q grid from 0, where it is
!-4 pi zion N (the Coulomb part's)
 nq=int(qend/dq)+1
 allocate(gg(nq),tt(nq),wq(nq),an(nq),av(nq))
!the sums of fn and fv times sin(q r) at every q: sin((k+1)dq r) from
!the recurrence 2 cos(dq r) sin(k dq r) - sin((k-1) dq r), whose error
!grows by ~eps a step (3e-13 at qend), in place of nq*nr sines
 an(:)=0.0d0
 av(:)=0.0d0
 do ii=2,nr
   cs=2.0d0*cos(dq*xx(ii))
   sp=0.0d0
   sq=sin(dq*xx(ii))
   do kk=2,nq
     an(kk)=an(kk)+fn(ii)*sq
     av(kk)=av(kk)+fv(ii)*sq
     sm=cs*sq-sp
     sp=sq
     sq=sm
   end do
 end do
 gg(1)=-4.0d0*pi*zion*nel
 do kk=2,nq
   qq=dq*dble(kk-1)
   sn=an(kk)
   vq=4.0d0*pi*(av(kk)-zion/qq)
   gg(kk)=sn*vq
!q^2 n(q)^2 and q^2 V(q)^2
   an(kk)=sn**2
   av(kk)=vq**2
 end do

!E_loc from the q grid
 call ug_weights(nq,dq,wq)
 eloc_q=sum(wq(:)*gg(:))/(2.0d0*pi**2)
 echk=(eloc_q-eloc_r)/nel

!the bound from dq (int q^2 V^2 diverges at 0)
 call ug_cumint_down(an(2:nq),nq-1,dq,tt(2:nq))
 call ug_cumint_down(av(2:nq),nq-1,dq,gg(2:nq))
 tt(2:nq)=sqrt(tt(2:nq)*gg(2:nq))/(2.0d0*pi**2*nel)

 err=0.01d0
 cvl(:,:)=0.0d0
 do ii=1,nerrlv
   call ug_cross_down(tt(2:nq),nq-1,dq,dq,err,qerr,found)
   if(found) then
     cvl(1,ii)=0.125d0*qerr**2
     cvl(2,ii)=err
   end if
   err=sqrt(0.1d0)*err
 end do

 deallocate(xx,wx,fn,fv,gg,tt,wq,an,av)
 return
 end subroutine vloc_cutoff
