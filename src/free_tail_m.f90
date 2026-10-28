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
!> The all-electron tail beyond rc that eresid transforms, on its linear
!> mesh rc + (i-1) hh (doc/smoothness_plan.md, phase 8). mesh_tail takes
!> the wave function on the log mesh to its last nonzero point, as for a
!> bound state. For a well state, the tail is the atom's free solution at
!> its energy, with no well, under a window that takes it to zero far from
!> rc,
!>
!>     t(r) = u_free(r) w(r),  w(r) = erfc((r - R0)/L)/2,  R0 = rc + 6 L,
!>
!> to R0 + x_e L, where it is below 1e-17 of its value at rc (window_end).
!> u_free solves the radial equation in the all-electron
!> potential at the state's energy outward from the state's u and u' at
!> rc, so it is the state itself up to the well's onset (>= rc), and the
!> tail no longer depends on the well: its onset, width and node count.
!> free_tail gives it on eresid's mesh, and free_tail_log on the log mesh
!> beyond the last node that main solves the state to (scatterstate, past
!> the reach of everything taken at or near rc), so that the state is the
!> free solution everywhere. The
!> window's spectrum falls as a Gaussian, and L (free_tail_width) puts its
!> share of E_resid above qcut below about 1e-13, for growing tails (ep <
!> 0) as well, smoothly in qcut and ep, while kappa L <= 12
!> (free_tail_defined).

 module free_tail_m
 use constants_m, only: dp, pi
 implicit none
 private
 public :: mesh_tail, free_tail, free_tail_log, free_tail_width
 public :: tail_kinetic, free_tail_k, free_tail_defined

!> the leakage bound: L solves (dq L)**2/4 - 9 kappa L = wbound
 real(dp), parameter :: wbound=30.0d0
!> the largest kappa L for which that bound holds (free_tail_width)
 real(dp), parameter :: klmax=12.0d0
!> the width in energy of the smooth split of ep into k and kappa (Ha)
 real(dp), parameter :: dsplit=1.0d-2
!> the window's centre from rc, in L
 real(dp), parameter :: wcen=6.0d0
!> the log of the tail's value at its end relative to that at rc
 real(dp), parameter :: wend=39.0d0
!> Runge-Kutta steps per step of eresid's linear mesh
 integer, parameter :: nsub=2
!> the longest Runge-Kutta step between the nodes of the log mesh (bohr)
 real(dp), parameter :: hlog=2.0d-3

 contains

 subroutine mesh_tail(rr,uu,mmax,rc,hh,nr,ulin,iout,ifail)

!! uu on the linear mesh rc + (i-1) hh to its last nonzero point (the
!! solvers zero it beyond their inward start; a node there is dropped or
!! added at a value at machine precision), interpolated by lm_interp
!!
!! @param rr  log radial mesh
!! @param uu  all-electron wave function (*rr)
!! @param mmax  size of the log radial mesh
!! @param rc  core radius (bohr)
!! @param hh  spacing of the linear mesh
!! @param nr  number of its nodes
!! @param ulin  uu on them
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use logmesh_m, only: lm_interp, lm_ncell

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),uu(mmax),rc,hh
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: nr
 real(dp), allocatable, intent(out) :: ulin(:)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: rend
 integer :: ii,nlast

 ifail=0
 nlast=0
 do ii=mmax,1,-1
   if(abs(uu(ii))>0.0d0) then
     nlast=ii
     exit
   end if
 end do
 rend=0.0d0
 if(nlast>0) rend=rr(nlast)
 nr=0
 if(rend>rc) nr=int((rend-rc)/hh)+1
 if(nr<2*lm_ncell) then
   write(iout,'(/2a,f10.6,a,f10.6)') &
&        'mesh_tail: ERROR wave function vanishes', &
&       ' too close to rc: rc=',rc,', last nonzero point',rend
   ifail=1
   return
 end if
 allocate(ulin(nr))
 do ii=1,nr
   call lm_interp(rr,uu,mmax,rc+dble(ii-1)*hh,ulin(ii))
 end do
 return
 end subroutine mesh_tail

 subroutine free_tail_width(qcut,ee,wl,iout,ifail)

!! the window's width L, the positive root of (dq L)**2/4 - 9 kappa L =
!! 30 with dq = qcut - k: the window's share of E_resid above qcut is
!! about exp(-(dq L)**2/4), and a tail growing as exp(kappa r) gains at
!! most exp(9 kappa L) over the window. k and kappa are sqrt(2 ep) and
!! sqrt(-2 ep) split smoothly at ep = 0, k = sqrt(sqrt(ep**2 + d**2) + ep)
!! and kappa = sqrt(sqrt(ep**2 + d**2) - ep), d = 0.01 Ha, so that L has
!! no kink there. qcut <= k is an input error: the state then oscillates
!! faster than the cutoff. So is kappa L > 12: the windowed tail
!! exp(kappa (r - rc)) erfc((r - R0)/L)/2 peaks at exp(6 kappa L +
!! (kappa L)**2/4), which exceeds the exp(9 kappa L) of the bound beyond
!! kappa L = 12 (an I d well state at -0.50 Ha with qcut 1.43: kappa L
!! 20.7, the tail 3e99 times its value at rc, E_resid 6e171 Ha)
!!
!! @param qcut  the cutoff wave vector (bohr^-1)
!! @param ee  the state's energy ep (Ha)
!! @param wl  the width L (bohr)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 real(dp), intent(in) :: qcut,ee
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: wl
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: ss,kk,kp,dq

 ifail=0
 ss=sqrt(ee**2+dsplit**2)
 kk=free_tail_k(ee)
 kp=sqrt(ss-ee)
 dq=qcut-kk
 if(dq<=0.0d0) then
   write(iout,'(/a,f10.6,a,f10.6/a)') &
&       'free_tail_width: ERROR qcut=',qcut,' is not above the wave vector', &
&       kk,'of the well state: raise qcut or lower its energy'
   ifail=1
   return
 end if
 wl=(9.0d0*kp+sqrt(81.0d0*kp**2+wbound*dq**2))/(0.5d0*dq**2)
 if(kp*wl>klmax) then
   write(iout,'(/a,f10.6,a,f10.6/a,f8.2,a/a)') &
&       'free_tail_width: ERROR the free tail of the well state at',ee, &
&       ' Ha with qcut=',qcut,'grows too fast for its window (kappa L =', &
&       kp*wl,' > 12): raise qcut or the state''s energy'
   ifail=1
   return
 end if
 return
 end subroutine free_tail_width

 logical function free_tail_defined(qcut,ee)
!! whether the free tail of a state at ee is defined at qcut: qcut above
!! its wave vector k, and kappa L <= 12 (free_tail_width)
!!
!! @param qcut  the cutoff wave vector (bohr^-1)
!! @param ee  the state's energy (Ha)

!Input variables
 real(dp), intent(in) :: qcut,ee

!Local variables
 real(dp) :: kp,dq,wl

 free_tail_defined=.false.
 dq=qcut-free_tail_k(ee)
 if(dq<=0.0d0) return
 kp=sqrt(sqrt(ee**2+dsplit**2)-ee)
 wl=(9.0d0*kp+sqrt(81.0d0*kp**2+wbound*dq**2))/(0.5d0*dq**2)
 free_tail_defined=kp*wl<=klmax
 end function free_tail_defined

 real(dp) function free_tail_k(ee)
!! the wave vector k of a state at ee, sqrt(2 ee) split smoothly at ee = 0
!! as in free_tail_width; the free tail is defined for qcut > k
!!
!! @param ee  the state's energy (Ha)

!Input variables
 real(dp), intent(in) :: ee

 free_tail_k=sqrt(sqrt(ee**2+dsplit**2)+ee)
 end function free_tail_k

 subroutine free_tail(irel,c_light,ll,kap,ee,rr,vv,uu,up,mmax,rc,qcut,hh,nr, &
&                     ulin,dlin,wscale,iout,ifail)

!! the windowed free tail of a well state on the linear mesh rc + (i-1) hh
!! to R0 + x_e L (above), nsub Runge-Kutta steps per mesh step
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2
!!              fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param ee  the state's energy (Ha)
!! @param rr  log radial mesh
!! @param vv  all-electron potential
!! @param uu  the well state (*rr; uu the large component for irel=2)
!! @param up  its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param rc  core radius (bohr)
!! @param qcut  the cutoff wave vector (bohr^-1)
!! @param hh  spacing of the linear mesh
!! @param nr  number of its nodes
!! @param ulin  the tail on them
!! @param dlin  its derivative
!! @param wscale  optional, a factor on L, for the tests (default 1)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ee,rc,qcut,hh
 real(dp), intent(in) :: rr(mmax),vv(mmax),uu(mmax),up(mmax)
 real(dp), intent(in), optional :: wscale
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: nr
 real(dp), allocatable, intent(out) :: ulin(:),dlin(:)
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: rlin(:)
 real(dp) :: wl
 integer :: ii

 ifail=0
 call free_tail_width(qcut,ee,wl,iout,ifail)
 if(ifail/=0) return
 if(present(wscale)) wl=wscale*wl
 nr=int((wcen+window_end(ee,wl))*wl/hh)+1
 allocate(rlin(nr),ulin(nr),dlin(nr))
 do ii=1,nr
   rlin(ii)=rc+dble(ii-1)*hh
 end do
 call free_solve(irel,c_light,ll,kap,ee,rr,vv,uu,up,mmax,rc,rc,wl,nr,rlin, &
&                hh/dble(nsub),ulin,dlin)
 deallocate(rlin)
 return
 end subroutine free_tail

 subroutine free_tail_log(irel,c_light,ll,kap,ee,rr,vv,uu,up,mmax,rc,qcut, &
&      istart,iout,ifail)

!! the state's values beyond node istart on the log mesh, where its solver
!! stopped, the windowed free tail (above) continued from that node, and
!! zero beyond R0 + x_e L; Runge-Kutta steps of at most hlog between the
!! nodes
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2
!!              fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param ee  the state's energy (Ha)
!! @param rr  log radial mesh
!! @param vv  all-electron potential
!! @param uu  the well state (*rr; uu the large component for irel=2)
!! @param up  its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param rc  core radius (bohr)
!! @param qcut  the cutoff wave vector (bohr^-1)
!! @param istart  the last node of the solver's values, beyond rc
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use logmesh_m, only: lm_step

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax,istart
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ee,rc,qcut
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: uu(mmax),up(mmax)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: ut(:),dt(:)
 real(dp) :: wl,al
 integer :: ii,i0,i1

 ifail=0
 call free_tail_width(qcut,ee,wl,iout,ifail)
 if(ifail/=0) return
 i0=istart+1
 i1=count(rr(:)<=rc+(wcen+window_end(ee,wl))*wl)
 if(i1<i0) return
 allocate(ut(i0:i1),dt(i0:i1))
 call free_solve(irel,c_light,ll,kap,ee,rr,vv,uu,up,mmax,rc,rr(istart),wl, &
&                i1-i0+1,rr(i0:i1),hlog,ut,dt)
 al=lm_step(rr,mmax)
 do ii=i0,i1
   uu(ii)=ut(ii)
   up(ii)=al*rr(ii)*dt(ii)
 end do
 uu(i1+1:mmax)=0.0d0
 up(i1+1:mmax)=0.0d0
 deallocate(ut,dt)
 return
 end subroutine free_tail_log

 real(dp) function window_end(ee,wl)

!! x_e, the end of the window R0 + x_e L in L from its centre: the tail
!! there, at most exp(kappa (r - rc)) erfc(x) of its value at rc, is below
!! exp(-39) = 1e-17 of it, x_e the root of x**2 - kappa L x - 6 kappa L =
!! 39 (erfc(x) < exp(-x**2)); 6.24 for kappa = 0. kappa as in
!! free_tail_width
!!
!! @param ee  the state's energy ep (Ha)
!! @param wl  the window's width L (bohr)

!Input variables
 real(dp), intent(in) :: ee,wl

!Local variables
 real(dp) :: kl

 kl=sqrt(sqrt(ee**2+dsplit**2)-ee)*wl
 window_end=0.5d0*(kl+sqrt(kl**2+4.0d0*(wend+wcen*kl)))
 end function window_end

 subroutine free_solve(irel,c_light,ll,kap,ee,rr,vv,uu,up,mmax,rc,rs,wl,nout, &
&                      rout,hmax,uout,dout)

!! the windowed free solution and its derivative at the radii rout (rs <=
!! rout(1) < rout(2) < ...), the window that of rc: the radial equation
!! in vv at ee, in the first-order form of ode_rc_der (the flux form for
!! irel <= 1, Dirac for irel = 2), by classical Runge-Kutta in equal steps
!! of at most hmax between consecutive radii, from uu and uu' = up/(al r)
!! at rs (lm_interp); vv interpolated by lm_interp, and taken as
!! vv(mmax) rr(mmax)/r beyond the log mesh

 use logmesh_m, only: lm_interp, lm_step

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax,nout
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ee,rc,rs,wl,hmax
 real(dp), intent(in) :: rr(mmax),vv(mmax),uu(mmax),up(mmax),rout(nout)

!Output variables
 real(dp), intent(out) :: uout(nout),dout(nout)

!Local variables
 real(dp) :: r0,al,u0,du0,rad,hs,vlo,vmid,vhi,xx
 real(dp) :: yy(2),k1(2),k2(2),k3(2),k4(2),dy(2)
 real(dp), allocatable :: dudr(:)
 integer :: ii,jj,ns

 r0=rc+wcen*wl
 al=lm_step(rr,mmax)
 allocate(dudr(mmax))
 dudr(:)=up(:)/(al*rr(:))
 call lm_interp(rr,uu,mmax,rs,u0)
 call lm_interp(rr,dudr,mmax,rs,du0)
 deallocate(dudr)

 rad=rs
 vlo=vat(rad)
 yy(1)=u0
 if(irel<=1) then
   yy(2)=(du0-u0/rs)/xmass(vlo)
 else
   yy(2)=(du0+dble(kap)*u0/rs)*c_light/(2.0d0*c_light**2+ee-vlo)
 end if

 do ii=1,nout
   ns=ceiling((rout(ii)-rad)/hmax)
   hs=0.0d0
   if(ns>0) hs=(rout(ii)-rad)/dble(ns)
   do jj=1,ns
     vmid=vat(rad+0.5d0*hs)
     vhi=vat(rad+hs)
     call deriv(rad,vlo,yy,k1)
     call deriv(rad+0.5d0*hs,vmid,yy+0.5d0*hs*k1,k2)
     call deriv(rad+0.5d0*hs,vmid,yy+0.5d0*hs*k2,k3)
     call deriv(rad+hs,vhi,yy+hs*k3,k4)
     yy(:)=yy(:)+hs*(k1(:)+2.0d0*k2(:)+2.0d0*k3(:)+k4(:))/6.0d0
     rad=rad+hs
     vlo=vhi
   end do
   rad=rout(ii)
   call deriv(rad,vlo,yy,dy)
   xx=(rad-r0)/wl
   uout(ii)=yy(1)*0.5d0*erfc(xx)
   dout(ii)=dy(1)*0.5d0*erfc(xx)-yy(1)*exp(-xx**2)/(wl*sqrt(pi))
 end do
 return

 contains

 real(dp) function vat(rad)
 real(dp), intent(in) :: rad
 if(rad<=rr(mmax)) then
   call lm_interp(rr,vv,mmax,rad,vat)
 else
   vat=vv(mmax)*rr(mmax)/rad
 end if
 end function vat

 real(dp) function xmass(vr)
 real(dp), intent(in) :: vr
 xmass=1.0d0
 if(irel==1) xmass=1.0d0+0.5d0*(ee-vr)/c_light**2
 end function xmass

 subroutine deriv(rad,vr,yv,dy)

!! y' of the first-order system at rad, the potential there vr

 real(dp), intent(in) :: rad,vr,yv(2)
 real(dp), intent(out) :: dy(2)
 real(dp) :: xm

 if(irel<=1) then
! uu'=uu/r + M qq, qq'=(l(l+1)/(M r**2) + 2(V-E)) uu - qq/r
   xm=xmass(vr)
   dy(1)=yv(1)/rad+xm*yv(2)
   dy(2)=(dble(ll*(ll+1))/(xm*rad**2)+2.0d0*(vr-ee))*yv(1)-yv(2)/rad
 else
! g'=-kap g/r + (2c**2 + E - V) f/c, f'=kap f/r + (V - E) g/c
   dy(1)=-dble(kap)*yv(1)/rad+(2.0d0*c_light**2+ee-vr)*yv(2)/c_light
   dy(2)=dble(kap)*yv(2)/rad+(vr-ee)*yv(1)/c_light
 end if
 return
 end subroutine deriv

 end subroutine free_solve

 real(dp) function tail_kinetic(ll,rc,hh,nr,ulin,dlin)

!! the kinetic energy of the tail, int_rc^inf t (-t''/2 + l(l+1) t/(2 r**2))
!! dr = t(rc) t'(rc)/2 + int_rc^inf (t'**2 + l(l+1) t**2/r**2)/2 dr (t
!! vanishes at the end), by the rule of ugrid_m
!!
!! @param ll  angular momentum
!! @param rc  core radius (bohr)
!! @param hh  spacing of the linear mesh rc + (i-1) hh
!! @param nr  number of its nodes
!! @param ulin  the tail on them
!! @param dlin  its derivative

 use ugrid_m, only: ug_weights

!Input variables
 integer, intent(in) :: ll,nr
 real(dp), intent(in) :: rc,hh,ulin(nr),dlin(nr)

!Local variables
 real(dp), allocatable :: wr(:)
 real(dp) :: rad
 integer :: ii

 allocate(wr(nr))
 call ug_weights(nr,hh,wr)
 tail_kinetic=0.5d0*ulin(1)*dlin(1)
 do ii=1,nr
   rad=rc+dble(ii-1)*hh
   tail_kinetic=tail_kinetic+0.5d0*wr(ii) &
&              *(dlin(ii)**2+dble(ll*(ll+1))*(ulin(ii)/rad)**2)
 end do
 deallocate(wr)
 return
 end function tail_kinetic

 end module free_tail_m
