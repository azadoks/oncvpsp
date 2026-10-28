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
!> The pieces of the model core charges at continuous radii
!> (doc/smoothness_plan.md, phase 5): the Teter function and its
!> derivatives with no cutoff (M-5), the crossover radius where the core
!> charge equals fcfact times the pseudo valence charge (M-1), the
!> polynomial of icmod 1 with its constant in closed form (M-3), and the
!> scale of the Teter function of icmod 2 (M-4).
!> The model core charges use the Teter function tapered to 0 between
!> x = 3 and 4 (teter_core_derivs) by a C-infinity switch: the bare
!> function decays only as x**-10, far slower than a core charge, and the
!> ionic tail would then differ from -zion/r out to 20 bohr and more.
!>
!> The Teter function is g = h**2 with h(x) = sin(2 pi x)/(2 pi x (1-4x**2)
!> (1-x**2)). The partial fractions of 1/(2 pi x (1-4x**2)(1-x**2)) have
!> their poles at p = 0, +-1/2, +-1, where sin(2 pi x) = cos(2 pi p)
!> sin(2 pi (x-p)), so h is a sum of five shifted sinc functions,
!>
!>   h(x) = S(u) + 2/3 [S(u-pi) + S(u+pi)] + 1/6 [S(u-2pi) + S(u+2pi)],
!>
!> S(y) = sin(y)/y, u = 2 pi x: entire, with no removable singularity left
!> to treat. S and its derivatives come from their power series for
!> |y| <= 2 and from the closed form beyond, both to a few eps there.
!>
!> The model core reaches the writers as pieces (core_pieces): a model
!> with a joint (icmod 1 at r*, 6 and 7 at r_s) is the fill inside it and
!> the core charge beyond, each a smooth function past the joint, so the
!> interpolants onto the output meshes do not span it (the psp8 file's
!> derivative columns, of which the 3rd and 4th jump there for icmod 6
!> and 7). The fill, and the Teter function of a model with no joint, are
!> closed forms (pw_form: poly_form, teter_form, the budget's fill_form),
!> evaluated at the output radii and at the origin; only the core charge
!> beyond a joint is interpolated.

 module modcore_m
 use constants_m, only: dp, pi, nrhomod
 implicit none
 private
 public :: teter_derivs, teter_core_derivs, core_crossover, poly_core_fit
 public :: teter_scale, core_pieces, core_pieces_smooth
 public :: poly_form, teter_form
 public :: d2exc_baseline, d2exc_report

!> the weights of the shifted sinc functions of h, at the shifts k pi
 real(dp), parameter :: wsinc(-2:2)=[1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0, &
&                                    2.0d0/3.0d0, 1.0d0/6.0d0]

!> the taper of the model's Teter function: 1 up to xtap1
 real(dp), parameter :: xtap1=3.0d0
!> the taper of the model's Teter function: 0 from xtap2
 real(dp), parameter :: xtap2=4.0d0

 contains

 subroutine teter_derivs(xx,nder,gd)

!! the Teter function g(x) = (sin(2 pi x)/(2 pi x (1-4x**2)(1-x**2)))**2
!! and its derivatives at x >= 0, for every x (no cutoff)
!!
!! @param xx  argument, x >= 0
!! @param nder  highest derivative wanted, 0<=nder<=4
!! @param gd  d^k g/dx^k in gd(k), k=0..nder

!Input variables
 real(dp), intent(in) :: xx
 integer, intent(in) :: nder

!Output variables
 real(dp), intent(out) :: gd(0:nder)

!Local variables
 real(dp) :: hd(0:4)
 integer :: kk,jj

 call teter_h(xx,nder,hd)
 do kk=0,nder
   gd(kk)=0.0d0
   do jj=0,kk
     gd(kk)=gd(kk) + binom(kk,jj)*hd(jj)*hd(kk-jj)
   end do
 end do
 end subroutine teter_derivs

 subroutine teter_core_derivs(xx,nder,gd)

!! the Teter function of the model core charges, g(x) s(x), and its
!! derivatives at x >= 0: s the C-infinity switch f(1-t)/(f(1-t)+f(t)),
!! f(t) = exp(-1/t), t = (x-xtap1)/(xtap2-xtap1), 1 up to xtap1 and 0 from
!! xtap2 (as cinf_taper); the derivatives of s by Taylor arithmetic
!!
!! @param xx  argument, x >= 0
!! @param nder  highest derivative wanted, 0<=nder<=4
!! @param gd  d^k (g s)/dx^k in gd(k), k=0..nder

!Input variables
 real(dp), intent(in) :: xx
 integer, intent(in) :: nder

!Output variables
 real(dp), intent(out) :: gd(0:nder)

!Local variables
 real(dp) :: gt(0:4),sw(0:4),tj(0:4),e0(0:4),e1(0:4),den(0:4),fac
 integer :: kk,jj

 if(xx<=xtap1) then
   call teter_derivs(xx,nder,gd)
   return
 end if
 gd(:)=0.0d0
 if(xx>=xtap2) return

! the switch's Taylor coefficients in x, s^(k)/k!
 fac=1.0d0/(xtap2-xtap1)
 tj(:)=0.0d0
 tj(0)=(xx-xtap1)*fac
 tj(1)=fac
 call jet_expinv(tj,e1)
 tj(0)=1.0d0-tj(0)
 tj(1)=-fac
 call jet_expinv(tj,e0)
 call jet_inv(e0+e1,den)
 call jet_mul(e0,den,sw)

! s^(k) = k! sw(k), and Leibniz
 fac=1.0d0
 do kk=1,4
   fac=fac*kk
   sw(kk)=sw(kk)*fac
 end do
 call teter_derivs(xx,nder,gt(0:nder))
 do kk=0,nder
   do jj=0,kk
     gd(kk)=gd(kk) + binom(kk,jj)*gt(kk-jj)*sw(jj)
   end do
 end do
 end subroutine teter_core_derivs

 subroutine jet_expinv(tj,ej)

!! Taylor coefficients of exp(-1/t) from those of t > 0, to order 4; 0 where
!! exp(-1/t) < 1e-30 (t < 1/ln(1e30)), so that the model core is never so
!! small that its value is written with no room for the E of a 3-digit
!! exponent (the step this leaves in g s is below 1e-37)
!!
!! @param tj  t^(k)/k!
!! @param ej  (exp(-1/t))^(k)/k!

!Input variables
 real(dp), intent(in) :: tj(0:4)

!Output variables
 real(dp), intent(out) :: ej(0:4)

!Local variables
 real(dp) :: aj(0:4)

 ej(:)=0.0d0
 if(tj(0)<1.0d0/(30.0d0*log(10.0d0))) return
 call jet_inv(tj,aj)
 aj(:)=-aj(:)
 call jet_exp(aj,ej)
 end subroutine jet_expinv

 subroutine jet_inv(aj,bj)

!! Taylor coefficients of 1/a from those of a, to order 4
!!
!! @param aj  a^(k)/k!, aj(0) /= 0
!! @param bj  (1/a)^(k)/k!

!Input variables
 real(dp), intent(in) :: aj(0:4)

!Output variables
 real(dp), intent(out) :: bj(0:4)

!Local variables
 integer :: kk,jj

 bj(0)=1.0d0/aj(0)
 do kk=1,4
   bj(kk)=0.0d0
   do jj=1,kk
     bj(kk)=bj(kk) - aj(jj)*bj(kk-jj)
   end do
   bj(kk)=bj(kk)/aj(0)
 end do
 end subroutine jet_inv

 subroutine jet_exp(aj,ej)

!! Taylor coefficients of exp(a) from those of a, to order 4
!!
!! @param aj  a^(k)/k!
!! @param ej  (exp(a))^(k)/k!

!Input variables
 real(dp), intent(in) :: aj(0:4)

!Output variables
 real(dp), intent(out) :: ej(0:4)

!Local variables
 integer :: kk,jj

 ej(0)=exp(aj(0))
 do kk=1,4
   ej(kk)=0.0d0
   do jj=1,kk
     ej(kk)=ej(kk) + jj*aj(jj)*ej(kk-jj)
   end do
   ej(kk)=ej(kk)/kk
 end do
 end subroutine jet_exp

 subroutine jet_mul(aj,bj,cj)

!! Taylor coefficients of a b from those of a and b, to order 4
!!
!! @param aj,bj  a^(k)/k!, b^(k)/k!
!! @param cj  (a b)^(k)/k!

!Input variables
 real(dp), intent(in) :: aj(0:4),bj(0:4)

!Output variables
 real(dp), intent(out) :: cj(0:4)

!Local variables
 integer :: kk,jj

 do kk=0,4
   cj(kk)=0.0d0
   do jj=0,kk
     cj(kk)=cj(kk) + aj(jj)*bj(kk-jj)
   end do
 end do
 end subroutine jet_mul

 subroutine teter_h(xx,nder,hd)

!! h = sqrt(g) of the Teter function and its derivatives, as the sum of
!! shifted sinc functions (module description)
!!
!! @param xx  argument
!! @param nder  highest derivative wanted, 0<=nder<=4
!! @param hd  d^k h/dx^k in hd(k), k=0..nder

!Input variables
 real(dp), intent(in) :: xx
 integer, intent(in) :: nder

!Output variables
 real(dp), intent(out) :: hd(0:4)

!Local variables
 real(dp) :: sd(0:4)
 integer :: ii

 hd(:)=0.0d0
 do ii=-2,2
   call sinc_derivs(2.0d0*pi*xx - ii*pi,nder,sd)
   hd(0:nder)=hd(0:nder) + wsinc(ii)*sd(0:nder)
 end do
 do ii=1,nder
   hd(ii)=hd(ii)*(2.0d0*pi)**ii
 end do
 end subroutine teter_h

 subroutine sinc_derivs(yy,nder,sd)

!! S(y) = sin(y)/y and its derivatives: for |y| <= 2 the power series
!! S^(k) = sum_n (-1)**n y**(2n-k)/((2n+1) (2n-k)!), else the closed form
!! S^(k) = sum_j C(k,j) sin^(k-j)(y) (-1)**j j!/y**(j+1)
!!
!! @param yy  argument
!! @param nder  highest derivative wanted, 0<=nder<=4
!! @param sd  d^k S/dy^k in sd(k), k=0..nder

!Input variables
 real(dp), intent(in) :: yy
 integer, intent(in) :: nder

!Output variables
 real(dp), intent(out) :: sd(0:4)

!Local variables
 integer, parameter :: nser=16
 real(dp) :: sc(0:3),term,fact
 integer :: kk,jj,nn

 sd(:)=0.0d0
 if(abs(yy)<=2.0d0) then
   do kk=0,nder
     do nn=nser,0,-1
       if(2*nn<kk) exit
! y**(2n-k)/(2n-k)!
       term=1.0d0
       do jj=1,2*nn-kk
         term=term*yy/jj
       end do
       sd(kk)=sd(kk) + (-1)**nn*term/(2*nn+1)
     end do
   end do
 else
! sin^(m)(y) = sin(y + m pi/2)
   sc=[sin(yy), cos(yy), -sin(yy), -cos(yy)]
   do kk=0,nder
     fact=1.0d0
     do jj=0,kk
       if(jj>0) fact=-fact*jj
       sd(kk)=sd(kk) + binom(kk,jj)*sc(mod(kk-jj,4))*fact/yy**(jj+1)
     end do
   end do
 end if
 end subroutine sinc_derivs

 real(dp) function binom(nn,kk)

!! the binomial coefficient C(nn,kk)
!!
!! @param nn,kk  0<=kk<=nn

!Input variables
 integer, intent(in) :: nn,kk

!Local variables
 integer :: jj

 binom=1.0d0
 do jj=1,kk
   binom=binom*(nn-jj+1)/jj
 end do
 end function binom

 subroutine core_crossover(rr,rhoc,rhotps,mmax,fcfact,rstar,ierr)

!! the crossover radius: the outermost radius where the core charge equals
!! fcfact times the pseudo valence charge, the root of
!! ln(rhoc) - ln(fcfact rhotps) on the log mesh, converged to round-off;
!! the mesh interval of the outermost change of sign is the bracket
!!
!! @param rr  log radial mesh (bohr)
!! @param rhoc  core charge density
!! @param rhotps  pseudo valence charge density
!! @param mmax  size of the log radial mesh
!! @param fcfact  factor of the pseudo valence charge
!! @param rstar  the crossover radius (bohr)
!! @param ierr  0, or 1 if the core charge never exceeds fcfact rhotps

 use logmesh_m, only: lm_root

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rhoc(mmax),rhotps(mmax),fcfact

!Output variables
 real(dp), intent(out) :: rstar
 integer, intent(out) :: ierr

!Local variables
 real(dp), allocatable :: lf(:)
 integer :: ii,ircc

 rstar=0.0d0
 ierr=1
 ircc=0
 do ii=mmax,1,-1
   if(rhoc(ii)>fcfact*rhotps(ii)) then
     ircc=ii
     exit
   end if
 end do
 if(ircc==0 .or. ircc==mmax) return

! the logarithms where both are positive; elsewhere the far nodes, which
! the interpolant near the root does not reach
 allocate(lf(mmax))
 do ii=1,mmax
   lf(ii)=log(max(rhoc(ii),tiny(1.0d0))) &
&        - log(max(fcfact*rhotps(ii),tiny(1.0d0)))
 end do
 call lm_root(rr,lf,mmax,rr(ircc),rr(ircc+1),rstar,ierr)
 deallocate(lf)
 end subroutine core_crossover

 subroutine poly_core_fit(rstar,fmatch,a0,aco,iout,ifail)

!! the polynomial model core of icmod 1, p(r) = a0 + sum_j aco(j) r**(j+2),
!! j=1..5, that matches fmatch at rstar and is barely monotonic on
!! [0, rstar]: a0 the least for which p' <= 0 there, in closed form.
!! The coefficients are affine in a0: p = pA + a0 pB, pA the match with
!! a0 = 0 and pB the polynomial of r**3..r**7 with pB(rstar) = -1 and its
!! first 4 derivatives 0 there, pB' = c r**2 (r-rstar)**4, c = -105/rstar**7.
!! With pA' = r**2 qA(r), p' <= 0 on (0, rstar) is a0 >= phi(r) =
!! qA(r) rstar**7/(105 (r-rstar)**4), so a0 = max phi on [0, rstar) (phi
!! falls to -infinity at rstar, as qA(rstar) = rho'/rstar**2 < 0), at 0 or
!! at a root of phi', of the cubic (r-rstar) qA' - 4 qA (companion matrix)
!!
!! @param rstar  matching radius (bohr)
!! @param fmatch  the core charge and its first 4 derivatives at rstar
!! @param a0  constant term
!! @param aco  coefficients of r**3..r**7
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use lapack_m, only: dgesv, dhseqr

!Input variables
 real(dp), intent(in) :: rstar,fmatch(5)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: a0,aco(5)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: amat(5,5),acoa(5),acob(5),qa(0:4),cub(0:3),hess(3,3)
 real(dp) :: wr(3),wi(3),zz(1,1),work(16),rq,fac
 integer :: ipvt(5),info,jj,kk,mm,ncub

 ifail=0
! pA: the match with a0 = 0
 do jj=1,5
   do kk=0,4
     fac=1.0d0
     do mm=0,kk-1
       fac=fac*(jj+2-mm)
     end do
     amat(kk+1,jj)=fac*rstar**(jj+2-kk)
   end do
 end do
 acoa(:)=fmatch(:)
 call dgesv(5,1,amat,5,ipvt,acoa,5,info)
 if(info/=0) then
   write(iout,'(/a,i0)') 'poly_core_fit: ERROR singular match, dgesv info =', &
&        info
   ifail=1
   return
 end if

! pB = c (r**7/7 - 2 rstar r**6/3 + 6 rstar**2 r**5/5 - rstar**3 r**4
!         + rstar**4 r**3/3)
 fac=-105.0d0/rstar**7
 acob=fac*[rstar**4/3.0d0, -rstar**3, 1.2d0*rstar**2, &
&          -2.0d0*rstar/3.0d0, 1.0d0/7.0d0]

! qA(r) = pA'(r)/r**2 = sum_j (j+2) acoa(j) r**(j-1)
 do jj=1,5
   qa(jj-1)=(jj+2)*acoa(jj)
 end do

! the cubic (r-rstar) qA' - 4 qA = sum_k ((k-4) qa(k) - rstar (k+1) qa(k+1))
! r**k, k=0..3
 do kk=0,3
   cub(kk)=(kk-4)*qa(kk) - rstar*(kk+1)*qa(kk+1)
 end do

 a0=phi(0.0d0)
 ncub=3
 do while(ncub>0)
   if(cub(ncub)/=0.0d0) exit
   ncub=ncub-1
 end do
 if(ncub>0) then
! companion matrix of the monic cubic (or lower)
   hess(:,:)=0.0d0
   do kk=1,ncub
     hess(1,kk)=-cub(ncub-kk)/cub(ncub)
     if(kk>1) hess(kk,kk-1)=1.0d0
   end do
   call dhseqr('E','N',ncub,1,ncub,hess,3,wr,wi,zz,1,work,size(work),info)
   if(info/=0) then
     write(iout,'(/a,i0)') 'poly_core_fit: ERROR dhseqr info =',info
     ifail=1
     return
   end if
! every point of [0, rstar) bounds a0 from below, so a nearly real root
! may count by its real part
   do kk=1,ncub
     rq=newton_cubic(wr(kk))
     if(rq>0.0d0 .and. rq<rstar) a0=max(a0,phi(rq))
   end do
 end if
 aco(:)=acoa(:) + a0*acob(:)

 contains

 real(dp) function phi(rq)

!! the least a0 for which p'(rq) <= 0
!!
!! @param rq  radius, 0 <= rq < rstar

!Input variables
 real(dp), intent(in) :: rq

!Local variables
 integer :: mq

 phi=qa(4)
 do mq=3,0,-1
   phi=phi*rq + qa(mq)
 end do
 phi=phi*rstar**7/(105.0d0*(rq-rstar)**4)
 end function phi

 real(dp) function newton_cubic(r0)

!! two Newton steps on the cubic from r0, if they stay finite
!!
!! @param r0  start

!Input variables
 real(dp), intent(in) :: r0

!Local variables
 real(dp) :: cv,cd
 integer :: it,mq

 newton_cubic=r0
 do it=1,2
   cv=cub(ncub)
   cd=0.0d0
   do mq=ncub-1,0,-1
     cd=cd*newton_cubic + cv
     cv=cv*newton_cubic + cub(mq)
   end do
   if(cd==0.0d0) exit
   newton_cubic=newton_cubic - cv/cd
 end do
 end function newton_cubic

 end subroutine poly_core_fit

 subroutine teter_scale(ymatch,xx,ierr,iout,ifail)

!! the argument x of the Teter function, 0 < x < 1.5, at which its
!! logarithmic slope x g'(x)/g(x) = 2 x h'/h is ymatch (< 0); it falls
!! monotonically from 0 to -infinity there. Newton's method kept in the
!! bracket (else bisection), ending with one more step after a step of a
!! few ulps
!!
!! @param ymatch  r rho'/rho of the core charge at the matching radius
!! @param xx  the argument
!! @param ierr  0, or 1 if ymatch is not in the range of x g'/g on (0, 1.5)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use root_m, only: root_tol

!Input variables
 real(dp), intent(in) :: ymatch
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: xx
 integer, intent(out) :: ierr
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: xa,xb,hd(0:4),ff,fd,dx,ll
 integer :: it,nsmall

 ifail=0
 xx=0.0d0
 ierr=1
 if(.not. ymatch<0.0d0) return
 ierr=0
 xa=0.0d0
 xb=1.5d0
 xx=0.75d0
 nsmall=0
 do it=1,200
   call teter_h(xx,2,hd)
   ll=hd(1)/hd(0)
   ff=2.0d0*xx*ll - ymatch
   fd=2.0d0*(ll + xx*hd(2)/hd(0) - xx*ll**2)
   if(ff>0.0d0) then
     xa=xx
   else
     xb=xx
   end if
   if(nsmall>0) exit
   dx=-ff/fd
   if(.not. (xx+dx>xa .and. xx+dx<xb)) dx=0.5d0*(xa+xb)-xx
   if(abs(dx)<=root_tol(xx,0.0d0)) nsmall=1
   xx=xx+dx
 end do
 if(it>200) then
   write(iout,'(/a)') 'teter_scale: ERROR not converged'
   ifail=1
   return
 end if
 end subroutine teter_scale

 subroutine core_pieces(rr,rhoc,mmax,rj,fill,par,corepw,rhomod)

!! the model core with a joint at rj as two pieces, each column of the
!! model and its first 4 derivatives a piecewise function: the fill on
!! [0, rj], its closed form, and the core charge's derivatives (lm_derivs)
!! beyond, each continued lm_fh_reach()+2 nodes past rj (the interpolant's
!! reach) and 0 farther, where the fill may overflow; and the model on the
!! mesh, the fill below rj
!!
!! @param rr  log radial mesh (bohr)
!! @param rhoc  core charge density
!! @param mmax  size of the log radial mesh
!! @param rj  the joint (bohr), inside the mesh
!! @param fill  the fill and its first 4 derivatives at a radius
!! @param par  the fill's parameters
!! @param corepw  the pieces of each column
!! @param rhomod  model core charge and its first 4 derivatives on the mesh

 use logmesh_m, only: lm_derivs, lm_fh_reach
 use pieces_m, only: piecewise_t, pw_alloc, pw_form, pw_set_form

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rhoc(mmax),rj,par(:)
 procedure(pw_form) :: fill

!Output variables
 type(piecewise_t), intent(out) :: corepw(nrhomod)
 real(dp), intent(out) :: rhomod(mmax,nrhomod)

!Local variables
 real(dp) :: dd(nrhomod)
 integer :: ii,jj,kk,next

 next=lm_fh_reach()+2
 jj=count(rr(:)<=rj)
 do kk=1,nrhomod
   call pw_alloc(corepw(kk),mmax,2,[0.0d0,rj,rr(mmax)])
 end do
 do ii=1,min(mmax,jj+next)
   call fill(par,rr(ii),dd)
   do kk=1,nrhomod
     corepw(kk)%vals(ii,1)=dd(kk)
   end do
 end do
 do kk=1,nrhomod
   call pw_set_form(corepw(kk),1,fill,par,kk)
 end do
 do ii=max(1,jj+1-next),mmax
   call lm_derivs(rr,rhoc,mmax,rr(ii),4,dd)
   do kk=1,nrhomod
     corepw(kk)%vals(ii,2)=dd(kk)
   end do
 end do
 do kk=1,nrhomod
   do ii=1,mmax
     if(rr(ii)<rj) then
       rhomod(ii,kk)=corepw(kk)%vals(ii,1)
     else
       rhomod(ii,kk)=corepw(kk)%vals(ii,2)
     end if
   end do
 end do
 return
 end subroutine core_pieces

 subroutine core_pieces_smooth(rr,rhomod,mmax,corepw,form,par)

!! a model core with no joint as pieces of one segment each, a column
!! ending at its last nonzero node (the end of a tapered Teter function),
!! with the model's closed form if it has one
!!
!! @param rr  log radial mesh (bohr)
!! @param rhomod  model core charge and its first 4 derivatives on the mesh
!! @param mmax  size of the log radial mesh
!! @param corepw  the pieces of each column
!! @param form  the model and its first 4 derivatives at a radius
!! @param par  its parameters

 use pieces_m, only: piecewise_t, pw_alloc, pw_form, pw_set_form

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rhomod(mmax,nrhomod)
 procedure(pw_form), optional :: form
 real(dp), intent(in), optional :: par(:)

!Output variables
 type(piecewise_t), intent(out) :: corepw(nrhomod)

!Local variables
 integer :: kk,nlast

 do kk=1,nrhomod
   nlast=mmax
   do while(nlast>1 .and. rhomod(nlast,kk)==0.0d0)
     nlast=nlast-1
   end do
   call pw_alloc(corepw(kk),mmax,1,[0.0d0,rr(nlast)])
   corepw(kk)%vals(:,1)=rhomod(:,kk)
   if(present(form)) call pw_set_form(corepw(kk),1,form,par,kk)
 end do
 return
 end subroutine core_pieces_smooth

 subroutine poly_form(par,rad,dd)

!! the polynomial a0 + sum_j aco(j) r**(j+2), j=1..5, of icmod 1
!! (poly_core_fit) and its first 4 derivatives
!!
!! @param par  a0, aco(1:5)
!! @param rad  radius (bohr)
!! @param dd  the polynomial and its derivatives at rad

!Input variables
 real(dp), intent(in) :: par(:),rad

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 real(dp) :: fac
 integer :: ii,jj,ll

 do ii=1,5
   if(ii==1) then
     dd(ii)=par(1)
   else
     dd(ii)=0.0d0
   end if
   do jj=1,5
! the term of r**(jj+2) differentiated ii-1 times, none past its degree
     if(jj+3-ii<0) cycle
     fac=1.0d0
     do ll=0,ii-2
       fac=fac*(jj+2-ll)
     end do
     dd(ii)=dd(ii)+fac*par(jj+1)*rad**(jj+3-ii)
   end do
 end do
 return
 end subroutine poly_form

 subroutine teter_form(par,rad,dd)

!! the model a g(r/b) s(r/b) of the Teter function (teter_core_derivs) and
!! its first 4 derivatives
!!
!! @param par  the amplitude a and the scale b (bohr)
!! @param rad  radius (bohr)
!! @param dd  the model and its derivatives at rad

!Input variables
 real(dp), intent(in) :: par(:),rad

!Output variables
 real(dp), intent(out) :: dd(5)

!Local variables
 real(dp) :: gd(0:4)
 integer :: kk

 call teter_core_derivs(rad/par(2),4,gd)
 do kk=0,4
   dd(kk+1)=par(1)*gd(kk)/par(2)**kk
 end do
 return
 end subroutine teter_form

 subroutine d2exc_baseline(rhotae,rhoc,rhoae,rhotps,rhops,rr,zion,iexc, &
&                          c_light,nc,nv,la,rcut,mmax,adams,rhomod,d2excae, &
&                          d2excps, &
&                          d2mdiff,iout)

!! the start of the model core analysis the model cores share: the
!! all-electron second derivatives of Exc in the occupations, and those of
!! the pseudo-atom with no core correction, written with their rms
!! difference
!!
!! @param rhotae  all-electron valence charge density
!! @param rhoc  core charge density
!! @param rhoae  all-electron valence charges
!! @param rhotps  pseudo valence charge density
!! @param rhops  pseudo valence charges
!! @param rr  log radial mesh (bohr)
!! @param zion  ion charge
!! @param iexc  exchange-correlation functional
!! @param c_light  speed of light (Hartree atomic units)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular momenta of the states
!! @param rcut  radius of the Exc integrals
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param rhomod  model core charge; zero on return
!! @param d2excae  all-electron second derivatives
!! @param d2excps  pseudo-atom second derivatives with no core correction
!! @param d2mdiff  their rms difference
!! @param iout  unit of the text output

 use der2exc_m, only: der2exc
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: iexc,nc,nv,mmax
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 integer, intent(in) :: la(nv+nc)
 real(dp), intent(in) :: rhotae(mmax),rhoc(mmax),rhoae(mmax,nv)
 real(dp), intent(in) :: rhotps(mmax),rhops(mmax,nv),rr(mmax)
 real(dp), intent(in) :: zion,rcut
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: d2excae(nv,nv),d2excps(nv,nv)

!Output variables
 real(dp), intent(out) :: rhomod(mmax,nrhomod)
 real(dp), intent(out) :: d2mdiff

!Local variables
 integer :: jj,kk

 write(iout,'(/a/a)') 'Model core correction analysis',&
&                  '  based on d2Exc/dn_idn_j'

! get derivatives of all-electron xc energy

 call der2exc(rhotae,rhoc,rhoae,rr,d2excae,d2excps,d2mdiff, &
&                   zion,iexc,c_light,nc,nv,la,rcut,mmax,adams)

 write(iout,'(/a/)') 'd2excae - all-electron derivatives'
 do kk=1,nv
   write(iout,'(1p,4d16.6)') (d2excae(kk,jj),jj=1,nv)
 end do

! compute d2excps with no core correction
 rhomod(:,:)=0.0d0

 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
&                   zion,iexc,c_light,nc,nv,la,rcut,mmax,adams)

 call d2exc_report('d2excps - pseudofunction derivatives with no core '// &
&                  'correction',d2excps,nv,d2mdiff,iout)

 return
 end subroutine d2exc_baseline

 subroutine d2exc_report(title,d2excps,nv,d2mdiff,iout)

!! writes the pseudo-atom's second derivatives of Exc under title, and
!! their rms difference from the all-electron ones
!!
!! @param title  the heading
!! @param d2excps  the pseudo-atom's second derivatives
!! @param nv  number of valence states
!! @param d2mdiff  their rms difference from the all-electron ones
!! @param iout  unit of the text output

!Input variables
 character(len=*), intent(in) :: title
 integer, intent(in) :: nv
 real(dp), intent(in) :: d2excps(nv,nv),d2mdiff
 integer, intent(in) :: iout

!Local variables
 integer :: jj,kk

 write(iout,'(/a/)') title
 do kk=1,nv
   write(iout,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
 end do
 write(iout,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

 return
 end subroutine d2exc_report

 end module modcore_m
