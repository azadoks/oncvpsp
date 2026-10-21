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

 module modcore_m
 use constants_m, only: dp, pi
 use logmesh_m, only: lm_root
 implicit none
 private
 public :: teter_derivs, teter_core_derivs, core_crossover, poly_core_fit
 public :: teter_scale

!> the weights of the shifted sinc functions of h, at the shifts k pi
 real(dp), parameter :: wsinc(-2:2)=[1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0, &
&                                    2.0d0/3.0d0, 1.0d0/6.0d0]

!> the taper of the model's Teter function: 1 up to xtap1, 0 from xtap2
 real(dp), parameter :: xtap1=3.0d0, xtap2=4.0d0

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
!! xtap2 (as ionic_switch); the derivatives of s by Taylor arithmetic
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
!! @param rr  log radial mesh
!! @param rhoc  core charge
!! @param rhotps  pseudo valence charge
!! @param mmax  size of the log radial mesh
!! @param fcfact  factor of the pseudo valence charge
!! @param rstar  the crossover radius (bohr)
!! @param ierr  0, or 1 if the core charge never exceeds fcfact rhotps

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

 subroutine poly_core_fit(rstar,fmatch,a0,aco)

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

!Input variables
 real(dp), intent(in) :: rstar,fmatch(5)

!Output variables
 real(dp), intent(out) :: a0,aco(5)

!Local variables
 real(dp) :: amat(5,5),acoa(5),acob(5),qa(0:4),cub(0:3),hess(3,3)
 real(dp) :: wr(3),wi(3),zz(1,1),work(16),rq,fac
 integer :: ipvt(5),info,jj,kk,mm,ncub
 external :: dgesv,dhseqr

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
   write(6,'(/a,i4)') 'poly_core_fit: ERROR singular match, dgesv info =', &
&        info
   error stop 1
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
     write(6,'(/a,i4)') 'poly_core_fit: ERROR dhseqr info =',info
     error stop 1
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

 subroutine teter_scale(ymatch,xx,ierr)

!! the argument x of the Teter function, 0 < x < 1.5, at which its
!! logarithmic slope x g'(x)/g(x) = 2 x h'/h is ymatch (< 0); it falls
!! monotonically from 0 to -infinity there. Newton's method kept in the
!! bracket (else bisection), ending with one more step after a step of a
!! few ulps
!!
!! @param ymatch  r rho'/rho of the core charge at the matching radius
!! @param xx  the argument
!! @param ierr  0, or 1 if ymatch is not in the range of x g'/g on (0, 1.5)

!Input variables
 real(dp), intent(in) :: ymatch

!Output variables
 real(dp), intent(out) :: xx
 integer, intent(out) :: ierr

!Local variables
 real(dp) :: xa,xb,hd(0:4),ff,fd,dx,ll
 integer :: it,nsmall

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
   if(abs(dx)<=4.0d0*spacing(xx)) nsmall=1
   xx=xx+dx
 end do
 if(it>200) then
   write(6,'(/a)') 'teter_scale: ERROR not converged'
   error stop 1
 end if
 end subroutine teter_scale

 end module modcore_m
