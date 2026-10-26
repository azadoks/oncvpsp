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

 module optimize_m

 implicit none
 private
 public :: optimize

 contains

! calculates convergence-optimized pseudo-wave-function coefficients
! in orthonormal and spherical Bessel function bases

 subroutine optimize(nnull,nbas,pswf0_sb,pswf0_or,nqout,qout,&
&                    eresid0,eresiddot,eresidmat,&
&                    pswfnull_sb,pswfnull_or,uunorm,ps0norm,eresidmin,&
&                    pswfopt_sb,pswfopt_or,ekin_anal,eresq,pswfopt_null, &
&                    wreg,optdg)

!! @param nnull  number of unconstrained basis vectors for residual minimization
!! @param nbas  number of sbf basis functions
!! @param pswf0_sb(nbas)  sbf basis  coefficients for constrained part of pseudo
!!                        wave function
!! @param pswf0_or(nbas)  or basis  coefficients
!! @param nqout  number of q values in [0,qmax] for which results are to be
!!               saved
!! @param qout(nqout)  q values in [0,qmax] for which residuals are to be
!!                     calculated
!! @param eresid0(nqout)  set of <pswf0| E_resid |pswf0> matrix elements
!! @param eresiddot(nnull,nqout)  set of <pswfnull| E_resid |pswf0> matrix
!!                                elements
!! @param eresidmat(nnull,nnull,nqout)  set of <pswfnull | E_resid |pswfnull'>
!! @param pswfnull_sb(nbas,nnull)  sbf coefficients of null-space eigenfunctions
!! @param pswfnull_or(nbas,nnull)  or basis coefficients of null-space
!!                                 eigenfunctions
!! @param uunorm  all-electron charge inside rc
!! @param ps0norm  ps0 charge inside rc
!! @param eresidmin  minimized E_resid for qout(1), the optimization cutoff
!! @param pswfopt_sb  optimized pseudowavefunction sbf coefficients
!! @param pswfopt_or  optimized pseudowavefunction or basis coefficients
!! @param ekin_anal  total pswf kinetic energy calculated analytically from
!!                   E_resid
!! @param eresq(nqout)  E_resid for optimized pswf as a function of cutoff
!! @param pswfopt_null  null-space coefficients of the optimized pswf,
!!                      pswfopt_sb = pswf0_sb
!!                                   + sum_i pswfopt_null(i)*pswfnull_sb(:,i)
!! @param wreg  the regularization: with wreg>0, wreg times E_resid at
!!              qout(3) is added to the E_resid at qout(1) that is
!!              minimized. Where qcut*rc is large, E_resid(qcut) is flat to
!!              round-off in the null space (e1 ~ 1e-18 Ha) and the sign of
!!              x_1, a pair of waves one of which has a node, is set by the
!!              quadrature's rounding; E_resid at qcut/4 is not flat, so the
!!              wave is determined and smooth in the inputs
!! @param optdg  E_resid(qout(1)) at the optimum, the regularization's
!!               term wreg E_resid(qout(3)) there, and the smallest
!!               eigenvalue e1 and |f1| of the problem at qout(1) alone

 use constants_m, only: dp
 use lapack_m, only: dsyev
 implicit none

!Input variables
 integer, intent(in) :: nnull,nbas,nqout
 real(dp), intent(in) :: pswf0_sb(nbas),pswfnull_sb(nbas,nnull)
 real(dp), intent(in) :: pswf0_or(nbas),pswfnull_or(nbas,nnull)
 real(dp), intent(in) :: eresid0(nqout),eresiddot(nnull,nqout)
 real(dp), intent(in) :: eresidmat(nnull,nnull,nqout)
 real(dp), intent(in) :: qout(nqout)
 real(dp), intent(in) :: uunorm,ps0norm

!Output variables
 real(dp), intent(out) :: eresq(nqout)
 real(dp), intent(out) :: eresidmin,ekin_anal
 real(dp), intent(out) :: pswfopt_sb(nbas),pswfopt_or(nbas)
 real(dp), intent(in) :: wreg
 real(dp), intent(out) :: pswfopt_null(nnull)
 real(dp), intent(out) :: optdg(4)

!Local variables
 real(dp) :: emin
 real(dp), allocatable :: work(:),wmat(:,:),wev(:),wvec(:)
 real(dp), allocatable :: eresidevec_nu(:,:)
 real(dp), allocatable :: pswfopt(:),pswfopt_nu(:)
 real(dp), allocatable :: eresiddot_ev(:),eresideval(:)
 real(dp), allocatable :: rmat(:,:),rev(:)
 integer :: ii,jj,kk,info

!find eigenvalues and eigenvectors of saved E_resid matrix for qout(1)
!which is meant to be the cutoff used for optimization

 allocate(work(5*nnull),wmat(nnull,nnull),wev(nnull),wvec(nnull))
 allocate(eresidevec_nu(nnull,nnull),eresideval(nnull))
 allocate(eresiddot_ev(nnull))

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

 wmat(:,:)=eresidmat(:,:,1)
 if(wreg>0.0d0) wmat(:,:)=wmat(:,:)+wreg*eresidmat(:,:,3)

 call dsyev( 'V', 'U', nnull, wmat, nnull, wev, work, 5*nnull, info )
 if(info /= 0) then
   write(6,'(a,i4)') 'optimize: residual-energy eigenvalue ERROR, info=',info
   error stop 1
 end if

!convert <ps0| E_resid |psnull> vector from the null-space representation
!to the eigenvector representation

 wvec(:)=eresiddot(:,1)
 if(wreg>0.0d0) wvec(:)=wvec(:)+wreg*eresiddot(:,3)
 eresiddot_ev(:)=0.0d0
 do jj=1,nnull
   do kk=1,nnull
     eresiddot_ev(jj)=eresiddot_ev(jj)+wvec(kk)*wmat(kk,jj)
   end do
 end do

!e1 and |f1| of the problem at qout(1) alone: those just found, or with
!the regularization those of E_resid(qout(1)) by itself. Both at
!round-off (1e-18 and 1e-16 Ha where it is flat, at least 9e-9 in the
!reference inputs) leave the wave undetermined
 if(wreg>0.0d0) then
   allocate(rmat(nnull,nnull),rev(nnull))
   rmat(:,:)=eresidmat(:,:,1)
   call dsyev( 'V', 'U', nnull, rmat, nnull, rev, work, 5*nnull, info )
   if(info /= 0) then
     write(6,'(a,i4)') 'optimize: residual-energy eigenvalue ERROR, info=', &
&          info
     error stop 1
   end if
   optdg(3)=rev(1)
   optdg(4)=abs(sum(eresiddot(:,1)*rmat(:,1)))
   deallocate(rmat,rev)
 else
   optdg(3)=wev(1)
   optdg(4)=abs(eresiddot_ev(1))
   if(optdg(3)<1.0d-12 .and. optdg(4)<1.0d-12) then
     write(6,'(/a,1p,2e10.2/a)') 'optimize: WARNING E_resid flat to '// &
&          'round-off at qcut, e1, |f1|=',optdg(3),optdg(4), &
&          '  the wave is not determined (see wave_regularization)'
   end if
 end if

 eresideval(:)=wev(:)
 eresidevec_nu(:,:)=wmat(:,:)

 deallocate(work,wmat,wev,wvec)

!find coefficients of E_resid eigenvectors which must be added to pswf0
!to minimize E_resid

!note that E_resid is an exteremely simple function in this representation,
!consisting of a constant plus as sum of linear and quadratic terms in each of
!these coefficients

 allocate(pswfopt(nnull),pswfopt_nu(nnull))
 pswfopt(:)=0.d0
 pswfopt_nu(:)=0.d0

! Interval-halving search to minimize E*r

! This is all based on Eqs. (14) - (16) in my paper as correccted in the
! erratum.  It is obvious from Eq.(14) that to minimize E^r, the signs of x_i
! must be the opposite of those of f_i.  Eq.(16) is derived by treating
! x_1 as depenent upon x_2...x_N-M through Eq.(15).  We can then turn
! around and treat those x's as functions all dependent on x_1.  Since the
! signs of x_1 and f_1 must be opposite at the minimum, and since e_1
! is the smallest eigenvalue, the denominator of Eq.(16) is always positive,
! and the magnitudes of the x_i are all monotonically increasing functions
! of x_1.  The limits of |x_1| are zero and D_norm (yy here), so a simple
! interval-halving search varying |x_1| to satisfy the norm constraint
! (Eq.(15)) should be robust.

! Translations between the notation in the paper and here are as follows:
!   pswfopt(ii) = x_i
!   eresiddot_ev(ii) = f_i
!   eresideval(ii) = e_i
!   D_norm == yy
!   E^r_00 = eresid0(1)
!   N-M = nnull


 call solve_norm_constraint(nnull,eresideval,eresiddot_ev,uunorm,ps0norm, &
&                            pswfopt)

 emin=eresid0(1)
 do ii=1,nnull
   emin=emin+(2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
&             +eresideval(ii)*pswfopt(ii)**2)
 end do
 eresidmin=emin

 ekin_anal=0.0d0
!find cutoff-dependence of E_resid for optimized pswf
!must transform pswfopt to the null basis representation to be
!consistent with the input dot product and
 allocate(work(nnull))
 do ii=1,nnull
   pswfopt_nu(ii)=0.0d0
   do jj=1,nnull
     pswfopt_nu(ii)=pswfopt_nu(ii)+eresidevec_nu(ii,jj)*pswfopt(jj)
   end do
 end do

 do kk=1,nqout
   emin=eresid0(kk)
   work(:)=0.0d0
   do ii=1,nnull
     emin=emin+2.0d0*eresiddot(ii,kk)*pswfopt_nu(ii)
     work(:)=work(:)+eresidmat(:,ii,kk)*pswfopt_nu(ii)
   end do
   do ii=1,nnull
     emin=emin+pswfopt_nu(ii)*work(ii)
   end do
   eresq(kk)=emin
   if(qout(kk)==0.0d0) ekin_anal=emin
 end do

!with the regularization, the E_resid at qout(1) of the wave found, and
!the regularization's term
 optdg(1)=eresq(1)
 if(wreg>0.0d0) then
   eresidmin=eresq(1)
   optdg(2)=wreg*eresq(3)
 else
   optdg(2)=0.0d0
 end if

!calculate optimized pswf in sbf representation
 pswfopt_sb(:)=pswf0_sb(:)
 pswfopt_or(:)=pswf0_or(:)
 do ii=1,nnull
   pswfopt_sb(:)=pswfopt_sb(:)+pswfopt_nu(ii)*pswfnull_sb(:,ii)
   pswfopt_or(:)=pswfopt_or(:)+pswfopt_nu(ii)*pswfnull_or(:,ii)
 end do

 pswfopt_null(:)=pswfopt_nu(:)
 deallocate(pswfopt,pswfopt_nu)
 return

 contains

 subroutine solve_norm_constraint(nnull,eresideval,eresiddot_ev,uunorm, &
&                                 ps0norm,pswfopt)
!! the x_i of Eq. (16) at the |x_1| that satisfies the norm constraint,
!! Eq. (15): by interval halving on |x_1| from 0 to D_norm, then Newton's
!! method; x_1 of the sign opposite to f_1's
!!
!! @param nnull  dimension of the null space
!! @param eresideval  eigenvalues e_i of E_resid in the null space
!! @param eresiddot_ev  linear terms f_i in the eigenvector representation
!! @param uunorm  all-electron norm inside rc
!! @param ps0norm  norm of the constrained wave alone
!! @param pswfopt  the x_i

!Input variables
 integer, intent(in) :: nnull
 real(dp), intent(in) :: eresideval(nnull),eresiddot_ev(nnull)
 real(dp), intent(in) :: uunorm,ps0norm

!Output variables
 real(dp), intent(inout) :: pswfopt(nnull)

!Local variables
 real(dp) :: yy,tt,rtsgn,x1,x1min,x1max,dt,dtold,dx,den
 real(dp), parameter :: eps=1.d-10
 integer :: ii,iter
 integer, parameter :: niter=100
 logical :: converged

 rtsgn=-sign(1.0d0, eresiddot_ev(1))
 if(ps0norm>uunorm) then
   write(6,'(/a)') 'optimize: ERROR ps0norm > uunorm, code will stop'
   error stop 1
 end if
 yy=sqrt((uunorm-ps0norm))

 x1min=0.0d0
 x1max=yy
 converged=.false.

 do iter=1,niter

   x1=0.5d0*(x1max+x1min)

   tt=x1**2-yy**2
   do ii=2,nnull
     pswfopt(ii)=-eresiddot_ev(ii)/(eresideval(ii)-eresideval(1) &
&                 +abs(eresiddot_ev(1))/x1)
     tt=tt+pswfopt(ii)**2
   end do

   if(x1max-x1min<eps) then
     pswfopt(1)=rtsgn*x1
     converged=.true.
     exit
   end if

   if(tt<0.0d0) then
     x1min=x1
   else
     x1max=x1
   end if
 end do

 if(converged .eqv. .false.) then
   write(6,'(/a)') 'optimize: ERROR interval-halving search not converged'
   error stop 1
 end if

! Newton's method on the constraint tt(x1) from the bracket's midpoint:
! tt is smooth and increasing in x1, so a step takes x1 from the bracket's
! eps to round-off, and the x_i do not depend on where the halving stopped
! (doc/smoothness_plan.md, Q-6). Converged when a step was within the
! round-off of tt (a few ulps of yy**2) or a few ulps of x1, the x_i then
! evaluated at the moved x1
 converged=.false.
 dx=0.0d0
 dtold=0.0d0
 do iter=1,6
   tt=x1**2-yy**2
   dt=2.0d0*x1
   do ii=2,nnull
     den=eresideval(ii)-eresideval(1)+abs(eresiddot_ev(1))/x1
     pswfopt(ii)=-eresiddot_ev(ii)/den
     tt=tt+pswfopt(ii)**2
     dt=dt+2.0d0*pswfopt(ii)**2*abs(eresiddot_ev(1))/(x1**2*den)
   end do
   if(iter>1 .and. (abs(dx)*dtold<=8.0d0*epsilon(1.0d0)*yy**2 .or. &
&                   abs(dx)<=4.0d0*spacing(x1))) then
     converged=.true.
     exit
   end if
   dx=-tt/dt
   dtold=dt
   x1=x1+dx
 end do
 if(converged .eqv. .false.) then
   write(6,'(/a,1p,e10.2)') &
&        'optimize: ERROR Newton polish not converged, last step',dx
   error stop 1
 end if
 pswfopt(1)=rtsgn*x1
 end subroutine solve_norm_constraint

 end subroutine optimize

 end module optimize_m
