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
!> The chirp-z transform X_m = sum_{n=0}^{N-1} x_n exp(i theta n m),
!> m = 0..M-1, for any theta, by Bluestein's convolution: with
!> n m = (n**2 + m**2 - (m-n)**2)/2,
!>
!>     X_m = w_m sum_n (x_n w_n) conj(w_{m-n}),  w_k = exp(i theta k**2/2),
!>
!> a convolution of length N + M - 1, done by radix-2 FFTs of the next
!> power of 2, L: O(L log L) for the O(N M) of the sums. eresid takes the
!> transforms of its tail on a uniform q grid this way, where the nodes of
!> both meshes are uniform. The phases theta k**2/2 reach 1e4 rad, so
!> they are reduced mod 2 pi in quad precision before their sine and
!> cosine (in double, 1e-12 of error at that size); the FFT's twiddles are
!> each evaluated directly. The error of X is then about 1e-16 of
!> sum |x_n| (log L of it), against the round-off of the sines of the
!> direct sums.

 module czt_m
 use constants_m, only: dp, qp, pi
 implicit none
 private
 public :: czt_plan_t, czt_plan, czt_apply, czt_free

!> 2 pi in quad precision
 real(qp), parameter :: twopiq= &
&  6.283185307179586476925286766559005768394_qp

!> a transform of N points to M, ready to apply to any x
 type :: czt_plan_t
   integer :: nn=0,mm=0,ll=0
!> the chirps w_n (n<N) and w_m (m<M), the FFT of the convolution's
!> kernel, and the FFT's twiddles exp(-2 pi i k/L)
   complex(dp), allocatable :: wn(:),wm(:),hf(:),tw(:)
 end type czt_plan_t

 contains

 subroutine czt_plan(nn,mm,theta,plan)

!! the plan of the transform of nn points to mm at the angle theta (above)
!!
!! @param nn  number of input points N
!! @param mm  number of outputs M
!! @param theta  the angle step
!! @param plan  the plan

!Input variables
 integer, intent(in) :: nn,mm
 real(dp), intent(in) :: theta

!Output variables
 type(czt_plan_t), intent(out) :: plan

!Local variables
 complex(dp), allocatable :: hh(:)
 integer :: kk,ll

 ll=1
 do while(ll<nn+mm-1)
   ll=2*ll
 end do
 plan%nn=nn
 plan%mm=mm
 plan%ll=ll
 allocate(plan%wn(0:nn-1),plan%wm(0:mm-1),plan%hf(0:ll-1))
 allocate(plan%tw(0:max(ll/2,1)-1))
 do kk=0,ll/2-1
   plan%tw(kk)=cmplx(cos(2.0d0*pi*dble(kk)/dble(ll)), &
&                    -sin(2.0d0*pi*dble(kk)/dble(ll)),dp)
 end do
 do kk=0,nn-1
   plan%wn(kk)=chirp(theta,kk)
 end do
 do kk=0,mm-1
   if(kk<nn) then
     plan%wm(kk)=plan%wn(kk)
   else
     plan%wm(kk)=chirp(theta,kk)
   end if
 end do

!the kernel conj(w_k), k = -(N-1)..M-1, wrapped onto 0..L-1
 allocate(hh(0:ll-1))
 hh(:)=(0.0d0,0.0d0)
 do kk=0,mm-1
   hh(kk)=conjg(plan%wm(kk))
 end do
 do kk=1,nn-1
   hh(ll-kk)=conjg(plan%wn(kk))
 end do
 call fft(hh,ll,plan%tw,.false.)
 plan%hf(:)=hh(:)
 deallocate(hh)
 return
 end subroutine czt_plan

 complex(dp) function chirp(theta,kk)

!! exp(i theta k**2/2), the phase reduced mod 2 pi in quad precision
!!
!! @param theta  the angle step
!! @param kk  the index k

!Input variables
 real(dp), intent(in) :: theta
 integer, intent(in) :: kk

!Local variables
 real(qp) :: ph
 real(dp) :: pd

 ph=real(theta,qp)*(real(kk,qp)*real(kk,qp)/2.0_qp)
 ph=ph-twopiq*anint(ph/twopiq)
 pd=real(ph,dp)
 chirp=cmplx(cos(pd),sin(pd),dp)
 return
 end function chirp

 subroutine czt_apply(plan,xx,yy)

!! X = the transform of x (above) by the plan
!!
!! @param plan  the plan
!! @param xx  the N inputs x_n
!! @param yy  the M outputs X_m

!Input variables
 type(czt_plan_t), intent(in) :: plan
 complex(dp), intent(in) :: xx(0:plan%nn-1)

!Output variables
 complex(dp), intent(out) :: yy(0:plan%mm-1)

!Local variables
 complex(dp), allocatable :: aa(:)
 integer :: ll

 ll=plan%ll
 allocate(aa(0:ll-1))
 aa(:)=(0.0d0,0.0d0)
 aa(0:plan%nn-1)=xx(:)*plan%wn(:)
 call fft(aa,ll,plan%tw,.false.)
 aa(:)=aa(:)*plan%hf(:)
 call fft(aa,ll,plan%tw,.true.)
 yy(:)=aa(0:plan%mm-1)*plan%wm(:)/dble(ll)
 deallocate(aa)
 return
 end subroutine czt_apply

 subroutine czt_free(plan)

!! releases the plan
!!
!! @param plan  the plan

!Input/output variables
 type(czt_plan_t), intent(inout) :: plan

 if(allocated(plan%wn)) deallocate(plan%wn,plan%wm,plan%hf,plan%tw)
 plan%nn=0
 plan%mm=0
 plan%ll=0
 return
 end subroutine czt_free

 subroutine fft(aa,nn,tw,inverse)

!! in place, sum_j a_j exp(-2 pi i j k/n) (inverse: exp(+...), not divided
!! by n), n a power of 2: radix-2 decimation in time
!!
!! @param aa  the values, then their transform
!! @param nn  their number n
!! @param tw  the twiddles exp(-2 pi i k/n), k<n/2
!! @param inverse  the inverse transform

!Input variables
 integer, intent(in) :: nn
 complex(dp), intent(in) :: tw(0:max(nn/2,1)-1)
 logical, intent(in) :: inverse

!Input/output variables
 complex(dp), intent(inout) :: aa(0:nn-1)

!Local variables
 complex(dp) :: tt,ww
 integer :: ii,jj,kk,mm,half,step

!bit-reversed order
 jj=0
 do ii=0,nn-2
   if(ii<jj) then
     tt=aa(ii)
     aa(ii)=aa(jj)
     aa(jj)=tt
   end if
   kk=nn/2
   do while(kk<=jj)
     jj=jj-kk
     kk=kk/2
   end do
   jj=jj+kk
 end do

!butterflies
 half=1
 do while(half<nn)
   step=nn/(2*half)
   do mm=0,half-1
     ww=tw(mm*step)
     if(inverse) ww=conjg(ww)
     do ii=mm,nn-1,2*half
       tt=ww*aa(ii+half)
       aa(ii+half)=aa(ii)-tt
       aa(ii)=aa(ii)+tt
     end do
   end do
   half=2*half
 end do
 return
 end subroutine fft

 end module czt_m
