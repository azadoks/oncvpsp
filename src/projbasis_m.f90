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
!> The steps of orthonormal projectors run_vkb and sr_so_r share: the
!> square roots S^(+-1/2) of an overlap matrix from its eigenvectors, the
!> congruence S^(1/2)^T A S^(1/2), linear combinations of piecewise
!> functions, the sort of the coefficients by magnitude, and the sign of a
!> projector, positive at its peak.
 module projbasis_m
 implicit none
 private
 public :: pb_half_powers, pb_congruence, pb_combine, pb_sort_abs
 public :: pb_sign_at_peak

 contains

 subroutine pb_half_powers(sv,sev,nn,ld,sphalf,smhalf,rfloor,ndrop)

!! S^(1/2) and S^(-1/2) in the eigenvectors sv of S, each column scaled by
!! the square root of its eigenvalue sev, or divided by it; with rfloor,
!! the columns whose eigenvalue is at or below rfloor times the largest
!! are zero in both (a pseudo-inverse: those directions are dropped)
!!
!! @param sv  the eigenvectors of S (dsyev)
!! @param sev  its eigenvalues
!! @param nn  the order of S
!! @param ld  the leading dimension
!! @param sphalf  S^(1/2) (output)
!! @param smhalf  S^(-1/2) (output)
!! @param rfloor  relative floor of the eigenvalues (optional)
!! @param ndrop  number of columns dropped (output, with rfloor)

 use constants_m, only: dp

!Input variables
 integer, intent(in) :: nn,ld
 real(dp), intent(in) :: sv(ld,ld),sev(ld)
 real(dp), intent(in), optional :: rfloor

!Input/output variables
 real(dp), intent(inout) :: sphalf(ld,ld),smhalf(ld,ld)

!Output variables
 integer, intent(out), optional :: ndrop

!Local variables
 real(dp) :: tt,smin
 integer :: ii,jj

 smin=-huge(1.0d0)
 if(present(rfloor)) smin=rfloor*maxval(sev(1:nn))
 if(present(ndrop)) ndrop=0
 do jj=1,nn
   if(sev(jj)<=smin) then
     sphalf(1:nn,jj)=0.0d0
     smhalf(1:nn,jj)=0.0d0
     if(present(ndrop)) ndrop=ndrop+1
     cycle
   end if
   tt=sqrt(sev(jj))
   do ii=1,nn
     sphalf(ii,jj)=tt*sv(ii,jj)
     smhalf(ii,jj)=sv(ii,jj)/tt
   end do
 end do
 return
 end subroutine pb_half_powers

 subroutine pb_congruence(aa,sh,nn,ld,ast)

!! A* = S^(1/2)^T A S^(1/2)
!!
!! @param aa  A
!! @param sh  S^(1/2)
!! @param nn  the order
!! @param ld  the leading dimension
!! @param ast  A* (output; zero beyond nn)

 use constants_m, only: dp

!Input variables
 integer, intent(in) :: nn,ld
 real(dp), intent(in) :: aa(ld,ld),sh(ld,ld)

!Output variables
 real(dp), intent(out) :: ast(ld,ld)

!Local variables
 real(dp) :: at(ld,ld)
 integer :: ii,jj,kk

 at(:,:)=0.0d0
 ast(:,:)=0.0d0
 do ii=1,nn
   do jj=1,nn
     do kk=1,nn
       at(ii,jj)=at(ii,jj)+aa(ii,kk)*sh(kk,jj)
     end do
   end do
 end do
 do ii=1,nn
   do jj=1,nn
     do kk=1,nn
       ast(ii,jj)=ast(ii,jj)+at(kk,jj)*sh(kk,ii)
     end do
   end do
 end do
 return
 end subroutine pb_congruence

 subroutine pb_combine(cc,ld,pin,nn,pout)

!! the combinations pout(j) = sum_i cc(i,j) pin(i) of piecewise functions
!! on the same pieces; each pout(j) keeps the pieces it has on input
!!
!! @param cc  the coefficients
!! @param ld  their leading dimension
!! @param pin  the functions
!! @param nn  how many
!! @param pout  the combinations (pieces on input, values on output)

 use constants_m, only: dp
 use pieces_m, only: piecewise_t

!Input variables
 integer, intent(in) :: ld,nn
 real(dp), intent(in) :: cc(ld,ld)
 type(piecewise_t), intent(in) :: pin(nn)

!Input/output variables
 type(piecewise_t), intent(inout) :: pout(nn)

!Local variables
 integer :: ii,jj

 do jj=1,nn
   pout(jj)%vals(:,:)=0.0d0
   do ii=1,nn
     pout(jj)%vals(:,:)=pout(jj)%vals(:,:) + cc(ii,jj)*pin(ii)%vals(:,:)
   end do
 end do
 return
 end subroutine pb_combine

 subroutine pb_sort_abs(ee,pw,nn)

!! the coefficients ee and their functions pw sorted by decreasing |ee|
!! (a bubble sort, of at most 100 passes)
!!
!! @param ee  the coefficients
!! @param pw  their functions
!! @param nn  how many

 use constants_m, only: dp
 use pieces_m, only: piecewise_t

!Input variables
 integer, intent(in) :: nn

!Input/output variables
 real(dp), intent(inout) :: ee(nn)
 type(piecewise_t), intent(inout) :: pw(nn)

!Local variables
 type(piecewise_t) :: pwswap
 real(dp) :: tt
 integer :: ii,jj
 logical :: sorted

 do ii=1,100
   sorted=.true.
   do jj=2,nn
     if(abs(ee(jj))>abs(ee(jj-1))) then
       tt=ee(jj-1)
       pwswap=pw(jj-1)
       ee(jj-1)=ee(jj)
       pw(jj-1)=pw(jj)
       ee(jj)=tt
       pw(jj)=pwswap
       sorted=.false.
     end if
   end do
   if(sorted) exit
 end do
 return
 end subroutine pb_sort_abs

 subroutine pb_sign_at_peak(rr,mmax,pw,vkb)

!! a projector on the mesh, vkb, from its pieces pw, both with the sign
!! (physically irrelevant) that makes it positive at its peak
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param pw  the projector's pieces
!! @param vkb  the projector on the mesh (output)

 use constants_m, only: dp
 use pieces_m, only: piecewise_t, pw_node_values

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax)

!Input/output variables
 type(piecewise_t), intent(inout) :: pw

!Output variables
 real(dp), intent(out) :: vkb(mmax)

!Local variables
 real(dp) :: apk
 integer :: ii,ipk

 call pw_node_values(rr,mmax,pw,vkb)
 apk=0.0d0
 ipk=1
 do ii=1,mmax
   if(abs(vkb(ii))>apk) then
     apk=abs(vkb(ii))
     ipk=ii
   end if
 end do
 if(vkb(ipk)<0.0d0) then
   vkb(:)=-vkb(:)
   pw%vals(:,:)=-pw%vals(:,:)
 end if
 return
 end subroutine pb_sign_at_peak

 end module projbasis_m
