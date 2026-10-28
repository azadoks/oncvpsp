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

 module sr_so_r_m

 implicit none
 private
 public :: sr_so_r

 contains

 subroutine sr_so_r(lmax,nproj,rr,mmax,evkb,vkbpw, &
&                       vsr,esr,vso,eso,vsrpw,vsopw,iout,ifail)

!! reformulates non-local potentials based on j = l +/- 1/2 to scalar-
!! relativistic and L dot S projectors
!! uses relationship <L dot S> = (J^2 - L^2 - S^2)/2
!! so L dot S = +/- l/2 for j = l +/- 1/2
!!
!! The projectors are functions of pieces (pieces_m): their overlaps are
!! taken per piece, and the combinations act on the pieces, so the
!! results are smooth in the core radii (doc/smoothness_plan.md, R-10)
!!
!! @param lmax  maximum angular momentum
!! @param nproj  number of projectors for each l
!! @param rr  log radial grid
!! @param mmax  size of radial grid
!! @param vkbpw  vkb projectors for each l and j=l+/-1/2, as functions of
!!               pieces
!! @param evkb  coefficients of VKB projectors for each l and j=l+/-1/2
!! @param vsr  normalized scalar-relativistic projectors on the mesh
!! @param esr  energy coefficients of vsr
!! @param vso  normalized spin-orbit projectors on the mesh
!! @param eso  energy coefficients of vso, to be multiplied by <L dot S>
!! @param vsrpw  vsr as functions of pieces
!! @param vsopw  vso as functions of pieces
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxl, mxlin, mxkap, mxprj
 use projbasis_m, only: pb_half_powers, pb_congruence, pb_combine, &
&                       pb_sort_abs, pb_sign_at_peak
 use pieces_m, only: piecewise_t, pw_product_integral, pw_node_values, &
&                    pw_free
 use lapack_m, only: dsyev

!Input variables
 integer, intent(in) :: lmax,mmax
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,mxkap)
 type(piecewise_t), intent(in) :: vkbpw(mxprj,mxl+1,mxkap)
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: vsr(mmax,2*mxprj,mxl+1),esr(2*mxprj,mxl+1)
 real(dp), intent(inout) :: vso(mmax,2*mxprj,mxl+1),eso(2*mxprj,mxl+1)

!Output variables
 type(piecewise_t), intent(out) :: vsrpw(2*mxprj,mxl+1)
 type(piecewise_t), intent(out) :: vsopw(2*mxprj,mxl+1)
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,jj,ik1,ik2,ip1,ip2,ll,l1,info,nn
 real(dp) :: sovl(2*mxprj,2*mxprj),sovlev(2*mxprj)
 real(dp) :: ascl(2*mxprj,2*mxprj),aso(2*mxprj,2*mxprj)
 real(dp) :: asclst(2*mxprj,2*mxprj),wsclst(2*mxprj)
 real(dp) :: asost(2*mxprj,2*mxprj),wsost(2*mxprj)
 real(dp) :: sphalf(2*mxprj,2*mxprj),smhalf(2*mxprj,2*mxprj)
 real(dp) :: fscl(mxprj),fso(mxprj),work(10*mxprj)
 type(piecewise_t) :: pwt(2*mxprj),pwst(2*mxprj)

 ifail=0
 do l1=1,lmax+1
   ll=l1-1

   if(ll==0) then
     vsr(:,:,l1)=0.0d0
     vso(:,:,l1)=0.0d0
     esr(:,l1)=0.0d0
     eso(:,l1)=0.0d0
     if(nproj(l1)>=1) then
       do ii=1,nproj(l1)
         vsrpw(ii,l1)=vkbpw(ii,l1,1)
         vsr(:,ii,l1)=0.0d0
         call pw_node_values(rr,mmax,vsrpw(ii,l1),vsr(:,ii,l1))
         esr(ii,l1)=evkb(ii,l1,1)
       end do
     end if
     cycle
   end if

   nn=2*nproj(l1)

   fscl(1)=(ll+1)/dble(2*ll+1)
   fscl(2)=ll/dble(2*ll+1)
   fso(1)=2/dble(2*ll+1)
   fso(2)=-2/dble(2*ll+1)

! construct overlap matrix and diagonal energy matrices

   sovl(:,:)=0.0d0
   ascl(:,:)=0.0d0
   aso(:,:)=0.0d0

   do ik1=1,2
     do ip1=1,nproj(l1)
       ii=ip1+(ik1-1)*nproj(l1)

       ascl(ii,ii)=fscl(ik1)*evkb(ip1,l1,ik1)
       aso(ii,ii)=fso(ik1)*evkb(ip1,l1,ik1)

       pwt(ii)=vkbpw(ip1,l1,ik1)

       do ik2=1,2
         do ip2=1,nproj(l1)
           jj=ip2+(ik2-1)*nproj(l1)

           call pw_product_integral(rr,mmax,2*l1,vkbpw(ip1,l1,ik1), &
&                                   vkbpw(ip2,l1,ik2),sovl(ii,jj))

         end do
       end do
     end do
   end do

   call dsyev( 'V', 'U', nn, sovl, 2*mxprj, sovlev, work, 10*mxprj, info )

   if(info /= 0) then
     write(iout,'(a,i4)') 'sr_so_r: S matrix eigenvalue ERROR, info=',info
     ifail=1
     return
   end if

! construct S^(-1/2) AND s^(1/2)

   call pb_half_powers(sovl,sovlev,nn,2*mxprj,sphalf,smhalf)

! take linear combinations to form orthonormal basis functions

   do jj=1,nn
     pwst(jj)=pwt(1)
   end do
   call pb_combine(smhalf,2*mxprj,pwt,nn,pwst)

! construct A^(-1)* = S^(1/2)^T A^(-1) S^(1/2)

   call pb_congruence(ascl,sphalf,nn,2*mxprj,asclst)
   call pb_congruence(aso,sphalf,nn,2*mxprj,asost)

! find eigenvalues and eigenvectors of the A* matrices

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   call dsyev( 'V', 'U', nn, asclst, 2*mxprj, wsclst, work, 10*mxprj, info )

   if(info /= 0) then
     write(iout,'(a,i4)') 'sr_so_r: A* matrix eigenvalue ERROR, info=',info
     ifail=1
     return
   end if

   call dsyev( 'V', 'U', nn,  asost, 2*mxprj,  wsost, work, 10*mxprj, info )


   if(info /= 0) then
     write(iout,'(a,i4)') 'sr_so_r: A* matrix eigenvalue ERROR, info=',info
     ifail=1
     return
   end if

! take linear combinations to form orthonormal projectors

   vsr(:,:,l1)=0.0d0
   vso(:,:,l1)=0.0d0
   esr(:,l1)=0.0d0
   eso(:,l1)=0.0d0

   do ii=1,nn
     esr(ii,l1)=wsclst(ii)
     eso(ii,l1)=  wsost(ii)
     vsrpw(ii,l1)=pwt(1)
     vsopw(ii,l1)=pwt(1)
   end do
   call pb_combine(asclst,2*mxprj,pwst,nn,vsrpw(1:nn,l1))
   call pb_combine(asost,2*mxprj,pwst,nn,vsopw(1:nn,l1))

! bubble-sort on coefficient magnitudes for scalar and then s-o
! (Yes, I know bubble-sort is the least-efficient sorting algorithm.)

   call pb_sort_abs(esr(1:nn,l1),vsrpw(1:nn,l1),nn)
   call pb_sort_abs(eso(1:nn,l1),vsopw(1:nn,l1),nn)

   write(iout,'(/a,i2)') &
&       ' Orthonormal scalar projector coefficients, l = ',ll
   write(iout,'(1p,6e12.4)') (esr(jj,l1),jj=1,nn)
   write(iout,'(/a,i2)') &
&       ' Orthonormal spin-orbit projector coefficients, l = ',ll
   write(iout,'(1p,6e12.4)') (eso(jj,l1),jj=1,nn)

! the projectors on the mesh, with their sign (physically irrelevant)
! positive at their peak (needed for compaisons apparently)

   do jj=1,nn
     call pb_sign_at_peak(rr,mmax,vsopw(jj,l1),vso(:,jj,l1))
     call pb_sign_at_peak(rr,mmax,vsrpw(jj,l1),vsr(:,jj,l1))
   end do

   do ii=1,nn
     call pw_free(pwt(ii))
     call pw_free(pwst(ii))
   end do

 end do ! l1

 return
 end subroutine sr_so_r

 end module sr_so_r_m
