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

 module run_vkb_m

 implicit none
 private
 public :: run_vkb

 contains

 subroutine run_vkb(lmax,lloc,lpopt,dvloc0,rc,irc,nproj,rr,mmax,pswf, &
&                   pswfx,chix,vfull,vp,evkb,vkb,nlim,dlta,vrloc,irel,nkap, &
&                   vkbpw,vlocpw,bbeig,iout,ifail)

!! computes Vanderbilt / Kleinman-Bylander non-local potentials. The
!! projectors are functions of pieces (pieces_m, doc/smoothness_plan.md
!! R-6, R-8, R-9): a projector of channel l is (eig - T - V_loc) psi inside
!! rc with a taper f Delta/2 on the last 5 % (D-4), zero beyond, and V_loc
!! itself is a polynomial (lloc=4) or a tapered semi-local potential
!! (lloc<4) inside its radius and the all-electron potential beyond, so a
!! projector has a kink at 0.95 rc, at rc, and where V_loc is matched; each
!! smooth piece is kept with its extension and every integral is taken per
!! piece, so the projectors and their coefficients are smooth in the radii
!!
!! @param lmax  maximum angular momentum
!! @param lloc  l for local potential (lloc==4 => use linear combination)
!! @param lpopt  choice of polynomial for lloc==4
!! @param dvloc0  amplitude at rr==0 to be smoothely added for lloc==4
!! @param rc  core radii (bohr)
!! @param irc  last mesh points at or inside the core radii
!! @param nproj  number of projectors for each l
!! @param rr  log radial grid
!! @param mmax  size of radial grid
!! @param pswf  pseudo wave functions
!! @param pswfx  their inside pieces (the Bessel sums) on the nodes to past
!!               rc (pspot)
!! @param chix  (eig - T) pswfx, the projectors' seeds, the same
!! @param vfull  all-electron potential
!! @param vp  semi-local pseudopotentials for first projectors
!!            (vp(:,5) is local potential if lloc=4)
!! @param vkb  VKB projectors on the mesh (output)
!! @param evkb  coefficients of VKB projectors (output)
!! @param nlim  index of maximum rc including that of vlocal (output)
!! @param dlta  the relativistic term Delta of each all-electron wave about
!!              its taper window (ae_delta)
!! @param vrloc  Delta/(2u) of the first wave of channel lloc about its
!!               window, for the taper of the local potential (lloc<4)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!!              (for pswf, vp, vkb, evkb and dlta the last index is 1 for
!!              j=l+1/2, 2 for j=l-1/2 when nkap=2)
!! @param vkbpw  the projectors as functions of pieces (output)
!! @param vlocpw  the screened local potential as a function of pieces
!!                (output): its inside pieces and the all-electron
!!                potential beyond rcloc, the last bound the mesh end
!! @param bbeig  the eigenvalues of each channel's symmetrized B matrix,
!!               ascending (B itself for one projector; 0 where there is
!!               none) (output)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxl, mxlin, rel_window, mxprj
 use kappa_m, only: kap_count, kap_of
 use projbasis_m, only: pb_half_powers, pb_congruence, pb_combine, &
&                       pb_sort_abs, pb_sign_at_peak
 use pieces_m, only: piecewise_t, pw_maxseg, pw_reach, pw_alloc, pw_free, &
&                    pw_integral, pw_product_integral, pw_node_values
 use vploc_m, only: vploc_t, vploc_fit, vploc_mesh, vploc_eval
 use lapack_m, only: dsyev

!Input variables
 integer, intent(in) :: lmax,lloc,lpopt,mmax,irel,nkap
 integer, intent(in) :: irc(mxlin+1),nproj(mxlin+1)
 real(dp), intent(in) :: rr(mmax),pswf(mmax,mxprj,mxl+1,nkap),vfull(mmax)
 real(dp), intent(in) :: pswfx(mmax,mxprj,mxl+1,nkap)
 real(dp), intent(in) :: chix(mmax,mxprj,mxl+1,nkap)
 real(dp), intent(in) :: dlta(mmax,mxprj,mxlin+1,nkap),vrloc(mmax,nkap)
 real(dp), intent(in) :: dvloc0,rc(mxlin+1)
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: vp(mmax,mxl+2,nkap)
 real(dp), intent(inout) :: vkb(mmax,mxprj,mxl+1,nkap)

!Output variables
 real(dp), intent(out) :: evkb(mxprj,mxl+1,nkap)
 integer, intent(out) :: nlim
 type(piecewise_t), intent(out) :: vkbpw(mxprj,mxl+1,nkap),vlocpw
 real(dp), intent(out) :: bbeig(mxprj,mxl+1,nkap)
 integer, intent(out) :: ifail

!Local variables
 integer :: lwork
 integer :: ii,jj,ikap,ll,l1,info,np,kap,mkap,nseg,iseg,nb,nbrk
 real(dp) :: sn,tt,xx,rcloc,rm,wj,bnd(0:pw_maxseg),brk(3)
 real(dp) :: bb(mxprj,mxprj),bbev(mxprj),bbi(mxprj,mxprj)
 real(dp) :: sovl(mxprj,mxprj),sovlev(mxprj),smhalf(mxprj,mxprj)
 real(dp) :: sphalf(mxprj,mxprj)
 real(dp) :: bbist(mxprj,mxprj),bbistev(mxprj)
 real(dp), allocatable :: vloc(:),vlin0(:),vlin1(:),vpiece(:),taper(:),work(:)
 type(piecewise_t) :: pw(mxprj),pwt(mxprj)
 type(vploc_t) :: cfloc

 ifail=0
 lwork=3*mxprj
 allocate(vloc(mmax),vlin0(mmax),vlin1(mmax),vpiece(mmax),taper(mmax))
 allocate(work(lwork))

! the local potential: its matching radius rcloc, its inside piece as
! smooth functions on the mesh (vlin0 below 0.95 rcloc, vlin1 on the taper
! window; the same polynomial for lloc=4), and its values on the nodes
! (vloc, and vp of the channel)
 if(lloc==4) then
   rcloc=rc(5)
   call vploc_fit(rr,vfull,mmax,rcloc,dvloc0,lpopt,cfloc)
   call vploc_mesh(cfloc,rr,vfull,mmax,vp(:,5,1))
   if(nkap==2) vp(:,5,2)=vp(:,5,1)
   vloc(:)=vp(:,5,1)
   do ii=1,mmax
     vlin0(ii)=vploc_eval(cfloc,rr(ii))
   end do
   vlin1(:)=vlin0(:)
 else
! the semi-local potential (eig - T) psi/psi of the channel's first wave,
! for nkap=2 the "scalar-relativistic" weighted average of the two j (there
! is then no SO for this l), with the taper f Delta/(2u) of S-1 as the
! polynomial of the window
   l1=lloc+1
   rcloc=rc(l1)
   nb=min(mmax,irc(l1)+pw_reach+16)
   mkap=kap_count(nkap,l1-1)
   vlin0(:)=0.0d0
   vlin1(:)=0.0d0
   do ikap=1,mkap
     if(mkap==1) then
       wj=1.0d0
     else
       wj=dble(l1-ikap+1)/dble(2*l1-1)
     end if
     do ii=1,nb
       if(pswfx(ii,1,l1,ikap)==0.0d0) then
         write(iout,'(2a,i2,a,f10.6)') 'run_vkb: ERROR node of the local', &
&              ' channel''s pseudo wave within reach of rc, l=',lloc, &
&              ' r=',rr(ii)
         ifail=1
         return
       end if
       tt=chix(ii,1,l1,ikap)/pswfx(ii,1,l1,ikap)
       xx=20.0d0*((rr(ii)/rcloc)-1.0d0)
       vlin0(ii)=vlin0(ii)+wj*tt
       vlin1(ii)=vlin1(ii)+wj*(tt-(1.0d0-xx**2)**2*vrloc(ii,ikap))
     end do
   end do
   do ii=1,mmax
     if(ii>irc(l1)) then
       vloc(ii)=vfull(ii)
     else if(rr(ii)<=rel_window*rcloc) then
       vloc(ii)=vlin0(ii)
     else
       vloc(ii)=vlin1(ii)
     end if
   end do
   vp(:,l1,1)=vloc(:)
   if(nkap==2) vp(:,l1,2)=vloc(:)
 end if

! the local potential's pieces, with the node values of vloc
 if(lloc==4) then
   bnd(0:2)=[0.0d0,rcloc,rr(mmax)]
   call pw_alloc(vlocpw,mmax,2,bnd(0:2))
   vlocpw%vals(:,1)=vlin0(:)
   vlocpw%vals(:,2)=vfull(:)
 else
   bnd(0:3)=[0.0d0,rel_window*rcloc,rcloc,rr(mmax)]
   call pw_alloc(vlocpw,mmax,3,bnd(0:3))
   vlocpw%vals(:,1)=vlin0(:)
   vlocpw%vals(:,2)=vlin1(:)
   vlocpw%vals(:,3)=vfull(:)
 end if

! Vanderbilt / Kleinman-Bylander projector construction
 nlim=irc(lloc+1)
 bbeig(:,:,:)=0.0d0
 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   ll=l1-1
   mkap=kap_count(nkap,l1-1)
   nlim=max(nlim,irc(l1))

! the segments of this channel's projectors: the kinks inside (0, rc) at
! 0.95 rc (the taper), and at the local potential's rcloc and 0.95 rcloc
! (lloc<4) if they fall inside, in order, then rc
   nbrk=1
   brk(1)=rel_window*rc(l1)
   if(rcloc<rc(l1)) then
     nbrk=nbrk+1
     brk(nbrk)=rcloc
   end if
   if(lloc<4 .and. rel_window*rcloc<rc(l1) .and. rel_window*rcloc/=brk(1)) then
     nbrk=nbrk+1
     brk(nbrk)=rel_window*rcloc
   end if
   do ii=1,nbrk-1
     do jj=ii+1,nbrk
       if(brk(jj)<brk(ii)) then
         tt=brk(ii)
         brk(ii)=brk(jj)
         brk(jj)=tt
       end if
     end do
   end do
   nseg=nbrk+1
   bnd(0)=0.0d0
   bnd(1:nbrk)=brk(1:nbrk)
   bnd(nseg)=rc(l1)

! loop on J = ll +/- 1/2
   do ikap=1,mkap
     kap=kap_of(ll,ikap)
     np=nproj(l1)

! the projectors' pieces: on each segment chi = (eig - T - V_loc) psi less
! the taper f Delta/2 on the window, with V_loc's piece of the segment
! and every term as its smooth extension
     do jj=1,np
       call pw_alloc(pw(jj),mmax,nseg,bnd(0:nseg))
       do iseg=1,nseg
         rm=0.5d0*(bnd(iseg-1)+bnd(iseg))
         if(rm>rcloc) then
           vpiece(:)=vfull(:)
         else if(lloc<4 .and. rm>rel_window*rcloc) then
           vpiece(:)=vlin1(:)
         else
           vpiece(:)=vlin0(:)
         end if
         taper(:)=0.0d0
         if(rm>rel_window*rc(l1)) then
           do ii=1,mmax
             xx=20.0d0*((rr(ii)/rc(l1))-1.0d0)
             taper(ii)=0.5d0*(1.0d0-xx**2)**2*dlta(ii,jj,l1,ikap)
           end do
         end if
         pw(jj)%vals(:,iseg)=chix(:,jj,l1,ikap)-vpiece(:)*pswfx(:,jj,l1,ikap) &
&                            -taper(:)
       end do
     end do

! Vanderbilt B-matrix construction
     bb(:,:)=0.0d0
     call pw_integral(rr,mmax,2*l1,pswfx(:,1,l1,ikap),pw(1),bb(1,1))
     evkb(1,l1,ikap)=1.0d0/bb(1,1)
     bbeig(1,l1,ikap)=bb(1,1)
! this is Kleinman-Bylander result if nproj==1

     if(nproj(l1)>=2) then

       do jj=1,np
         do ii=1,np
           call pw_integral(rr,mmax,2*l1,pswfx(:,ii,l1,ikap),pw(jj),bb(ii,jj))
         end do
       end do

! symmetrize exactly
       if(irel==2) then
         write(iout,'(/a,i3,a,i3)') 'B matrix Hermiticity error, ll=',l1-1, &
&              '  kap=',kap
       else
         write(iout,'(/a,i4)') 'B matrix Hermiticity error, ll=',l1-1
       end if
       do jj=2,np
         do ii=1,jj-1
           if(irel==2) then
             write(iout,'(2i4,1p,d14.4)') ii,jj,bb(ii,jj)-bb(jj,ii)
           else
             write(iout,'(2i4,1p,e14.4)') ii,jj,bb(ii,jj)-bb(jj,ii)
           end if
           tt=0.5d0*(bb(ii,jj)+bb(jj,ii))
           bb(ii,jj)=tt
           bb(jj,ii)=tt
         end do
       end do

! find eigenvalues and eigenvectors of the B matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

       call dsyev( 'V', 'U', np, bb, mxprj, bbev, work, lwork, info )
       if(info /= 0) then
         write(iout,'(a,i4)') 'run_vkb: B matrix eigenvalue ERROR, info=',info
         ifail=1
         return
       end if
       bbeig(1:np,l1,ikap)=bbev(1:np)

! take linear combinations to form diagonal projectors

       do jj=1,np
         pwt(jj)=pw(jj)
       end do

       do jj=1,np
         pw(jj)%vals(:,:)=0.0d0
         evkb(jj,l1,ikap)=1.0d0/bbev(jj)
         do ii=1,np
           pw(jj)%vals(:,:)=pw(jj)%vals(:,:)+bb(ii,jj)*pwt(ii)%vals(:,:)
         end do
       end do

     end if !nproj>=2

! normalize projectors
     do jj=1,np
       call pw_product_integral(rr,mmax,2*l1,pw(jj),pw(jj),sn)
       call pw_node_values(rr,mmax,pw(jj),vkb(:,jj,l1,ikap))
       if(vkb(irc(l1)-10,jj,l1,ikap)>=0.0d0) then
         pw(jj)%vals(:,:)=pw(jj)%vals(:,:)/sqrt(sn)
       else
         pw(jj)%vals(:,:)=-pw(jj)%vals(:,:)/sqrt(sn)
       end if
       evkb(jj,l1,ikap)=sn*evkb(jj,l1,ikap)
     end do


! calculate overlap
     if(nproj(l1)>=2) then


       bbi(:,:)=0.0d0
       do jj=1,np
         sovl(jj,jj)=1.0d0
         bbi(jj,jj)=evkb(jj,l1,ikap)
         do ii=jj+1,np
           call pw_product_integral(rr,mmax,2*l1,pw(ii),pw(jj),sn)
           sovl(ii,jj)=sn
           sovl(jj,ii)=sn
         end do
       end do

! find eigenvalues and eigenvectors of the overlap matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

       call dsyev( 'V', 'U', np, sovl, mxprj, sovlev, work, lwork, info )
       if(info /= 0) then
         write(iout,'(a,i4)') 'run_vkb: sovl matrix eigenvalue ERROR, info=', &
&              info
         ifail=1
         return
       end if

! construct S^(-1/2) AND s^(1/2)

       call pb_half_powers(sovl,sovlev,np,mxprj,sphalf,smhalf)

! construct B^(-1)* = S^(1/2)^T B^(-1) S^(1/2)

       call pb_congruence(bbi,sphalf,np,mxprj,bbist)

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

       call dsyev( 'V', 'U', np , bbist, mxprj, bbistev, work, lwork, info )
       if(info /= 0) then
         write(iout,'(a,i4)') &
&              'run_vkb: B^(-1)* matrix eigenvalue ERROR, info=',info
         ifail=1
         return
       end if


! take linear combinations to form orthonormal basis functions

       do jj=1,np
         pwt(jj)=pw(jj)
       end do
       call pb_combine(smhalf,mxprj,pwt,np,pw)


! take linear combinations to form orthonormal projectors

       do jj=1,np
         pwt(jj)=pw(jj)
         evkb(jj,l1,ikap)=bbistev(jj)
       end do
       call pb_combine(bbist,mxprj,pwt,np,pw)

! re-order if necessary so first projector has largest magnitude coefficient
! which is consistent with relativistic sr_so_r.f90 output

       call pb_sort_abs(evkb(1:np,l1,ikap),pw(1:np),np)

       write(iout,'(/a,1p,5e12.4)') '  Orthonormal projector coefficients',&
&          (evkb(jj,l1,ikap),jj=1,np)

     end if !nproj(l1)>=2

! the projectors on the mesh; their sign (physically irrelevant) positive
! at their peak (needed for compaisons apparently)

     do jj=1,np
       call pb_sign_at_peak(rr,mmax,pw(jj),vkb(:,jj,l1,ikap))
       vkbpw(jj,l1,ikap)=pw(jj)
       call pw_free(pw(jj))
       call pw_free(pwt(jj))
     end do

     if(nproj(l1)<mxprj) then
       do jj=nproj(l1)+1, mxprj
         evkb(jj,l1,ikap)=0.0d0
         vkb(:,jj,l1,ikap)=0.0d0
       end do
     end if
   end do !ikap
 end do !l1

 deallocate(vloc,vlin0,vlin1,vpiece,taper,work)

 return
 end subroutine run_vkb

 end module run_vkb_m
