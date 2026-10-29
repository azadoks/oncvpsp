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
 use constants_m, only: dp
 use pieces_m, only: piecewise_t
 implicit none
 private
 public :: run_vkb, local_potential, channel_vkb, local_t

!> the local potential (local_potential): its pieces inside its radius
!> and on the mesh, for the projectors of every channel (channel_vkb)
 type :: local_t
!> its matching radius (bohr)
   real(dp) :: rcloc=0.0d0
!> its inside piece below 0.95 rcloc, as a smooth function on the mesh
   real(dp), allocatable :: vlin0(:)
!> its piece on the taper window, the same as vlin0 for lloc=4
   real(dp), allocatable :: vlin1(:)
!> the screened local potential on the mesh (Ha)
   real(dp), allocatable :: vloc(:)
!> its pieces: the inside pieces and the all-electron potential beyond
!> rcloc, the last bound the mesh end
   type(piecewise_t) :: vlocpw
 end type local_t

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
!! piece, so the projectors and their coefficients are smooth in the radii.
!! The local potential (local_potential) and each channel's projectors
!! (channel_vkb) are stages of their own
!!
!! @param lmax  maximum angular momentum
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param lpopt  choice of polynomial for lloc==4
!! @param dvloc0  amplitude at rr==0 to be smoothely added for lloc==4
!! @param rc  core radius (bohr) of each l
!! @param irc  last mesh point at or inside rc, of each l
!! @param nproj  number of projectors for each l
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param pswf  pseudo wave functions
!! @param pswfx  their inside pieces (the Bessel sums) on the nodes to past
!!               rc (pspot)
!! @param chix  (eig - T) pswfx, the projectors' seeds, the same
!! @param vfull  all-electron potential (Ha)
!! @param vp  semi-local pseudopotentials for first projectors
!!            (vp(:,5) is local potential if lloc=4)
!! @param evkb  coefficients of the projectors (Ha) (output)
!! @param vkb  VKB projectors on the mesh (output)
!! @param nlim  index of maximum rc including that of vlocal (output)
!! @param dlta  the relativistic term Delta of each all-electron wave about
!!              its taper window (ae_delta)
!! @param vrloc  Delta/(2u) of the first wave of channel lloc about its
!!               window, for the taper of the local potential (lloc<4)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!!              (for pswf, vp, vkb, evkb and dlta the last index is 1 for
!!              j=l+1/2, 2 for j=l-1/2 when nkap=2)
!! @param vkbpw  the projectors as pieces (output)
!! @param vlocpw  the screened local potential as pieces (output): its
!!                inside pieces and the all-electron potential beyond
!!                rcloc, the last bound the mesh end
!! @param bbeig  the eigenvalues of each channel's symmetrized B matrix,
!!               ascending (B itself for one projector; 0 where there is
!!               none) (output)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: mxl, mxlin, mxprj

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
 type(local_t) :: lp
 integer :: l1,lc

! the local potential, from the local channel's first wave for lloc<4
 lc=min(lloc+1,mxl+1)
 call local_potential(lloc,lpopt,dvloc0,rc,irc,rr,mmax,vfull,nkap, &
&                     pswfx(:,:,lc,:),chix(:,:,lc,:),vrloc,lp,iout,ifail)
 if(ifail/=0) return
 if(lloc==4) then
   vp(:,5,1)=lp%vloc(:)
   if(nkap==2) vp(:,5,2)=vp(:,5,1)
 else
   vp(:,lloc+1,1)=lp%vloc(:)
   if(nkap==2) vp(:,lloc+1,2)=lp%vloc(:)
 end if
 vlocpw=lp%vlocpw

! Vanderbilt / Kleinman-Bylander projector construction
 nlim=irc(lloc+1)
 bbeig(:,:,:)=0.0d0
 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   nlim=max(nlim,irc(l1))
   call channel_vkb(l1,lloc,lp,rc(l1),irc(l1),nproj(l1),rr,mmax,vfull, &
&                   pswfx(:,:,l1,:),chix(:,:,l1,:),dlta(:,:,l1,:),irel,nkap, &
&                   evkb(:,l1,:),vkb(:,:,l1,:),vkbpw(:,l1,:),bbeig(:,l1,:), &
&                   iout,ifail)
   if(ifail/=0) return
 end do

 return
 end subroutine run_vkb

 subroutine local_potential(lloc,lpopt,dvloc0,rc,irc,rr,mmax,vfull,nkap, &
&                           pswfx,chix,vrloc,lp,iout,ifail)

!! the local potential: a polynomial (lloc=4, vploc_m) inside rc(5), or
!! the semi-local potential (eig - T) psi/psi of the first wave of channel
!! lloc, for nkap=2 the "scalar-relativistic" weighted average of the two
!! j, with the taper f Delta/(2u) on its window; the all-electron
!! potential beyond
!!
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param lpopt  choice of polynomial for lloc==4
!! @param dvloc0  amplitude at rr==0 to be smoothely added for lloc==4
!! @param rc  core radius (bohr) of each l
!! @param irc  last mesh point at or inside rc, of each l
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param vfull  all-electron potential (Ha)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param pswfx  the local channel's pseudo-waves' inside pieces (lloc<4)
!! @param chix  (eig - T) pswfx, the same
!! @param vrloc  Delta/(2u) of the first wave of channel lloc about its
!!               window, for the taper of the local potential (lloc<4)
!! @param lp  the local potential
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_construction, failed, msg_vkb_local_node
 use constants_m, only: mxlin, rel_window, mxprj
 use kappa_m, only: kap_count
 use pieces_m, only: pw_reach, pw_maxseg, pw_alloc
 use vploc_m, only: vploc_t, vploc_fit, vploc_mesh, vploc_eval

!Input variables
 integer, intent(in) :: lloc,lpopt,mmax,nkap
 integer, intent(in) :: irc(mxlin+1)
 real(dp), intent(in) :: dvloc0,rc(mxlin+1),rr(mmax),vfull(mmax)
 real(dp), intent(in) :: pswfx(mmax,mxprj,nkap),chix(mmax,mxprj,nkap)
 real(dp), intent(in) :: vrloc(mmax,nkap)
 integer, intent(in) :: iout

!Output variables
 type(local_t), intent(out) :: lp
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,ikap,l1,mkap,nb
 real(dp) :: tt,xx,wj,bnd(0:pw_maxseg)
 type(vploc_t) :: cfloc

 ifail=0
 allocate(lp%vlin0(mmax),lp%vlin1(mmax),lp%vloc(mmax))
! the local potential: its matching radius rcloc, its inside piece as
! smooth functions on the mesh (vlin0 below 0.95 rcloc, vlin1 on the taper
! window; the same polynomial for lloc=4), and its values on the nodes
! (vloc, and vp of the channel)
 if(lloc==4) then
   lp%rcloc=rc(5)
   call vploc_fit(rr,vfull,mmax,lp%rcloc,dvloc0,lpopt,cfloc,ifail)
   if(ifail/=0) return
   call vploc_mesh(cfloc,rr,vfull,mmax,lp%vloc)
   do ii=1,mmax
     lp%vlin0(ii)=vploc_eval(cfloc,rr(ii))
   end do
   lp%vlin1(:)=lp%vlin0(:)
 else
! the semi-local potential (eig - T) psi/psi of the channel's first wave,
! for nkap=2 the "scalar-relativistic" weighted average of the two j (there
! is then no SO for this l), with the taper f Delta/(2u) of S-1 as the
! polynomial of the window
   l1=lloc+1
   lp%rcloc=rc(l1)
   nb=min(mmax,irc(l1)+pw_reach+16)
   mkap=kap_count(nkap,l1-1)
   lp%vlin0(:)=0.0d0
   lp%vlin1(:)=0.0d0
   do ikap=1,mkap
     if(mkap==1) then
       wj=1.0d0
     else
       wj=dble(l1-ikap+1)/dble(2*l1-1)
     end if
     do ii=1,nb
       if(pswfx(ii,1,ikap)==0.0d0) then
         write(iout,'(/2a,i0,a,f12.6)') 'run_vkb: ERROR node of the local', &
&              ' channel''s pseudo wave within reach of rc, l=',lloc, &
&              ' r=',rr(ii)
         ifail=failed(fail_construction,msg_vkb_local_node)
         return
       end if
       tt=chix(ii,1,ikap)/pswfx(ii,1,ikap)
       xx=20.0d0*((rr(ii)/lp%rcloc)-1.0d0)
       lp%vlin0(ii)=lp%vlin0(ii)+wj*tt
       lp%vlin1(ii)=lp%vlin1(ii)+wj*(tt-(1.0d0-xx**2)**2*vrloc(ii,ikap))
     end do
   end do
   do ii=1,mmax
     if(ii>irc(l1)) then
       lp%vloc(ii)=vfull(ii)
     else if(rr(ii)<=rel_window*lp%rcloc) then
       lp%vloc(ii)=lp%vlin0(ii)
     else
       lp%vloc(ii)=lp%vlin1(ii)
     end if
   end do
 end if

! the local potential's pieces, with the node values of vloc
 if(lloc==4) then
   bnd(0:2)=[0.0d0,lp%rcloc,rr(mmax)]
   call pw_alloc(lp%vlocpw,mmax,2,bnd(0:2),ifail)
   if(ifail/=0) return
   lp%vlocpw%vals(:,1)=lp%vlin0(:)
   lp%vlocpw%vals(:,2)=vfull(:)
 else
   bnd(0:3)=[0.0d0,rel_window*lp%rcloc,lp%rcloc,rr(mmax)]
   call pw_alloc(lp%vlocpw,mmax,3,bnd(0:3),ifail)
   if(ifail/=0) return
   lp%vlocpw%vals(:,1)=lp%vlin0(:)
   lp%vlocpw%vals(:,2)=lp%vlin1(:)
   lp%vlocpw%vals(:,3)=vfull(:)
 end if

 return
 end subroutine local_potential

 subroutine channel_vkb(l1,lloc,lp,rcl,ircl,nprj,rr,mmax,vfull,pswfx,chix, &
&                       dlta,irel,nkap,evkb,vkb,vkbpw,bbeig,iout,ifail)

!! the projectors of channel l = l1-1 (each kappa): on each piece
!! (eig - T - V_loc) psi less the taper f Delta/2, the B matrix, and for
!! two or more projectors their orthonormal combinations, normalized, with
!! their sign positive at their peak
!!
!! @param l1  l+1
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param lp  the local potential (local_potential)
!! @param rcl  core radius (bohr)
!! @param ircl  last mesh point at or inside rcl
!! @param nprj  its number of projectors
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param vfull  all-electron potential (Ha)
!! @param pswfx  the pseudo-waves' inside pieces on the nodes to past rc
!! @param chix  (eig - T) pswfx, the projectors' seeds, the same
!! @param dlta  the relativistic term Delta of each all-electron wave about
!!              its taper window (ae_delta)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param evkb  coefficients of the projectors (Ha)
!! @param vkb  the projectors on the mesh
!! @param vkbpw  the projectors as pieces
!! @param bbeig  the eigenvalues of the symmetrized B matrix, ascending (B
!!               itself for one projector; 0 where there is none)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_numerical, failed, msg_vkb_b_dsyev, &
&                    msg_vkb_binv_dsyev, msg_vkb_sovl_dsyev
 use constants_m, only: rel_window, mxprj
 use kappa_m, only: kap_count, kap_of
 use projbasis_m, only: pb_half_powers, pb_congruence, pb_combine, &
&                       pb_sort_abs, pb_sign_at_peak
 use pieces_m, only: pw_maxseg, pw_alloc, pw_free, pw_integral, &
&                    pw_product_integral, pw_node_values
 use lapack_m, only: dsyev

!Input variables
 integer, intent(in) :: l1,lloc,ircl,nprj,mmax,irel,nkap
 type(local_t), intent(in) :: lp
 real(dp), intent(in) :: rcl,rr(mmax),vfull(mmax)
 real(dp), intent(in) :: pswfx(mmax,mxprj,nkap),chix(mmax,mxprj,nkap)
 real(dp), intent(in) :: dlta(mmax,mxprj,nkap)
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: evkb(mxprj,nkap),vkb(mmax,mxprj,nkap)
 type(piecewise_t), intent(inout) :: vkbpw(mxprj,nkap)
 real(dp), intent(inout) :: bbeig(mxprj,nkap)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: lwork
 integer :: ii,jj,ikap,ll,info,kap,mkap,nseg,iseg,nbrk
 real(dp) :: sn,tt,xx,rm,bnd(0:pw_maxseg),brk(3)
 real(dp) :: bb(mxprj,mxprj),bbev(mxprj),bbi(mxprj,mxprj)
 real(dp) :: sovl(mxprj,mxprj),sovlev(mxprj),smhalf(mxprj,mxprj)
 real(dp) :: sphalf(mxprj,mxprj)
 real(dp) :: bbist(mxprj,mxprj),bbistev(mxprj)
 real(dp), allocatable :: vpiece(:),taper(:),work(:)
 type(piecewise_t) :: pw(mxprj),pwt(mxprj)

 ifail=0
 lwork=3*mxprj
 allocate(vpiece(mmax),taper(mmax),work(lwork))
 bbeig(:,:)=0.0d0
 ll=l1-1
 mkap=kap_count(nkap,l1-1)

! the segments of this channel's projectors: the kinks inside (0, rc) at
! 0.95 rc (the taper), and at the local potential's rcloc and 0.95 rcloc
! (lloc<4) if they fall inside, in order, then rc
 nbrk=1
 brk(1)=rel_window*rcl
 if(lp%rcloc<rcl) then
   nbrk=nbrk+1
   brk(nbrk)=lp%rcloc
 end if
 if(lloc<4 .and. rel_window*lp%rcloc<rcl .and. &
&   rel_window*lp%rcloc/=brk(1)) then
   nbrk=nbrk+1
   brk(nbrk)=rel_window*lp%rcloc
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
 bnd(nseg)=rcl

! loop on J = ll +/- 1/2
 do ikap=1,mkap
   kap=kap_of(ll,ikap)

! the projectors' pieces: on each segment chi = (eig - T - V_loc) psi less
! the taper f Delta/2 on the window, with V_loc's piece of the segment
! and every term as its smooth extension
   do jj=1,nprj
     call pw_alloc(pw(jj),mmax,nseg,bnd(0:nseg),ifail)
     if(ifail/=0) return
     do iseg=1,nseg
       rm=0.5d0*(bnd(iseg-1)+bnd(iseg))
       if(rm>lp%rcloc) then
         vpiece(:)=vfull(:)
       else if(lloc<4 .and. rm>rel_window*lp%rcloc) then
         vpiece(:)=lp%vlin1(:)
       else
         vpiece(:)=lp%vlin0(:)
       end if
       taper(:)=0.0d0
       if(rm>rel_window*rcl) then
         do ii=1,mmax
           xx=20.0d0*((rr(ii)/rcl)-1.0d0)
           taper(ii)=0.5d0*(1.0d0-xx**2)**2*dlta(ii,jj,ikap)
         end do
       end if
       pw(jj)%vals(:,iseg)=chix(:,jj,ikap)-vpiece(:)*pswfx(:,jj,ikap) &
&                            -taper(:)
     end do
   end do

! Vanderbilt B-matrix construction
   bb(:,:)=0.0d0
   call pw_integral(rr,mmax,2*l1,pswfx(:,1,ikap),pw(1),bb(1,1))
   evkb(1,ikap)=1.0d0/bb(1,1)
   bbeig(1,ikap)=bb(1,1)
! this is Kleinman-Bylander result if nproj==1

   if(nprj>=2) then

     do jj=1,nprj
       do ii=1,nprj
         call pw_integral(rr,mmax,2*l1,pswfx(:,ii,ikap),pw(jj),bb(ii,jj))
       end do
     end do

! symmetrize exactly
     if(irel==2) then
       write(iout,'(/a,i3,a,i3)') 'B matrix Hermiticity error, ll=',l1-1, &
&              '  kap=',kap
     else
       write(iout,'(/a,i4)') 'B matrix Hermiticity error, ll=',l1-1
     end if
     do jj=2,nprj
       do ii=1,jj-1
         if(irel==2) then
           write(iout,'(2i4,1p,e14.4)') ii,jj,bb(ii,jj)-bb(jj,ii)
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

     call dsyev( 'V', 'U', nprj, bb, mxprj, bbev, work, lwork, info )
     if(info /= 0) then
       write(iout,'(a,i4)') 'run_vkb: B matrix eigenvalue ERROR, info=',info
       ifail=failed(fail_numerical,msg_vkb_b_dsyev)
       return
     end if
     bbeig(1:nprj,ikap)=bbev(1:nprj)

! take linear combinations to form diagonal projectors

     do jj=1,nprj
       pwt(jj)=pw(jj)
     end do

     do jj=1,nprj
       pw(jj)%vals(:,:)=0.0d0
       evkb(jj,ikap)=1.0d0/bbev(jj)
       do ii=1,nprj
         pw(jj)%vals(:,:)=pw(jj)%vals(:,:)+bb(ii,jj)*pwt(ii)%vals(:,:)
       end do
     end do

   end if !nproj>=2

! normalize projectors
   do jj=1,nprj
     call pw_product_integral(rr,mmax,2*l1,pw(jj),pw(jj),sn,ifail)
     if(ifail/=0) return
     call pw_node_values(rr,mmax,pw(jj),vkb(:,jj,ikap))
     if(vkb(ircl-10,jj,ikap)>=0.0d0) then
       pw(jj)%vals(:,:)=pw(jj)%vals(:,:)/sqrt(sn)
     else
       pw(jj)%vals(:,:)=-pw(jj)%vals(:,:)/sqrt(sn)
     end if
     evkb(jj,ikap)=sn*evkb(jj,ikap)
   end do


! calculate overlap
   if(nprj>=2) then


     bbi(:,:)=0.0d0
     do jj=1,nprj
       sovl(jj,jj)=1.0d0
       bbi(jj,jj)=evkb(jj,ikap)
       do ii=jj+1,nprj
         call pw_product_integral(rr,mmax,2*l1,pw(ii),pw(jj),sn,ifail)
         if(ifail/=0) return
         sovl(ii,jj)=sn
         sovl(jj,ii)=sn
       end do
     end do

! find eigenvalues and eigenvectors of the overlap matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

     call dsyev( 'V', 'U', nprj, sovl, mxprj, sovlev, work, lwork, info )
     if(info /= 0) then
       write(iout,'(a,i4)') 'run_vkb: sovl matrix eigenvalue ERROR, info=', &
&              info
       ifail=failed(fail_numerical,msg_vkb_sovl_dsyev)
       return
     end if

! construct S^(-1/2) AND s^(1/2)

     call pb_half_powers(sovl,sovlev,nprj,mxprj,sphalf,smhalf)

! construct B^(-1)* = S^(1/2)^T B^(-1) S^(1/2)

     call pb_congruence(bbi,sphalf,nprj,mxprj,bbist)

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

     call dsyev( 'V', 'U', nprj , bbist, mxprj, bbistev, work, lwork, info )
     if(info /= 0) then
       write(iout,'(a,i4)') &
&              'run_vkb: B^(-1)* matrix eigenvalue ERROR, info=',info
       ifail=failed(fail_numerical,msg_vkb_binv_dsyev)
       return
     end if


! take linear combinations to form orthonormal basis functions

     do jj=1,nprj
       pwt(jj)=pw(jj)
     end do
     call pb_combine(smhalf,mxprj,pwt,nprj,pw)


! take linear combinations to form orthonormal projectors

     do jj=1,nprj
       pwt(jj)=pw(jj)
       evkb(jj,ikap)=bbistev(jj)
     end do
     call pb_combine(bbist,mxprj,pwt,nprj,pw)

! re-order if necessary so first projector has largest magnitude coefficient
! which is consistent with relativistic sr_so_r.f90 output

     call pb_sort_abs(evkb(1:nprj,ikap),pw(1:nprj),nprj)

     write(iout,'(/a,1p,5e12.4)') '  Orthonormal projector coefficients',&
&          (evkb(jj,ikap),jj=1,nprj)

   end if !nprj>=2

! the projectors on the mesh; their sign (physically irrelevant) positive
! at their peak (needed for compaisons apparently)

   do jj=1,nprj
     call pb_sign_at_peak(rr,mmax,pw(jj),vkb(:,jj,ikap))
     vkbpw(jj,ikap)=pw(jj)
     call pw_free(pw(jj))
     call pw_free(pwt(jj))
   end do

   if(nprj<mxprj) then
     do jj=nprj+1, mxprj
       evkb(jj,ikap)=0.0d0
       vkb(:,jj,ikap)=0.0d0
     end do
   end if
 end do !ikap

 deallocate(vpiece,taper,work)
 return
 end subroutine channel_vkb

 end module run_vkb_m
