!
! The ionic tail of the pseudopotentials, after PseudoPotentialGenerators.jl
! (src/unscreening.jl: ionic_tail, ionic_switch, ionic_tail_window,
! with_ionic_tail, ionic_screening).
!
! Beyond the largest core radius rcmax the screened potentials vp are V_AE,
! screened by the all-electron valence density n_v.  Unscreening them with
! the pseudo valence density rho, which differs from n_v there (the
! scalar- or fully relativistic reference is solved with a different
! equation than the pseudo-atom), leaves V_Hxc[n_v] - V_Hxc[rho] in vpuns:
! a long-ranged term that is not ionic, and that plane-wave codes see
! wherever they cut the table.  Here the screened potentials are switched
! over a window (ra,rb) to
!
!   V_ion + V_Hxc[rho + rhomod],
!   V_ion = -Z/r + V_H[n_c] + V_xc[n_c + n_v] - V_xc[n_v + rhomod],
!
! the frozen-core ionic potential (-zion/r plus the core penetration and,
! without a model core, the core-valence xc nonlinearity), so that vpuns
! is vp - V_Hxc[rho + rhomod] up to ra and V_ion beyond rb.  The pseudo
! valence states see the switched local potential, so rho is iterated
! with the tail to a fixed point: the reference configuration stays a
! fixed point of the pseudo-atom SCF.  The projectors vanish beyond rc
! and are unchanged.
!

! the frozen-core ionic potential V_ion of the all-electron atom
 subroutine ionic_tail(zz,zion,rhoc,rhotae,rhomod,rr,mmax,iexc,vion)
!zz  nuclear charge
!zion  ionic charge (zz minus the core occupation)
!rhoc  all-electron core charge
!rhotae  all-electron valence charge
!rhomod  model core charge
!vion  output V_ion
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 integer :: mmax,iexc
 real(dp) :: zz,zion
 real(dp) :: rhoc(mmax),rhotae(mmax),rhomod(mmax),rr(mmax),vion(mmax)
 real(dp) :: eeel,eexc
 real(dp), allocatable :: vt(:),vxct(:)

 allocate(vt(mmax),vxct(mmax))

! V_H[n_c], its charge fixed at the end of the mesh to zz-zion so that
! -Z/r + V_H[n_c] is -zion/r there
 call vout(0,rhoc,rhoc,vt,vxct,zz-zion,eeel,eexc,rr,mmax,iexc)
 vion(:)=-zz/rr(:)+vt(:)-vxct(:)

 call vout(1,rhotae,rhoc,vt,vxct,zion,eeel,eexc,rr,mmax,iexc)
 vion(:)=vion(:)+vxct(:)
 call vout(1,rhotae,rhomod,vt,vxct,zion,eeel,eexc,rr,mmax,iexc)
 vion(:)=vion(:)-vxct(:)

 deallocate(vt,vxct)
 return
 end subroutine ionic_tail

! the window (ra,rb) of the switch: ra is the largest core radius, beyond
! which the screened potentials are V_AE, and rb is where V_ion has become
! -zion/r to tol, but at least 2*ra (a switch over less than ra would be
! harder than the pseudization)
 subroutine ionic_window(rr,mmax,zion,vion,rcmax,tol,ra,rb)
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 integer :: mmax
 real(dp) :: rr(mmax),vion(mmax)
 real(dp) :: zion,rcmax,tol,ra,rb
 real(dp) :: rcore
 integer :: ii

 rcore=0.0d0
 do ii=mmax,1,-1
   if(abs(vion(ii)+zion/rr(ii))>=tol) then
     rcore=rr(min(ii+1,mmax))
     exit
   end if
 end do
 ra=rcmax
 rb=max(rcore,2.0d0*ra)
 return
 end subroutine ionic_window

! C-infinity switch from 1 (r <= ra) to 0 (r >= rb): f(1-x)/(f(1-x)+f(x)),
! f(t) = exp(-1/t), x = (r-ra)/(rb-ra); its transform decays faster than any
! power of q, so it adds no plane-wave tail of its own
 function ionic_switch(r,ra,rb)
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp) :: ionic_switch
 real(dp) :: r,ra,rb
 real(dp) :: x,f0,f1

 if(r<=ra) then
   ionic_switch=1.0d0
 else if(r>=rb) then
   ionic_switch=0.0d0
 else
   x=(r-ra)/(rb-ra)
   f0=exp(-1.0d0/(1.0d0-x))
   f1=exp(-1.0d0/x)
   ionic_switch=f0/(f0+f1)
 end if
 return
 end function ionic_switch

! the occupied pseudo valence states of the local potential vp(:,lloc+1,:)
! and the VKB projectors: rho, the state densities rhops, the orbitals
! uupsa and the eigenvalue sum eeig.  nk=1 for the scalar- and
! non-relativistic drivers, nk=2 (both kappa) for the fully relativistic one.
 subroutine ionic_valence(nk,nc,nv,la,fa,ea,lloc,nproj,rr,vp,vkb,evkb, &
&                         mmax,mxprj,rho,rhops,uupsa,eeig)
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 integer :: nk,nc,nv,lloc,mmax,mxprj
 integer :: la(30),nproj(6)
 real(dp) :: fa(30),ea(30,nk)
 real(dp) :: rr(mmax),vp(mmax,5,nk),vkb(mmax,mxprj,4,nk),evkb(mxprj,4,nk)
 real(dp) :: rho(mmax),rhops(mmax,nv),uupsa(mmax,nk,nv),eeig
 real(dp), allocatable :: uu(:),up(:)
 real(dp) :: et,emin,emax,fj,de
 integer :: kk,ll,l1,ikap,mkap,mch,ierr
 integer :: nodes(4)

 allocate(uu(mmax),up(mmax))
 if(nk==1) then
   de=0.25d0
 else
   de=0.10d0
 end if
 rho(:)=0.0d0
 rhops(:,:)=0.0d0
 eeig=0.0d0
 nodes(:)=0
 do kk=1,nv
   ll=la(nc+kk)
   l1=ll+1
   mkap=nk
   if(ll==0) mkap=1
   do ikap=1,mkap
     if(nk==1) then
       fj=1.0d0
     else if(ikap==1) then
       fj=dble(ll+1)/dble(2*ll+1)
     else
       fj=dble(ll)/dble(2*ll+1)
     end if
     et=ea(nc+kk,ikap)
     emax=(1.0d0-de)*et
     emin=(1.0d0+de)*et
     call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
&                  rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                  uu,up,mmax,mch)
     if(ierr/=0) then
       write(6,'(a,3i4)') 'ionic_valence: lschvkbb ERROR',ll+nodes(l1)+1,ll,ierr
       flush(6)
       stop
     end if
     uupsa(:,ikap,kk)=uu(:)
     rhops(:,kk)=rhops(:,kk)+fj*(uu(:)/rr(:))**2
     eeig=eeig+fa(nc+kk)*fj*et
   end do
   rho(:)=rho(:)+fa(nc+kk)*rhops(:,kk)
   nodes(l1)=nodes(l1)+1
 end do
 deallocate(uu,up)
 return
 end subroutine ionic_valence

! switch the screened potentials vp to V_ion + V_Hxc[rho + rhomod] over
! the ionic window, iterating the pseudo valence states (rho, rhops, uupsa,
! eeig) and their screening (vo, vxc, eeel, eexc) to a fixed point
 subroutine ionic_screening(nk,nc,nv,la,fa,ea,lmax,lloc,nproj,rr,vp,vkb, &
&                           evkb,mmax,mxprj,zz,zion,zval,iexc,rcmax,rhoc, &
&                           rhotae,rhomod,rho,rhops,uupsa,eeig,vo,vxc,eeel,eexc)
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), parameter :: tol=1.0d-10
 integer, parameter :: maxiter=50
 integer :: nk,nc,nv,lmax,lloc,mmax,mxprj,iexc
 integer :: la(30),nproj(6)
 real(dp) :: fa(30),ea(30,nk)
 real(dp) :: rr(mmax),vp(mmax,5,nk),vkb(mmax,mxprj,4,nk),evkb(mxprj,4,nk)
 real(dp) :: zz,zion,zval,rcmax,eeig,eeel,eexc
 real(dp) :: rhoc(mmax),rhotae(mmax),rhomod(mmax)
 real(dp) :: rho(mmax),rhops(mmax,nv),uupsa(mmax,nk,nv),vo(mmax),vxc(mmax)
 real(dp), allocatable :: vion(:),ww(:),vp0(:,:,:),vnew(:),vnext(:,:),vsw(:)
 real(dp) :: ra,rb,change,sn,sd
 real(dp), external :: ionic_switch
 integer :: ii,l1,ik,iter

 allocate(vion(mmax),ww(mmax),vp0(mmax,5,nk),vnew(mmax),vnext(mmax,nk),vsw(mmax))

 call ionic_tail(zz,zion,rhoc,rhotae,rhomod,rr,mmax,iexc,vion)
 call ionic_window(rr,mmax,zion,vion,rcmax,tol,ra,rb)
 do ii=1,mmax
   ww(ii)=ionic_switch(rr(ii),ra,rb)
 end do
 vp0(:,:,:)=vp(:,:,:)

 do ik=1,nk
   vp(:,lloc+1,ik)=ww(:)*vp0(:,lloc+1,ik)+(1.0d0-ww(:))*(vion(:)+vo(:))
 end do
 do iter=1,maxiter
   call ionic_valence(nk,nc,nv,la,fa,ea,lloc,nproj,rr,vp,vkb,evkb, &
&                     mmax,mxprj,rho,rhops,uupsa,eeig)
   call vout(1,rho,rhomod,vo,vxc,zval,eeel,eexc,rr,mmax,iexc)
! change of the local potential where the valence electrons are; once it
! is below tol, vp is kept as the potential whose states rho, eeig and
! (through rho) vo are, so that the pseudo-atom SCF of the reference
! configuration starts at its own fixed point
   sn=0.0d0; sd=0.0d0
   do ik=1,nk
     vnew(:)=ww(:)*vp0(:,lloc+1,ik)+(1.0d0-ww(:))*(vion(:)+vo(:))
     sn=sn+sum(rho(:)*rr(:)**3*(vnew(:)-vp(:,lloc+1,ik))**2)
     sd=sd+sum(rho(:)*rr(:)**3)
     vnext(:,ik)=vnew(:)
   end do
   change=sqrt(sn/sd)
   if(change<=tol) exit
   if(iter==maxiter) then
     write(6,'(/a,i4,a,1p,e10.2)') 'ionic_screening: ERROR the ionic tail &
&      &did not converge in',maxiter,' steps, change =',change
     stop
   end if
   vp(:,lloc+1,:)=vnext(:,:)
 end do

! the semi-local potentials, switched with the screening their local
! potential was switched with
 vsw(:)=vp(:,lloc+1,1)-ww(:)*vp0(:,lloc+1,1)
 do ik=1,nk
   do l1=1,max(lmax+1,lloc+1)
     if(l1==lloc+1) cycle
     vp(:,l1,ik)=ww(:)*vp0(:,l1,ik)+vsw(:)
   end do
 end do

 write(6,'(/a,f8.4,a,f8.4,a,i3,a)') 'Ionic tail switched over r =',ra, &
&      ' to',rb,' a.u., converged in',iter,' steps'

 deallocate(vion,ww,vp0,vnew,vnext,vsw)
 return
 end subroutine ionic_screening
