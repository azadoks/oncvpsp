! Unit tests of the physics of the all-electron atom (aeatom): relations the
! exact Kohn-Sham atom satisfies, each checked over a table of atoms (He, O,
! Si, Cu, Au neutral; O+, Si+ and Cu(0.5+) ions), nr, sr and fr, the VWN LDA
! (iexc 6) and PBE (iexc 4), on the program's default mesh (step 1 of
! doc/physics_tests_plan.md). PZ (iexc 3) is not among them: its potential
! steps at rs = 1, so it is not the derivative of its energy there and
! Janak's theorem and the virial miss by ~1e-6 Ha. The tolerances are those
! of the mesh (amesh 1.006, 4th-order Adams steps: ~1e-9 relative in the
! energies, as test_adams measures it) and of the self-consistency
! (ae_tolres 1e-11 in the residual; the total energy is stationary, so its
! error is second order in that):
!   charge     int r**2 rho dr = sum f, and int r**2 rhoc dr = the core's
!              sum f; 1e-10 relative
!   gauss      the electrostatic part of the potential, r (vfull - vxc[rho]),
!              is -Z at the nucleus and -(Z - N) where the density has died
!              off: the Hartree potential of the density's total charge.
!              At rr(1) the Hartree term adds rr(1) V_H(0) >= 0, at most
!              rr(1) int rho r dr, to within 1e-7 Z: vfull is the last
!              input potential, which at the first points is up to ~1e-8 Z
!              from the output (most for sr and fr, whose mixing filters
!              the residual near the nucleus); at rmax, 1e-8
!   virial     (nr) the scaling n -> g**3 n(g r) leaves E stationary at
!              g = 1: 2 T + E_ne + E_H + dE_xc[n_g]/dg = 0, the derivative
!              by 4th-order differences in g = exp(k al), k = -2..2, which
!              shift the density by k points of the log mesh; 1e-7 of T
!   janak      dE/df_i = eps_i: the 4th-order difference of etot in the
!              occupation of the highest valence state and of the first
!              (steps of 0.01); for fr the shares of f between j = l +- 1/2
!              are fixed (kap_j_share), so dE/df is their average of the two
!              eigenvalues; 1e-7 Ha + 1e-8 |eps|
!   hellmann   (VWN) dE/dlambda = int r**2 rho W dr for the external potential
!              lambda W (vext), W = -exp(-r**2) Ha, by the same differences
!              in lambda (steps of 0.01); 1e-7 of |int r**2 rho W dr|
!              For sr, both theorems hold only to O((Z/c)**2) relative: the
!              scalar-relativistic equation is not Dirac's, and its total
!              energy is that of the large components' density, normalized
!              alone, so it is not the variational energy whose derivatives
!              are the eigenvalues (~2e-5 Ha for O's 2s and 2p). For sr the
!              tolerance is (Z/c)**2 relative instead
!   nr_limit   (VWN) at c 10 and 20 times the physical value, sr and fr give
!              the nr eigenvalues and total energy, and the difference
!              falls by 4 (O(1/c**2)) when c doubles: the ratio within
!              4 (Z/c)**2 of 4 (the next order), and the eigenvalues'
!              noise, 1e-9 Ha, over the difference
!   so_order   (fr, VWN) for l > 0, the j = l - 1/2 level is below j = l + 1/2,
!              and their splitting falls by 4 when c doubles, to the same
!              tolerance
!   aufbau     within each l the eigenvalues rise with n
! Stops with a non-zero code if any check failed.
 program test_ae_physics
 use constants_m, only: dp, mxstate
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check, check_done
 use aeatom_m, only: aeatom
 use adams_m, only: adams_t, adams_init
 use main_setup_m, only: make_log_mesh
 use quad_m, only: quad_gregory_sum
 use logmesh_m, only: lm_head
 use xc_dispatch_m, only: exc_dispatch
 use kappa_m, only: kap_j_share
 use failure_m, only: warnings_t
 use relativity_m, only: c_light_default
 implicit none
 integer, parameter :: iout=6

!> an atom of the table: its configuration, the first nc states the core
 type :: atom_t
   character(len=6) :: name=' '
   real(dp) :: zz=0.0_dp
   integer :: nc=0,ncv=0
   integer :: na(mxstate)=0,la(mxstate)=0
   real(dp) :: fa(mxstate)=0.0_dp
 end type atom_t

!> a solution of aeatom
 type :: sol_t
   real(dp) :: etot=0.0_dp,eterm(4)=0.0_dp
   real(dp) :: ea(mxstate,2)=0.0_dp
   real(dp), allocatable :: rho(:),rhoc(:),vi(:)
 end type sol_t

 integer, parameter :: natom=8
 real(dp), parameter :: cfac=10.0_dp
 integer, parameter :: iexcs(2)=[6,4]
 type(atom_t) :: atoms(natom)
 type(adams_t) :: adams
 type(sol_t) :: s0
 real(dp), allocatable :: rr(:),vzero(:)
 real(dp) :: al,log_mesh(4)
 integer :: mmax,ia,irel,ix

 call table(atoms)
 adams=adams_init(4)

 do ia=1,natom
   call make_log_mesh(0,atoms(ia)%zz,-1.0_dp,-1.0_dp,-1.0_dp,4,al,mmax,rr, &
&                     log_mesh)
   allocate(vzero(mmax))
   vzero(:)=0.0_dp
   do ix=1,size(iexcs)
     do irel=0,2
       s0=solve(atoms(ia),irel,iexcs(ix),c_light_default,atoms(ia)%fa,vzero)
       call check_charge(atoms(ia),irel,iexcs(ix),s0)
       call check_gauss(atoms(ia),irel,iexcs(ix),s0)
       if(irel==0) call check_virial(atoms(ia),iexcs(ix),s0)
       call check_janak(atoms(ia),irel,iexcs(ix),s0)
       call check_aufbau(atoms(ia),irel,iexcs(ix),s0)
! the limits in c and the external potential need one functional
       if(irel==2 .and. ix==1) call check_so(atoms(ia),iexcs(ix),s0)
     end do
     if(ix==1) call check_hellmann(atoms(ia),iexcs(ix))
     if(ix==1) call check_nr_limit(atoms(ia),iexcs(ix))
   end do
   deallocate(vzero,rr)
 end do

 call check_done('ae_physics')

 contains

 subroutine table(atoms)
!! the atoms of the tests, with their reference configurations
 type(atom_t), intent(out) :: atoms(natom)
 atoms(1)=atom('He',2.0_dp,0,[1],[0],[2.0_dp])
 atoms(2)=atom('O',8.0_dp,1,[1,2,2],[0,0,1],[2.0_dp,2.0_dp,4.0_dp])
 atoms(3)=atom('O+',8.0_dp,1,[1,2,2],[0,0,1],[2.0_dp,2.0_dp,3.0_dp])
 atoms(4)=atom('Si',14.0_dp,3,[1,2,2,3,3],[0,0,1,0,1], &
&              [2.0_dp,2.0_dp,6.0_dp,2.0_dp,2.0_dp])
 atoms(5)=atom('Si+',14.0_dp,3,[1,2,2,3,3],[0,0,1,0,1], &
&              [2.0_dp,2.0_dp,6.0_dp,2.0_dp,1.0_dp])
 atoms(6)=atom('Cu',29.0_dp,3,[1,2,2,3,3,3,4],[0,0,1,0,1,2,0], &
&              [2.0_dp,2.0_dp,6.0_dp,2.0_dp,6.0_dp,10.0_dp,1.0_dp])
 atoms(7)=atom('Cu.5+',29.0_dp,3,[1,2,2,3,3,3,4],[0,0,1,0,1,2,0], &
&              [2.0_dp,2.0_dp,6.0_dp,2.0_dp,6.0_dp,10.0_dp,0.5_dp])
 atoms(8)=atom('Au',79.0_dp,10, &
&              [1,2,2,3,3,3,4,4,4,4,5,5,5,6],[0,0,1,0,1,2,0,1,2,3,0,1,2,0], &
&              [2.0_dp,2.0_dp,6.0_dp,2.0_dp,6.0_dp,10.0_dp,2.0_dp,6.0_dp, &
&               10.0_dp,14.0_dp,2.0_dp,6.0_dp,10.0_dp,1.0_dp])
 end subroutine table

 type(atom_t) function atom(name,zz,nc,na,la,fa)
 character(len=*), intent(in) :: name
 real(dp), intent(in) :: zz,fa(:)
 integer, intent(in) :: nc,na(:),la(:)
 atom%name=name
 atom%zz=zz
 atom%nc=nc
 atom%ncv=size(na)
 atom%na(1:atom%ncv)=na
 atom%la(1:atom%ncv)=la
 atom%fa(1:atom%ncv)=fa
 end function atom

 type(sol_t) function solve(at,irel,iexc,c_light,fa,vext)
!! the self-consistent atom at of occupations fa in vext, on the mesh rr
 type(atom_t), intent(in) :: at
 integer, intent(in) :: irel,iexc
 real(dp), intent(in) :: c_light,fa(mxstate),vext(mmax)
 type(warnings_t) :: warns
 integer :: nkap,it,ierr,ifail
 real(dp) :: res
 nkap=1
 if(irel==2) nkap=2
 allocate(solve%rho(mmax),solve%rhoc(mmax),solve%vi(mmax))
 solve%vi(:)=0.0_dp
 ierr=0
 call aeatom(irel,c_light,nkap,at%na,at%la,solve%ea(:,1:nkap), &
&            fa(1:at%ncv),at%nc,at%ncv,it,solve%rhoc,solve%rho,rr,solve%vi, &
&            at%zz,mmax,adams,iexc,solve%etot,solve%eterm,res,ierr,vext, &
&            .false.,warns,iout,ifail)
 if(ifail/=0 .or. ierr/=0) then
   write(stderr,'(3a,3i3)') 'test_ae_physics: aeatom failed for ', &
&    trim(at%name),' irel, iexc, ierr',irel,iexc,ierr
   error stop 1
 end if
 if(irel/=2) then
   solve%ea(:,2)=solve%ea(:,1)
 else
   where(at%la(1:at%ncv)==0) solve%ea(1:at%ncv,2)=solve%ea(1:at%ncv,1)
 end if
 end function solve

 real(dp) function radint(ff,rho,pow)
!! int_0^rmax ff dr on the log mesh, ff = rho r**pow; below rr(1) the
!! head of ff going as r**pow times the density's power (lm_head)
 real(dp), intent(in) :: ff(mmax),rho(mmax),pow
 radint=al*quad_gregory_sum(ff*rr,mmax) &
&       +lm_head(ff(1)*rr(1),rho(1),rho(2),al,pow+1.0_dp)
 end function radint

 character(len=40) function tag(at,irel,iexc)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: irel,iexc
 write(tag,'(a,1x,a,i0,a,i0)') trim(at%name),'irel=',irel,' iexc=',iexc
 end function tag

 real(dp) function srel(irel,at)
!! (Z/c)**2 for irel=1, else 0: the size of the difference between the
!! scalar-relativistic equation and Dirac's, by which the sr total energy
!! (of the large components' density) misses Janak's and Hellmann-Feynman's
!! theorems
!! TODO: sratom's density is the large components' alone; with the small
!! components in it (as the fr density has) the sr energy might satisfy both
!! theorems, and this bound could go. That would change every sr output
!! (doc/physics_tests_plan.md, Step 1)
 integer, intent(in) :: irel
 type(atom_t), intent(in) :: at
 srel=0.0_dp
 if(irel==1) srel=(at%zz/c_light_default)**2
 end function srel

 real(dp) function c_ratio_tol(at,dd)
!! the tolerance of a ratio d(c)/d(2c) = 4 at c = cfac c_light: the
!! O((Z/c)**2) next term, relative, and the eigenvalues' noise 1e-9 Ha
!! over the difference dd at 2c
 type(atom_t), intent(in) :: at
 real(dp), intent(in) :: dd
 c_ratio_tol=4.0_dp*(at%zz/(cfac*c_light_default))**2 &
&           +4.0_dp*1.0e-9_dp/abs(dd)
 end function c_ratio_tol

 character(len=2) function state(at,is)
!! the label of state is of at, as 3d
 type(atom_t), intent(in) :: at
 integer, intent(in) :: is
 character(len=4), parameter :: lsym='spdf'
 state=char(48+at%na(is))//lsym(at%la(is)+1:at%la(is)+1)
 end function state

 subroutine check_charge(at,irel,iexc,s)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: irel,iexc
 type(sol_t), intent(in) :: s
 real(dp) :: nn,nc
 nn=sum(at%fa(1:at%ncv))
 nc=sum(at%fa(1:at%nc))
 call check('charge: total '//tag(at,irel,iexc), &
&           abs(radint(rr**2*s%rho,s%rho,2.0_dp)-nn)/nn,1.0e-10_dp)
 if(at%nc>0) call check('charge: core '//tag(at,irel,iexc), &
&           abs(radint(rr**2*s%rhoc,s%rhoc,2.0_dp)-nc)/nc,1.0e-10_dp)
 end subroutine check_charge

 subroutine check_gauss(at,irel,iexc,s)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: irel,iexc
 type(sol_t), intent(in) :: s
 real(dp) :: vxc(mmax),exc(mmax),ww(mmax),vh0,nn
 integer :: ifail
 call exc_dispatch(iexc,c_light_default,al,s%rho,vxc,exc,rr,mmax, &
&                  'test_ae_physics',ifail)
 ww(:)=rr*(s%vi-vxc)
 vh0=radint(rr*s%rho,s%rho,1.0_dp)
 nn=sum(at%fa(1:at%ncv))
 call check('gauss: nucleus '//tag(at,irel,iexc), &
&           max(-(ww(1)+at%zz),ww(1)+at%zz-rr(1)*vh0),1.0e-7_dp*at%zz)
 call check('gauss: tail '//tag(at,irel,iexc), &
&           abs(ww(mmax)+(at%zz-nn)),1.0e-8_dp)
 end subroutine check_gauss

 subroutine check_virial(at,iexc,s)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: iexc
 type(sol_t), intent(in) :: s
 real(dp) :: ex(-2:2),dex,lhs
 integer :: kk
 do kk=-2,2
   if(kk/=0) ex(kk)=exc_scaled(iexc,s%rho,kk)
 end do
! g = exp(k al): dE/dg at g = 1 from the 4th-order difference in k al
 dex=(ex(-2)-8.0_dp*ex(-1)+8.0_dp*ex(1)-ex(2))/(12.0_dp*al)
 lhs=2.0_dp*s%eterm(1)+s%eterm(2)+s%eterm(3)+dex
 call check('virial: '//tag(at,0,iexc),abs(lhs)/s%eterm(1),1.0e-7_dp)
 end subroutine check_virial

 real(dp) function exc_scaled(iexc,rho,kk)
!! E_xc of g**3 rho(g r), g = exp(kk al): the density shifted by kk points
!! of the log mesh, held at its first value below rr(1) and zero past rmax
 integer, intent(in) :: iexc,kk
 real(dp), intent(in) :: rho(mmax)
 real(dp) :: rhog(mmax),vxc(mmax),exc(mmax)
 integer :: ii,jj,ifail
 do ii=1,mmax
   jj=ii+kk
   if(jj<1) then
     rhog(ii)=rho(1)
   else if(jj>mmax) then
     rhog(ii)=0.0_dp
   else
     rhog(ii)=rho(jj)
   end if
 end do
 rhog(:)=exp(3.0_dp*kk*al)*rhog
 call exc_dispatch(iexc,c_light_default,al,rhog,vxc,exc,rr,mmax, &
&                  'test_ae_physics',ifail)
 exc_scaled=radint(rr**2*rhog*exc,rhog,2.0_dp)
 end function exc_scaled

 subroutine check_janak(at,irel,iexc,s)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: irel,iexc
 type(sol_t), intent(in) :: s
 real(dp), parameter :: df=0.01_dp
 real(dp) :: fa(mxstate),ee(-2:2),dedf,eps
 type(sol_t) :: sf
 integer :: is,kk,ik,iv
 do iv=1,2
   is=at%ncv
   if(iv==2) is=at%nc+1
   if(iv==2 .and. is==at%ncv) cycle
   do kk=-2,2
     if(kk==0) cycle
     fa(:)=at%fa
     fa(is)=fa(is)+kk*df
     sf=solve(at,irel,iexc,c_light_default,fa,vzero)
     ee(kk)=sf%etot
   end do
   dedf=(ee(-2)-8.0_dp*ee(-1)+8.0_dp*ee(1)-ee(2))/(12.0_dp*df)
   if(irel==2) then
     eps=0.0_dp
     do ik=1,2
       eps=eps+kap_j_share(at%la(is),ik,1.0_dp)*s%ea(is,ik)
     end do
   else
     eps=s%ea(is,1)
   end if
   call check('janak: '//state(at,is)//' '//tag(at,irel,iexc),abs(dedf-eps), &
&             1.0e-7_dp+max(1.0e-8_dp,srel(irel,at))*abs(eps))
 end do
 end subroutine check_janak

 subroutine check_hellmann(at,iexc)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: iexc
 real(dp), parameter :: dl=0.01_dp
 real(dp) :: ww(mmax),ee(-2:2),dedl,expect,scale
 type(sol_t) :: s
 integer :: irel,kk
 ww(:)=-exp(-rr**2)
 do irel=0,2
   do kk=-2,2
     s=solve(at,irel,iexc,c_light_default,at%fa,kk*dl*ww)
     ee(kk)=s%etot
     if(kk==0) then
       expect=radint(rr**2*s%rho*ww,s%rho,2.0_dp)
       scale=abs(expect)
     end if
   end do
   dedl=(ee(-2)-8.0_dp*ee(-1)+8.0_dp*ee(1)-ee(2))/(12.0_dp*dl)
   call check('hellmann: '//tag(at,irel,iexc),abs(dedl-expect)/scale, &
&             max(1.0e-7_dp,srel(irel,at)))
 end do
 end subroutine check_hellmann

 subroutine check_nr_limit(at,iexc)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: iexc
 type(sol_t) :: snr,s1,s2
 real(dp) :: d1,d2
 integer :: irel,is,ik,nkap
 snr=solve(at,0,iexc,c_light_default,at%fa,vzero)
 do irel=1,2
   nkap=1
   if(irel==2) nkap=2
   s1=solve(at,irel,iexc,cfac*c_light_default,at%fa,vzero)
   s2=solve(at,irel,iexc,2.0_dp*cfac*c_light_default,at%fa,vzero)
   d1=s1%etot-snr%etot
   d2=s2%etot-snr%etot
   call check('nr_limit: etot '//tag(at,irel,iexc),abs(d1/d2-4.0_dp), &
&             c_ratio_tol(at,d2))
   do is=1,at%ncv
     do ik=1,nkap
       if(ik==2 .and. at%la(is)==0) cycle
       d1=s1%ea(is,ik)-snr%ea(is,1)
       d2=s2%ea(is,ik)-snr%ea(is,1)
       call check('nr_limit: eps '//state(at,is)//char(48+ik)//' '// &
&                 tag(at,irel,iexc),abs(d1/d2-4.0_dp),c_ratio_tol(at,d2))
     end do
   end do
 end do
 end subroutine check_nr_limit

 subroutine check_so(at,iexc,s)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: iexc
 type(sol_t), intent(in) :: s
 type(sol_t) :: s1,s2
 integer :: is
 s1=solve(at,2,iexc,cfac*c_light_default,at%fa,vzero)
 s2=solve(at,2,iexc,2.0_dp*cfac*c_light_default,at%fa,vzero)
 do is=1,at%ncv
   if(at%la(is)==0) cycle
   call check('so_order: '//state(at,is)//' '//tag(at,2,iexc), &
&             s%ea(is,2)<s%ea(is,1))
   call check('so_scaling: '//state(at,is)//' '//tag(at,2,iexc), &
&             abs((s1%ea(is,1)-s1%ea(is,2))/(s2%ea(is,1)-s2%ea(is,2)) &
&                 -4.0_dp),c_ratio_tol(at,s2%ea(is,1)-s2%ea(is,2)))
 end do
 end subroutine check_so

 subroutine check_aufbau(at,irel,iexc,s)
 type(atom_t), intent(in) :: at
 integer, intent(in) :: irel,iexc
 type(sol_t), intent(in) :: s
 integer :: is,js,ik
 do is=1,at%ncv
   do js=1,at%ncv
     if(at%la(js)/=at%la(is) .or. at%na(js)/=at%na(is)+1) cycle
     do ik=1,2
       call check('aufbau: '//tag(at,irel,iexc),s%ea(is,ik)<s%ea(js,ik))
     end do
   end do
 end do
 end subroutine check_aufbau

 end program test_ae_physics
