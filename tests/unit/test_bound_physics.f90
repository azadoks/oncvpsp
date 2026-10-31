! Unit tests of the physics of the bound-state solvers (lschfb, ldiracfb) on
! potentials whose spectra are known (step 2 of doc/physics_tests_plan.md).
! The hydrogen-like atoms are zz = 1, 6, 26 and 79, every state to n = 5,
! each on the program's mesh scaled by 1/zz (amesh 1.006, rr1 5e-5/zz, rmax
! 150/zz, far enough out for 5s's tail); test_solver_path checks only
! hydrogen's 1s, 2p and 3d, for the independence of the iteration path.
! Each quantity is computed on that mesh and on the mesh of half its step
! (amesh 1.003): the 4th-order Adams steps make the difference from the
! exact value fall by 2**4 (test_adams measures the order), so it is the
! mesh's alone, and the Richardson value (16 q(1.003) - q(1.006))/15 leaves
! its next order, a few 1e-3 of it (the difference falls by 15.6 to 15.9):
! 1e-10 relative. The invariants are checked on that value, and the
! program's mesh within 1e-7 (the highest states, with the most nodes, are
! off by up to 5e-8):
!   hydrogen   (nr) E = -zz**2/(2 n**2) for every l < n; 1e-10 relative
!   dirac      (fr) for every kappa, E = c**2 (1/sqrt(1 + t**2) - 1),
!              t = zz/c/(n - |kappa| + sqrt(kappa**2 - (zz/c)**2)), the
!              point nucleus' Dirac energy, evaluated without the
!              cancellation as -c**2 t**2/(s (1 + s)), s = sqrt(1 + t**2);
!              1e-10 relative
!   j_average  (sr) the scalar-relativistic equation leaves out spin-orbit
!              only, so for l > 0 its energy is the Dirac energies' average
!              over j = l +- 1/2 weighted by 2j + 1 to first order in
!              (zz/c)**2: they differ by O((zz/c)**4 zz**2/n**3), the
!              coefficient below 1. For s it is Dirac's equation for the
!              large component of kappa = -1 (l (l+1) = kappa (kappa+1) = 0,
!              and the spin-orbit term, ~ 1 + kappa, vanishes), so the
!              energies are equal: 1e-10 relative
!   scaling    (nr) on meshes scaled by 1/zz the equation is the same with
!              r -> r/zz, so E(zz) = zz**2 E(1) to round-off: 1e-12
!              relative (eig_tol, the solvers' convergence)
!   moments    (nr) <r> = (3 n**2 - l (l+1))/(2 zz) and <1/r> = zz/n**2;
!              1e-10 relative
!   nodes      u (nr) and the large component g (fr) have n - l - 1 nodes
!   orthogonal two states of one l (nr), or of one kappa (fr), are
!              orthogonal: int u_n u_n' dr, int (g_n g_n' + f_n f_n') dr
!              within 1e-9 of their norms, on the program's mesh
!   oscillator (nr) V = omega**2 r**2/2, omega 0.25, 1 and 4 (zz = 0):
!              E = omega (2 n_r + l + 3/2), n_r = 0..2, l = 0..3, on the
!              program's mesh unscaled (rr1 5e-5, rmax 45); 1e-10 relative
!   hellmann   (nr, fr) dE/dlambda = <W> for V = -zz/r + lambda W, W =
!              -zz**2 exp(-(zz r)**2): the 4th-order difference in lambda
!              (steps of 0.01) against int u**2 W dr, or int (g**2 + f**2)
!              W dr, for 1s, 2p and 3d (and 2p's j = 1/2), on the
!              program's mesh; 1e-8 relative. sr is left out: lschfb
!              returns only the large component, whose norm is not the one
!              the sr eigenvalue is stationary in
! Stops with a non-zero code if any check failed.
 program test_bound_physics
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check, check_done
 use adams_m, only: adams_t, adams_init
 use lschfb_m, only: lschfb
 use ldiracfb_m, only: ldiracfb
 use quad_m, only: quad_log_norm
 use relativity_m, only: c_light => c_light_default
 implicit none
 integer, parameter :: iout=6

!> a log mesh
 type :: mesh_t
   integer :: mmax=0
   real(dp), allocatable :: rr(:)
 end type mesh_t

!Local variables
 integer, parameter :: nz=4, nmax=5
 real(dp), parameter :: zs(nz)=[1.0_dp,6.0_dp,26.0_dp,79.0_dp]
!> the program's mesh ratio, and half its step
 real(dp), parameter :: ameshes(2)=[1.006_dp,1.003_dp]
 real(dp), parameter :: tol_mesh=1.0e-7_dp, tol_e=1.0e-10_dp
 real(dp), parameter :: tol_scale=1.0e-12_dp, tol_ovl=1.0e-9_dp
 real(dp), parameter :: tol_hf=1.0e-8_dp
 type(adams_t) :: adams
 type(mesh_t) :: meshes(2)
 real(dp) :: e1(nmax,0:nmax-1)
 integer :: iz

 adams=adams_init(4)

 do iz=1,nz
   call make_meshes(5.0e-5_dp/zs(iz),150.0_dp/zs(iz))
   call check_hydrogen(zs(iz))
   call check_dirac(zs(iz))
   call check_j_average(zs(iz))
   call check_hellmann(zs(iz))
 end do

 call make_meshes(5.0e-5_dp,45.0_dp)
 call check_oscillator()

 call check_done('bound_physics')

 contains

 subroutine make_meshes(rr1,rmax)
!! the log meshes of ratios ameshes from rr1 to rmax
 real(dp), intent(in) :: rr1,rmax
 real(dp) :: al
 integer :: im,ii
 do im=1,2
   al=log(ameshes(im))
   meshes(im)%mmax=int(log(rmax/rr1)/al)
   if(allocated(meshes(im)%rr)) deallocate(meshes(im)%rr)
   allocate(meshes(im)%rr(meshes(im)%mmax))
   do ii=1,meshes(im)%mmax
     meshes(im)%rr(ii)=rr1*exp(al*(ii-1))
   end do
 end do
 end subroutine make_meshes

 real(dp) function richardson(qq)
!! the value without the 4th-order mesh error, from qq on the two meshes
 real(dp), intent(in) :: qq(2)
 richardson=(16.0_dp*qq(2)-qq(1))/15.0_dp
 end function richardson

 real(dp) function integral(mesh,ff,gam)
!! int ff dr over the mesh, ff going as r**(2 gam) at the origin
 type(mesh_t), intent(in) :: mesh
 real(dp), intent(in) :: ff(:),gam
 integral=quad_log_norm(mesh%rr,ff,mesh%mmax,0,gam)
 end function integral

 integer function nodes(uu)
!! the sign changes of uu, ignoring values below 1e-10 of its largest
 real(dp), intent(in) :: uu(:)
 real(dp) :: umax,slast
 integer :: ii
 umax=maxval(abs(uu))
 slast=0.0_dp
 nodes=0
 do ii=1,size(uu)
   if(abs(uu(ii))<1.0e-10_dp*umax) cycle
   if(slast/=0.0_dp .and. sign(1.0_dp,uu(ii))/=slast) nodes=nodes+1
   slast=sign(1.0_dp,uu(ii))
 end do
 end function nodes

 function label(what,zz,nn,ll,kap)
!! a check's label: what, zz, the state and kappa (fr only)
 character(len=*), intent(in) :: what
 real(dp), intent(in) :: zz
 integer, intent(in) :: nn,ll
 integer, intent(in), optional :: kap
 character(len=64) :: label
 write(label,'(a,a,i0,a,i0,a1)') what,' zz ',nint(zz),' ',nn, &
&  'spdf'(ll+1:ll+1)
 if(present(kap)) write(label,'(2a,i0)') trim(label),' kap ',kap
 end function label

 subroutine solve_sch(mesh,irel,zz,nn,ll,vv,ee,uu)
!! the bound state nn, ll of vv by lschfb, from the guess ee
 type(mesh_t), intent(in) :: mesh
 integer, intent(in) :: irel,nn,ll
 real(dp), intent(in) :: zz,vv(:)
 real(dp), intent(inout) :: ee
 real(dp), intent(out) :: uu(:)
 real(dp) :: up(mesh%mmax)
 integer :: ierr,mch,ifail
 mch=0
 uu(:)=0.0_dp
 up(:)=0.0_dp
 call lschfb(irel,c_light,nn,ll,ierr,ee,mesh%rr,zz,vv,uu,up,mesh%mmax, &
&            adams,mch,iout,ifail)
 if(ifail/=0 .or. ierr/=0) then
   write(stderr,'(a,i0,a,i0,a,i0)') 'lschfb failed: ierr ',ierr,' n ',nn, &
&    ' l ',ll
   error stop 1
 end if
 end subroutine solve_sch

 subroutine solve_dirac(mesh,zz,nn,ll,kap,vv,ee,uu)
!! the bound state nn, kap of vv by ldiracfb, from the guess ee
 type(mesh_t), intent(in) :: mesh
 integer, intent(in) :: nn,ll,kap
 real(dp), intent(in) :: zz,vv(:)
 real(dp), intent(inout) :: ee
 real(dp), intent(out) :: uu(:,:)
 real(dp) :: up(mesh%mmax,2)
 integer :: ierr,mch,ifail
 mch=0
 uu(:,:)=0.0_dp
 up(:,:)=0.0_dp
 call ldiracfb(c_light,nn,ll,kap,ierr,ee,mesh%rr,zz,vv,uu,up,mesh%mmax, &
&              adams,mch,iout,ifail)
 if(ifail/=0 .or. ierr/=0) then
   write(stderr,'(a,i0,a,i0,a,i0)') 'ldiracfb failed: ierr ',ierr,' n ', &
&    nn,' kap ',kap
   error stop 1
 end if
 end subroutine solve_dirac

 subroutine check_hydrogen(zz)
!! hydrogen, scaling, moments, nodes and orthogonal (nr)
 real(dp), intent(in) :: zz
 real(dp), allocatable :: uu(:),u1(:,:,:)
 real(dp) :: ee(2),rm(2),rmi(2),ex,norm,ovl
 integer :: nn,ll,mm,im
 allocate(u1(meshes(1)%mmax,nmax,0:nmax-1))
 do ll=0,nmax-1
   do nn=ll+1,nmax
     ex=-zz**2/(2.0_dp*nn**2)
     do im=1,2
       associate(rr=>meshes(im)%rr)
         allocate(uu(meshes(im)%mmax))
         ee(im)=1.05_dp*ex
         call solve_sch(meshes(im),0,zz,nn,ll,-zz/rr,ee(im),uu)
         norm=integral(meshes(im),uu**2,ll+1.0_dp)
         rm(im)=integral(meshes(im),uu**2*rr,ll+1.5_dp)/norm
         rmi(im)=integral(meshes(im),uu**2/rr,ll+0.5_dp)/norm
         if(im==1) u1(:,nn,ll)=uu
         deallocate(uu)
       end associate
     end do
     call check(label('hydrogen mesh',zz,nn,ll),abs(ee(1)/ex-1.0_dp), &
&               tol_mesh)
     call check(label('hydrogen',zz,nn,ll),abs(richardson(ee)/ex-1.0_dp), &
&               tol_e)
     if(zz==zs(1)) e1(nn,ll)=ee(1)
     call check(label('scaling',zz,nn,ll), &
&               abs(ee(1)/(zz**2*e1(nn,ll))-1.0_dp),tol_scale)
     call check(label('moments <r>',zz,nn,ll), &
&               abs(richardson(rm)*2.0_dp*zz/(3*nn**2-ll*(ll+1))-1.0_dp), &
&               tol_e)
     call check(label('moments <1/r>',zz,nn,ll), &
&               abs(richardson(rmi)*nn**2/zz-1.0_dp),tol_e)
     call check(label('nodes',zz,nn,ll),nodes(u1(:,nn,ll))==nn-ll-1)
     norm=integral(meshes(1),u1(:,nn,ll)**2,ll+1.0_dp)
     do mm=ll+1,nn-1
       ovl=integral(meshes(1),u1(:,nn,ll)*u1(:,mm,ll),ll+1.0_dp)
       call check(label('orthogonal',zz,nn,ll),abs(ovl) &
&                 /sqrt(norm*integral(meshes(1),u1(:,mm,ll)**2,ll+1.0_dp)), &
&                 tol_ovl)
     end do
   end do
 end do
 end subroutine check_hydrogen

 real(dp) function e_dirac(zz,nn,kap)
!! the point nucleus' Dirac energy, without the rest energy
 real(dp), intent(in) :: zz
 integer, intent(in) :: nn,kap
 real(dp) :: tt,ss
 tt=(zz/c_light)/(nn-abs(kap)+sqrt(kap**2-(zz/c_light)**2))
 ss=sqrt(1.0_dp+tt**2)
 e_dirac=-c_light**2*tt**2/(ss*(1.0_dp+ss))
 end function e_dirac

 subroutine check_dirac(zz)
!! dirac, nodes and orthogonal (fr)
 real(dp), intent(in) :: zz
 real(dp), allocatable :: uu(:,:),u1(:,:,:)
 real(dp) :: ee(2),ex,gam,norm(nmax),ovl
 integer :: nn,ll,kap,mm,is,im
 allocate(u1(meshes(1)%mmax,2,nmax))
 do ll=0,nmax-1
   do is=1,2
! kap = -(l+1) (j = l+1/2), then l (j = l-1/2)
     kap=merge(-(ll+1),ll,is==1)
     if(kap==0) cycle
     gam=sqrt(kap**2-(zz/c_light)**2)
     do nn=ll+1,nmax
       ex=e_dirac(zz,nn,kap)
       do im=1,2
         allocate(uu(meshes(im)%mmax,2))
         ee(im)=1.05_dp*ex
         call solve_dirac(meshes(im),zz,nn,ll,kap,-zz/meshes(im)%rr,ee(im),uu)
         if(im==1) u1(:,:,nn)=uu
         deallocate(uu)
       end do
       call check(label('dirac mesh',zz,nn,ll,kap),abs(ee(1)/ex-1.0_dp), &
&                 tol_mesh)
       call check(label('dirac',zz,nn,ll,kap),abs(richardson(ee)/ex-1.0_dp), &
&                 tol_e)
       call check(label('nodes',zz,nn,ll,kap),nodes(u1(:,1,nn))==nn-ll-1)
       norm(nn)=integral(meshes(1),u1(:,1,nn)**2+u1(:,2,nn)**2,gam)
       do mm=ll+1,nn-1
         ovl=integral(meshes(1),u1(:,1,nn)*u1(:,1,mm)+u1(:,2,nn)*u1(:,2,mm), &
&                     gam)
         call check(label('orthogonal',zz,nn,ll,kap), &
&                   abs(ovl)/sqrt(norm(nn)*norm(mm)),tol_ovl)
       end do
     end do
   end do
 end do
 end subroutine check_dirac

 subroutine check_j_average(zz)
!! j_average (sr): the Dirac energies' (2j+1)-weighted average over j
 real(dp), intent(in) :: zz
 real(dp), allocatable :: uu(:)
 real(dp) :: ee(2),ex,eo
 integer :: nn,ll,im
 do ll=0,nmax-1
   do nn=ll+1,nmax
     if(ll==0) then
       ex=e_dirac(zz,nn,-1)
     else
       ex=((ll+1)*e_dirac(zz,nn,-(ll+1))+ll*e_dirac(zz,nn,ll))/(2*ll+1)
     end if
     do im=1,2
       allocate(uu(meshes(im)%mmax))
       ee(im)=1.05_dp*ex
       call solve_sch(meshes(im),1,zz,nn,ll,-zz/meshes(im)%rr,ee(im),uu)
       deallocate(uu)
     end do
     if(ll==0) then
       call check(label('j_average',zz,nn,ll),abs(richardson(ee)/ex-1.0_dp), &
&                 tol_e)
     else
! the next order, relative to the first
       eo=(zz/c_light)**4*zz**2/nn**3
       call check(label('j_average',zz,nn,ll),abs(richardson(ee)-ex)/eo, &
&                 1.0_dp)
     end if
   end do
 end do
 end subroutine check_j_average

 subroutine check_oscillator()
!! oscillator (nr, zz = 0)
 real(dp), parameter :: omegas(3)=[0.25_dp,1.0_dp,4.0_dp]
 real(dp), allocatable :: uu(:)
 real(dp) :: ee(2),ex
 integer :: io,nr,ll,im
 character(len=64) :: lab
 do io=1,size(omegas)
   do ll=0,3
     do nr=0,2
       ex=omegas(io)*(2*nr+ll+1.5_dp)
       do im=1,2
         allocate(uu(meshes(im)%mmax))
         ee(im)=0.95_dp*ex
         call solve_sch(meshes(im),0,0.0_dp,nr+ll+1,ll, &
&                       0.5_dp*omegas(io)**2*meshes(im)%rr**2,ee(im),uu)
         deallocate(uu)
       end do
       write(lab,'(a,f0.2,a,i0,a,i0)') 'oscillator omega ',omegas(io), &
&        ' n_r ',nr,' l ',ll
       call check(trim(lab)//' mesh',abs(ee(1)/ex-1.0_dp),tol_mesh)
       call check(lab,abs(richardson(ee)/ex-1.0_dp),tol_e)
     end do
   end do
 end do
 end subroutine check_oscillator

 subroutine check_hellmann(zz)
!! hellmann (nr, fr): 1s, 2p and 3d, and fr's 2p j = 1/2
 real(dp), intent(in) :: zz
 real(dp), parameter :: hh=0.01_dp
 integer, parameter :: nst=4
 integer, parameter :: nns(nst)=[1,2,3,2],lls(nst)=[0,1,2,1]
 integer, parameter :: kaps(nst)=[-1,-2,-3,1]
! lambda = 0 last, for its function
 integer, parameter :: korder(5)=[-2,-1,1,2,0]
 real(dp), allocatable :: ww(:),uu(:,:),dens(:)
 real(dp) :: ee(-2:2),ex,dedl,wexp,gam
 integer :: is,kk,ik,irel
 associate(mesh=>meshes(1),rr=>meshes(1)%rr)
   allocate(ww(mesh%mmax),uu(mesh%mmax,2),dens(mesh%mmax))
   ww(:)=-zz**2*exp(-(zz*rr)**2)
   do is=1,nst
     do irel=0,2,2
       if(irel==0 .and. kaps(is)>0) cycle
       ex=-zz**2/(2.0_dp*nns(is)**2)
       do ik=1,5
         kk=korder(ik)
         ee(kk)=ex
         if(irel==0) then
           call solve_sch(mesh,0,zz,nns(is),lls(is),-zz/rr+kk*hh*ww,ee(kk), &
&                         uu(:,1))
         else
           call solve_dirac(mesh,zz,nns(is),lls(is),kaps(is),-zz/rr+kk*hh*ww, &
&                           ee(kk),uu)
         end if
       end do
       dedl=(ee(-2)-8.0_dp*ee(-1)+8.0_dp*ee(1)-ee(2))/(12.0_dp*hh)
       if(irel==0) then
         dens(:)=uu(:,1)**2
         gam=lls(is)+1.0_dp
       else
         dens(:)=uu(:,1)**2+uu(:,2)**2
         gam=sqrt(kaps(is)**2-(zz/c_light)**2)
       end if
       wexp=integral(mesh,dens*ww,gam)/integral(mesh,dens,gam)
       if(irel==0) then
         call check(label('hellmann nr',zz,nns(is),lls(is)), &
&                   abs(dedl/wexp-1.0_dp),tol_hf)
       else
         call check(label('hellmann fr',zz,nns(is),lls(is),kaps(is)), &
&                   abs(dedl/wexp-1.0_dp),tol_hf)
       end if
     end do
   end do
 end associate
 end subroutine check_hellmann

 end program test_bound_physics
