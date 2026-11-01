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
!> Unit test of the physics of the pseudopotential: a whole run
!> (oncvpsp_run) and its channels (optimize_channel) on the input's mesh,
!> each check an invariant of the construction:
!>   norm       int_0^rc phi_i phi_j = int_0^rc u_i u_j for the pseudo-waves
!>              phi and all-electron states u of each channel's projectors
!>              (generalized norm conservation; the large components for
!>              irel 2, as fpovlp); 2e-10, the error of the mesh quadrature
!>              against the optimizer's integrals of its basis (it falls 32
!>              times with the log step halved: O 2.4e-12, 8e-14)
!>   match      phi = u beyond rc; 0, the construction copies u
!>   charge     int r**2 rho_ps dr = zion; 1e-12
!>   tail       r V_loc(r) + zion is in [-Q_c(r) - d, d] beyond max(2 rcmax,
!>              r_core) (rcmax the largest core radius, beyond r_core the
!>              model core the core): the unscreened local potential is the
!>              ionic one, -Z/r + V_H[n_c], whose Hartree term is at most
!>              (Z - zion)/r and at least that less the core charge Q_c(r)
!>              outside r over r; d = |int r**2 n_c dr - (Z - zion)|, the
!>              core density's charge error on the mesh (5e-10 for W sr),
!>              as the Hartree integrals start from Z - zion (vout);
!>              1e-10, the round-off of the Hartree integrals.
!>              For irel 1, 2 only beyond where Q_c(r)/r < 1e-10, as
!>              ionic_tail_m switches to the ionic potential only there
!>   local      the projectors vanish beyond their rc, on the mesh and as
!>              pieces; 0
!>   core       the model core is the core beyond r_core, and is >= 0; 0
!>   curves     each convergence curve E_resid(q), an integral from q of
!>              q**4 |phi(q)|**2 >= 0, does not rise; 1e-9 of its start,
!>              the round-off of its quadratic form, evaluated at each q
!>              apart, whose terms cancel to the curve's value
!>   so_kernel  (irel 2) the scalar-relativistic and spin-orbit projectors
!>              give back the kernels of the two j:
!>              sum_m esr_m |vsr_m><vsr_m| + <L.S> sum_m eso_m |vso_m><vso_m|
!>              = sum_i e_ij |beta_ij><beta_ij|, <L.S> = l/2, -(l+1)/2; 1e-13
!>              of the largest kernel value (test_sr_so_r on a real run)
!> and, for irel 0, the invariants that hold to the mesh's error, from the
!> residuals r1 on the input's mesh and r2 on the mesh of half its log
!> step: they are the mesh's error if |r2| <= 2 |r1 - r2|/15 + floor (the
!> error falls 16 times at the Adams order; twice that estimate, as the
!> heavy atoms at amesh 1.006 are not yet asymptotic and fall faster):
!>   eigenvalues  the pseudo-atom's reference states have the all-electron
!>                eigenvalues; floor 1e-10 Ha (the solvers' 1e-12)
!>   logder       atan(r phi'/phi) at r = rc 1.006**2 and each projector's
!>                energy is the all-electron one; floor 1e-9
!>   dlogder      and so is its energy derivative, -2 int_0^r u**2/u(r)**2
!>                times d atan/dx, equal by norm conservation; by 4th-order
!>                differences in steps of 1e-4 Ha, relative; floor 1e-8
!>   density      rho_ps is the all-electron valence density beyond rcmax,
!>                1e-9 of its largest value (the pseudo states' tails
!>                decay as their eigenvalues, equal to ~1e-8 Ha)
!> For irel 1 and 2 these four hold only to O(1/c**2): the pseudo-atom's
!> states miss the relativistic terms beyond rc (ionic_tail_m), 1e-7 Ha
!> (Si) to 1e-5 Ha (W) in the eigenvalues, falling 4 times as c doubles
!> (doc/physics_tests_plan.md, Step 4).
!> Stops with a non-zero code if any check failed.
!>
!> usage: test_ps_physics IREL INPUT.dat

 program test_ps_physics
 use constants_m, only: dp, mxcnf, mxlin, mxl, mxprj
 use unit_check_m, only: check, check_done
 use input_m, only: input_t, config_t
 use input_text_m, only: read_input_text
 use oncvpsp_run_m, only: oncvpsp_run, complete_input
 use failure_m, only: failure_t
 use ae_state_m, only: ae_state_t, solve_ae
 use channel_m, only: channel_t, optimize_channel
 use main_setup_m, only: core_radius_points
 use results_m, only: results_t
 use logmesh_m, only: lm_integral
 use phsft_m, only: phsft_pair
 use pieces_m, only: pw_node_values
 use kappa_m, only: kap_of
 implicit none

!> a run on one mesh: the pseudopotential and the stages it was made from
 type :: build_t
   type(results_t) :: psp
   type(ae_state_t) :: ae
   type(channel_t) :: chs(mxlin+1)
   integer :: mmax=0
 end type build_t

!> the residuals of the mesh-limited invariants of a build
 type :: resid_t
   real(dp), allocatable :: eig(:),ld(:),dld(:)
 end type resid_t

 real(dp), parameter :: tol_norm=2.0d-10,tol_charge=1.0d-12
 real(dp), parameter :: tol_tail=1.0d-10,tol_curve=1.0d-9
 real(dp), parameter :: tol_so=1.0d-13,tol_density=1.0d-9
 real(dp), parameter :: floor_eig=1.0d-10,floor_ld=1.0d-9
 real(dp), parameter :: floor_dld=1.0d-8
 type(build_t) :: b1,b2
 type(resid_t) :: r1,r2
 type(input_t) :: inp
 integer :: irel,nkap,iout
 character(len=1024) :: arg,input_filename

 call get_command_argument(1,arg)
 read(arg,*) irel
 call get_command_argument(2,input_filename)
 nkap=1
 if(irel==2) nkap=2
 open(newunit=iout,file='/dev/null',action='write')

 call build(-1.0d0,b1)
 call check_norm(b1)
 call check_matching(b1)
 call check_charge(b1)
 call check_tail(b1)
 call check_local(b1)
 call check_core(b1)
 call check_curves(b1)
 if(irel==2) call check_so(b1)
 if(irel==0) then
   call check_density(b1)
   call build(sqrt(b1%ae%rr(2)/b1%ae%rr(1)),b2)
   call residuals(b1,r1)
   call residuals(b2,r2)
   call check_mesh('eigenvalues',r1%eig,r2%eig,floor_eig)
   call check_mesh('logder',r1%ld,r2%ld,floor_ld)
   call check_mesh('dlogder',r1%dld,r2%dld,floor_dld)
 end if
 close(iout)
 call check_done('ps_physics')

 contains

 subroutine read_input(inp)
 type(input_t), intent(out) :: inp
 integer :: unit,inline,ifail
 open(newunit=unit,file=trim(input_filename),status='old',action='read')
 call read_input_text(unit,inline,inp,iout,ifail)
 close(unit)
 if(ifail/=0) error stop 'test_ps_physics: read_input_text failed'
 end subroutine read_input

 subroutine build(amesh,bb)
!! the run and its stages on the mesh of ratio amesh (the input's if <= 0)
 real(dp), intent(in) :: amesh
 type(build_t), intent(out) :: bb
 type(failure_t) :: at
 type(config_t) :: cnf(mxcnf)
 integer :: irc(mxlin+1),l1,nprj,ierr,ifail
 real(dp) :: rcmax

 call read_input(inp)
 if(amesh>0.0d0) inp%amesh=amesh
 call oncvpsp_run(irel,0,inp,'','',.false.,bb%psp,at,iout,ifail)
 if(ifail/=0) error stop 'test_ps_physics: oncvpsp_run failed'
 call read_input(inp)
 if(amesh>0.0d0) inp%amesh=amesh
 call complete_input(inp,cnf,iout,ifail)
 if(ifail/=0) error stop 'test_ps_physics: complete_input failed'
 call solve_ae(irel,nkap,0,inp,cnf,'',bb%ae,ierr,iout,ifail)
 if(ifail/=0 .or. ierr/=0) error stop 'test_ps_physics: solve_ae failed'
 call core_radius_points(inp%lmax,inp%lloc,inp%rc,bb%ae%rr,bb%ae%mmax,irc, &
&                        rcmax,iout,ifail)
 if(ifail/=0) error stop 'test_ps_physics: core_radius_points failed'
 bb%mmax=bb%ae%mmax
 do l1=1,inp%lmax+1
   nprj=inp%nproj(l1)
   if(l1==inp%lloc+1) nprj=1
   call optimize_channel(bb%ae,inp,l1,nprj,irc(l1),bb%chs(l1),iout,ifail)
   if(ifail/=0) error stop 'test_ps_physics: optimize_channel failed'
 end do
 end subroutine build

 function tag(what,l1,ikap,ii) result(str)
!! the label of a check of l = l1-1, kappa ikap and projector ii (0 for all)
 character(len=*), intent(in) :: what
 integer, intent(in) :: l1,ikap,ii
 character(len=64) :: str
 if(ii>0) then
   write(str,'(2a,i0,a,i0,a,i0)') what,' l=',l1-1,' ikap=',ikap,' i=',ii
 else
   write(str,'(2a,i0,a,i0)') what,' l=',l1-1,' ikap=',ikap
 end if
 end function tag

 subroutine check_norm(bb)
 type(build_t), intent(in) :: bb
 integer :: l1,ikap,ii,jj
 real(dp) :: sps,sae,err
 do l1=1,inp%lmax+1
   associate(ch=>bb%chs(l1))
     do ikap=1,ch%mkap
       err=0.0d0
       do ii=1,ch%nproj
         do jj=ii,ch%nproj
           call lm_integral(bb%ae%rr,ch%pswf(:,ii,ikap)*ch%pswf(:,jj,ikap), &
&                           bb%mmax,2*l1,0.0d0,inp%rc(l1),sps)
           call lm_integral(bb%ae%rr,ch%uua(:,ii,ikap)*ch%uua(:,jj,ikap), &
&                           bb%mmax,2*l1,0.0d0,inp%rc(l1),sae)
           err=max(err,abs(sps-sae))
         end do
       end do
       call check(tag('norm',l1,ikap,0),err,tol_norm)
     end do
   end associate
 end do
 end subroutine check_norm

 subroutine check_matching(bb)
 type(build_t), intent(in) :: bb
 integer :: l1,ikap,ii
 real(dp) :: err
 do l1=1,inp%lmax+1
   associate(ch=>bb%chs(l1))
     do ikap=1,ch%mkap
       err=0.0d0
       do ii=1,ch%nproj
         err=max(err,maxval(abs(ch%pswf(:,ii,ikap)-ch%uua(:,ii,ikap)), &
&                           mask=bb%ae%rr>=inp%rc(l1)))
       end do
       call check(tag('match',l1,ikap,0),err,0.0d0)
     end do
   end associate
 end do
 end subroutine check_matching

 subroutine check_charge(bb)
 type(build_t), intent(in) :: bb
 real(dp) :: qq
 associate(pot=>bb%psp%pot)
   call lm_integral(pot%rr,pot%rho*pot%rr**2,bb%mmax,2,0.0d0,pot%rr(bb%mmax),qq)
   call check('charge',abs(qq-pot%zion),tol_charge)
 end associate
 end subroutine check_charge

 subroutine check_density(bb)
 type(build_t), intent(in) :: bb
 associate(pot=>bb%psp%pot)
   call check('density beyond rcmax', &
&             maxval(abs(pot%rho-bb%psp%rhotae),mask=pot%rr>=pot%rcmax) &
&             /maxval(bb%psp%rhotae),tol_density)
 end associate
 end subroutine check_density

 subroutine check_tail(bb)
 type(build_t), intent(in) :: bb
 integer :: kk
 real(dp) :: qc,qout,rv,dq,err,r0
 associate(pot=>bb%psp%pot,rr=>bb%psp%pot%rr,mmax=>bb%mmax)
   call lm_integral(rr,bb%ae%rhoc*rr**2,mmax,2,0.0d0,rr(mmax),qc)
   dq=abs(qc-(pot%zz-pot%zion))
   r0=max(2.0d0*pot%rcmax,bb%psp%rcore)
   err=0.0d0
   do kk=mmax,1,-1
     if(rr(kk)<r0) exit
     call lm_integral(rr,bb%ae%rhoc*rr**2,mmax,2,rr(kk),rr(mmax),qout)
! beyond the switch of ionic_tail_m for irel 1, 2
     if(irel/=0 .and. qout>1.0d-10*rr(kk)) exit
     rv=rr(kk)*pot%vpuns(kk,inp%lloc+1)+pot%zion
     err=max(err,rv-dq,-dq-qout-rv)
   end do
   call check('ionic tail',err,tol_tail)
 end associate
 end subroutine check_tail

 subroutine check_local(bb)
 type(build_t), intent(in) :: bb
 integer :: l1,ikap,ii
 real(dp) :: err,ff(bb%mmax)
 associate(pot=>bb%psp%pot)
   do l1=1,min(inp%lmax+1,mxl+1)
     if(l1==inp%lloc+1) cycle
     do ikap=1,bb%chs(l1)%mkap
       err=0.0d0
       do ii=1,inp%nproj(l1)
         call pw_node_values(pot%rr,bb%mmax,pot%vkbpw(ii,l1,ikap),ff)
         err=max(err,maxval(abs(ff),mask=pot%rr>inp%rc(l1)), &
&                maxval(abs(pot%vkb(:,ii,l1,ikap)),mask=pot%rr>inp%rc(l1)))
       end do
       call check(tag('vkb beyond rc',l1,ikap,0),err,0.0d0)
     end do
   end do
 end associate
 end subroutine check_local

 subroutine check_core(bb)
 type(build_t), intent(in) :: bb
 if(inp%icmod==0) return
 associate(pot=>bb%psp%pot)
   call check('core beyond r_core',maxval(abs(pot%rhomod(:,1)-bb%ae%rhoc), &
&             mask=pot%rr>=bb%psp%rcore),0.0d0)
   call check('core >= 0',max(0.0d0,-minval(pot%rhomod(:,1))),0.0d0)
 end associate
 end subroutine check_core

 subroutine check_curves(bb)
 type(build_t), intent(in) :: bb
 integer :: l1,ikap,ii,kk
 real(dp) :: err
 do l1=1,inp%lmax+1
   do ikap=1,bb%chs(l1)%mkap
     do ii=1,bb%chs(l1)%nproj
       associate(cc=>bb%psp%tests%convergence%curve(:,:,ii,l1,ikap))
         err=0.0d0
         do kk=2,size(cc,2)
           if(cc(1,kk)<=0.0d0) exit
           err=max(err,cc(2,kk)-cc(2,kk-1))
         end do
         call check(tag('curve',l1,ikap,ii),err/cc(2,1),tol_curve)
       end associate
     end do
   end do
 end do
 end subroutine check_curves

 subroutine check_so(bb)
 type(build_t), intent(in) :: bb
 integer :: l1,ikap,ii,kk,ll,nn
 real(dp) :: fso
 real(dp), allocatable :: kj(:,:),ks(:,:),ff(:)
 associate(pot=>bb%psp%pot,mmax=>bb%mmax)
   allocate(ff(mmax))
   do l1=2,min(inp%lmax+1,mxl+1)
     if(l1==inp%lloc+1) cycle
     ll=l1-1
! the kernels on the mesh as far as the projectors reach
     nn=count(pot%rr<=inp%rc(l1))+2
     allocate(kj(nn,nn),ks(nn,nn))
     do ikap=1,2
       fso=0.5d0*ll
       if(ikap==2) fso=-0.5d0*(ll+1)
       kj(:,:)=0.0d0
       do ii=1,inp%nproj(l1)
         call pw_node_values(pot%rr,mmax,pot%vkbpw(ii,l1,ikap),ff)
         do kk=1,nn
           kj(:,kk)=kj(:,kk)+pot%evkb(ii,l1,ikap)*ff(1:nn)*ff(kk)
         end do
       end do
       ks(:,:)=0.0d0
       do ii=1,2*mxprj
         if(pot%esr(ii,l1)/=0.0d0) then
           call pw_node_values(pot%rr,mmax,pot%vsrpw(ii,l1),ff)
           do kk=1,nn
             ks(:,kk)=ks(:,kk)+pot%esr(ii,l1)*ff(1:nn)*ff(kk)
           end do
         end if
         if(pot%eso(ii,l1)/=0.0d0) then
           call pw_node_values(pot%rr,mmax,pot%vsopw(ii,l1),ff)
           do kk=1,nn
             ks(:,kk)=ks(:,kk)+fso*pot%eso(ii,l1)*ff(1:nn)*ff(kk)
           end do
         end if
       end do
       call check(tag('so kernel',l1,ikap,0), &
&                 maxval(abs(kj-ks))/maxval(abs(kj)),tol_so)
     end do
     deallocate(kj,ks)
   end do
 end associate
 end subroutine check_so

 subroutine residuals(bb,rs)
!! the signed residuals of the mesh-limited invariants, in a fixed order
 type(build_t), intent(in) :: bb
 type(resid_t), intent(out) :: rs
 integer, parameter :: np=5
 real(dp), parameter :: de=1.0d-4
 integer :: l1,ikap,ii,ll,nvkb,ifail,kk
 real(dp) :: pae(np,1),pps(np,1),rq,dae,dps

 associate(cf=>bb%psp%tests%configs)
   rs%eig=[(cf%ee(3,kk,1)-cf%ee(2,kk,1),kk=1,cf%nrow(1))]
 end associate
 allocate(rs%ld(0),rs%dld(0))
 do l1=1,inp%lmax+1
   ll=l1-1
! a radius of both meshes' interpolants, the same for both
   rq=inp%rc(l1)*1.006d0**2
   nvkb=inp%nproj(l1)
   if(ll==inp%lloc) nvkb=0
   associate(ch=>bb%chs(l1),pot=>bb%psp%pot)
     do ikap=1,ch%mkap
       do ii=1,ch%nproj
! energies epa+2de, epa+de, epa, epa-de, epa-2de
         call phsft_pair(irel,bb%ae%c_light,ll,kap_of(ll,ikap),nvkb,bb%ae%rr, &
&                        bb%ae%vfull,bb%ae%zz,pot%vloc(:,ikap),pot%vlocpw, &
&                        pot%vkbpw(1:nvkb,l1,ikap),pot%evkb(1:nvkb,l1,ikap), &
&                        bb%mmax,bb%ae%adams,1,[rq],np,ch%epa(ii,ikap)+2*de, &
&                        de,ch%epa(ii,ikap),pae,pps,iout,ifail)
         if(ifail/=0) error stop 'test_ps_physics: phsft_pair failed'
         dae=(-pae(1,1)+8*pae(2,1)-8*pae(4,1)+pae(5,1))/(12*de)
         dps=(-pps(1,1)+8*pps(2,1)-8*pps(4,1)+pps(5,1))/(12*de)
         rs%ld=[rs%ld,pps(3,1)-pae(3,1)]
         rs%dld=[rs%dld,(dps-dae)/dae]
       end do
     end do
   end associate
 end do
 end subroutine residuals

 subroutine check_mesh(what,res1,res2,floor)
!! the residuals res2 within the mesh error estimated from res1 and res2
 character(len=*), intent(in) :: what
 real(dp), intent(in) :: res1(:),res2(:),floor
 real(dp) :: err
 err=maxval(abs(res2)-2.0d0*abs(res1-res2)/15.0d0)
 call check(what//' (less the mesh error)',max(err,0.0d0),floor)
 end subroutine check_mesh

 end program test_ps_physics
