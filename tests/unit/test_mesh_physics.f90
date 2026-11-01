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
!> Unit test of the mesh convergence of the self-consistent all-electron
!> atom (aeatom; test_adams has the single states): its total energy and
!> eigenvalues q converge at the Adams order, 4, as the log step al -> 0.
!> On the meshes of amesh 1.012, 1.006 and 1.003 (the built-in first point
!> and end), the ratio of the successive differences
!>   rho = (q1 - q2)/(q2 - q3)
!> is (al1**4 - al2**4)/(al2**4 - al3**4) = 15.80 for an error C al**4;
!> checked: |log2(rho/15.80)| <= 0.5, the order to 0.5 (a 3rd- or 5th-order
!> error is 1 off). The next order changes rho by a relative O(al) with a
!> coefficient of ~10 for deep states (Au 2s: rho 13.7-14.4 at al 0.012),
!> and more where the mesh error is small, a near cancellation of terms of
!> both signs (O 2s, Si 3s, Cu 3p: q2 - q3 down to 7e-11 Ha, rho up to 20,
!> 0.36 in the order); the noise of the eigenvalues, ~1e-12 Ha, is under 2%
!> of the smallest q2 - q3.
!> O, Si, Cu and Au, nr, sr and fr, VWN and PBE. Not PZ (iexc 3): its
!> energy steps at r_s = 1 (rounded coefficients), so its total energies
!> jump by up to ~5e-7 Ha as the mesh moves the nodes across the radius of
!> r_s = 1, and the differences have no order (rho from -7 to 17)
!> (doc/physics_tests_plan.md, Step 5).
!> Stops with a non-zero code if any check failed.

 program test_mesh_physics
 use constants_m, only: dp, mxstate
 use unit_check_m, only: check, check_done
 use aeatom_m, only: aeatom
 use adams_m, only: adams_t, adams_init
 use main_setup_m, only: make_log_mesh
 use failure_m, only: warnings_t
 use relativity_m, only: c_light_default
 implicit none

!> an atom of the table: its configuration, the first nc states the core
 type :: atom_t
   character(len=6) :: name=' '
   real(dp) :: zz=0.0_dp
   integer :: nc=0,ncv=0
   integer :: na(mxstate)=0,la(mxstate)=0
   real(dp) :: fa(mxstate)=0.0_dp
 end type atom_t

 integer, parameter :: natom=4
 real(dp), parameter :: ameshes(3)=[1.012_dp,1.006_dp,1.003_dp]
 real(dp), parameter :: tol_order=0.5_dp
 integer, parameter :: iexcs(2)=[6,4]
 type(atom_t) :: atoms(natom)
 type(adams_t) :: adams
 real(dp) :: etot(3),ea(mxstate,2,3),al(3),rho4
 integer :: iout,ia,irel,ix,im,ii,ik
 character(len=64) :: lab

 call table(atoms)
 adams=adams_init(4)
 open(newunit=iout,file='/dev/null',action='write')
 do ia=1,natom
   do ix=1,size(iexcs)
     do irel=0,2
       do im=1,3
         call solve(atoms(ia),irel,iexcs(ix),ameshes(im),al(im),etot(im), &
&                   ea(:,:,im))
       end do
       rho4=(al(1)**4-al(2)**4)/(al(2)**4-al(3)**4)
       write(lab,'(a,1x,a,i0,a,i0)') trim(atoms(ia)%name),'irel=',irel, &
&        ' iexc=',iexcs(ix)
       call check_order(trim(lab)//' etot',etot)
       do ii=1,atoms(ia)%ncv
         do ik=1,merge(2,1,irel==2 .and. atoms(ia)%la(ii)>0)
           write(lab,'(a,1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
&            trim(atoms(ia)%name),'irel=',irel,' iexc=',iexcs(ix),' n=', &
&            atoms(ia)%na(ii),' l=',atoms(ia)%la(ii),' ikap=',ik
           call check_order(trim(lab),ea(ii,ik,:))
         end do
       end do
     end do
   end do
 end do
 close(iout)
 call check_done('mesh_physics')

 contains

 subroutine table(atoms)
!! the atoms of the tests, with their reference configurations
 type(atom_t), intent(out) :: atoms(natom)
 atoms(1)=atom('O',8.0_dp,1,[1,2,2],[0,0,1],[2.0_dp,2.0_dp,4.0_dp])
 atoms(2)=atom('Si',14.0_dp,3,[1,2,2,3,3],[0,0,1,0,1], &
&              [2.0_dp,2.0_dp,6.0_dp,2.0_dp,2.0_dp])
 atoms(3)=atom('Cu',29.0_dp,3,[1,2,2,3,3,3,4],[0,0,1,0,1,2,0], &
&              [2.0_dp,2.0_dp,6.0_dp,2.0_dp,6.0_dp,10.0_dp,1.0_dp])
 atoms(4)=atom('Au',79.0_dp,10, &
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

 subroutine solve(at,irel,iexc,amesh,al,etot,ea)
!! the self-consistent atom at on the built-in mesh of ratio amesh: its
!! log step, total energy and eigenvalues
 type(atom_t), intent(in) :: at
 integer, intent(in) :: irel,iexc
 real(dp), intent(in) :: amesh
 real(dp), intent(out) :: al,etot,ea(mxstate,2)
 type(warnings_t) :: warns
 real(dp), allocatable :: rr(:),rho(:),rhoc(:),vi(:),vext(:)
 real(dp) :: log_mesh(4),eterm(4),res
 integer :: nkap,it,ierr,ifail,mmax
 call make_log_mesh(irel,at%zz,amesh,-1.0_dp,-1.0_dp,4,al,mmax,rr,log_mesh)
 nkap=1
 if(irel==2) nkap=2
 allocate(rho(mmax),rhoc(mmax),vi(mmax),vext(mmax))
 vi(:)=0.0_dp
 vext(:)=0.0_dp
 ea(:,:)=0.0_dp
 ierr=0
 call aeatom(irel,c_light_default,nkap,at%na,at%la,ea(:,1:nkap), &
&            at%fa(1:at%ncv),at%nc,at%ncv,it,rhoc,rho,rr,vi, &
&            at%zz,mmax,adams,iexc,etot,eterm,res,ierr,vext, &
&            .false.,warns,iout,ifail)
 if(ifail/=0 .or. ierr/=0) error stop 'test_mesh_physics: aeatom failed'
 end subroutine solve

 subroutine check_order(what,qq)
!! the order of the differences of qq on the three meshes, against 4
 character(len=*), intent(in) :: what
 real(dp), intent(in) :: qq(3)
 real(dp) :: rr,err
 rr=(qq(1)-qq(2))/(qq(2)-qq(3))
 err=huge(1.0_dp)
 if(rr>0.0_dp) err=abs(log(rr/rho4)/log(2.0_dp))
 call check(what//' order - 4',err,tol_order)
 end subroutine check_order

 end program test_mesh_physics
