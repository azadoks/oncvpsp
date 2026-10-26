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
!> The atom compressed of oncvpsp_confine.x (whose header describes the
!> wall lines): the self-consistent all-electron atom (aeatom) and
!> pseudo-atom (psatom) of the reference configuration in a smooth wall
!> beyond r0, for each strength from the free atom.
 module confine_wall_m
 use constants_m, only: dp, mxstate, mxprj, mxl, mxlin
 use ps_file_m, only: ps_file_t
 use confine_summary_m, only: confine_summary_t
 implicit none
 private
 public :: confine_walls

!> the onset of the wall, x**2 exp(-dwall/x): smooth at r0, so that outwf_m
!> takes it as part of the local potential (bohr)
 real(dp), parameter :: dwall=0.5d0

 contains

 subroutine confine_walls(pf,irel,zz,mmax,rr,vfull,nc,nv,iexc,na,la,fa,ea, &
&                         nk,kwall,r0,sm)
!! the wall lines of the strengths kwall(1:nk) beyond r0, kwall(0)=0 the
!! free atom, and the walls' part of the summary
 use aeatom_m, only: aeatom
 use psatom_m, only: psatom
 type(ps_file_t), intent(in) :: pf
 integer, intent(in) :: irel,mmax,nc,nv,iexc,na(mxstate),la(mxstate),nk
 real(dp), intent(in) :: zz,rr(mmax),vfull(mmax),fa(mxstate)
 real(dp), intent(in) :: ea(mxstate,2),r0
 real(dp), intent(inout) :: kwall(0:nk)
 type(confine_summary_t), intent(inout) :: sm
 integer :: na_ps(mxstate),la_ps(mxstate),nprj(mxlin+1)
 integer :: ik,ii,ierrae,ierrps,itae,itps,nkap
 real(dp) :: fa_ps(mxstate,2),etae(0:nk),etps(0:nk)
 real(dp) :: eaw(mxstate,2),rpk(mxstate,2),eterm(4),res,ea_ps(mxstate,2)
 real(dp), allocatable :: vext(:),rhoc(:),rhoa(:),vi(:),vpuns(:,:),rhops(:)
 real(dp), allocatable :: eigae(:,:,:),eigps(:,:,:)

 kwall(0)=0.0d0
 nkap=pf%nkap
 allocate(vext(mmax),rhoc(mmax),rhoa(mmax),vi(mmax),rhops(mmax))
 allocate(vpuns(mmax,mxl+2),eigae(mxstate,2,0:nk),eigps(mxstate,2,0:nk))
 eigae(:,:,:)=0.0d0
 eigps(:,:,:)=0.0d0
 vpuns(:,:)=0.0d0
 vpuns(:,pf%lloc+1)=pf%vpuns(:)
 nprj(:)=0
 do ii=1,4
   nprj(ii)=pf%nvkb(ii,1)
 end do
 call valence_states(nc,nv,la,fa,na_ps,la_ps,fa_ps)
 write(6,'(a)') '# wall_r0 r0(bohr)'
 write(6,'(2a)') '# wall k(Ha/bohr**2) etot_ae(Ha) etot_ps(Ha) de_ae(Ha)', &
&     ' de_ps(Ha)'
 write(6,'(a)') '# wall_eig k n l ikap e_ae(Ha) e_ps(Ha)'
 write(6,'(a)') '# wall_failed k atom'
 write(6,'(a,1p,e24.16)') 'wall_r0',r0
 sm%nwfail=0
 ierrae=0
 ierrps=0
 do ik=0,nk
   do ii=1,mmax
     vext(ii)=wall(kwall(ik),rr(ii)-r0)
   end do
! each atom from the wall before (the free atom from the data file's
! potential and the run's pseudo-density): the walls ramp the atom in
   if(ik==0) then
     eaw(:,:)=ea(:,:)
     vi(:)=vfull(:)
   else if(ierrae/=0) then
     eaw(:,:)=eigae(:,:,0)
     vi(:)=vfull(:)
   end if
   rpk(:,:)=0.0d0
   ierrae=0
   call aeatom(irel,nkap,na,la,eaw(:,1:nkap),fa,rpk(:,1:nkap),nc,nc+nv, &
&              itae,rhoc,rhoa,rr,vi,zz,mmax,iexc,etae(ik),eterm,res,ierrae, &
&              vext,.true.)
   eigae(:,:,ik)=eaw(:,:)
   if(ik==0 .or. ierrps/=0) then
     rhops(:)=pf%rhov(:)
     ea_ps(1:nv,:)=ea(nc+1:nc+nv,:)
   end if
   ierrps=0
   call psatom(na_ps,la_ps,ea_ps(:,1:nkap),fa_ps,nv,itps,pf%rhomod,rhops, &
&              rr,pf%rcmax,mmax,mxprj,iexc,etps(ik),nprj,vpuns,pf%lloc, &
&              pf%vlocpw,pf%vkbpw(:,:,1:nkap),pf%evkb(:,:,1:nkap),ierrps, &
&              nkap,vext)
   eigps(1:nv,:,ik)=ea_ps(1:nv,:)
   if(ierrae/=0) write(6,'(a,1p,e24.16,a)') 'wall_failed',kwall(ik),' ae'
   if(ierrps/=0) write(6,'(a,1p,e24.16,a)') 'wall_failed',kwall(ik),' ps'
   if(ik==0 .and. (ierrae/=0 .or. ierrps/=0)) then
     write(6,'(a)') 'oncvpsp_confine: ERROR the free atom failed'
     error stop 1
   end if
   if(ierrae/=0 .or. ierrps/=0) then
     sm%nwfail=sm%nwfail+1
     cycle
   end if
   write(6,'(a,1p,5e24.16)') 'wall',kwall(ik),etae(ik),etps(ik), &
&       etae(ik)-etae(0),etps(ik)-etps(0)
   if(ik>0) sm%sdde=max(sm%sdde,abs(etps(ik)-etps(0)-etae(ik)+etae(0)))
   call write_wall_eig(nc,nv,na,la,nkap,kwall(ik),eigae(:,:,ik), &
&                      eigps(:,:,ik))
 end do
 deallocate(vext,rhoc,rhoa,vi,rhops,vpuns,eigae,eigps)
 end subroutine confine_walls

 subroutine valence_states(nc,nv,la,fa,na_ps,la_ps,fa_ps)
!! the valence states as run_config_ps numbers them, n from l+1
 integer, intent(in) :: nc,nv,la(mxstate)
 real(dp), intent(in) :: fa(mxstate)
 integer, intent(out) :: na_ps(mxstate),la_ps(mxstate)
 real(dp), intent(out) :: fa_ps(mxstate,2)
 integer :: nav(4),kk,l1
 do l1=1,4
   nav(l1)=l1-1
 end do
 do kk=1,nv
   la_ps(kk)=la(nc+kk)
   l1=la_ps(kk)+1
   nav(l1)=nav(l1)+1
   na_ps(kk)=nav(l1)
   fa_ps(kk,:)=fa(nc+kk)
 end do
 end subroutine valence_states

 subroutine write_wall_eig(nc,nv,na,la,nkap,kk,eigae,eigps)
!! the wall_eig lines of the strength kk: the all-electron valence
!! eigenvalues eigae(nc+1:nc+nv) and the pseudo-atom's eigps(1:nv)
 integer, intent(in) :: nc,nv,na(mxstate),la(mxstate),nkap
 real(dp), intent(in) :: kk,eigae(mxstate,2),eigps(mxstate,2)
 integer :: jj,ll,ikap
 do jj=1,nv
   ll=la(nc+jj)
   do ikap=1,merge(1,nkap,ll==0)
     write(6,'(a,1p,e24.16,3i4,2e24.16)') 'wall_eig',kk,na(nc+jj),ll, &
&         ikap,eigae(nc+jj,ikap),eigps(jj,ikap)
   end do
 end do
 end subroutine write_wall_eig

 real(dp) function wall(kk,xx)
!! the wall, 1/2 kk xx**2 exp(-dwall/xx) for xx > 0, 0 inside
 real(dp), intent(in) :: kk,xx
 wall=0.0d0
 if(xx>0.0d0) wall=0.5d0*kk*xx**2*exp(-dwall/xx)
 end function wall

 end module confine_wall_m
