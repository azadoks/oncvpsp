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
!> The all-electron states of a channel (l, and kappa when fully
!> relativistic) that oncvpsp_main builds its projectors from: the bound
!> valence states, then the well states, each with its energy epa and n
!> npa; the relativistic term Delta about the taper window; their overlap
!> inside rc; and the weights of their pseudo-waves in the cutoff.
 module ae_projectors_m
 use constants_m, only: dp, mxprj, mxstate, mxkap
 use pieces_m, only: pw_reach
 use free_tail_m, only: free_tail_log
 use well_mesh_m, only: well_mesh_t
 use logmesh_m, only: lm_fh_reach
 implicit none
 private
 public :: ae_bound_projectors, ae_well_projectors, ae_delta_window
 public :: ae_overlap, ecut_weights

 contains

 subroutine ae_bound_projectors(irel,c_light,ll,kap,ikap,nc,nv,na,la,ea,nproj, &
&                               it,rr,zz,vfull,mmax,uu,up,iprj,epa,npa,uua, &
&                               upa,iout,ifail)
!! the bound valence states of l (at most nproj), the first iprj
!! projectors' states; and npa(1), the n above the highest core state of l
 use boundstate_m, only: boundstate
 use renorm_r_m, only: renorm_r
 integer, intent(in) :: irel,ll,kap,ikap,nc,nv,na(mxstate),la(mxstate)
 real(dp), intent(in) :: c_light
 integer, intent(in) :: nproj,it,mmax
 real(dp), intent(in) :: ea(mxstate,mxkap),rr(mmax),zz,vfull(mmax)
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 integer, intent(out) :: iprj
 real(dp), intent(inout) :: epa(mxprj),uua(mmax,mxprj),upa(mmax,mxprj)
 integer, intent(inout) :: npa(mxprj)
 real(dp) :: et,cnorm
 integer :: kk,ierr,mch
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

 ifail=0
 iprj=0
!get principal quantum number for the highest core state for this l
 npa(1)=ll+1
 do kk=1,nc
   if(la(kk)==ll) npa(1)=na(kk)+1
 end do !kk

!get all-electron bound states for projectors
 if(nv/=0) then
   do kk=nc+1,nc+nv
     if(la(kk)==ll) then
       iprj=iprj+1
       et=ea(kk,ikap)
       call boundstate(irel,c_light,na(kk),ll,kap,ierr,et, &
&                      rr,zz,vfull,uu,up,mmax,mch,iout,ifail)
       if(ifail/=0) return
       if(ierr /= 0) then
         if(irel==2) then
           write(iout,'(a,i2,a,i2,a,i2)') &
&                'oncvpsp_r 409: ldiracfb convergence ERROR, n=',&
&                 na(kk),'  l=',ll, 'kap-',kap
         else
           write(iout,'(/a,3i4)') &
&                'oncvpsp-387: lschfb convergence ERROR n,l,iter=', &
&                na(kk),la(kk),it
         end if
         ifail=1
         return
       end if

       epa(iprj)=ea(kk,ikap)
       npa(iprj)=na(kk)
       if(irel==2) then
         call renorm_r(uu,rr,ll,kap,c_light,zz,mmax,cnorm,iout,ifail)
         if(ifail/=0) return
       else
         cnorm=1.0d0
       end if
       uua(:,iprj)=uu(:,1)
       upa(:,iprj)=cnorm*up(:,1)
     end if !la(kk)==ll
     if(iprj==nproj) exit
   end do !kk
 end if !nv/=0
 end subroutine ae_bound_projectors

 subroutine ae_well_projectors(irel,c_light,ll,kap,ikap,nproj,ep,debl,well_n, &
&                              ep_as_given,well_onset,well_tail, &
&                              well_norm_rc,rc,irc,qcut,rr,zz,vfull,mmax, &
&                              uu,up,nnwell,iprj,epa,npa,uua,upa,vwell, &
&                              inwell,wmesh,iout,ifail)
!! the well states of the projectors iprj+1 to nproj: the first at ep if
!! there are no bound ones, the next ones debl above; their wells are
!! those of ikap=1, kept for ikap=2 (vwell, wmesh). nnwell is the n of the
!! last well state of ikap=1 that was not fully relativistic
 use fpovlp_m, only: fpovlp
 use boundstate_m, only: boundstate
 use renorm_r_m, only: renorm_r
 use wellstate_m, only: wellstate
 integer, intent(in) :: irel,ll,kap,ikap,nproj,well_n,well_tail,irc,mmax
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ep,debl,well_onset,rc,qcut,rr(mmax),zz
 real(dp), intent(in) :: vfull(mmax)
 logical, intent(in) :: ep_as_given,well_norm_rc
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 integer, intent(inout) :: nnwell,iprj,npa(mxprj)
 real(dp), intent(inout) :: epa(mxprj),uua(mmax,mxprj),upa(mmax,mxprj)
 real(dp), intent(inout) :: vwell(mmax,mxprj)
 logical, intent(inout) :: inwell(mxprj)
 type(well_mesh_t), intent(inout) :: wmesh(mxprj)
 real(dp) :: et,cnorm,eptgt
 integer :: kk,ierr,mch,mchx
 logical :: well_tail_free
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

 ifail=0
 well_tail_free=well_tail==1
!get all-electron well states for projectors
!if there were no valence states, use ep from input data for 1st well state
!otherwise shift up by input debl

 if(iprj==0) epa(1)=ep
 if(iprj>=nproj) return
 do kk=1,nproj-iprj
   iprj=iprj+1
   if(iprj>1 .and. debl<=0.0d0) then
     if(irel==2) then
       write(iout,'(a,f8.3,a,i2/a/a)') 'oncvpsp: ERROR debl =',debl, &
&            ' for l=',ll, &
&            ' not allowed with 2 or more scattering states', &
&            'program will stop'
     else
       write(iout,'(a,f8.3,a,i2/a/a)') 'oncvpsp: ERROR debl =',debl, &
&            ' for l=',ll, &
&            ' ERROR not allowed with 2 or more scattering states', &
&            'program will stop'
     end if
     ifail=1
     return
   end if

! npa is not reset for ikap=2 to retain same node count for consistency with
! keeping same well potential as for ikap=1.
   if(iprj>1) then
     epa(iprj)=epa(iprj-1)+debl
     if(ikap==1) then
       npa(iprj)=npa(iprj-1)+1
     end if
   end if
!with the confined well, the n of the channel's first well state is the
!input's well_n, if set; the next ones follow by one each
   if(well_tail==2 .and. kk==1 .and. ikap==1 .and. well_n>0) then
     if(well_n<npa(iprj)) then
       write(iout,'(/a,i2,a,i3,a,i3)') 'oncvpsp: ERROR well_n for l=',ll, &
&           ' is',well_n,', below the first well state''s n,', &
&           npa(iprj)
       ifail=1
       return
     end if
     npa(iprj)=well_n
   end if

! the well potential is calculated for the less-bound ikap=1 state and kept
! for ikap=2, which is more deeply bound and should always have a bound
! state for this combination of npa and vwell.  For irel=2 npa takes the n
! of the well state found, which can have more nodes than asked for; else
! npa keeps the input n.
   if(ikap==1 .and. irel==2) then
     call wellstate(npa(iprj),ll,kap,rc,epa(iprj), &
&                   ep_as_given,well_onset,well_tail==2,.false.,rr, &
&                   vfull, &
&                   vwell(1,iprj),uu,up,zz,mmax,mch,irel,c_light, &
&                   inwell(iprj),eptgt,wmesh(iprj),iout,ifail)
     if(ifail/=0) return
   else if(ikap==1) then
     nnwell=npa(iprj)
     call wellstate(nnwell,ll,kap,rc,epa(iprj),ep_as_given, &
&                   well_onset,well_tail==2,well_tail_free,rr,vfull, &
&                   vwell(1,iprj),uu,up,zz,mmax,mch,irel,c_light, &
&                   inwell(iprj),eptgt,wmesh(iprj),iout,ifail)
     if(ifail/=0) return
   else
     et=epa(iprj)
     if(wmesh(iprj)%mx>0) then
!the j = l - 1/2 state in the well of j = l + 1/2, on its continued mesh
       wmesh(iprj)%uu(:,:)=0.0d0
       wmesh(iprj)%up(:,:)=0.0d0
       call boundstate(irel,c_light,npa(iprj),ll,kap,ierr, &
&                      epa(iprj),wmesh(iprj)%rr,zz, &
&                      wmesh(iprj)%vv,wmesh(iprj)%uu,wmesh(iprj)%up, &
&                      wmesh(iprj)%mx,mchx,iout,ifail)
       if(ifail/=0) return
       uu(:,:)=wmesh(iprj)%uu(1:mmax,:)
       up(:,:)=wmesh(iprj)%up(1:mmax,:)
     else
       call boundstate(irel,c_light,npa(iprj),ll,kap,ierr, &
&                      epa(iprj),rr,zz,vwell(1,iprj),uu,up, &
&                      mmax,mch,iout,ifail)
       if(ifail/=0) return
     end if

     if(ierr /= 0) then
       write(iout,'(a,i2,a,i2,a,i3,a,f12.6)') &
&            'oncvpsp_r 474: ldiracfb convergence ERROR, n=',&
&             npa(iprj),'  l=',ll, 'kap=',kap, ' epa=', &
&             epa(iprj)

       ifail=1
       return
     end if

     write(iout,'(/a,i3,a,i3,a,i3/a,f8.4)') &
&             '   Wellstate for l =',ll,'  n =',npa(iprj), &
&             ' kap=',kap, &
&             '     eigenvalue = ',et
   end if

   if(irel==2 .and. wmesh(iprj)%mx>0) then
     call renorm_r(wmesh(iprj)%uu,wmesh(iprj)%rr,ll,kap,c_light,zz, &
&                  wmesh(iprj)%mx,cnorm,iout,ifail)
     if(ifail/=0) return
     uu(:,:)=wmesh(iprj)%uu(1:mmax,:)
   else if(irel==2) then
     call renorm_r(uu,rr,ll,kap,c_light,zz,mmax,cnorm,iout,ifail)
     if(ifail/=0) return
   else
     cnorm=1.0d0
   end if
   uua(:,iprj)=uu(:,1)
   upa(:,iprj)=cnorm*up(:,1)
!a state the well confines is normalized inside rc with well_norm "rc",
!and by default with the confined well, as with the free tail: its norm
!over all space is the well's, and would scale its E_resid with the width
!of the well; else over all space, as the solvers give it (the kernels
!see neither)
   if(well_norm_rc .and. wmesh(iprj)%mx>0) then
     call fpovlp(uua(1,iprj),uua(1,iprj),mmax,rc,ll, &
&                kap,zz,cnorm,rr,irel,c_light)
     cnorm=1.0d0/sqrt(cnorm)
     uua(:,iprj)=cnorm*uua(:,iprj)
     upa(:,iprj)=cnorm*upa(:,iprj)
     wmesh(iprj)%uu(:,:)=cnorm*wmesh(iprj)%uu(:,:)
     wmesh(iprj)%up(:,:)=cnorm*wmesh(iprj)%up(:,:)
   end if
   if(well_tail_free .and. inwell(iprj)) then
     if(ikap==1) epa(iprj)=eptgt
     call free_tail_state(irel,c_light,ll,kap,ikap,nnwell,npa(iprj),epa(iprj), &
&                         rc,irc,qcut,rr,zz,vfull,mmax,uu,up,uua(:,iprj), &
&                         upa(:,iprj),iout,ifail)
     if(ifail/=0) return
   end if
 end do !kk
 end subroutine ae_well_projectors

 subroutine free_tail_state(irel,c_light,ll,kap,ikap,nnwell,nn,ee,rc,irc,qcut, &
&                           rr,zz,vfull,mmax,uu,up,ua,pa,iout,ifail)
!! with the free tail, a well state is the atom's solution at exactly the
!! energy asked for (the well's eigenvalue only for j = l - 1/2, D-7c),
!! integrated outward past the reach of everything taken at or near rc
!! (ae_delta's nodes below, and the interpolants'); normalized inside rc,
!! as its norm over all space is the well's; and the windowed free tail
!! beyond, so that nothing sees the well but that energy
!! (doc/smoothness_plan.md, phase 8). Only for ep >= 0 (D-7, 2026-10-23):
!! below 0 the tail grows, and the round-off of E_resid, computed from the
!! transform of the windowed tail, grows with its peak squared, from kappa
!! L ~ 1.5 (free_tail_m)
 use fpovlp_m, only: fpovlp
 use scatterstate_m, only: scatterstate
 integer, intent(in) :: irel,ll,kap,ikap,nn,irc,mmax
 real(dp), intent(in) :: c_light
 integer, intent(inout) :: nnwell
 real(dp), intent(in) :: ee,rc,qcut,rr(mmax),zz,vfull(mmax)
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 real(dp), intent(out) :: ua(mmax),pa(mmax)
 real(dp) :: cnorm
 integer :: mfree,ierr
 integer, intent(in) :: iout
 integer, intent(out) :: ifail
 ifail=0
 if(ikap==1 .and. ee<0.0d0) then
   write(iout,'(/a/a,i2,a,i3,a,f12.6,a/a)') &
&       'oncvpsp: ERROR well_tail = "free" needs ep >= 0:', &
&       'the well state of l=',ll,', n=',nn,' at ep=', &
&       ee,' Ha has a growing tail;', &
&       'raise its ep (debl), or use well_tail = "confined"'// &
&       ' or "well"'
   ifail=1
   return
 end if
 mfree=min(mmax,irc+pw_reach+16+lm_fh_reach()+2)
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 call scatterstate(irel,c_light,nnwell,ll,kap,ierr,ee,rr, &
&                  zz,vfull,uu,up,mmax,mfree,iout)
 if(ierr/=0) then
   write(iout,'(a,i2,a,f12.8)') &
&        'oncvpsp: ERROR scattering state failed for l=',ll, &
&        ' at energy',ee
   ifail=1
   return
 end if
 ua(:)=uu(:,1)
 pa(:)=up(:,1)
 call fpovlp(ua,ua,mmax,rc,ll, &
&            kap,zz,cnorm,rr,irel,c_light)
 ua(:)=ua(:)/sqrt(cnorm)
 pa(:)=pa(:)/sqrt(cnorm)
 call free_tail_log(irel,c_light,ll,kap,ee,rr,vfull, &
&                   ua,pa,mmax, &
&                   rc,qcut,mfree,iout,ifail)
 if(ifail/=0) return
 end subroutine free_tail_state

 subroutine ae_delta_window(irel,c_light,ll,kap,nproj,epa,rr,vfull,uua,upa, &
&                           mmax,al,rc,irc,local,dlta,vrloc,iout,ifail)
!! the relativistic term Delta of the all-electron wave about the taper
!! window [0.95 rc, rc] (ae_delta, D-4), for the projectors; the local
!! potential of lloc<4 (local) takes Delta/(2u) on the window (S-1)
 use ae_delta_m, only: ae_delta
 integer, intent(in) :: irel,ll,kap,nproj,mmax,irc
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: epa(mxprj),rr(mmax),vfull(mmax)
 real(dp), intent(in) :: uua(mmax,mxprj),upa(mmax,mxprj),al,rc
 logical, intent(in) :: local
 real(dp), intent(inout) :: dlta(mmax,mxprj),vrloc(mmax)
 integer :: iprj,ii,jj,kk
 integer, intent(in) :: iout
 integer, intent(out) :: ifail
 ifail=0
 do iprj=1,nproj
   ii=max(1,1+int(log(0.95d0*rc/rr(1))/al)-pw_reach-8)
   jj=min(mmax,irc+pw_reach+16)
   call ae_delta(irel,c_light,ll,kap,epa(iprj),rr,vfull,uua(1,iprj), &
&                upa(1,iprj),mmax,ii,jj,dlta(1,iprj))
   if(local .and. iprj==1) then
     do kk=ii,jj
       if(uua(kk,1)==0.0d0) then
         write(iout,'(a,i2,a,f10.6)') &
&              'oncvpsp: ERROR node of the local channel l=',ll, &
&              ' all-electron wave near rc, at r=',rr(kk)
         ifail=1
         return
       end if
       vrloc(kk)=0.5d0*dlta(kk,1)/uua(kk,1)
     end do
   end if
 end do
 end subroutine ae_delta_window

 subroutine ae_overlap(irel,c_light,ll,kap,nproj,uua,rr,zz,mmax,rc,qq)
!! the overlap matrix of the projectors' all-electron states inside rc
 use fpovlp_m, only: fpovlp
 integer, intent(in) :: irel,ll,kap,nproj,mmax
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: uua(mmax,mxprj),rr(mmax),zz,rc
 real(dp), intent(inout) :: qq(mxprj,mxprj)
 integer :: ii,jj
 do jj=1,nproj
   do ii=1,jj
     call fpovlp(uua(1,ii),uua(1,jj),mmax,rc,ll,kap, &
&                zz,qq(ii,jj),rr,irel,c_light)
     qq(jj,ii)=qq(ii,jj)
   end do
 end do
 end subroutine ae_overlap

 subroutine ecut_weights(ll,nc,nv,la,fa,fj,nproj,ecut_wt,wcurve,iout,ifail)
!! the weights of the pseudo-waves of l in the cutoff of the pseudo-waves
!! (the weighted mean, and which waves the largest is over): the
!! occupation of each valence state (its share fj of J), or for the first
!! wave of a channel with none occupied ecut_wt (ecut_weight), by default
!! 2(2l+1), as if the channel were full; in a solid it holds about the
!! charge that hybridization puts inside rc, an input as the atom cannot
!! give it
 integer, intent(in) :: ll,nc,nv,la(mxstate),nproj
 real(dp), intent(in) :: fa(mxstate),fj,ecut_wt
 real(dp), intent(inout) :: wcurve(mxprj)
 integer :: jj,kk
 integer, intent(in) :: iout
 integer, intent(out) :: ifail
 ifail=0
 jj=0
 do kk=nc+1,nc+nv
   if(la(kk)==ll .and. jj<nproj) then
     jj=jj+1
     wcurve(jj)=fj*fa(kk)
   end if
 end do
 if(.not. any(wcurve(:)>0.0d0)) then
   if(ecut_wt>=0.0d0) then
     wcurve(1)=fj*ecut_wt
   else
     wcurve(1)=fj*dble(2*(2*ll+1))
   end if
 else if(ecut_wt>=0.0d0) then
   write(iout,'(/a,i2,a)') 'oncvpsp: ERROR ecut_weight for l=',ll, &
&       ' is for a channel with no occupied state; this one has one'
   ifail=1
   return
 end if
 end subroutine ecut_weights

 end module ae_projectors_m
