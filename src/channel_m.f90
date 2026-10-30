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
!> The second stage of a run (oncvpsp_run): one angular momentum's
!> channel, from the all-electron atom (ae_state_m) and the channel's own
!> inputs alone. Its projectors' all-electron states (bound and well
!> states), their windowed Delta, and for each the optimized pseudo-wave,
!> its semi-local potential and its convergence curve; both kappas of an
!> l when fully relativistic, as the second kappa's well states take the
!> first's well. No channel depends on another, so a channel can be redone
!> alone for new rc, qcut, ncon, nbas, ep or debl.

 module channel_m
 use constants_m, only: dp, mxprj, mxkap, nerrdg
 use failure_m, only: warnings_t
 implicit none
 private
 public :: channel_t, optimize_channel, channel_key, ncurve, nkey

!> the points of a convergence curve: run_optimize's q grid,
!> 1 + 24 rc max(2 q of the last basis function, 2 qcut, 20), ~2400 for
!> nbas=14 and less than 3000 for qcut up to 15 at rc up to 4; a buffer,
!> with room for qcut rc up to 160
 integer, parameter :: ncurve=8000

!> the size of a channel's key (channel_key)
 integer, parameter :: nkey=28

!> one channel's results; the arrays are per projector and kappa
 type :: channel_t
!> l+1
   integer :: l1=0
!> number of kappas: 1, or 2 for irel=2 and l>0
   integer :: mkap=0
!> number of projectors
   integer :: nproj=0
!> the last log mesh point at or inside rc (run_optimize may move it)
   integer :: irc=0
!> the principal quantum number of each projector's state
   integer :: npa(mxprj)=0
!> whether each projector is a well state
   logical :: wellprj(mxprj)=.false.
!> the q of the last basis function (run_optimize)
   real(dp) :: qmsbf=0.0d0
!> the projectors' energies (Ha)
   real(dp) :: epa(mxprj,mxkap)=0.0d0
!> the all-electron states
   real(dp), allocatable :: uua(:,:,:)
!> their derivatives in the mesh index
   real(dp), allocatable :: upa(:,:,:)
!> Delta of each state about the window (ae_delta_window)
   real(dp), allocatable :: dlta(:,:,:)
!> Delta/(2u) of the first state, for the local channel only
   real(dp), allocatable :: vrloc(:,:)
!> the optimized pseudo-waves
   real(dp), allocatable :: pswf(:,:,:)
!> their continuations past rc
   real(dp), allocatable :: pswfx(:,:,:)
!> (E - T) applied to the pseudo-waves' continuations
   real(dp), allocatable :: chix(:,:,:)
!> the screened semi-local potential (Ha), per kappa
   real(dp), allocatable :: vp(:,:)
!> the projectors before run_vkb
   real(dp), allocatable :: vkb(:,:,:)
!> the cutoffs of the error levels (run_optimize)
   real(dp), allocatable :: cvgplt(:,:,:,:)
!> the convergence curves (run_optimize)
   real(dp), allocatable :: ecurve(:,:,:,:)
!> the curves' tails and checks (run_optimize)
   real(dp), allocatable :: ektail(:,:,:)
!> the optimization's diagnostics (run_optimize)
   real(dp), allocatable :: optdg(:,:,:)
!> the waves' shares (run_optimize)
   real(dp), allocatable :: wshare(:,:,:)
!> the weights of the pseudo-waves in the weighted cutoff (ecut_weights)
   real(dp), allocatable :: wcurve(:,:)
!> the warnings of optimize_channel
   type(warnings_t) :: warn
!> what it was made from (channel_key), so that a run assembled from
!> channels made before takes only those of its own input
   real(dp) :: key(nkey)=0.0d0
 end type channel_t

 contains

 subroutine optimize_channel(ae,inp,l1,nprj,irc,ch,iout,ifail,node_only)

!! the channel of l = l1-1 with nprj projectors (1 for the local
!! channel, whose pseudo-wave the local potential needs): for each kappa,
!! the bound and well states of the projectors, their Delta, the
!! overlaps inside rc, the optimized pseudo-waves (run_optimize) and
!! their weights in the weighted cutoff
!!
!! @param ae  the all-electron atom (ae_state_m)
!! @param inp  the input values (input_m), completed
!! @param l1  l+1
!! @param nprj  its number of projectors
!! @param irc  last mesh point at or inside rc
!! @param ch  the channel
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
!! @param node_only  to stop each kappa after the node test of its first
!!                   projector (pspot), the channel's results unset
!!                   (default false)

 use ae_state_m, only: ae_state_t
 use input_m, only: input_t
 use kappa_m, only: kap_count, kap_of, kap_j_share
 use well_mesh_m, only: well_mesh_t
 use ae_projectors_m, only: ae_bound_projectors, ae_well_projectors, &
&                           ae_delta_window, ae_overlap, ecut_weights
 use failure_m, only: stage_channel
 use run_optimize_m, only: run_optimize

!Input variables
 type(ae_state_t), intent(in) :: ae
 type(input_t), intent(in) :: inp
 integer, intent(in) :: l1,nprj,irc
 integer, intent(in) :: iout
 logical, intent(in), optional :: node_only

!Output variables
 type(channel_t), intent(out) :: ch
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: uu(:,:),up(:,:),vwell(:,:)
 real(dp) :: qq(mxprj,mxprj),fj
 type(well_mesh_t) :: wmesh(mxprj)
 logical :: inwell(mxprj),well_tail_free,well_norm_rc
 integer :: ll,ikap,kap,iprj,nnwell,mmax,nkap

 ifail=0
 mmax=ae%mmax
 nkap=ae%nkap
 ll=l1-1
 well_tail_free=inp%well_tail==1
 well_norm_rc=inp%well_norm==2 .or. (inp%well_norm==0 .and. inp%well_tail/=0)

 ch%l1=l1
 ch%warn%stage=stage_channel
 ch%warn%ll=ll
 ch%mkap=kap_count(nkap,ll)
 ch%nproj=nprj
 ch%irc=irc
 ch%key(:)=channel_key(ae,inp,l1,nprj)
 allocate(ch%uua(mmax,mxprj,nkap),ch%upa(mmax,mxprj,nkap))
 allocate(ch%dlta(mmax,mxprj,nkap),ch%vrloc(mmax,nkap))
 allocate(ch%pswf(mmax,mxprj,nkap),ch%pswfx(mmax,mxprj,nkap))
 allocate(ch%chix(mmax,mxprj,nkap),ch%vp(mmax,nkap),ch%vkb(mmax,mxprj,nkap))
 allocate(ch%cvgplt(2,nerrdg,mxprj,nkap),ch%ecurve(2,ncurve,mxprj,nkap))
 allocate(ch%ektail(2,mxprj,nkap),ch%optdg(4,mxprj,nkap))
 allocate(ch%wshare(3,mxprj,nkap),ch%wcurve(mxprj,nkap))
 ch%uua(:,:,:)=0.0d0
 ch%upa(:,:,:)=0.0d0
 ch%dlta(:,:,:)=0.0d0
 ch%vrloc(:,:)=0.0d0
 ch%pswf(:,:,:)=0.0d0
 ch%pswfx(:,:,:)=0.0d0
 ch%chix(:,:,:)=0.0d0
 ch%vp(:,:)=0.0d0
 ch%vkb(:,:,:)=0.0d0
 ch%cvgplt(:,:,:,:)=0.0d0
 ch%ecurve(:,:,:,:)=0.0d0
 ch%ektail(:,:,:)=0.0d0
 ch%optdg(:,:,:)=0.0d0
 ch%wshare(:,:,:)=0.0d0
 ch%wcurve(:,:)=0.0d0

 allocate(uu(mmax,2),up(mmax,2),vwell(mmax,mxprj))
 up(:,:)=0.0d0
 vwell(:,:)=0.0d0
 inwell(:)=.false.

! loop on J = ll +/- 1/2 when nkap=2; fj is the weight of each J
 do ikap=1,ch%mkap
   kap=kap_of(ll,ikap)
   ch%warn%ikap=ikap
   fj=1.0d0
   if(nkap==2) fj=kap_j_share(ll,ikap,1.0d0)

   uu(:,:)=0.0d0; qq(:,:)=0.0d0
   if(ikap==1) inwell(:)=.false.

   call ae_bound_projectors(ae%irel,ae%c_light,ll,kap,ikap,inp%nc,inp%nv, &
&                           inp%na,inp%la,ae%ea,nprj,ae%it,ae%rr,inp%zz, &
&                           ae%vfull,mmax,ae%adams,uu,up,iprj, &
&                           ch%epa(:,ikap),ch%npa,ch%uua(:,:,ikap), &
&                           ch%upa(:,:,ikap),iout,ifail)
   if(ifail/=0) return
   call ae_well_projectors(ae%irel,ae%c_light,ll,kap,ikap,nprj, &
&                          inp%ep(l1,ikap),inp%debl(l1),inp%well_n(l1), &
&                          inp%ep_as_given,inp%well_onset,inp%well_tail, &
&                          well_norm_rc,inp%rc(l1),ch%irc,inp%qcut(l1),ae%rr, &
&                          inp%zz,ae%vfull,mmax,ae%adams,uu,up,nnwell,iprj, &
&                          ch%epa(:,ikap),ch%npa,ch%uua(:,:,ikap), &
&                          ch%upa(:,:,ikap),vwell,inwell,wmesh,ch%warn,iout, &
&                          ifail)
   if(ifail/=0) return
   call ae_delta_window(ae%irel,ae%c_light,ll,kap,nprj,ch%epa(:,ikap), &
&                       ae%rr,ae%vfull,ch%uua(:,:,ikap),ch%upa(:,:,ikap), &
&                       mmax,inp%rc(l1),ch%irc,l1==inp%lloc+1, &
&                       ch%dlta(:,:,ikap),ch%vrloc(:,ikap),iout,ifail)
   if(ifail/=0) return
   call ae_overlap(ae%irel,ae%c_light,ll,kap,nprj,ch%uua(:,:,ikap),ae%rr, &
&                  inp%zz,mmax,inp%rc(l1),qq)
   call run_optimize(ch%epa(1,ikap),ll,kap,ae%irel,ae%c_light,mmax,ae%rr, &
&                    ch%uua(1,1,ikap),ch%upa(1,1,ikap),qq, &
&                    well_tail_free .and. inwell(:),wmesh,inp%rc(l1), &
&                    ch%irc, &
&                    inp%qcut(l1),ch%qmsbf,inp%ncon(l1),inp%nbas(l1), &
&                    nprj, &
&                    ch%pswf(1,1,ikap),ch%vp(1,ikap),ch%vkb(1,1,ikap), &
&                    ae%vfull,ch%cvgplt(1,1,1,ikap),ncurve, &
&                    ch%ecurve(1,1,1,ikap), &
&                    ch%pswfx(1,1,ikap),ch%chix(1,1,ikap), &
&                    ch%ektail(1,1,ikap),inp%wave_regularization, &
&                    ch%optdg(1,1,ikap), &
&                    ch%wshare(1,1,ikap),ch%warn,iout,ifail,node_only)
   if(ifail/=0) return

!the weights of its pseudo-waves in the cutoff of the pseudo-waves
   call ecut_weights(ll,inp%nc,inp%nv,inp%la,inp%fa,fj,nprj, &
&                    inp%ecut_weight(l1), &
&                    ch%wcurve(:,ikap),iout,ifail)
   if(ifail/=0) return

 end do !ikap
 ch%wellprj(:)=inwell(:)
 return
 end subroutine optimize_channel

 function channel_key(ae,inp,l1,nprj) result(key)

!! what the channel l1 with nprj projectors depends on: the inputs of
!! optimize_channel (its l's and the well's and regularization's) and the
!! all-electron atom, by its numbers that tell one atom from another;
!! two channels of the same key are the same channel
!!
!! @param ae  the all-electron atom (ae_state_m)
!! @param inp  the input values (input_m)
!! @param l1  l+1
!! @param nprj  number of projectors

 use ae_state_m, only: ae_state_t
 use input_m, only: input_t

!Input variables
 type(ae_state_t), intent(in) :: ae
 type(input_t), intent(in) :: inp
 integer, intent(in) :: l1,nprj

!Output variables
 real(dp) :: key(nkey)

 key(:)=[real(l1,dp),real(nprj,dp),inp%rc(l1),inp%qcut(l1), &
&        real(inp%ncon(l1),dp),real(inp%nbas(l1),dp),inp%ep(l1,1), &
&        inp%ep(l1,mxkap),inp%debl(l1),real(inp%well_n(l1),dp), &
&        inp%ecut_weight(l1),merge(1.0d0,0.0d0,inp%ep_as_given), &
&        inp%well_onset,real(inp%well_tail,dp),real(inp%well_norm,dp), &
&        inp%wave_regularization,real(inp%lloc,dp),ae%zz, &
&        real(ae%irel,dp),real(ae%iexc,dp),real(ae%nc,dp),real(ae%nv,dp), &
&        real(ae%mmax,dp),ae%al,ae%rr(1),ae%c_light,ae%etot,ae%res]
 return
 end function channel_key

 end module channel_m
