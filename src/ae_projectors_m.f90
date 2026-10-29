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
 implicit none
 private
 public :: ae_bound_projectors, ae_well_projectors, ae_delta_window
 public :: ae_overlap, ecut_weights

 contains

 subroutine ae_bound_projectors(irel,c_light,ll,kap,ikap,nc,nv,na,la,ea,nproj, &
&                               it,rr,zz,vfull,mmax,adams,uu,up,iprj,epa,npa, &
&                               uua, &
&                               upa,iout,ifail)
!! the bound valence states of l (at most nproj), the first iprj
!! projectors' states; and npa(1), the n above the highest core state of l
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2
!! @param ikap  the component of l, 1 (j=l+1/2) or 2 (j=l-1/2)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param na  principal quantum numbers
!! @param la  angular momenta
!! @param ea  the all-electron eigenvalues (Ha) of each state and kappa
!! @param nproj  number of projectors for each l, here that of l
!! @param it  the iterations of the reference atom, for the error message
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param vfull  all-electron potential (Ha)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param uu  work space for the radial wave function (times r)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr);
!!            work space as uu
!! @param iprj  the number of bound states found
!! @param epa  the energies of the projectors' states (Ha), set for 1..iprj
!! @param npa  their principal quantum numbers; npa(1) also when iprj=0
!! @param uua  the projectors' all-electron states (times r), the large
!!             component for irel=2
!! @param upa  their derivatives in the mesh index (al*r*d(uua)/dr),
!!             renormalized with uua for irel=2
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use failure_m, only: fail_convergence, failed, msg_ae_proj_bound_conv
 use boundstate_m, only: boundstate
 use renorm_r_m, only: renorm_r
 use constants_m, only: dp, mxprj, mxstate, mxkap
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: irel,ll,kap,ikap,nc,nv,na(mxstate),la(mxstate)
 real(dp), intent(in) :: c_light
 integer, intent(in) :: nproj,it,mmax
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: ea(mxstate,mxkap),rr(mmax),zz,vfull(mmax)
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 real(dp), intent(inout) :: epa(mxprj),uua(mmax,mxprj),upa(mmax,mxprj)
 integer, intent(inout) :: npa(mxprj)

!Output variables
 integer, intent(out) :: iprj
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: et,cnorm
 integer :: kk,ierr,mch

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
&                      rr,zz,vfull,uu,up,mmax,adams,mch,iout,ifail)
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
         ifail=failed(fail_convergence,msg_ae_proj_bound_conv)
         return
       end if

       epa(iprj)=ea(kk,ikap)
       npa(iprj)=na(kk)
       if(irel==2) then
         call renorm_r(uu,rr,ll,kap,c_light,zz,mmax,cnorm,iout=iout,ifail=ifail)
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
&                              adams, &
&                              uu,up,nnwell,iprj,epa,npa,uua,upa,vwell, &
&                              inwell,wmesh,iout,ifail)
!! the well states of the projectors iprj+1 to nproj: the first at ep if
!! there are no bound ones, the next ones debl above; their wells are
!! those of ikap=1, kept for ikap=2 (vwell, wmesh). nnwell is the n of the
!! last well state of ikap=1 that was not fully relativistic
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2
!! @param ikap  the component of l, 1 (j=l+1/2) or 2 (j=l-1/2)
!! @param nproj  number of projectors for each l, here that of l
!! @param ep  the input energy of the first well state (Ha)
!! @param debl  the energy step between well states (Ha), > 0 for two or
!!              more
!! @param well_n  the n of the first well state with well_tail 2, or 0
!! @param ep_as_given  take ep as given (wellstate, D-1)
!! @param well_onset  onset of the well as a multiple of rc, 0 or less
!!                    for the default
!! @param well_tail  0 Hamann's well, 1 the free tail, 2 the confined well
!! @param well_norm_rc  normalize a confined state inside rc
!! @param rc  core radius (bohr)
!! @param irc  last mesh point at or inside rc
!! @param qcut  the cutoff wave vector (bohr**-1), for the free tail
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param vfull  all-electron potential (Ha)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param uu  work space for the radial wave function (times r)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr);
!!            work space as uu
!! @param nnwell  the n of the last well state of ikap=1 that was not fully
!!               relativistic; set for ikap=1
!! @param iprj  the number of the projectors' states, nproj on return
!! @param epa  the energies of the projectors' states (Ha)
!! @param npa  their principal quantum numbers
!! @param uua  the projectors' all-electron states (times r)
!! @param upa  their derivatives in the mesh index (al*r*d(uua)/dr)
!! @param vwell  the potential with the well of each well state (Ha), set
!!               for ikap=1 and kept for ikap=2
!! @param inwell  whether each well state is confined by its well
!! @param wmesh  each confined state on its continued mesh (well_mesh_m)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use failure_m, only: fail_construction, fail_convergence, fail_input, failed, &
&                    msg_ae_proj_debl_scattering, msg_ae_proj_well_conv, &
&                    msg_ae_proj_well_n
 use fpovlp_m, only: fpovlp
 use boundstate_m, only: boundstate
 use renorm_r_m, only: renorm_r
 use wellstate_m, only: wellstate
 use constants_m, only: dp, mxprj
 use well_mesh_m, only: well_mesh_t
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: irel,ll,kap,ikap,nproj,well_n,well_tail,irc,mmax
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ep,debl,well_onset,rc,qcut,rr(mmax),zz
 real(dp), intent(in) :: vfull(mmax)
 logical, intent(in) :: ep_as_given,well_norm_rc
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 integer, intent(inout) :: nnwell,iprj,npa(mxprj)
 real(dp), intent(inout) :: epa(mxprj),uua(mmax,mxprj),upa(mmax,mxprj)
 real(dp), intent(inout) :: vwell(mmax,mxprj)
 logical, intent(inout) :: inwell(mxprj)
 type(well_mesh_t), intent(inout) :: wmesh(mxprj)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: et,cnorm,eptgt
 integer :: kk,ierr,mch,mchx
 logical :: well_tail_free

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
     ifail=failed(fail_construction,msg_ae_proj_debl_scattering)
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
       write(iout,'(/a,i0,a,i0,a,i0)') &
&           'ae_well_projectors: ERROR well_n for l=',ll, &
&           ' is',well_n,', below the first well state''s n,', &
&           npa(iprj)
       ifail=failed(fail_input,msg_ae_proj_well_n)
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
&                   vwell(1,iprj),uu,up,zz,mmax,adams,mch,irel,c_light, &
&                   inwell(iprj),eptgt,wmesh(iprj),iout,ifail)
     if(ifail/=0) return
   else if(ikap==1) then
     nnwell=npa(iprj)
     call wellstate(nnwell,ll,kap,rc,epa(iprj),ep_as_given, &
&                   well_onset,well_tail==2,well_tail_free,rr,vfull, &
&                   vwell(1,iprj),uu,up,zz,mmax,adams,mch,irel,c_light, &
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
&                      wmesh(iprj)%mx,adams,mchx,iout,ifail)
       if(ifail/=0) return
       uu(:,:)=wmesh(iprj)%uu(1:mmax,:)
       up(:,:)=wmesh(iprj)%up(1:mmax,:)
     else
       call boundstate(irel,c_light,npa(iprj),ll,kap,ierr, &
&                      epa(iprj),rr,zz,vwell(1,iprj),uu,up, &
&                      mmax,adams,mch,iout,ifail)
       if(ifail/=0) return
     end if

     if(ierr /= 0) then
       write(iout,'(a,i2,a,i2,a,i3,a,f12.6)') &
&            'oncvpsp_r 474: ldiracfb convergence ERROR, n=',&
&             npa(iprj),'  l=',ll, 'kap=',kap, ' epa=', &
&             epa(iprj)

       ifail=failed(fail_convergence,msg_ae_proj_well_conv)
       return
     end if

     write(iout,'(/a,i3,a,i3,a,i3/a,f8.4)') &
&             '   Wellstate for l =',ll,'  n =',npa(iprj), &
&             ' kap=',kap, &
&             '     eigenvalue = ',et
   end if

   if(irel==2 .and. wmesh(iprj)%mx>0) then
     call renorm_r(wmesh(iprj)%uu,wmesh(iprj)%rr,ll,kap,c_light,zz, &
&                  wmesh(iprj)%mx,cnorm,iout=iout,ifail=ifail)
     if(ifail/=0) return
     uu(:,:)=wmesh(iprj)%uu(1:mmax,:)
   else if(irel==2) then
     call renorm_r(uu,rr,ll,kap,c_light,zz,mmax,cnorm,iout=iout,ifail=ifail)
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
&                         rc,irc,qcut,rr,zz,vfull,mmax,adams,uu,up,uua(:, &
&                         iprj), &
&                         upa(:,iprj),iout,ifail)
     if(ifail/=0) return
   end if
 end do !kk
 end subroutine ae_well_projectors

 subroutine free_tail_state(irel,c_light,ll,kap,ikap,nnwell,nn,ee,rc,irc,qcut, &
&                           rr,zz,vfull,mmax,adams,uu,up,ua,pa,iout,ifail)
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
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2
!! @param ikap  the component of l, 1 (j=l+1/2) or 2 (j=l-1/2)
!! @param nnwell  the n of the last well state of ikap=1, the node count
!!               scatterstate is given and sets
!! @param nn  the n of the well state, for the error message
!! @param ee  the energy of the state (Ha), >= 0 for ikap=1
!! @param rc  core radius (bohr)
!! @param irc  last mesh point at or inside rc
!! @param qcut  the cutoff wave vector (bohr**-1)
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param vfull  all-electron potential (Ha)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param uu  work space for the radial wave function (times r)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr);
!!            work space as uu
!! @param ua  the state (times r), normalized inside rc, with the free
!!            tail
!! @param pa  its derivative in the mesh index (al*r*d(ua)/dr)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use failure_m, only: fail_input, failed, msg_free_tail_ep
 use fpovlp_m, only: fpovlp
 use scatterstate_m, only: scatterstate
 use constants_m, only: dp
 use free_tail_m, only: free_tail_log
 use logmesh_m, only: lm_fh_reach
 use pieces_m, only: pw_reach
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: irel,ll,kap,ikap,nn,irc,mmax
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ee,rc,qcut,rr(mmax),zz,vfull(mmax)
 integer, intent(in) :: iout

!Input/output variables
 integer, intent(inout) :: nnwell
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)

!Output variables
 real(dp), intent(out) :: ua(mmax),pa(mmax)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: cnorm
 integer :: mfree
 ifail=0
 if(ikap==1 .and. ee<0.0d0) then
   write(iout,'(/a/a,i0,a,i0,a,f12.6,a/a)') &
&       'free_tail_state: ERROR well_tail = "free" needs ep >= 0:', &
&       'the well state of l=',ll,', n=',nn,' at ep=', &
&       ee,' Ha has a growing tail;', &
&       'raise its ep (debl), or use well_tail = "confined"'// &
&       ' or "well"'
   ifail=failed(fail_input,msg_free_tail_ep)
   return
 end if
 mfree=min(mmax,irc+pw_reach+16+lm_fh_reach()+2)
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 call scatterstate(irel,c_light,nnwell,ll,kap,ee,rr, &
&                  zz,vfull,uu,up,mmax,adams,mfree,ifail)
 if(ifail/=0) return
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
&                           mmax,rc,irc,local,dlta,vrloc,iout,ifail)
!! the relativistic term Delta of the all-electron wave about the taper
!! window [0.95 rc, rc] (ae_delta, D-4), for the projectors; the local
!! potential of lloc<4 (local) takes Delta/(2u) on the window (S-1)
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2
!! @param nproj  number of projectors for each l, here that of l
!! @param epa  the energies of the projectors' states (Ha)
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param uua  the projectors' all-electron states (times r)
!! @param upa  their derivatives in the mesh index (al*r*d(uua)/dr)
!! @param mmax  size of the log radial mesh
!! @param rc  core radius (bohr)
!! @param irc  last mesh point at or inside rc
!! @param local  whether l is that of the local potential
!! @param dlta  Delta of each state (ae_delta), set about the window
!! @param vrloc  Delta/(2u) of the first state about the window when local
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use failure_m, only: fail_construction, failed, msg_delta_window_node
 use ae_delta_m, only: ae_delta
 use constants_m, only: dp, mxprj, rel_window
 use logmesh_m, only: lm_index
 use pieces_m, only: pw_reach

!Input variables
 integer, intent(in) :: irel,ll,kap,nproj,mmax,irc
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: epa(mxprj),rr(mmax),vfull(mmax)
 real(dp), intent(in) :: uua(mmax,mxprj),upa(mmax,mxprj),rc
 logical, intent(in) :: local
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: dlta(mmax,mxprj),vrloc(mmax)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: iprj,ii,jj,kk
 ifail=0
 do iprj=1,nproj
   ii=max(1,lm_index(rr,mmax,rel_window*rc)-pw_reach-8)
   jj=min(mmax,irc+pw_reach+16)
   call ae_delta(irel,c_light,ll,kap,epa(iprj),rr,vfull,uua(1,iprj), &
&                upa(1,iprj),mmax,ii,jj,dlta(1,iprj))
   if(local .and. iprj==1) then
     do kk=ii,jj
       if(uua(kk,1)==0.0d0) then
         write(iout,'(/a,i0,a,f12.6)') &
&              'ae_delta_window: ERROR node of the local channel l=',ll, &
&              ' all-electron wave near rc, at r=',rr(kk)
         ifail=failed(fail_construction,msg_delta_window_node)
         return
       end if
       vrloc(kk)=0.5d0*dlta(kk,1)/uua(kk,1)
     end do
   end if
 end do
 end subroutine ae_delta_window

 subroutine ae_overlap(irel,c_light,ll,kap,nproj,uua,rr,zz,mmax,rc,qq)
!! the overlap matrix of the projectors' all-electron states inside rc
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2
!! @param nproj  number of projectors for each l, here that of l
!! @param uua  the projectors' all-electron states (times r)
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param rc  core radius (bohr)
!! @param qq  the overlaps inside rc (fpovlp), set for 1..nproj
 use fpovlp_m, only: fpovlp
 use constants_m, only: dp, mxprj

!Input variables
 integer, intent(in) :: irel,ll,kap,nproj,mmax
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: uua(mmax,mxprj),rr(mmax),zz,rc

!Input/output variables
 real(dp), intent(inout) :: qq(mxprj,mxprj)

!Local variables
 integer :: ii,jj
 do jj=1,nproj
   do ii=1,jj
     call fpovlp(uua(1,ii),uua(1,jj),mmax,rc,ll,kap, &
&                zz,qq(ii,jj),rr,irel,c_light)
     qq(jj,ii)=qq(ii,jj)
   end do
 end do
 end subroutine ae_overlap

 subroutine ecut_weights(ll,nc,nv,la,fa,fj,nproj,ecut_weight,wcurve,iout,ifail)
!! the weights of the pseudo-waves of l in the cutoff of the pseudo-waves
!! (the weighted mean, and which waves the largest is over): the
!! occupation of each valence state (its share fj of J), or for the first
!! wave of a channel with none occupied ecut_weight, by default
!! 2(2l+1), as if the channel were full; in a solid it holds about the
!! charge that hybridization puts inside rc, an input as the atom cannot
!! give it
!!
!! @param ll  angular momentum
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular momenta
!! @param fa  occupancies
!! @param fj  the share of the occupation of l in this kappa (kap_j_share), 1
!!            when nkap=1
!! @param nproj  number of projectors for each l, here that of l
!! @param ecut_weight  the input weight of the first wave of a channel with no
!!                     occupied state, negative for the default 2(2l+1)
!! @param wcurve  the weights of the waves; zero on entry
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use failure_m, only: fail_input, failed, msg_ecut_weight_unoccupied
 use constants_m, only: dp, mxprj, mxstate

!Input variables
 integer, intent(in) :: ll,nc,nv,la(mxstate),nproj
 real(dp), intent(in) :: fa(mxstate),fj,ecut_weight
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: wcurve(mxprj)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: jj,kk
 ifail=0
 jj=0
 do kk=nc+1,nc+nv
   if(la(kk)==ll .and. jj<nproj) then
     jj=jj+1
     wcurve(jj)=fj*fa(kk)
   end if
 end do
 if(.not. any(wcurve(:)>0.0d0)) then
   if(ecut_weight>=0.0d0) then
     wcurve(1)=fj*ecut_weight
   else
     wcurve(1)=fj*dble(2*(2*ll+1))
   end if
 else if(ecut_weight>=0.0d0) then
   write(iout,'(/a,i0,a)') 'ecut_weights: ERROR ecut_weight for l=',ll, &
&       ' is for a channel with no occupied state; this one has one'
   ifail=failed(fail_input,msg_ecut_weight_unoccupied)
   return
 end if
 end subroutine ecut_weights

 end module ae_projectors_m
