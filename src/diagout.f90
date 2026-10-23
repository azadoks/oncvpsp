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
! The machine-readable block of the diagnostics, for programs that read the
! output (an optimizer taking finite differences of them): the values that
! the tables above print to 6 or 2 decimals, to 17 significant digits.

 subroutine diagout(lmax,nproj,mxprj,nkap,cvgplt,npsh,epsh2,depsh,rphs, &
&                   pshae,pshps,ngh,ghn,ghe,ghp,wcut,errw,ecutw,wcurve, &
&                   ektail,icmod,rr,mmax,rcmax,rhoc,corepw,lloc,evkb, &
&                   bbeig,ncfg,ntcf,itcf,etcf,dtcf,cvloc,vlchk,optdg, &
&                   wshare)

!! writes the block "Begin DIAGNOSTICS" ... "End DIAGNOSTICS", one value
!! per line, each line tagged with its quantity:
!!
!! ecut l ikap iproj err ecut: the cutoff (Ha) at which the residual
!! kinetic energy per electron of the optimized wave of projector iproj
!! falls to err (Ha), for the levels err that it reaches, from 1e-2 down
!! by half decades to 1e-7, below the tables' 1e-5 (run_optimize), on the
!! curve with its tail: the q integral to q_infinity plus ekin_tail
!!
!! ekin_tail l ikap iproj t: the curve's tail, the residual kinetic energy
!! of that wave from q_infinity to three times that (Ha), up to 1.4e-7 Ha
!! in the reference outputs
!!
!! ekin_check l ikap iproj t: the r-space kinetic energy of that wave less
!! the Fourier one to three times q_infinity (Ha): by Parseval the energy
!! above that (below 1e-10 Ha in the reference outputs), plus the errors
!! of the q grid's quadrature below q_infinity (up to 8e-9 Ha) and of the
!! r-space integral (up to 5e-9 Ha at the default mesh), of either sign
!!
!! ld_r l ikap r: the radius (bohr) of the log derivatives of l, ikap
!!
!! ld l ikap e ae ps: atan(r u'/u) at that radius and energy e (Ha), of
!! the all-electron and the pseudo wave (run_phsft)
!!
!! ghost l ikap n e_ref e_basis: a bound state of the nonlocal
!! Schroedinger equation (Ha; 0 for the vacuum level when l has none) and
!! the eigenvalue of the barrier-basis diagonalization matched to it
!! (run_ghosts); a lower eigenvalue is a GHOST(-), reported below
!! min(e_ref - 1e-4, e_ref*(1 + 1e-4)), so e_basis - e_ref is the margin
!!
!! ghost_r l ikap e r: of the basis eigenstates more than 0.05 Ha above
!! the reference, the one of least <r>/rc: its energy and <r>/rc, a
!! GHOST(+) below 1
!!
!! ecut_w err ecut: the cutoff (Ha) at which the weighted mean of the
!! residual kinetic energies per electron falls to err (Ha), the cutoff of
!! the model core's budget with budget_cutoff = "weighted"
!! (modcore_budget_table); only then written, with
!!
!! ecut_w_weight l ikap iproj w: the weight of each pseudo-wave in the mean
!! (those with w > 0)
!!
!! core_rc r rho_ae rho_model: the radius rcmax and the all-electron and
!! model core charges there (4 pi rho, as rhoc); with a model core
!! (icmod>0) only, as the next three
!!
!! core_extent frac r_ae r_model: the radius at which each falls to frac
!! of its value at rcmax, outwards from it, for frac = 1e-2 and 1e-4
!!
!! core_model_max r rho, core_model_min r rho: the largest and smallest
!! value of the model core on (0, rcmax] and where it is, an end of the
!! interval or a stationary point (a hollow core's maximum is inside)
!!
!! vkb_b l ikap i b: the eigenvalues of the channel's symmetrized B matrix
!! <pswf_i|chi_j>, ascending, and vkb_cond l ikap c their largest over their
!! smallest magnitude, the condition number of B (1 for one projector)
!!
!! vkb_e l ikap i e: the coefficients of the orthonormal projectors (Ha),
!! the eigenvalues of the nonlocal kernel on its span, largest magnitude
!! first, as printed
!!
!! test_config icfg de_ae de_ps err: of test configuration icfg (0 the
!! reference), the all-electron and pseudo total-energy differences from
!! the reference configuration and the excitation error de_ae - de_ps
!! (Ha); test_config_failed icfg in place of these lines where either atom
!! failed
!!
!! test_eig icfg n l ikap f eae eps: its valence states as in its table,
!! the occupation and the all-electron and pseudo eigenvalues (Ha)
!!
!! ecut_vloc err ecut: the wave cutoff (Ha) at which the local energy of
!! the reference configuration's pseudo valence density, summed over the
!! density's sphere (2 sqrt(2 ecut)), is converged to err per electron
!! (Ha), bounded by the tail of |n(q) V_loc(q)|, for the levels 1e-2 down
!! by half decades to 1e-5 that it reaches (vloc_cutoff)
!!
!! vloc_check e: the local energy per electron from the q grid less that
!! from r space (Ha), the error of the transforms
!!
!! opt l ikap iproj e reg e1 f1: of the optimized wave of projector iproj,
!! E_resid at qcut (Ha), the regularization's term (wave_regularization
!! times E_resid at qcut/4, 0 without it) and the smallest eigenvalue e1
!! and |f1| of the minimization at qcut (optimize): both at round-off
!! leave the wave undetermined
!!
!! well_share l ikap iproj e_well e_free: of the optimized wave of a well
!! state confined by its well (well_tail = "well"), E_resid at qcut with the
!! well's tail and with the free tail at its energy (Ha), the same wave
!! inside rc: e_well - e_free is the well's share; well_share_undefined l
!! ikap iproj in its place where the free tail is not defined
!! (free_tail_defined): qcut not above the state's wave vector, which
!! makes its cutoffs meaningless, or kappa L > 12
!!
!! ikap is 1 for j=l+1/2 (and for nkap=1), 2 for j=l-1/2
!!
!! @param lmax  maximum angular momentum
!! @param nproj  number of projectors for each l
!! @param mxprj  dimension of the number of projectors
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!! @param cvgplt  the crossings, (Ecut, err) per level, projector, l, ikap
!! @param npsh  number of energies of the log derivatives
!! @param epsh2  the first (highest) energy
!! @param depsh  the energy step (downwards)
!! @param rphs  the radius of the log derivatives of each l
!! @param pshae  the all-electron log derivatives
!! @param pshps  the pseudopotential log derivatives
!! @param ngh  rows of the bound-state ghost test per l and kap
!! @param ghn  their principal quantum numbers
!! @param ghe  their energies, e_ref and e_basis
!! @param ghp  the least-localized-ghost margin per l and kap, e and <r>/rc
!!             (0 if no state qualifies)
!! @param wcut  whether to write the weighted cutoff
!! @param errw  its error per electron (Ha)
!! @param ecutw  the weighted cutoff (Ha)
!! @param wcurve  the weights of the pseudo-waves
!! @param ektail  of each wave: (1,:) the curve's tail, its E_resid from
!!                q_infinity to three times that; (2,:) the r-space less
!!                the Fourier kinetic energy to three times q_infinity
!! @param icmod  the model core's kind, 0 for none
!! @param rr  log radial mesh
!! @param mmax  size of the log radial mesh
!! @param rcmax  the largest core radius (bohr)
!! @param rhoc  the all-electron core charge
!! @param corepw  the model core charge and its derivatives as pieces
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param evkb  the projector coefficients
!! @param bbeig  the eigenvalues of the B matrices
!! @param ncfg  number of configurations, the reference and the tests
!! @param ntcf  rows of each one's eigenvalue table, -1 if it failed
!! @param itcf  n, l, ikap of each row
!! @param etcf  f, eae, eps of each row
!! @param dtcf  the AE and PS energy differences and the excitation error
!! @param cvloc  the local potential's (Ecut, err) per level, Ecut 0 for
!!               a level not reached
!! @param vlchk  the local energy from the q grid less that from r space
!! @param optdg  of each optimized wave: E_resid at qcut, the
!!               regularization's term, e1 and |f1|
!! @param wshare  of each optimized wave: E_resid at qcut with the well's
!!                tail and with the free tail, and 1 for a confined well
!!                state, -1 where the free tail is undefined, 0 otherwise

 use constants_m, only: dp, mxl, mxlin, nerrdg, nerrlv, mxstate, mxkap, &
&                       nrhomod, mxcnf
 use pieces_m, only: piecewise_t, pw_value, pw_node_values, pw_segment
 use logmesh_m, only: lm_interp, lm_root, lm_derivs
 use ugrid_m, only: ug_brent
 use input_echo_m, only: fe
 implicit none

!Input variables
 integer, intent(in) :: lmax,mxprj,nkap,npsh
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: cvgplt(2,nerrdg,mxprj,mxl+1,nkap)
 real(dp), intent(in) :: epsh2,depsh
 real(dp), intent(in) :: rphs(4),pshae(npsh,4,nkap),pshps(npsh,4,nkap)
 integer, intent(in) :: ngh(mxl+1,mxkap),ghn(mxstate,mxl+1,mxkap)
 real(dp), intent(in) :: ghe(2,mxstate,mxl+1,mxkap),ghp(2,mxl+1,mxkap)
 logical, intent(in) :: wcut
 real(dp), intent(in) :: errw,ecutw,wcurve(mxprj,mxl+1,nkap)
 integer, intent(in) :: icmod,mmax
 real(dp), intent(in) :: rr(mmax),rcmax,rhoc(mmax)
 type(piecewise_t), intent(in) :: corepw(nrhomod)
 real(dp), intent(in) :: ektail(2,mxprj,mxl+1,nkap)
 integer, intent(in) :: lloc,ncfg
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap),bbeig(mxprj,mxl+1,nkap)
 integer, intent(in) :: ntcf(mxcnf),itcf(3,2*mxstate,mxcnf)
 real(dp), intent(in) :: etcf(3,2*mxstate,mxcnf),dtcf(3,mxcnf)
 real(dp), intent(in) :: cvloc(2,nerrlv),vlchk
 real(dp), intent(in) :: optdg(4,mxprj,mxl+1,nkap)
 real(dp), intent(in) :: wshare(3,mxprj,mxl+1,nkap)

!Local variables
 integer :: ii,jj,l1,ikap,mkap,iprj,kseg,np
 real(dp) :: aerc,mdrc,frac,rae,rmd,rext(2),fext(2)
 real(dp), allocatable :: fmd(:)

 write(6,'(/a)') 'Begin DIAGNOSTICS'
 write(6,'(a)') '# ecut l ikap iproj err(Ha) ecut(Ha)'
 do l1=1,lmax+1
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     do iprj=1,nproj(l1)
       do jj=1,nerrdg
         if(cvgplt(1,jj,iprj,l1,ikap)/=0.0d0) then
           write(6,'(a,3i4,1p,2e24.16)') 'ecut',l1-1,ikap,iprj, &
&                fe(cvgplt(2,jj,iprj,l1,ikap)),cvgplt(1,jj,iprj,l1,ikap)
         end if
       end do
     end do
   end do
 end do

 write(6,'(a)') '# opt l ikap iproj e(Ha) reg(Ha) e1(Ha) f1(Ha)'
 do l1=1,lmax+1
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     do iprj=1,nproj(l1)
       write(6,'(a,3i4,1p,4e24.16)') 'opt',l1-1,ikap,iprj, &
&            fe(optdg(:,iprj,l1,ikap))
     end do
   end do
 end do

 write(6,'(a)') '# well_share l ikap iproj e_well(Ha) e_free(Ha)'
 do l1=1,lmax+1
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     do iprj=1,nproj(l1)
       if(wshare(3,iprj,l1,ikap)>0.0d0) then
         write(6,'(a,3i4,1p,2e24.16)') 'well_share',l1-1,ikap,iprj, &
&              fe(wshare(1:2,iprj,l1,ikap))
       else if(wshare(3,iprj,l1,ikap)<0.0d0) then
         write(6,'(a,3i4)') 'well_share_undefined',l1-1,ikap,iprj
       end if
     end do
   end do
 end do

 write(6,'(a)') '# ecut_vloc err(Ha) ecut(Ha)'
 do jj=1,nerrlv
   if(cvloc(1,jj)/=0.0d0) then
     write(6,'(a,1p,2e24.16)') 'ecut_vloc',fe(cvloc(2,jj)),cvloc(1,jj)
   end if
 end do
 write(6,'(a)') '# vloc_check e(Ha)'
 write(6,'(a,1p,e24.16)') 'vloc_check',fe(vlchk)

 write(6,'(a)') '# ekin_tail l ikap iproj t(Ha)'
 do l1=1,lmax+1
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     do iprj=1,nproj(l1)
       write(6,'(a,3i4,1p,e24.16)') 'ekin_tail',l1-1,ikap,iprj, &
&            fe(ektail(1,iprj,l1,ikap))
     end do
   end do
 end do

 write(6,'(a)') '# ekin_check l ikap iproj t(Ha)'
 do l1=1,lmax+1
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     do iprj=1,nproj(l1)
       write(6,'(a,3i4,1p,e24.16)') 'ekin_check',l1-1,ikap,iprj, &
&            fe(ektail(2,iprj,l1,ikap))
     end do
   end do
 end do

 if(wcut) then
   write(6,'(a)') '# ecut_w err(Ha) ecut(Ha)'
   write(6,'(a,1p,2e24.16)') 'ecut_w',fe(errw),ecutw
   write(6,'(a)') '# ecut_w_weight l ikap iproj w'
   do l1=1,lmax+1
     mkap=nkap
     if(l1==1) mkap=1
     do ikap=1,mkap
       do iprj=1,mxprj
         if(wcurve(iprj,l1,ikap)>0.0d0) then
           write(6,'(a,3i4,1p,e24.16)') 'ecut_w_weight',l1-1,ikap,iprj, &
&                fe(wcurve(iprj,l1,ikap))
         end if
       end do
     end do
   end do
 end if

!none with --no-log-derivatives (npsh 0)
 if(npsh>0) then
   write(6,'(a)') '# ld_r l ikap r(bohr)'
   write(6,'(a)') '# ld l ikap e(Ha) atan(r*u''/u)_ae atan(r*u''/u)_ps'
 end if
 do l1=1,4
   if(npsh==0) exit
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     write(6,'(a,2i4,1p,e24.16)') 'ld_r',l1-1,ikap,rphs(l1)
     do ii=1,npsh
       write(6,'(a,2i4,1p,3e24.16)') 'ld',l1-1,ikap, &
&            epsh2-depsh*real(ii-1,dp),pshae(ii,l1,ikap),pshps(ii,l1,ikap)
     end do
   end do
 end do
 write(6,'(a)') '# ghost l ikap n e_ref(Ha) e_basis(Ha)'
 write(6,'(a)') '# ghost_r l ikap e(Ha) <r>/rc'
 do l1=1,lmax+1
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     do jj=1,ngh(l1,ikap)
       write(6,'(a,3i4,1p,2e24.16)') 'ghost',l1-1,ikap,ghn(jj,l1,ikap), &
&            fe(ghe(:,jj,l1,ikap))
     end do
     if(ghp(2,l1,ikap)>0.0d0) then
       write(6,'(a,2i4,1p,2e24.16)') 'ghost_r',l1-1,ikap,ghp(:,l1,ikap)
     end if
   end do
 end do
 write(6,'(a)') '# vkb_b l ikap i b'
 write(6,'(a)') '# vkb_cond l ikap cond(B)'
 write(6,'(a)') '# vkb_e l ikap i e(Ha)'
 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     np=nproj(l1)
     if(np==0) cycle
     do iprj=1,np
       write(6,'(a,3i4,1p,e24.16)') 'vkb_b',l1-1,ikap,iprj, &
&            fe(bbeig(iprj,l1,ikap))
     end do
     write(6,'(a,2i4,1p,e24.16)') 'vkb_cond',l1-1,ikap, &
&          maxval(abs(bbeig(1:np,l1,ikap)))/minval(abs(bbeig(1:np,l1,ikap)))
     do iprj=1,np
       write(6,'(a,3i4,1p,e24.16)') 'vkb_e',l1-1,ikap,iprj, &
&            fe(evkb(iprj,l1,ikap))
     end do
   end do
 end do

 write(6,'(a)') '# test_config icfg de_ae(Ha) de_ps(Ha) err(Ha)'
 write(6,'(a)') '# test_eig icfg n l ikap f eae(Ha) eps(Ha)'
 do jj=1,ncfg
   if(ntcf(jj)<0) then
     write(6,'(a,i4)') 'test_config_failed',jj-1
     cycle
   end if
   write(6,'(a,i4,1p,3e24.16)') 'test_config',jj-1,fe(dtcf(:,jj))
   do ii=1,ntcf(jj)
     write(6,'(a,4i4,1p,3e24.16)') 'test_eig',jj-1,itcf(:,ii,jj), &
&          fe(etcf(:,ii,jj))
   end do
 end do

 if(icmod>0) then
   allocate(fmd(mmax))
   call pw_node_values(rr,mmax,corepw(1),fmd)
   call lm_interp(rr,rhoc,mmax,rcmax,aerc)
   mdrc=pw_value(rr,mmax,corepw(1),rcmax)
   write(6,'(a)') '# core_rc r(bohr) rho_ae rho_model'
   write(6,'(a,1p,3e24.16)') 'core_rc',rcmax,fe(aerc),fe(mdrc)
   write(6,'(a)') '# core_extent frac r_ae(bohr) r_model(bohr)'
   frac=1.0d-2
   do jj=1,2
     rae=fall_to(rhoc,frac*aerc,0)
     rmd=fall_to(fmd,frac*mdrc,1)
     write(6,'(a,1p,3e24.16)') 'core_extent',frac,rae,rmd
     frac=1.0d-2*frac
   end do
   call model_extrema(rext,fext)
   write(6,'(a)') '# core_model_max r(bohr) rho'
   write(6,'(a,1p,2e24.16)') 'core_model_max',rext(1),fe(fext(1))
   write(6,'(a)') '# core_model_min r(bohr) rho'
   write(6,'(a,1p,2e24.16)') 'core_model_min',rext(2),fe(fext(2))
   deallocate(fmd)
 end if
 write(6,'(a)') 'End DIAGNOSTICS'

 return

 contains

 real(dp) function fall_to(ff,level,model)

!! the radius beyond rcmax at which the core charge ff on the mesh first
!! falls to level: the root of its interpolant in the cell of the first
!! node below it (of the model core's piece of that node, model=1); 0 if
!! it never does
 real(dp), intent(in) :: ff(mmax),level
 integer, intent(in) :: model
 real(dp), allocatable :: gg(:)
 real(dp) :: ra
 integer :: kk,ierr

 fall_to=0.0d0
 do kk=1,mmax
   if(rr(kk)>rcmax .and. ff(kk)<level) exit
 end do
 if(kk>mmax) return
 allocate(gg(mmax))
 if(model==1 .and. pw_segment(corepw(1),rr(kk))>0) then
   gg(:)=corepw(1)%vals(:,pw_segment(corepw(1),rr(kk)))-level
 else
   gg(:)=ff(:)-level
 end if
 ra=max(rcmax,rr(kk-1))
 call lm_root(rr,gg,mmax,ra,rr(kk),fall_to,ierr)
 if(ierr/=0) then
   write(6,'(a)') 'diagout: ERROR the core charge''s bracket lost its sign'
   error stop 1
 end if
 deallocate(gg)
 end function fall_to

 subroutine model_extrema(rq,fq)

!! the largest (rq(1), fq(1)) and smallest (rq(2), fq(2)) value of the
!! model core on (0, rcmax]: from its nodes, at an end of the interval
!! (rr(1) for the origin) or at the root of its derivative between the
!! neighbours of the extreme node
 real(dp), intent(out) :: rq(2),fq(2)
 integer :: kk,nrc,iext,ierr
 real(dp) :: sg

 nrc=1
 do kk=1,mmax
   if(rr(kk)>rcmax) exit
   nrc=kk
 end do
 do kk=1,2
   sg=1.0d0
   if(kk==2) sg=-1.0d0
   iext=maxloc(sg*fmd(1:nrc),1)
   if(iext==1) then
     rq(kk)=rr(1)
     fq(kk)=fmd(1)
   else if(iext==nrc .and. sg*mdrc>=sg*fmd(nrc)) then
     rq(kk)=rcmax
     fq(kk)=mdrc
   else
     kseg=pw_segment(corepw(1),rr(iext))
     call ug_brent(slope,rr(iext-1),min(rr(iext+1),rcmax),rq(kk),ierr)
     if(ierr/=0) then
       rq(kk)=rr(iext)
       fq(kk)=fmd(iext)
     else
       fq(kk)=pw_value(rr,mmax,corepw(1),rq(kk))
     end if
   end if
 end do
 end subroutine model_extrema

 real(dp) function slope(rq)
!! d/dr of the model core's piece kseg at rq
 real(dp), intent(in) :: rq
 real(dp) :: dd(2)
 call lm_derivs(rr,corepw(1)%vals(:,kseg),mmax,rq,1,dd)
 slope=dd(2)
 end function slope

 end subroutine diagout
