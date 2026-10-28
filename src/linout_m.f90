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

 module linout_m

 implicit none
 private
 public :: linout

 contains

 subroutine linout(inp,c_light,vkbpw,evkb,vsrpw,esr,vsopw,eso,rr,vpuns,vlocpw, &
&                  rho,corepw,rhotae,rhoc,zion,mmax,nrl,soscale,rc_in, &
&                  fcfact_in,fcgrid,rcgrid,log_mesh,nkap,iout,ifail)

!! the pseudopotential interpolated to the linear output mesh in the
!! psp8 format of Abinit (pspcod=8), with the spin-orbit projectors when
!! fully relativistic, written to iout with the echo of the input
!!
!! @param inp  the input values (input_m)
!! @param c_light  speed of light (Hartree atomic units)
!! @param vkbpw  the projectors as pieces, written when nkap=1
!! @param evkb  coefficients of the projectors (Ha)
!! @param vsrpw  normalized scalar-relativistic projectors' pieces, written
!!               with vsopw when nkap=2
!! @param esr  energy coefficients of vsrpw
!! @param vsopw  normalized spin-orbit projectors' pieces
!! @param eso  energy coefficients of vsopw
!! @param rr  log radial mesh (bohr)
!! @param vpuns  unscreened semi-local pseudopotentials (vpuns(:,5) is local
!!               potential if linear combination is used)
!! @param vlocpw  the screened local potential as pieces; vpuns(:,lloc+1)
!!                less their node values is smooth
!! @param rho  valence pseudocharge
!! @param corepw  the model core charge and its first 4 derivatives, as
!!                pieces split at its joint (modcore_m, core_pieces)
!! @param rhotae  all-electron valence charge density
!! @param rhoc  core charge density
!! @param zion  ion charge
!! @param mmax  size of the log radial mesh
!! @param nrl  size of the output mesh
!! @param soscale  possible factor to boost or cut spin-orbit strength, not in
!!                 use yet
!! @param rc_in  input core radii
!! @param fcfact_in  input fcfact, for the echo
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param log_mesh  amesh, rr1, rmax and the Adams order of the log mesh
!!                  line, echoed when the input sets it (else 0)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxl, mxlin, nrhomod, mxprj
 use input_echo_m, only: write_input_echo, fe, date_stamp
 use outmesh_m, only: om_values, om_piece_values, om_local_values
 use pieces_m, only: piecewise_t
 use input_m, only: input_t

!Input variables
 type(input_t), intent(in) :: inp
 real(dp), intent(in) :: c_light
 integer, intent(in) :: mmax,nrl,nkap
 real(dp), intent(in) :: fcfact_in,zion,soscale
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)
 real(dp), intent(in) :: rr(mmax),vpuns(mmax,mxl+2),rho(mmax)
 real(dp), intent(in) :: rhotae(mmax),rhoc(mmax)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 real(dp), intent(in) :: esr(2*mxprj,mxl+1),eso(2*mxprj,mxl+1)
 type(piecewise_t), intent(in) :: vkbpw(mxprj,mxl+1,nkap),vlocpw
 type(piecewise_t), intent(in) :: vsrpw(2*mxprj,mxl+1),vsopw(2*mxprj,mxl+1)
 type(piecewise_t), intent(in) :: corepw(nrhomod)

!the core radii as read (rc_in), for the echo
 real(dp), intent(in) :: rc_in(mxlin+1)
 real(dp), intent(in) :: log_mesh(4)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,jj,ll,l1,ixc_abinit
 integer :: npr_sr(mxl+1),npr_so(mxl+1)
 real(dp) :: fchrg
 real(dp), allocatable :: rhomodl(:,:)
 real(dp),allocatable :: rhol(:),rl(:),vpl(:,:)
 real(dp),allocatable :: rhotael(:),rhocl(:)
 real(dp), allocatable :: vkbl(:,:,:),vsrl(:,:,:),vsol(:,:,:)
 character(len=2) :: pspd(3)
 character(len=48) :: version

 ifail=0
 allocate(rhol(nrl),rl(nrl),vpl(nrl,mxl+2),rhomodl(nrl,nrhomod))
 allocate(rhotael(nrl),rhocl(nrl))
 allocate(vkbl(nrl,mxprj,mxl+1))
 allocate(vsrl(nrl,2*mxprj,mxl+1),vsol(nrl,2*mxprj,mxl+1))

! every term of the sr_so decomposition (O-2): 2 nproj per l, nproj
! scalar and none spin-orbit for l = 0
 npr_sr(:)=0 ; npr_so(:)=0
 if(nkap==2) then
   do l1=1,inp%lmax+1
     if(l1==1) then
       npr_sr(l1)=inp%nproj(l1)
     else
       npr_sr(l1)=2*inp%nproj(l1)
       npr_so(l1)=2*inp%nproj(l1)
     end if
     write(iout,'(a,3i4)') 'l1-1,npr_sr,npr_so',l1-1,npr_sr(l1),npr_so(l1)
   end do
 end if

! everything on the linear output mesh, the projectors and the local
! potential from their pieces (outmesh_m)
!
 do  ii=1,nrl
   rl(ii)=inp%drl*dble(ii-1)
 end do
!
 vpl(:,:)=0.0d0
 l1=inp%lloc+1
 call om_local_values(rr,vpuns(1,l1),mmax,vlocpw,rl,vpl(1,l1),nrl,iout,ifail)
 if(ifail/=0) return

 do l1=1,inp%lmax+1
   if(l1 /= inp%lloc+1) then
     do jj=1,inp%nproj(l1)
       if(nkap==1) then
         call om_piece_values(rr,mmax,vkbpw(jj,l1,1),rl,vkbl(1,jj,l1),nrl, &
&                             iout,ifail)
         if(ifail/=0) return
       end if
     end do
     do jj=1,npr_sr(l1)
       call om_piece_values(rr,mmax,vsrpw(jj,l1),rl,vsrl(1,jj,l1),nrl, &
&            iout,ifail)
       if(ifail/=0) return
     end do
     do jj=1,npr_so(l1)
       call om_piece_values(rr,mmax,vsopw(jj,l1),rl,vsol(1,jj,l1),nrl, &
&            iout,ifail)
       if(ifail/=0) return
     end do
   end if
 end do


 call om_values(rr,rho,mmax,rl,rhol,nrl,iout,ifail)
 if(ifail/=0) return
 call om_values(rr,rhotae,mmax,rl,rhotael,nrl,iout,ifail)
 if(ifail/=0) return
 call om_values(rr,rhoc,mmax,rl,rhocl,nrl,iout,ifail)
 if(ifail/=0) return

 do jj=1,5
   call om_piece_values(rr,mmax,corepw(jj),rl,rhomodl(1,jj),nrl,iout,ifail)
   if(ifail/=0) return
 end do

! Output for Abinit input using pspcod=8

 if(inp%iexc==1) then
   ixc_abinit=4
 else if(inp%iexc==2) then
   ixc_abinit=5
 else if(inp%iexc==3) then
   ixc_abinit=2
 else if(inp%iexc==4) then
   ixc_abinit=11
!Abinit has no built-in VWN; the Libxc codes (lda_x_rel has its own c)
 else if(inp%iexc==6) then
   ixc_abinit=-001007
 else if(inp%iexc==7) then
   ixc_abinit=-532007
 else if(inp%iexc<0) then
   ixc_abinit=inp%iexc
 end if

 call date_stamp(pspd)

!fchrg flags the model core for Abinit: fcfact, or 1 when fcfact<=0
!(icmod=6 with the cutoff taken from the pseudo-waves)
 fchrg=0.0d0
 if(inp%icmod>=1) then
   fchrg=inp%fcfact
   if(fchrg<=0.0d0) fchrg=1.0d0
 end if

 write(iout,'(/a)') 'Begin PSPCODE8'
 write(iout,'(3a,4f10.5)') inp%atsym,'    ONCVPSP-4.0.1' &
&  ,'  r_core=',(inp%rc(l1),l1=1,inp%lmax+1)
 write(iout,'(2f12.4, 5a)') inp%zz,zion, '      ', pspd,  &
&  '    zatom,zion,pspd'
 write(iout,'(i6,i8,i4,3i6, a)') 8, ixc_abinit,inp%lmax,inp%lloc, &
&  nrl, 0, '    pspcod,pspxc,lmax,lloc,mmax,r2well'
 write(iout,'(3f12.8, a)') rl(nrl),fchrg, 0.0d0,  &
&  '    rchrg fchrg qchrg'
 if(nkap==1) then
   write(iout,'(5i6, a)') (inp%nproj(l1),l1=1,5),  &
&    '    nproj'
   write(iout,'(2i6, a)') 1,1, &
&    '           extension_switch'
 else
   write(iout,'(4i6, a)') npr_sr(1:4),'    nproj'

! this is where abinit is informed that we have spin-orbit
   write(iout,'(2i6, a)') 3,1, &
&    '           extension_switch'
   if(soscale==1.0d0) then
     write(iout,'(3i6, a)') npr_so(2:4),'    nprojso'
   else
     write(iout,'(3i6, a,f6.2)') npr_so(2:4),'    nprojso  X,soscale', soscale
   end if
 end if

! write the VKB projectors, or the scalar-relativistic ones, and the local
! potential if it is one of the semi-local potentials

 do l1=1,inp%lmax+1
   ll=l1-1
   if(ll==inp%lloc) then
     write(iout,'(i4)') ll
     do ii = 1,nrl
       write(iout,'(i6,1p,2e21.13)') ii,rl(ii),fe(vpl(ii,l1))
     end do
   else if(nkap==1 .and. inp%nproj(l1)>0) then
     call proj_block(ll,inp%rc(l1),inp%nproj(l1),evkb(:,l1,1),vkbl(:,:,l1), &
&                    mxprj, &
&                    .false.)
   else if(nkap==2 .and. npr_sr(l1)>0) then
     call proj_block(ll,inp%rc(l1),npr_sr(l1),esr(:,l1),vsrl(:,:,l1),2*mxprj, &
&                    .true.)
   end if
 end do

! write general local potential if called for
 if(inp%lloc>inp%lmax) then
   write(iout,'(i4)') inp%lloc
   do ii = 1,nrl
     write(iout,'(i6,1p,2e21.13)') ii,rl(ii),fe(vpl(ii,inp%lloc+1))
   end do
 end if

! write spin-orbit projectors (none when nkap=1)
 do l1=2,inp%lmax+1
   ll=l1-1
   if(npr_so(l1)>0) then
     call proj_block(ll,inp%rc(l1),npr_so(l1),soscale*eso(:,l1),vsol(:,:,l1), &
&                    2*mxprj,.true.)
   end if
 end do

! write the model core charge if called for
 if(fchrg>0.0d0) then
   do ii=1,nrl
     write(iout,'(i6,1p,6e21.13)') ii,rl(ii),  &
&      (fe(rhomodl(ii,jj)),jj=1,5)
   end do
 end if

! write pseudo valence, all-electron valence and all-electron core charges
 do ii=1,nrl
   write(iout,'(i6,1p,4e21.13)') ii,rl(ii), &
&        fe(rhol(ii)),fe(rhotael(ii)),fe(rhocl(ii))
 end do

 if(nkap==1) then
   version='#scalar-relativistic version 4.0.1 0r/27/2018'
 else
   version='#relativistic version 4.0.1 03/03/2109'
 end if
 write(iout,'(a)') '<INPUT>'
 write(iout,'(a/a/a/a)') &
&    '#', &
&    '#ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&    trim(version), &
&    '#'

 write(iout,'(a/a/a/a)') &
&    '#While it is not required under the terms of the GNU GPL, it is',&
&    '#suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&    '#in any publication utilizing these pseudopotentials.', &
&    '#'

 if(nkap==2) write(iout,'(a)') '#Echo of input data for oncvpsp-4.0.1'
 call write_input_echo(nkap,1,inp%atsym,inp%zz,inp%nc,inp%nv,inp%iexc,c_light, &
&                      inp%psfile, &
&                       inp%na,inp%la,inp%fa,inp%lmax,rc_in,inp%ep(:,1), &
&                       inp%ncon,inp%nbas,inp%qcut,inp%ep_as_given, &
&                       inp%well_onset,inp%well_tail,inp%wave_regularization, &
&                       inp%lloc, &
&                       inp%lpopt,inp%dvloc0, &
&                       inp%nproj,inp%debl,inp%well_n,inp%well_norm, &
&                       inp%ecut_weight, &
&                       inp%icmod,fcfact_in,inp%rcfact,inp%slack, &
&                       inp%budget_cutoff, &
&                       fcgrid,rcgrid, &
&                       inp%epsh1,inp%epsh2,inp%depsh,inp%rxpsh,inp%rlmax, &
&                       inp%drl, &
&                       inp%ncnf,inp%cnf, &
&                       log_mesh,inp%upf_mesh,iout)

! write termination signal
 write(iout,'(a)') '</INPUT>'
 write(iout,'(/a)') 'END_PSP'

 deallocate(rhol,rl,vpl,rhomodl,rhotael,rhocl)
 deallocate(vkbl,vsrl,vsol)

 return

 contains

 subroutine proj_block(ll,rcl,np,ee,vkb,ldv,zero_tail)
!! the block of the np projectors vkb of l: l and their coefficients ee, then
!! each point of the linear mesh, their values out to rc+2 drl and beyond
!! 0, or with zero_tail their values in the same short field
!!
!! @param ll  angular momentum
!! @param rcl  core radius of l
!! @param np  number of projectors
!! @param ee  coefficients of the projectors (Ha)
!! @param vkb  projectors on the linear mesh
!! @param ldv  second dimension of ee and vkb
!! @param zero_tail  0 beyond rc+2 drl, else the values
 integer, intent(in) :: ll,np,ldv
 real(dp), intent(in) :: rcl,ee(ldv),vkb(nrl,ldv)
 logical, intent(in) :: zero_tail
 real(dp) :: zero(20)
 integer :: ii,jj
 zero(:)=0.0d0
 write(iout,'(i4,23x,1p,20e21.13)') ll,(fe(ee(jj)),jj=1,np)
 do ii = 1,nrl
   if(rl(ii)<=rcl+2.0d0*inp%drl) then
     write(iout,'(i6,1p,21e21.13)') ii,rl(ii), &
&          (fe(vkb(ii,jj)),jj=1,np)
   else if(zero_tail) then
     write(iout,'(i6,f6.2,20f4.0)') ii,rl(ii),(zero(jj),jj=1,np)
   else
     write(iout,'(i6,f6.2,20f4.0)') ii,rl(ii),(vkb(ii,jj),jj=1,np)
   end if
 end do
 end subroutine proj_block

 end subroutine linout

 end module linout_m
