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
! interpolates various arrays onto linear radial mesh to create file
! for PWSCF input using the UPF file format

subroutine upfout(unit,lmax,lloc,rc,vkb,evkb,nproj,rr,drr,vpuns,rho,rhomod, &
   zz,zion,mmax,mxprj,iexc,icmod,drl,atsym,epstot, &
   na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
   fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
   epsh1,epsh2,depsh,rlmax,psfile,uupsa,ea,upfgrid)
   ! lmax  maximum angular momentum
   ! lloc  l for local potential
   ! rc  core radii
   ! vkb  VKB projectors
   ! evkb  coefficients of VKB projectors
   ! nproj  number of vkb projectors for each l
   ! rr  log radial grid
   ! vpuns  unscreened semi-local pseudopotentials (vp(:,5) is local potential
   !   if linear combination is used)
   ! rho  valence pseudocharge
   ! rhomod  model core charge
   ! zz  atomic number
   ! zion  at this point, total valence charge (becomes psuedoion charge)
   ! mmax  size of log radial grid
   ! mxprj  dimension of number of projectors
   ! iexc  type of exchange-correlation
   ! icmod  1 if model core charge is used, otherwise 0
   ! drl spacing of linear radial grid
   ! atsym  atomic symbol
   ! epstot  pseudoatom total energy
   ! psfile  should be 'upf'
   ! remaining input variables to be echoed:
   !   na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf
   !   fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact
   !   epsh1,epsh2,depsh,rlmax,psfile
   ! uupsa  pseudo-atomic orbital array
   ! ea  psuedo-orbital eigenvalues
   use, intrinsic :: iso_fortran_env, only: stdout => output_unit
   implicit none
   integer, parameter :: dp=kind(1.0d0)
   real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

   ! Input variables
   integer, intent(in) :: unit
   integer, intent(in) :: lmax,lloc,iexc,mmax,mxprj,icmod
   integer, intent(in) :: nproj(6)
   real(dp), intent(in) :: drl,fcfact,rcfact,zz,zion,epstot
   real(dp), intent(in) :: rr(mmax),drr(mmax),vpuns(mmax,5),rho(mmax),vkb(mmax,mxprj,4)
   real(dp), intent(in) :: rhomod(mmax,5)
   real(dp), intent(in) :: rc(6),evkb(mxprj,4)
   character*2, intent(in) :: atsym
   real(dp), intent(in) :: uupsa(mmax,nv)
   character*16, intent(in) :: upfgrid

   ! additional input for upf output to echo input file, all as defined
   ! in the main progam
   integer, intent(in) :: na(30),la(30),ncon(6),nbas(6)
   integer, intent(in) :: nvcnf(5),nacnf(30,5),lacnf(30,5)
   integer, intent(in) :: nc,nv,lpopt,ncnf
   real(dp), intent(in) :: fa(30),rc0(6),ep(6),qcut(6),debl(6),facnf(30,5),ea(30)
   real(dp), intent(in) :: dvloc0,epsh1,epsh2,depsh,rlmax
   character*4, intent(in) :: psfile

   ! Output variables - printing only

   ! Local variables
   integer :: ii, jj, l1, iproj, ntotproj, nrproj(mxprj)
   integer :: dtime(8)
   real(dp), allocatable :: dmat(:,:)
   character*5 :: lnames
   character*2 :: pspd(3)
   real(dp) :: xmin, dx

   lnames = "SPDFG"

   call date_and_time(VALUES=dtime)
   ii=dtime(1)-2000
   if(ii<10) then
      write(pspd(1),'(a,i1)') '0',ii
   else
      write(pspd(1),'(i2)') ii
   end if
   ii=dtime(2)
   if(ii<10) then
      write(pspd(2),'(a,i1)') '0',ii
   else
      write(pspd(2),'(i2)') ii
   end if
   ii=dtime(3)
   if(ii<10) then
      write(pspd(3),'(a,i1)') '0',ii
   else
      write(pspd(3),'(i2)') ii
   end if

   ! calculate total number of projectors and find cutoff indices
   ! for each angular momentum
   ntotproj = 0
   nrproj(:) = 0
   do l1 = 1, lmax + 1
      if (l1 /= lloc + 1) then
         ntotproj = ntotproj + nproj(l1)
         do ii = 1, mmax
            if (rr(ii) >= rc(l1)) then
               nrproj(l1) = ii
               exit
            end if
         end do
      end if
   end do

   allocate(dmat(ntotproj,ntotproj))

   ! section for upf output for pwscf
   write(stdout,'(/a)') 'Begin PSP_UPF'

   write(unit,'(a)') '<UPF version="2.0.1">'

   write(unit,'(a)') '  <PP_INFO>'

   write(unit, '(/t2,a)') 'This pseudopotential file has been produced using the code'
   write(unit, '(t2,a)') 'ONCVPSP (Optimized Norm-Conserving Vanderbilt PSeudopotential)'
   write(unit, '(t2,a)') 'scalar-relativistic version 4.0.1 06/20/2107 by D. R. Hamann'
   write(unit, '(t2,a)') 'Mat-Sim Research LLC and Rutgers University'
   write(unit, '(/t2,a)') 'Documentation with the package provides a full discription of the'
   write(unit, '(t2,a/)') 'input data below.'
   write(unit,'(/t2,a)') 'While it is not required under the terms of the GNU GPL, it is'
   write(unit,'(t2,a)') 'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)'
   write(unit,'(t2,a/)') 'in any publication using these pseudopotentials.'

   write(unit,'(/a)') '    <PP_INPUTFILE>'
   ! output printing (echos input data)
   write(unit,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
   write(unit,'(a)') '# atsym  z   nc   nv     iexc    psfile'
   write(unit,'(a,a,f6.2,2i5,i8,2a)') '  ', trim(atsym), zz, nc, nv, iexc, '      ', trim(psfile)
   write(unit,'(a)') '#'
   write(unit,'(a)') '#   n    l    f        energy (Ha)'
   do ii=1,nc+nv
      write(unit,'(2i5,f8.2)') na(ii), la(ii), fa(ii)
   end do
   write(unit, '(a)') '#'
   write(unit, '(a)') '# PSEUDOPOTENTIAL AND OPTIMIZATION'
   write(unit, '(a)') '# lmax'
   write(unit, '(i5)') lmax
   write(unit, '(a)') '#'
   write(unit, '(a)') '#   l,   rc,      ep,       ncon, nbas, qcut'
   do l1 = 1, lmax + 1
      write(unit,'(i5,2f10.5,2i5,f10.5)') l1 - 1, rc0(l1), ep(l1), ncon(l1), nbas(l1), qcut(l1)
   end do
   write(unit,'(a)') '#'
   write(unit,'(a)') '# LOCAL POTENTIAL'
   write(unit,'(a)') '# lloc, lpopt,  rc(5),   dvloc0'
   write(unit,'(2i5,f10.5,a,f10.5)') lloc, lpopt, rc0(5), '   ', dvloc0
   write(unit,'(a)') '#'
   write(unit,'(a)') '# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs'
   write(unit,'(a)') '# l, nproj, debl'
   do l1 = 1, lmax + 1
      write(unit,'(2i5,f10.5)') l1 - 1, nproj(l1), debl(l1)
   end do
   write(unit,'(a)') '#'
   write(unit,'(a)') '# MODEL CORE CHARGE'
   write(unit,'(a)') '# icmod, fcfact, rcfact'
   write(unit,'(i5,2f10.5)') icmod,fcfact,rcfact
   write(unit,'(a)') '#'
   write(unit,'(a)') '# LOG DERIVATIVE ANALYSIS'
   write(unit,'(a)') '# epsh1, epsh2, depsh'
   write(unit,'(3f8.2)') epsh1,epsh2,depsh
   write(unit,'(a)') '#'
   write(unit,'(a)') '# OUTPUT GRID'
   write(unit,'(a)') '# rlmax, drl'
   write(unit,'(2f8.2)') rlmax,drl
   write(unit,'(a)') '#'
   write(unit,'(a)') '# TEST CONFIGURATIONS'
   write(unit,'(a)') '# ncnf'
   write(unit,'(i5)') ncnf
   write(unit,'(a)') '# nvcnf'
   write(unit,'(a)') '#   n    l    f'
   do jj = 2, ncnf + 1
      write(unit,'(i5)') nvcnf(jj)
      do ii = nc + 1, nc + nvcnf(jj)
         write(unit,'(2i5,f8.2)') nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
      end do
      write(unit,'(a)') '#'
   end do
   write(unit,'(t4,a)') '</PP_INPUTFILE>'

   write(unit, '(t2,a)') '</PP_INFO>'

   write(unit, '(t2,a)') '<!--                               -->'
   write(unit, '(t2,a)') '<!-- END OF HUMAN READABLE SECTION -->'
   write(unit, '(t2,a)') '<!--                               -->'

   write(unit,'(t4,a)') '<PP_HEADER'
   write(unit,'(t8,a)') 'generated="Generated using ONCVPSP code by D. R. Hamann"'
   write(unit,'(t8,a)') 'author="anonymous"'
   write(unit,'(t8,5a)') 'date="', pspd(:), '"'
   write(unit,'(t8,a)') 'comment=""'
   write(unit,'(t8,a,a2,a)') 'element="', atsym, '"'
   write(unit,'(t8,a)') 'pseudo_type="NC"'
   write(unit,'(t8,a)') 'relativistic="scalar"'
   write(unit,'(t8,a)') 'is_ultrasoft="F"'
   write(unit,'(t8,a)') 'is_paw="F"'
   write(unit,'(t8,a)') 'is_coulomb="F"'
   write(unit,'(t8,a)') 'has_so="F"'
   write(unit,'(t8,a)') 'has_wfc="F"'
   write(unit,'(t8,a)') 'has_gipaw="F"'
   if(icmod>=1) then
      write(unit,'(t8,a)') 'core_correction="T"'
   else
      write(unit,'(t8,a)') 'core_correction="F"'
   end if
   if(iexc==2 .or. iexc==-000004) then
      write(unit,'(t8,a)') 'functional="NOX+HL"'
   else if(iexc==3 .or. iexc==-001009) then
      write(unit,'(t8,a)') 'functional="PZ"'
   else if(iexc==4 .or. iexc==-101130) then
      write(unit,'(t8,a)') 'functional="PBE"'
   else if(iexc==-001012) then
      write(unit,'(t8,a)') 'functional="SLA  PW   NOGX NOGC"'
   else if(iexc==-109134) then
      write(unit,'(t8,a)') 'functional="PW91"'
   else if(iexc==-116133) then
      write(unit,'(t8,a)') 'functional="PBESOL"'
   else if(iexc==-102130) then
      write(unit,'(t8,a)') 'functional="REVPBE"'
   else if(iexc==-106132) then
      write(unit,'(t8,a)') 'functional="BP"'
   else if(iexc==-106131) then
      write(unit,'(t8,a)') 'functional="BLYP"'
   else if(iexc==-118130) then
      write(unit,'(t8,a)') 'functional="WC"'
   else
      write(unit,'(t8,a,i0,a)') 'upfout: ERROR iexc = ', iexc, ' is presently unsupported for UPF output'
      stop
   end if
   write(unit,'(t8,a,f8.2,a)') 'z_valence="', zion, '"'
   write(unit,'(t8,a,1p,e20.11,a)') 'total_psenergy="', 2 * epstot, '"'
   write(unit,'(t8,a,1p,e20.11,a)') 'rho_cutoff="', rr(mmax), '"'
   write(unit,'(t8,a,i1,a)') 'l_max="', lmax, '"'
   if(lloc==4) then
      write(unit,'(t8,a)') 'l_local="-1"'
   else
      write(unit,'(t8,a,i1,a)') 'l_local="',lloc,'"'
   end if
   write(unit,'(t8,a,i6,a)') 'mesh_size="',mmax,'"'
   write(unit,'(t8,a,i1,a)') 'number_of_wfc="',nv,'"'
   write(unit,'(t8,a,i1,a)') 'number_of_proj="',ntotproj,'"/>' !end of PP_HEADER

   if (trim(upfgrid) == "log") then
      ! From UPF docs:
      !    r(i) = exp(xmin + i * dx / zmesh) = exp(xmin) * exp(i * dx / zmesh)
      ! r(2) / r(1) = exp(xmin) * exp(dx / zmesh)^2 / (exp(xmin) * exp(dx / zmesh)^1)
      !    => dx = log(r(2) / r(1)) * zmesh
      ! r(1) = exp(xmin) * exp(dx / zmesh)
      !    => xmin = log(r(1)) - dx / zmesh
      ! In ld1.x, zmesh = atomic number; we do the same here
      dx = log(rr(2) / rr(1)) * zz
      xmin = log(rr(1)) - dx / zz
      write(unit,'(t2,a,a,e20.10,a,i0,a,e20.10,a,e20.10,a,e20.10,a)') &
      '<PP_MESH dx="', dx, '" mesh="', mmax, '" xmin="', xmin, '" rmax="', rr(mmax), '" zmesh="', zz, '">'
   else
      write(unit,'(t2,a)') '<PP_MESH>'
   end if
   write(unit,'(t4,a,i4,a)') '<PP_R type="real"  size="', mmax, '" columns="8">'
   write(unit,'(8f10.4)') (rr(ii), ii=1, mmax)
   write(unit,'(t4,a)') '</PP_R>'
   write(unit,'(t4,a,i4,a)') '<PP_RAB type="real"  size="', mmax, '" columns="8">'
   write(unit,'(8f10.4)') (drr, ii=1, mmax)
   write(unit,'(t4,a)') '</PP_RAB>'
   write(unit,'(t2,a)') '</PP_MESH>'

   ! write local potential with factor of 2 for Rydberg units
   write(unit,'(t2,a,i4,a)') '<PP_LOCAL type="real"  size="', mmax, '" columns="4">'
   write(unit,'(1p,4e20.10)') (2.0_dp * vpuns(ii, lloc + 1), ii=1, mmax)
   write(unit,'(t2,a)') '</PP_LOCAL>'

   ! loop on angular mommentum for projector outputs
   write(unit,'(t2,a)') '<PP_NONLOCAL>'
   dmat(:, :) = 0.0_dp
   iproj = 0
   do l1 = 1, lmax + 1
      if (l1 == lloc + 1) cycle
      do jj = 1, nproj(l1)
         iproj = iproj + 1
         dmat(iproj, iproj) = 2.0_dp * evkb(jj, l1) !2 for Rydbergs
         if (iproj <= 9) then
            write(unit,'(t4,a,i1)') '<PP_BETA.', iproj
         else
            write(unit,'(t4,a,i2)') '<PP_BETA.', iproj
         end if
         write(unit,'(t8,a)') 'type="real"'
         write(unit,'(t8,a,i4,a)') 'size="', mmax, '"'
         write(unit,'(t8,a)') 'columns="4"'
         write(unit,'(t8,a,i1,a)') 'index="', iproj, '"'
         write(unit,'(t8,a,i1,a)') 'angular_momentum="', l1 - 1, '"'
         write(unit,'(t8,a,i4,a)') 'cutoff_radius_index="', nrproj(l1), '"'
         write(unit,'(t8,a,1p,e20.10,a)') 'cutoff_radius="', rc(l1), '" >'
         write(unit,'(1p,4e20.10)') (vkb(ii, jj, l1), ii = 1, mmax)
         if(iproj<=9) then
            write(unit,'(t4,a,i1,a)') '</PP_BETA.', iproj, '>'
         else
            write(unit,'(t4,a,i2,a)') '</PP_BETA.', iproj, '>'
         end if
      end do
   end do

   write(unit,'(t4,a,i4,a)') '<PP_DIJ type="real"  size="', ntotproj**2, '" columns="4">'
   write(unit,'(1p,4e20.10)') ((dmat(ii, jj), ii = 1, ntotproj), jj = 1, ntotproj)
   write(unit,'(t4,a)') '</PP_DIJ>'
   write(unit,'(t2,a)') '</PP_NONLOCAL>'
   write(unit,'(t2,a)') '<PP_PSWFC>'
   do ii=1,nv
      l1 = la(nc+ii)
      if(ii <= 9) then
         write(unit,'(t4,a,i1)') '<PP_CHI.', ii
      else
         write(unit,'(t4,a,i2)') '<PP_CHI.', ii
      end if
      write(unit,'(t8,a)') 'type="real"'
      write(unit,'(t8,a,i4,a)') 'size="', mmax, '"'
      write(unit,'(t8,a)') 'columns="4"'
      if (ii <= 9) then
         write(unit,'(t8,a,i1,a)') 'index="', ii, '"'
      else
         write(unit,'(t8,a,i2,a)') 'index="', ii, '"'
      end if
      write(unit,'(t8,a,f6.3,a)') 'occupation="', fa(nc + ii), '"'
      write(unit,'(t8,a,e20.10,a)') 'pseudo_energy="', 2 * ea(nc + ii), '"'
      write(unit,'(t8,a,i1,a,a)') 'label="', na(nc + ii), lnames(l1 + 1:l1 + 1), '"'
      write(unit,'(t8,a,i1,a)') 'l="', l1, '" >'
      write(unit,'(1p,4e20.10)') (uupsa(jj, ii), jj = 1, mmax)
      if (ii <= 9) then
         write(unit,'(t4,a,i1,a)') '</PP_CHI.', ii, '>'
      else
         write(unit,'(t4,a,i2,a)') '</PP_CHI.', ii, '>'
      end if
   end do
   write(unit,'(t2,a)') '</PP_PSWFC>'

   if (icmod >= 1) then
      write(unit,'(t2,a,i4,a)') '<PP_NLCC type="real"  size="', mmax, '" columns="4">'
      write(unit,'(1p,4e20.10)') (rhomod(ii, 1) / (4.0_dp * pi), ii = 1, mmax)
      write(unit,'(t2,a)') '</PP_NLCC>'
   end if

   write(unit,'(t2,a,i4,a)') '<PP_RHOATOM type="real"  size="', mmax, '" columns="4">'
   write(unit,'(1p,4e20.10)') ((rr(ii)**2) * rho(ii), ii = 1, mmax)
   write(unit,'(t2,a)') '</PP_RHOATOM>'

   write(unit,'(a)') '</UPF>'
   deallocate(dmat)
   ! write termination flag
   write(stdout,'(/a)') 'END_PSP'

   return
end subroutine upfout
