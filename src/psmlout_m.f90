!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! Copyright (c) 2015 by Alberto Garcia, ICMAB-CSIC for xmlf90-wxml calls
! in support of the PSML format
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
! Support for PSML file creation
!
 module psmlout_m
 implicit none

 public :: psmlout
 private

 character(len=1), dimension(0:4), parameter :: lsymb = ['s','p','d','f','g']

 contains

 subroutine psmlout(inp,vkb,evkb,rr,vpuns,rho,rhomod,rcore,zion,mmax, &
&                   nrl,epstot,rc0,fcfact,fcgrid,rcgrid,log_mesh,irel,c_light, &
&                   nkap,iout)

!! writes the pseudopotential to the output in PSML format, between
!! 'Begin PSML' and 'END_PSP' markers
!!
!! @param inp  the input values (input_m)
!! @param vkb  VKB projectors, last index 1 for j=l+1/2, 2 for j=l-1/2 when
!!             nkap=2
!! @param evkb  coefficients of VKB projectors, last index as for vkb
!! @param rr  log radial grid
!! @param vpuns  unscreened semi-local pseudopotentials (vp(:,5,1) is local
!!               potential if linear combination is used), last index as for vkb
!! @param rho  valence pseudocharge
!! @param rhomod  model core charge and its first four derivatives
!! @param rcore  radius beyond which the model core charge is the core
!!               charge: its joint (doc/smoothness_plan.md, M-14)
!! @param zion  at this point, total valence charge (becomes pseudoion charge)
!! @param mmax  size of log radial grid
!! @param nrl  size of linear radial grid (not used)
!! @param epstot  pseudoatom total energy (not used)
!! @param rc0  core radii as read from the input (bohr), rc0(5) for lloc=4
!! @param fcfact  prefactor for model core charge amplitude, as read from the
!!                input
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param log_mesh  amesh, rr1, rmax of the log mesh, echoed when the input sets
!!                  it (else 0)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!! @param iout  unit of the text output
!!
!! Alberto Garcia, February 1, 2015

 use xmlf90_wxml, only: xmlf_t, str, xml_OpenUnit, xml_Close, &
&                       xml_AddXMLDeclaration, xml_NewElement, &
&                       xml_EndElement, xml_AddArray
 use input_echo_m, only: fe
 use m_libxc_list, only: libxc_t, XC_NOT_IMPL, XC_LDA_X, XC_LDA_C_PZ, &
&                        XC_LDA_C_WIGNER, XC_LDA_C_HL, XC_GGA_X_PBE, &
&                        XC_GGA_C_PBE, XC_LDA_C_VWN, XC_LDA_X_REL
 use constants_m, only: dp, mxl, mxlin, nrhomod, mxprj

 use constants_m, only: pi
 use input_m, only: input_t


!Input variables
 type(input_t), intent(in) :: inp
 integer, intent(in) :: mmax,nrl,irel,nkap
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: fcfact,zion,epstot,rcore
 real(dp), target, intent(in) :: rr(mmax),vpuns(mmax,mxl+2,nkap),rho(mmax)
 real(dp), target, intent(in) :: vkb(mmax,mxprj,mxl+1,nkap)
 real(dp), target, intent(in) :: rhomod(mmax,nrhomod)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)

!the core radii as read (rc0), for the echo
 real(dp), intent(in) :: rc0(mxlin+1)
 real(dp), intent(in) :: fcgrid(3),rcgrid(3),log_mesh(3)
 integer, intent(in) :: iout

!Local variables
 type(xmlf_t) :: xf

 integer   :: npts,nr,ncp,norbs,npots
 real(dp)  :: total_valence_charge
 real(dp), pointer :: rgrid(:)
 real(dp), allocatable :: r0(:), f0(:)

 character(len=40):: psflavor

!nsh, lsh, fsh  n, l and occupation of the reference configuration states
!nn, ll, ff  n, l and occupation that label each semi-local potential
 integer, allocatable  :: nsh(:), lsh(:)
 integer, allocatable  :: nn(:), ll(:)
 real(dp), allocatable :: fsh(:), ff(:)

!---

! the document goes to the output, between markers, as do the psp8 and
! upf files
 write(iout,'(/a)') 'Begin PSML'
 call xml_OpenUnit(6,xf, indent=.false.)

 call xml_AddXMLDeclaration(xf,"UTF-8")

 call xml_NewElement(xf,"psml")
 call my_add_attribute(xf,"version","0.8")
 call my_add_attribute(xf,"energy_unit","hartree")
 call my_add_attribute(xf,"length_unit","bohr")

 call write_provenance()
 call label_potentials()

 psflavor ="Hamann's oncvpsp"

 call write_header()

!AG -- decide how to handle the case of Vlocal as one of the sl pots.


!AG: decide whether to use a single mesh (Hamann's own) or have the
! projectors use another one (linear, shorter)

 npts = mmax
 rgrid => rr(:)
 nr = npts + 1
 call write_grid()

 allocate(f0(nr))
 call write_semilocal()
 call write_charges()

 call xml_NewElement(xf,"pseudopotential-operator")
 call my_add_attribute(xf,"version","0.1")
 call my_add_attribute(xf,"energy_unit","hartree")
 call my_add_attribute(xf,"length_unit","bohr")

 call xml_NewElement(xf,"provenance")
 call my_add_attribute(xf,"creator","oncvpsp 2.1.2+psml")
 call xml_EndElement(xf,"provenance")

 call write_local()
 call write_projectors()

 call xml_EndElement(xf,"pseudopotential-operator")

 call xml_EndElement(xf,"psml")


 call xml_Close(xf)
 write(iout,'(/a)') 'END_PSP'

 deallocate(f0,r0)
 return

 contains

 subroutine write_provenance()
!! the provenance element: the code, the date and the input file
 integer :: dtime(8)
 character(len=10)     :: datestr

 call date_and_time(VALUES=dtime)
 write(datestr,"(i4,'-',i2.2,'-',i2.2)") dtime(1:3)

 call xml_NewElement(xf,"provenance")
 call my_add_attribute(xf,"creator","ONCVPSP-2.1.2+psml")
 call my_add_attribute(xf,"date",datestr)
 call xml_NewElement(xf,"input-file")
 call my_add_attribute(xf,"name","oncvpsp-input")
!
 call echo_input(xf,inp%atsym,inp%zz,inp%nc,inp%nv,inp%iexc,c_light, &
&                inp%psfile,inp%na,inp%la,inp%fa,inp%lmax,rc0,inp%ep,inp%ncon, &
&                inp%nbas,inp%qcut,inp%lloc,inp%lpopt,inp%dvloc0,inp%nproj, &
&                inp%debl,inp%icmod,fcfact,inp%rcfact,inp%slack, &
&                inp%budget_cut,fcgrid,rcgrid,inp%epsh1,inp%epsh2,inp%depsh, &
&                inp%rxpsh,inp%rlmax,inp%drl,inp%ncnf,inp%cnf,log_mesh, &
&                inp%upf_mesh,inp%ep_as_given,inp%well_onset,inp%well_tail, &
&                inp%wreg,inp%well_n,inp%well_norm,inp%ecut_wt)

 call xml_EndElement(xf,"input-file")
!
 call xml_EndElement(xf,"provenance")
 end subroutine write_provenance

 subroutine label_potentials()
!! the shells of the reference configuration and the valence charge; and
!! the n, l and occupation that label each semi-local potential: those of
!! the valence state of its l, else n above the core states of l and no
!! occupation
 integer :: ii,jj
 logical :: found
 norbs = inp%nc + inp%nv
 allocate(nsh(norbs), lsh(norbs), fsh(norbs))

 total_valence_charge = 0.0d0
 ncp = inp%nc + 1
 do ii = 1, norbs
   nsh(ii) = inp%na(ii)
   lsh(ii) = inp%la(ii)
   fsh(ii) = inp%fa(ii)
   if(ii > inp%nc) then
     total_valence_charge =   total_valence_charge + fsh(ii)
   end if
 end do

 npots = inp%lmax + 1
 allocate(ll(npots), nn(npots), ff(npots))
 do ii = 1, npots
   ll(ii) = ii - 1
   found = .false.
! look for the appropriate shell in the valence
   do jj = ncp, norbs
     if(lsh(jj) == ll(ii)) then
       found = .true.
       nn(ii) = nsh(jj)
       ff(ii) = fsh(jj)
       exit
     end if
   end do
   if(.not. found) then
! generate the appropriate effective n
     nn(ii) = ll(ii) + 1
     do jj = 1, inp%nc
       if(lsh(jj) == ll(ii)) then
         nn(ii) = nn(ii) + 1
       end if
       ff(ii) = 0.0d0
     end do
   end if
 end do
 end subroutine label_potentials

 subroutine write_header()
!! the header element: the atom, the relativity, the model core, the
!! exchange-correlation functional (its libxc functionals) and the valence
!! configuration
 type(libxc_t) :: libxc_id(2)
 character(len=30) :: xcfuntype, xcfunparam
 character(len=4)      :: coreattrib
 integer :: ii

 if(inp%icmod >= 1) then
   coreattrib = "yes"
 else
   coreattrib = "no"
 end if

! XC name handling

 select case(inp%iexc)

 case(1)
   xcfuntype    = 'LDA'
   xcfunparam   = 'Wigner'
   libxc_id = [ XC_LDA_X, XC_LDA_C_WIGNER ]
 case(2)
   xcfuntype    = 'LDA'
   xcfunparam   = 'Hedin-Lundqvist'
   libxc_id = [ XC_LDA_X, XC_LDA_C_HL ]
 case(3)
   xcfuntype    = 'LDA'
   xcfunparam   = 'Ceperley-Alder PZ'
   libxc_id = [ XC_LDA_X, XC_LDA_C_PZ ]

 case(4)
   xcfuntype    = 'GGA'
   xcfunparam   = 'Perdew-Burke-Ernzerhof'
   libxc_id = [ XC_GGA_X_PBE, XC_GGA_C_PBE ]

 case(6)
   xcfuntype    = 'LDA'
   xcfunparam   = 'Vosko-Wilk-Nusair'
   libxc_id = [ XC_LDA_X, XC_LDA_C_VWN ]
!Libxc's lda_x_rel has its own speed of light
 case(7)
   xcfuntype    = 'LDA'
   xcfunparam   = 'MacDonald-Vosko + VWN'
   libxc_id = [ XC_LDA_X_REL, XC_LDA_C_VWN ]

 case  default

   xcfuntype    = '-----'
   xcfunparam   = '-----'
   libxc_id = [  XC_NOT_IMPL, XC_NOT_IMPL ]   ! ???

 end select
!
!
 call xml_NewElement(xf,"header")
 call my_add_attribute(xf,"atomic-label",inp%atsym)
 call my_add_attribute(xf,"atomic-number",str(inp%zz))
 call my_add_attribute(xf,"z-pseudo",str(zion))
 call my_add_attribute(xf,"flavor",psflavor)
 if(irel==2) then
   call my_add_attribute(xf,"relativity","dirac")
 else if(irel==1) then
   call my_add_attribute(xf,"relativity","scalar")
 else
   call my_add_attribute(xf,"relativity","no")
 end if
 call my_add_attribute(xf,"polarized","no")
 call my_add_attribute(xf,"core-corrections",coreattrib)

 call xml_NewElement(xf,"exchange-correlation")
 call xml_NewElement(xf,"annotation")
 call my_add_attribute(xf,"oncvpsp-xc-code",str(inp%iexc))
 call my_add_attribute(xf,"oncvpsp-xc-type",trim(xcfuntype))
 call my_add_attribute(xf,"oncvpsp-xc-authors",trim(xcfunparam))
 call xml_EndElement(xf,"annotation")

 call xml_NewElement(xf,"libxc-info")
 call my_add_attribute(xf,"number-of-functionals","2")
 do ii = 1, 2
   call xml_NewElement(xf,"functional")
   call my_add_attribute(xf,"name",trim(libxc_id(ii)%name))
   call my_add_attribute(xf,"type",trim(libxc_id(ii)%xc_kind%str))
   call my_add_attribute(xf,"id",str(libxc_id(ii)%code))
   call xml_EndElement(xf,"functional")
 end do
 call xml_EndElement(xf,"libxc-info")
 call xml_EndElement(xf,"exchange-correlation")
!

!
 call do_configuration()
 call xml_EndElement(xf,"header")
 end subroutine write_header

 subroutine do_configuration()

!! writes the valence-configuration element of the header (not
!! polarized)
 integer :: ii

 call xml_NewElement(xf,"valence-configuration")
 call my_add_attribute(xf,"total-valence-charge", str(total_valence_charge))
 do ii = ncp, norbs
   if(fsh(ii) < 1.0d-10) cycle
   call xml_NewElement(xf,"shell")
   call my_add_attribute(xf,"n",str(nsh(ii)))
   call my_add_attribute(xf,"l",lsymb(lsh(ii)))
   call my_add_attribute(xf,"occupation",str(fsh(ii)))
   call xml_EndElement(xf,"shell")
 end do
 call xml_EndElement(xf,"valence-configuration")
 return
 end subroutine do_configuration

 subroutine write_grid()
!! the grid element: the log mesh with r=0 before it (r0)
 allocate(r0(nr))
 call add_zero_r(rgrid,rgrid,r0)
 r0(1) = 0.0d0

 call xml_NewElement(xf,"grid")
 call my_add_attribute(xf,"npts",str(nr))

 call xml_NewElement(xf,"annotation")
 call my_add_attribute(xf,"type","oncvpsp plus r=0")
 call my_add_attribute(xf,"oncvpsp-npts",str(npts))
 call my_add_attribute(xf,"oncvpsp-r1",trim(str(rgrid(1))))
 call my_add_attribute(xf,"oncvpsp-factor",trim(str(rgrid(2)/rgrid(1))))
 call xml_EndElement(xf,"annotation")

 call xml_NewElement(xf,"grid-data")
 call xml_AddArray(xf,r0(1:nr))
 call xml_EndElement(xf,"grid-data")

 call xml_EndElement(xf,"grid")
 end subroutine write_grid

 subroutine write_semilocal()
!! the semilocal-potentials element: of each l and j when fully
!! relativistic (last index of vpuns 1 for j=l+1/2, 2 for j=l-1/2, l=0,
!! j=0 in 1), else of each l
 integer :: ii,jj,l1
 real(dp)  :: jval
 real(dp), allocatable :: vps(:)
 allocate(vps(npts))
 if(irel==2) then
   call xml_NewElement(xf,"semilocal-potentials")
   call my_add_attribute(xf,"set","lj")

   do ii = 1, npots
     l1 = ii
     if(ii == 1) then
       jval = 0.0
       vps(:) = vpuns(:,l1,1)
       call add_zero_r(vps(1:npts),rgrid,f0)
       call write_psml_item(xf, class="slps", &
&                     nn=nn(ii), ll=ll(ii), jj=jval, &
&                     rc=inp%rc(ii),flavor=psflavor,&
&                     ff=f0)
     else
       do jj=1,2
         jval = ll(ii) + (3-2*jj)*0.5d0  ! convert (1,2) to (1,-1)*1/2
         vps(:) = vpuns(:,l1,jj)
         call add_zero_r(vps(1:npts),rgrid,f0)
         call write_psml_item(xf, class="slps", &
&                    nn=nn(ii), ll=ll(ii), jj=jval, &
&                    rc=inp%rc(ii),flavor=psflavor,&
&                    ff=f0)
       end do
     end if
   end do
   call xml_EndElement(xf,"semilocal-potentials")
 else
   call xml_NewElement(xf,"semilocal-potentials")
   if(irel==1) then
     call my_add_attribute(xf,"set","scalar_relativistic")
   else
     call my_add_attribute(xf,"set","non_relativistic")
   end if
   vpsd: do ii = 1, npots
     vps(:) = vpuns(:,ii,1)
     call add_zero_r(vps(1:npts),rgrid,f0)
     call write_psml_item(xf, class="slps", &
&                         nn=nn(ii), ll=ll(ii), &
&                         rc=inp%rc(ii),flavor=psflavor,&
&                         ff=f0)
   end do vpsd
   call xml_EndElement(xf,"semilocal-potentials")
 end if
 deallocate(vps)
 end subroutine write_semilocal

 subroutine write_charges()
!! the valence-charge element, and with a model core the
!! pseudocore-charge element
 call xml_NewElement(xf,"valence-charge")
 call my_add_attribute(xf,"total-charge",  &
&                    str(total_valence_charge))
 call xml_NewElement(xf,"radfunc")

 call xml_NewElement(xf,"data")
 call add_zero_r(rho(1:npts),rgrid,f0)
 call xml_AddArray(xf,fe(f0(1:nr)))
 call xml_EndElement(xf,"data")
 call xml_EndElement(xf,"radfunc")
 call xml_EndElement(xf,"valence-charge")


 if(inp%icmod >= 1) then
   call xml_NewElement(xf,"pseudocore-charge")
   call my_add_attribute(xf,"matching-radius",str(rcore))
! icmod 6 and 7 match the core charge to C2 at r_s
   if(inp%icmod==6 .or. inp%icmod==7) then
     call my_add_attribute(xf,"number-of-continuous-derivatives",str(2))
   else
     call my_add_attribute(xf,"number-of-continuous-derivatives",str(4))
   end if
   if(inp%icmod==1) then
     call my_add_attribute(xf,"annotation",  &
&                 "Monotonic 8th-order polynomial with no linear term")
   else if(inp%icmod==6 .or. inp%icmod==7) then
     call my_add_attribute(xf,"annotation",  &
&                 "Budgeted fill inside the matching radius")
   else
     call my_add_attribute(xf,"annotation","Teter function")
   end if
   call xml_NewElement(xf,"radfunc")

   call xml_NewElement(xf,"data")
   call add_zero_r(rhomod(1:npts,1),rgrid,f0)
   call xml_AddArray(xf,fe(f0(1:nr)))
   call xml_EndElement(xf,"data")
   call xml_EndElement(xf,"radfunc")
   call xml_EndElement(xf,"pseudocore-charge")
 end if
 end subroutine write_charges

 subroutine write_local()
!! the local-potential element
 call xml_NewElement(xf,"local-potential")
 if(inp%lloc > inp%lmax) then
   call my_add_attribute(xf,"type","oncv-fit")
 else
   call my_add_attribute(xf,"type","l="//str(inp%lloc))
 end if
 call xml_NewElement(xf,"radfunc")
 call xml_NewElement(xf,"data")
 call add_zero_r(vpuns(1:npts,inp%lloc+1,1),rgrid,f0)
 call xml_AddArray(xf,fe(f0(1:nr)))
 call xml_EndElement(xf,"data")
 call xml_EndElement(xf,"radfunc")
 call xml_EndElement(xf,"local-potential")
 end subroutine write_local

 subroutine write_projectors()
!! the projectors element: of each l and j when fully relativistic (l=0
!! with j=0), else of each l
 integer :: jj,l1,jk
 real(dp)  :: jval
 if(irel==2) then
   call xml_NewElement(xf,"projectors")
   call my_add_attribute(xf,"set","lj")

   do l1=1,inp%lmax+1
     if(l1==inp%lloc+1) cycle

     if(l1 == 1) then
! l=0, only one j=0
       jval = 0.0

       do jj=1,inp%nproj(l1)
         call add_zero_r(vkb(:,jj,l1,1),rgrid,f0)
         call write_psml_item(xf, class="proj", &
&                          seq=jj, ll=l1-1, jj=jval, &
&                          ekb=evkb(jj,l1,1), &
&                          type="oncv", ff=f0)
       end do

     else  ! l1 /=1  (l/=0)

! two j values
       do jk=1,2
         jval = l1-1 + (3-2*jk)*0.5d0  ! convert (1,2) to (1,-1)*1/2
         do jj=1,inp%nproj(l1)
           call add_zero_r(vkb(:,jj,l1,jk),rgrid,f0)
           call write_psml_item(xf, class="proj", &
&                          seq=jj, ll=l1-1, jj=jval, &
&                          ekb=evkb(jj,l1,jk), &
&                          type="oncv", ff=f0)
         end do
       end do

     end if  ! l1 == 1

   end do ! over l shells
   call xml_EndElement(xf,"projectors")
 else
   call xml_NewElement(xf,"projectors")
   if(irel==1) then
     call my_add_attribute(xf,"set","scalar_relativistic")
   else
     call my_add_attribute(xf,"set","non_relativistic")
   end if

   do l1=1,inp%lmax+1
     if(l1==inp%lloc+1) cycle
     do jj=1,inp%nproj(l1)
       call add_zero_r(vkb(:,jj,l1,1),rgrid,f0)
       call write_psml_item(xf, class="proj", &
&                           seq=jj, ll=l1-1, &
&                           ekb=evkb(jj,l1,1), &
&                           type="oncv", ff=f0)
     end do
   end do

   call xml_EndElement(xf,"projectors")
 end if
 end subroutine write_projectors

 end subroutine psmlout

 subroutine echo_input(xf,atsym,zz,nc,nv,iexc,c_light,psfile,na,la,fa,lmax, &
&                      rc0,ep,ncon,nbas,qcut,lloc,lpopt,dvloc0,nproj,debl, &
&                      icmod,fcfact,rcfact,slack,budget_cut,fcgrid,rcgrid, &
&                      epsh1,epsh2,depsh,rxpsh,rlmax,drl,ncnf,cnf,log_mesh, &
&                      upf_mesh,ep_as_given,well_onset,well_tail,wreg,well_n, &
&                      well_norm,ecut_wt)

!! echoes the input data into the input-file element of the provenance, a
!! line at a time, in the form of the psp8 echo (linout), so that it can be
!! run again
!!
!! @param xf  xml file handle
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param iexc  type of exchange-correlation
!! @param c_light  speed of light (Hartree atomic units)
!! @param lmax  maximum angular momentum
!! @param lloc  l for local potential
!! @param nproj  number of vkb projectors for each l
!! @param icmod  model core charge option
!! @param drl  spacing of linear radial grid
!! @param na  principal quantum numbers of the reference configuration states
!! @param la  angular momenta of the reference configuration states
!! @param ncon  number of matching constraints for each l
!! @param nbas  number of basis functions for each l
!! @param cnf  the reference and test configurations (input_m)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param lpopt  type of polynomial continuation for lloc=4
!! @param ncnf  number of test configurations
!! @param fa  occupations of the reference configuration states
!! @param rc0  core radii as read from the input (bohr), rc0(5) for lloc=4
!! @param ep  energies at which the pseudopotentials are generated (Ha)
!! @param qcut  residual-energy cutoff wave vectors for each l
!! @param debl  energy shifts for the second projectors for each l (Ha)
!! @param dvloc0  r=0 shift of the lloc=4 local potential (Ha)
!! @param fcfact  prefactor for model core charge amplitude, as read from the
!!                input
!! @param rcfact  prefactor for model core charge scale (icmod=3), or error
!!                tolerance (icmod=6)
!! @param slack  share of the cutoff the model core must fit (icmod=6, 7)
!! @param budget_cut  the cutoff of the pseudo-waves (icmod=6, 7): 0 the
!!                    largest of the ground-state waves, 1 their weighted
!!                    mean, 2 the largest of all waves (budget_cutoff)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param epsh1  lower limit of the log derivative analysis energy range (Ha)
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha)
!! @param depsh  energy step of the log derivative analysis (Ha)
!! @param rxpsh  radius of the log derivative analysis (bohr), negative for
!!               default
!! @param rlmax  maximum radius of the linear output mesh (bohr)
!! @param psfile  pseudopotential file format as read from the input
!! @param log_mesh  amesh, rr1, rmax of the log mesh, echoed when the input sets
!!                  it (else 0)
!! @param upf_mesh  most points of the log mesh of the UPF file, 0 for the
!!                  linear output mesh (echoed when not 0)
!! @param ep_as_given, well_onset, well_tail, wreg  echoed on the lmax
!!        line
!!                  when one is not its default (lmax_lines)
!! @param well_n  echoed on the projector lines when set (vkb_line)
!! @param well_norm  echoed on the lmax line when not 0 (lmax_lines)
!! @param ecut_wt  echoed on the projector lines when set (vkb_line)

 use xmlf90_wxml, only: xmlf_t, xml_AddPcData
 use input_echo_m, only: input_echo_lines, echo_len
 use constants_m, only: dp, mxstate, mxcnf, mxlin
 use input_m, only: config_t

!Input variables
 character(len=2), intent(in) :: atsym
 character(len=4), intent(in) :: psfile
 integer, intent(in) :: nc,nv,iexc,lmax,lloc,lpopt,icmod,ncnf,upf_mesh
 real(dp), intent(in) :: c_light
 integer, intent(in) :: na(mxstate),la(mxstate),ncon(mxlin+1),nbas(mxlin+1)
 integer, intent(in) :: nproj(mxlin+1)
 type(config_t), intent(in) :: cnf(mxcnf)
 real(dp), intent(in) :: zz,fa(mxstate),rc0(mxlin+1),ep(mxlin+1),qcut(mxlin+1)
 real(dp), intent(in) :: debl(mxlin+1)
 real(dp), intent(in) :: dvloc0,fcfact,epsh1,epsh2,depsh,rlmax,drl
 real(dp), intent(in) :: rcfact,fcgrid(3),rcgrid(3),rxpsh,log_mesh(3)
 real(dp), intent(in) :: slack
 integer, intent(in) :: budget_cut
 logical, intent(in) :: ep_as_given
 integer, intent(in) :: well_tail,well_n(mxlin+1),well_norm
 real(dp), intent(in) :: ecut_wt(mxlin+1)
 real(dp), intent(in) :: wreg
 real(dp), intent(in) :: well_onset

!Input/output variables
 type(xmlf_t), intent(inout) :: xf

!Local variables
 character(len=echo_len), allocatable :: text(:)
 integer :: ii,ntext

! the psp8 echo of the scalar-relativistic code
 call input_echo_lines(1,1,atsym,zz,nc,nv,iexc,c_light,psfile, &
&                      na,la,fa,lmax,rc0,ep,ncon,nbas,qcut,ep_as_given, &
&                      well_onset,well_tail,wreg,lloc,lpopt,dvloc0, &
&                      nproj,debl,well_n,well_norm,ecut_wt, &
&                      icmod,fcfact,rcfact,slack,budget_cut,fcgrid,rcgrid, &
&                      epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                      ncnf,cnf,log_mesh,upf_mesh, &
&                      text,ntext)
 do ii=1,ntext
   call add(text(ii))
 end do
 return

 contains

 subroutine add(text)

!! adds a line of text to the input-file element
!!
!! @param text  the line

!Input variables
 character(len=*), intent(in) :: text
 call xml_AddPcData(xf,trim(text),line_feed=.true.)
 return
 end subroutine add

 end subroutine echo_input

 subroutine my_add_attribute(xf,name,value)

!! adds an attribute to the current xml element, with trailing blanks
!! removed from its value
!!
!! @param xf  xml file handle
!! @param name  attribute name
!! @param value  attribute value

 use xmlf90_wxml, only: xmlf_t, xml_AddAttribute

!Input variables
 character(len=*), intent(in)  :: name
 character(len=*), intent(in)  :: value

!Input/output variables
 type(xmlf_t), intent(inout)   :: xf

 call xml_AddAttribute(xf,name,trim(value))
 return
 end subroutine my_add_attribute

 subroutine add_zero_r(ff,rr,f0)

!! Adds an r=0 element to a grid function f0: the interpolant in r of
!! the log-mesh values at r = 0 (lm_origin)
!!
!! @param ff  function on the log radial grid
!! @param rr  log radial grid
!! @param f0  ff with an extrapolated value at r=0 prepended

 use constants_m, only: dp
 use logmesh_m, only: lm_origin

!Input variables
 real(dp), intent(in)  :: ff(:), rr(:)

!Output variables
 real(dp), intent(out) :: f0(:)

!Local variables
 integer :: ii, npts, ierr

 npts = size(ff)
 if(size(f0) /= npts +1) error stop "nr /= npts + 1 in add_zero_r"
 do ii = 1, npts
   f0(ii+1) = ff(ii)
 end do
 call lm_origin(rr(1:npts),ff,npts,f0(1),ierr)
 if(ierr /= 0) error stop "add_zero_r: log mesh too short for lm_origin"
 return

 end subroutine add_zero_r


 subroutine write_psml_item(xf,class, &
&                           nn, ll, jj, ss, &
&                           seq, &
&                           rc, ekb, &
&                           flavor, type, set, &
&                           ff)

!! writes one radial function (a semi-local potential or a projector) as a
!! PSML element with its attributes
!!
!! @param xf  xml file handle
!! @param class  element name, "slps" or "proj"
!! @param nn  principal quantum number (optional)
!! @param ll  angular momentum (optional)
!! @param jj  total angular momentum (optional)
!! @param ss  spin, +0.5 or -0.5 (optional)
!! @param seq  index of the projector for this l (optional)
!! @param rc  core radius (bohr, optional)
!! @param ekb  projector energy coefficient (Ha, optional)
!! @param flavor  generation method (optional)
!! @param type  projector type (optional)
!! @param set  set the function belongs to (optional)
!! @param ff  the function on the grid, including r=0

 use xmlf90_wxml, only: xmlf_t, str, xml_NewElement, xml_EndElement, &
&                       xml_AddArray
 use input_echo_m, only: fe
 use constants_m, only: dp

!Input variables
 character(len=*), intent(in)  :: class

 integer, intent(in), optional  :: nn
 integer, intent(in), optional  :: ll
 real(dp), intent(in), optional  :: jj
 real(dp), intent(in), optional  :: ss

! for sl potentials
 real(dp), intent(in), optional  :: rc
 character(len=*), intent(in), optional  :: flavor

! for projectors
 integer, intent(in), optional   :: seq
 real(dp), intent(in), optional  :: ekb
 character(len=*), intent(in), optional  :: type

 character(len=*), intent(in), optional  :: set
 real(dp), intent(in)  :: ff(:)

!Input/output variables
 type(xmlf_t), intent(inout)   :: xf

 call xml_NewElement(xf,trim(class))
 if(present(set))  call my_add_attribute(xf,"set",set)

! we might want to check input values
 if(present(nn))  call my_add_attribute(xf,"n",str(nn))
 if(present(ll))  call my_add_attribute(xf,"l",lsymb(ll))
 if(present(jj))  call my_add_attribute(xf,"j", &
&                                        str(jj,format="(f3.1)"))
! spin: +0.5 or -0.5
 if(present(ss))  call my_add_attribute(xf,"s", &
&                                        str(ss,format="(f4.1)"))

 if(present(seq)) call my_add_attribute(xf,"seq",str(seq))

 if(present(rc))  call my_add_attribute(xf,"rc",str(rc))
 if(present(ekb))  call my_add_attribute(xf,"ekb",str(ekb))

 if(present(flavor))  call my_add_attribute(xf,"flavor",flavor)
 if(present(type))  call my_add_attribute(xf,"type",type)

 call xml_NewElement(xf,"radfunc")
 call xml_NewElement(xf,"data")
 call xml_AddArray(xf,fe(ff(:)))
 call xml_EndElement(xf,"data")
 call xml_EndElement(xf,"radfunc")

 call xml_EndElement(xf,trim(class))
 return

 end subroutine write_psml_item

 end module psmlout_m
